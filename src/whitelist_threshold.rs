use std::collections::BTreeMap;
use std::ffi::OsString;
use std::fs::{self, File};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use clap::{Args, Parser, ValueEnum};
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use serde::Serialize;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, ValueEnum)]
#[serde(rename_all = "snake_case")]
pub enum WhitelistMethod {
    #[default]
    Legacy,
    Steepest,
}

#[derive(Clone, Debug, Args, Serialize)]
pub struct ThresholdOptions {
    #[arg(long, value_enum, default_value = "legacy")]
    pub whitelist_method: WhitelistMethod,
    #[arg(long, default_value_t = 100)]
    pub whitelist_slope_window: usize,
    #[arg(long, default_value_t = 100)]
    pub whitelist_slope_min_count: usize,
    #[arg(long, default_value_t = 0.10)]
    pub whitelist_slope_max_rank_spread: f64,
}

impl Default for ThresholdOptions {
    fn default() -> Self {
        Self {
            whitelist_method: WhitelistMethod::Legacy,
            whitelist_slope_window: 100,
            whitelist_slope_min_count: 100,
            whitelist_slope_max_rank_spread: 0.10,
        }
    }
}

impl ThresholdOptions {
    pub fn validate(&self) -> Result<()> {
        if self.whitelist_slope_window == 0 || self.whitelist_slope_window.checked_mul(10).is_none()
        {
            bail!("--whitelist-slope-window must be positive and small enough to form complete windows");
        }
        let spread = self.whitelist_slope_max_rank_spread;
        if !spread.is_finite() || spread <= 0.0 || spread >= 1.0 {
            bail!("--whitelist-slope-max-rank-spread must be finite and strictly between 0 and 1");
        }
        Ok(())
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct SlopePoint {
    pub center_rank: usize,
    pub count: usize,
    pub slope: f64,
    pub window_count_drop: f64,
    pub left_rank: usize,
    pub right_rank: usize,
    pub left_count: usize,
    pub right_count: usize,
}

#[derive(Debug, Serialize)]
pub struct WindowDiagnostic {
    pub half_window: usize,
    pub status: &'static str,
    #[serde(flatten)]
    pub peak: Option<SlopePoint>,
    pub peak_to_median_fold: Option<f64>,
}

#[derive(Debug, Serialize)]
pub struct SteepestDiagnostic {
    pub threshold: Option<usize>,
    pub center_rank: Option<usize>,
    pub reason: &'static str,
    pub eligible_barcodes: usize,
    pub relative_rank_spread: Option<f64>,
    pub retained_including_ties: Option<usize>,
    pub windows: Vec<WindowDiagnostic>,
}

pub fn legacy_threshold(counts_desc: &[usize], exp_cells: usize) -> f64 {
    let n = counts_desc.len().min(exp_cells.max(1));
    if n == 0 {
        return 0.0;
    }
    // Same ascending top-N interpolation and operation order as main 291de65.
    let pos = (n - 1) as f64 * 0.95;
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    let value = if lo == hi {
        counts_desc[n - 1 - lo] as f64
    } else {
        let frac = pos - lo as f64;
        counts_desc[n - 1 - lo] as f64 * (1.0 - frac) + counts_desc[n - 1 - hi] as f64 * frac
    };
    value / 20.0
}

fn slope_profile(counts: &[usize], half: usize, min_count: usize) -> Vec<SlopePoint> {
    let n = counts
        .iter()
        .take_while(|&&count| count > min_count)
        .count();
    if half >= n || half >= n - half {
        return Vec::new();
    }
    let logs: Vec<_> = counts[..n]
        .iter()
        .map(|&count| (count as f64).ln())
        .collect();
    let ranks: Vec<_> = (1..=n).map(|rank| (rank as f64).ln()).collect();
    (half..n - half)
        .map(|center| {
            let (left, right) = (center - half, center + half);
            SlopePoint {
                center_rank: center + 1,
                count: counts[center],
                slope: (logs[right] - logs[left]) / (ranks[right] - ranks[left]),
                window_count_drop: counts[left] as f64 / counts[right] as f64,
                left_rank: left + 1,
                right_rank: right + 1,
                left_count: counts[left],
                right_count: counts[right],
            }
        })
        .collect()
}

fn median(values: &mut [f64]) -> f64 {
    values.sort_unstable_by(f64::total_cmp);
    let n = values.len();
    if n % 2 == 1 {
        values[n / 2]
    } else {
        (values[n / 2 - 1] + values[n / 2]) / 2.0
    }
}

pub fn steepest(counts_desc: &[usize], options: &ThresholdOptions) -> SteepestDiagnostic {
    let window = options.whitelist_slope_window;
    let floor = options.whitelist_slope_min_count;
    let mut out = SteepestDiagnostic {
        threshold: None,
        center_rank: None,
        reason: "insufficient points for the requested main slope window",
        eligible_barcodes: counts_desc.iter().filter(|&&count| count > floor).count(),
        relative_rank_spread: None,
        retained_including_ties: None,
        windows: Vec::new(),
    };
    let mut halves = vec![(window / 2).max(1), window, window * 2, window * 5];
    halves.sort_unstable();
    halves.dedup();
    for half in halves {
        let profile = slope_profile(counts_desc, half, floor);
        let Some(point) = profile.iter().min_by(|a, b| {
            a.slope
                .total_cmp(&b.slope)
                .then(a.center_rank.cmp(&b.center_rank))
        }) else {
            out.windows.push(WindowDiagnostic {
                half_window: half,
                status: "insufficient_points",
                peak: None,
                peak_to_median_fold: None,
            });
            continue;
        };
        let mut magnitudes: Vec<_> = profile.iter().map(|point| -point.slope).collect();
        let fold = -point.slope / median(&mut magnitudes).max(1e-12);
        out.windows.push(WindowDiagnostic {
            half_window: half,
            status: "computed",
            peak: Some(point.clone()),
            peak_to_median_fold: Some(fold),
        });
    }
    let Some(main) = out
        .windows
        .iter()
        .find(|item| item.half_window == window && item.peak.is_some())
    else {
        return out;
    };
    let peak = main.peak.as_ref().unwrap();
    out.center_rank = Some(peak.center_rank);
    let ranks: Vec<_> = out
        .windows
        .iter()
        .filter_map(|item| item.peak.as_ref().map(|p| p.center_rank))
        .collect();
    let reason = if ranks.len() < 2 {
        Some("fewer than two slope windows can be compared")
    } else {
        let spread = (*ranks.iter().max().unwrap() - *ranks.iter().min().unwrap()) as f64
            / peak.center_rank as f64;
        out.relative_rank_spread = Some(spread);
        if spread > options.whitelist_slope_max_rank_spread {
            Some("slope peak location is inconsistent across windows")
        } else if peak.slope >= 0.0 || main.peak_to_median_fold.unwrap() < 3.0 {
            Some("local descent is not sufficiently stronger than the curve background")
        } else if peak.window_count_drop < 1.5 {
            Some("support changes too little across the main slope window")
        } else {
            let retained = counts_desc
                .iter()
                .filter(|&&count| count >= peak.count)
                .count();
            out.retained_including_ties = Some(retained);
            let tolerance = window.max(
                (options.whitelist_slope_max_rank_spread * peak.center_rank as f64).ceil() as usize,
            );
            if retained - peak.center_rank > tolerance {
                Some(
                    "cutoff falls on a broad tied plateau; automatic barcode boundary is ambiguous",
                )
            } else {
                None
            }
        }
    };
    if let Some(reason) = reason {
        out.reason = reason;
    } else {
        out.threshold = Some(peak.count);
        out.reason = "stable local steepest descent found";
    }
    out
}

#[derive(Debug, Serialize)]
pub struct SideReport {
    pub status: &'static str,
    pub observed_theoretical_barcodes: usize,
    pub legacy_threshold: f64,
    pub legacy_selected: usize,
    pub steepest_selected: Option<usize>,
    pub applied_threshold: Option<f64>,
    pub comparison: &'static str,
    pub selected: Option<usize>,
    pub steepest: SteepestDiagnostic,
}

#[derive(Debug, Serialize)]
pub struct ThresholdReport {
    pub schema_version: usize,
    pub method: WhitelistMethod,
    pub exp_cells: usize,
    pub options: ThresholdOptions,
    pub count_definition: &'static str,
    pub interpretation: &'static str,
    pub sides: BTreeMap<String, SideReport>,
}

pub fn evaluate_side(
    counts: &HashMap<String, usize>,
    full: &HashSet<String>,
    exp_cells: usize,
    options: &ThresholdOptions,
) -> SideReport {
    let mut values: Vec<_> = counts
        .iter()
        .filter(|(bc, _)| full.contains(*bc))
        .map(|(_, &count)| count)
        .collect();
    values.sort_unstable_by(|a, b| b.cmp(a));
    let legacy = legacy_threshold(&values, exp_cells);
    let diagnostic = steepest(&values, options);
    let steepest_selected = diagnostic
        .threshold
        .map(|threshold| values.iter().filter(|&&count| count >= threshold).count());
    let (applied, comparison) = match options.whitelist_method {
        WhitelistMethod::Legacy => (Some(legacy), ">"),
        WhitelistMethod::Steepest => (diagnostic.threshold.map(|count| count as f64), ">="),
    };
    let legacy_selected = values
        .iter()
        .filter(|&&count| count as f64 > legacy)
        .count();
    SideReport {
        status: if applied.is_some() {
            "ready"
        } else {
            "needs_manual_review"
        },
        observed_theoretical_barcodes: values.len(),
        legacy_threshold: legacy,
        legacy_selected,
        steepest_selected,
        applied_threshold: applied,
        comparison,
        selected: if options.whitelist_method == WhitelistMethod::Legacy {
            Some(legacy_selected)
        } else {
            steepest_selected
        },
        steepest: diagnostic,
    }
}

pub fn selected_barcodes(
    counts: &HashMap<String, usize>,
    full: &HashSet<String>,
    report: &SideReport,
) -> Vec<String> {
    let mut selected: Vec<_> = counts
        .iter()
        .filter(|(bc, count)| {
            full.contains(*bc)
                && report.applied_threshold.is_some_and(|threshold| {
                    if report.comparison == ">=" {
                        **count as f64 >= threshold
                    } else {
                        **count as f64 > threshold
                    }
                })
        })
        .map(|(bc, _)| bc.clone())
        .collect();
    selected.sort();
    selected
}

impl ThresholdReport {
    pub fn new(
        exp_cells: usize,
        options: &ThresholdOptions,
        side3: SideReport,
        side5: SideReport,
    ) -> Self {
        Self {
            schema_version: 1, method: options.whitelist_method, exp_cells, options: options.clone(),
            count_definition: "High-quality putative barcode read counts intersected with the oriented theoretical whitelist; pre-correction, pre-merge; 3p/5p independent",
            interpretation: "Numerical whitelist-selection heuristic, not final cell calls, cell probabilities or FDR",
            sides: BTreeMap::from([("3p".into(), side3), ("5p".into(), side5)]),
        }
    }

    pub fn ensure_ready(&self) -> Result<()> {
        let failed: Vec<_> = self
            .sides
            .iter()
            .filter(|(_, side)| side.applied_threshold.is_none())
            .map(|(name, side)| format!("{name}: {}", side.steepest.reason))
            .collect();
        if !failed.is_empty() {
            bail!("no reliable steepest whitelist threshold ({}); inspect barcode_thresholds.json; rerun explicitly with --whitelist-method legacy or review slope parameters in a NEW output directory", failed.join("; "));
        }
        Ok(())
    }

    pub fn write(&self, directory: &Path) -> Result<()> {
        serde_json::to_writer_pretty(
            File::create(directory.join("barcode_thresholds.json"))?,
            self,
        )?;
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(directory.join("barcode_thresholds.tsv"))?;
        writer.write_record([
            "side",
            "method",
            "status",
            "legacy_threshold",
            "legacy_selected",
            "steepest_threshold",
            "steepest_selected",
            "applied_threshold",
            "comparison",
            "selected",
            "reason",
        ])?;
        for (name, side) in &self.sides {
            writer.write_record([
                name.clone(),
                format!("{:?}", self.method).to_lowercase(),
                side.status.into(),
                side.legacy_threshold.to_string(),
                side.legacy_selected.to_string(),
                side.steepest
                    .threshold
                    .map(|n| n.to_string())
                    .unwrap_or_default(),
                side.steepest_selected
                    .map(|n| n.to_string())
                    .unwrap_or_default(),
                side.applied_threshold
                    .map(|n| n.to_string())
                    .unwrap_or_default(),
                side.comparison.into(),
                side.selected.map(|n| n.to_string()).unwrap_or_default(),
                side.steepest.reason.into(),
            ])?;
        }
        writer.flush()?;
        Ok(())
    }
}

#[derive(Parser)]
#[command(
    about = "Preview pre-merge barcode whitelist thresholds from existing count tables; no FASTQ or alignment is read"
)]
struct PreviewCli {
    #[arg(long)]
    counts_3p: PathBuf,
    #[arg(long)]
    counts_5p: PathBuf,
    #[arg(long)]
    full_whitelist_3p: PathBuf,
    #[arg(long)]
    full_whitelist_5p: PathBuf,
    #[arg(long)]
    no_revcomp_whitelist: bool,
    #[arg(long, default_value_t = 5000)]
    exp_cells: usize,
    #[arg(long)]
    out_dir: PathBuf,
    #[command(flatten)]
    options: ThresholdOptions,
}

fn read_counts(path: &Path) -> Result<HashMap<String, usize>> {
    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').from_path(path)?;
    if reader.headers()?.iter().collect::<Vec<_>>() != ["barcode", "count"] {
        bail!("{}: expected barcode<TAB>count header", path.display());
    }
    let mut counts = HashMap::default();
    for record in reader.records() {
        let record = record?;
        let barcode = record[0].trim();
        if barcode.is_empty() {
            bail!("{}: empty barcode", path.display());
        }
        let count: usize = record[1]
            .parse()
            .with_context(|| format!("{}: invalid count for {barcode}", path.display()))?;
        if counts.insert(barcode.into(), count).is_some() {
            bail!("{}: duplicate barcode {barcode}", path.display());
        }
    }
    Ok(counts)
}

pub fn run_preview(args: impl IntoIterator<Item = OsString>) -> Result<()> {
    let cli = PreviewCli::parse_from(args);
    cli.options.validate()?;
    if cli.out_dir.exists() {
        bail!("preview --out-dir must not exist; use a new directory");
    }
    let counts3 = read_counts(&cli.counts_3p)?;
    let counts5 = read_counts(&cli.counts_5p)?;
    let full3 = crate::pipeline::read_whitelist(&cli.full_whitelist_3p, !cli.no_revcomp_whitelist)?;
    let full5 = crate::pipeline::read_whitelist(&cli.full_whitelist_5p, !cli.no_revcomp_whitelist)?;
    let report = ThresholdReport::new(
        cli.exp_cells,
        &cli.options,
        evaluate_side(&counts3, &full3, cli.exp_cells, &cli.options),
        evaluate_side(&counts5, &full5, cli.exp_cells, &cli.options),
    );
    if let Some(parent) = cli
        .out_dir
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)?;
    }
    fs::create_dir(&cli.out_dir)?;
    report.write(&cli.out_dir)?;
    report.ensure_ready()?;
    for (side, counts, full) in [("3p", &counts3, &full3), ("5p", &counts5, &full5)] {
        let selected = selected_barcodes(counts, full, &report.sides[side]);
        fs::write(
            cli.out_dir.join(format!("whitelist_{side}.csv")),
            selected
                .iter()
                .map(|bc| format!("{bc}\n"))
                .collect::<String>(),
        )?;
        println!(
            "[barcode-thresholds] {side}: threshold={:?}, comparison={}, selected={}",
            report.sides[side].applied_threshold,
            report.sides[side].comparison,
            selected.len()
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn signal() -> Vec<usize> {
        (1..=10000)
            .map(|i| {
                (100000.0
                    * (-0.1 * (i as f64).ln()
                        - 4.0 / (1.0 + (-((i as f64 - 3000.0) / 100.0)).exp()))
                    .exp()) as usize
            })
            .collect()
    }

    #[test]
    fn signal_is_stable_and_uses_center_count() {
        let values = signal();
        let result = steepest(&values, &ThresholdOptions::default());
        let rank = result.center_rank.unwrap();
        assert!((2900..=3200).contains(&rank));
        assert_eq!(result.threshold, Some(values[rank - 1]));
        assert!(result.relative_rank_spread.unwrap() <= 0.10);
    }

    #[test]
    fn flat_short_zero_powerlaw_and_tied_plateaus_are_rejected() {
        for values in [
            vec![],
            vec![0; 1000],
            vec![1000; 1000],
            vec![1000, 900, 800],
            (1..=10000).map(|i| 10000000 / i).collect(),
            [vec![10000; 10000], vec![200; 20000]].concat(),
        ] {
            let diagnostic = steepest(&values, &ThresholdOptions::default());
            assert_eq!(diagnostic.threshold, None);
            serde_json::to_string(&diagnostic).unwrap();
        }
    }

    #[test]
    fn low_count_tail_cannot_shift_cutoff() {
        let mut values = signal();
        let before = steepest(&values, &ThresholdOptions::default());
        values.extend(vec![100; 1000]);
        values.extend((1..100).rev());
        let after = steepest(&values, &ThresholdOptions::default());
        assert_eq!(
            (before.threshold, before.center_rank),
            (after.threshold, after.center_rank)
        );
    }

    #[test]
    fn inconsistent_transition_locations_are_rejected() {
        let values: Vec<_> = (1..=10000)
            .map(|i| {
                let drop = 1.0 / (1.0 + (-(i as f64 - 2500.0) / 25.0).exp())
                    + 1.0 / (1.0 + (-(i as f64 - 6000.0) / 100.0).exp());
                (100000.0 * (-0.1 * (i as f64).ln() - drop).exp()) as usize
            })
            .collect();
        assert_eq!(
            steepest(&values, &ThresholdOptions::default()).threshold,
            None
        );
    }
}
