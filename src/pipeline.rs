use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Context, Result};
use csv::WriterBuilder;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use rayon::prelude::*;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use serde::Serialize;
use zip::ZipArchive;

use crate::barcode::{
    correct_one_side_indexed, extract_putative, revcomp_upper, reverse_complement, strip_fixed_3p,
    strip_fixed_5p, umi_a_ratio, BarcodeIndex, CorrectedRead, PutativeRow,
};
use crate::fastq::{for_each_fastq_batch, write_fastq_record};

#[derive(Clone, Debug)]
pub struct PipelineConfig {
    pub fastq_fns: Vec<PathBuf>,
    pub full_bc_whitelist_3p: PathBuf,
    pub full_bc_whitelist_5p: PathBuf,
    pub revcomp_whitelist: bool,
    pub out_dir: PathBuf,
    pub batch_size: usize,
    pub assign_batchsize: usize,
    pub min_q: i32,
    pub exp_cells: usize,
    pub max_ed: usize,
    pub pair_min: Option<usize>,
    pub auto_pair_min_floor: usize,
    pub auto_pair_min_quantile: f64,
    pub dominance_min: f64,
    pub drop_umi_a_ratio_gt: f64,
    pub top1_alpha: f64,
    pub top1_alpha_umi: f64,
    pub require_pass_both_ends: bool,
    pub include_other_components: bool,
    pub max_other_component_barcodes: usize,
    pub absorb_unassigned_paired: bool,
    pub min_reads_per_cell: usize,
    pub bc_fixed_3p: String,
    pub umi_fixed_3p: String,
    pub bc_fixed_5p: String,
    pub umi_fixed_5p: String,
    pub skip_matched_fastq: bool,
    pub skip_unmatched_fastq: bool,
    pub skip_cell_fastq: bool,
    pub save_intermediate: bool,
}

#[derive(Debug)]
pub struct PipelineSummary {
    pub fastq_files: usize,
    pub reads_total: usize,
    pub reads_demultiplexed: usize,
    pub barcode_validity_stats: BarcodeValidityStats,
    pub putative_rows_3p: usize,
    pub putative_rows_5p: usize,
    pub merged_rows: usize,
    pub filtered_rows: usize,
    pub trimmed_barcode_uniques: TrimmedBarcodeUniques,
    pub clean_reads_rows: usize,
    pub pair_counts_rows: usize,
    pub assigned_rows: usize,
    pub core_cells: usize,
    pub core_barcodes: usize,
    pub pair_stats: PairStats,
    pub assign_stats: AssignStats,
    pub cell_read_stats_head: Vec<CellReadStat>,
}

#[derive(Clone, Debug, Serialize)]
pub struct BarcodeValidityStats {
    pub total_rows: usize,
    pub valid_any_n: usize,
    pub valid_any_ratio: f64,
    pub only_3p_n: usize,
    pub only_3p_ratio: f64,
    pub only_5p_n: usize,
    pub only_5p_ratio: f64,
    pub both_n: usize,
    pub both_ratio: f64,
    pub neither_n: usize,
    pub neither_ratio: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct TrimmedBarcodeUniques {
    pub unique_bc3_20bp_rc: usize,
    pub unique_bc5_20bp: usize,
    pub unique_union_bc3_bc5: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct PairStats {
    #[serde(rename = "PAIR_MIN")]
    pub pair_min: usize,
    #[serde(rename = "PAIR_MIN_mode")]
    pub pair_min_mode: String,
    #[serde(rename = "TOP1_ALPHA")]
    pub top1_alpha: f64,
    #[serde(rename = "TOP1_ALPHA_UMI")]
    pub top1_alpha_umi: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct AssignStats {
    pub n_total_reads: usize,
    pub n_assigned_reads: usize,
    pub n_cells_before_min_reads: usize,
    pub n_cells_after_min_reads: usize,
    pub min_reads_per_cell: usize,
    pub n_core_barcodes: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct CellReadStat {
    pub cell_id: String,
    pub n_reads: usize,
}

#[derive(Clone, Debug, Serialize)]
struct CorrectionMapRow {
    read_id: String,
    side: String,
    putative_barcode: String,
    corrected_barcode: String,
    status: String,
    final_corrected_20bp: String,
}

#[derive(Clone, Debug, Serialize)]
struct CleanRead {
    read_id: String,
    putative_umi: String,
    putative_umi_5p: String,
    #[serde(rename = "BC5n")]
    bc5n: String,
    #[serde(rename = "BC3n")]
    bc3n: String,
}

#[derive(Clone, Debug, Serialize)]
struct PairCount {
    #[serde(rename = "BC5n")]
    bc5n: String,
    #[serde(rename = "BC3n")]
    bc3n: String,
    support_reads: usize,
    support_umis: usize,
}

#[derive(Clone, Debug)]
struct PairPruneStats {
    small_components_seen: usize,
    small_components_motif_kept: usize,
    edges_kept_by_motif: usize,
    edges_kept_by_fallback: usize,
}

#[derive(Clone, Debug)]
struct Top1Dominance {
    top1_partner: String,
    dominance: f64,
}

#[derive(Clone, Debug)]
struct EdgeRecord {
    u: String,
    v: String,
    support_reads: usize,
    support_umis: usize,
}

#[derive(Clone, Debug)]
struct GraphComponent {
    nodes: Vec<String>,
    edges: Vec<EdgeRecord>,
}

#[derive(Clone, Debug, Serialize)]
struct AssignedRead {
    read_id: String,
    putative_umi: String,
    putative_umi_5p: String,
    #[serde(rename = "BC5n")]
    bc5n: String,
    #[serde(rename = "BC3n")]
    bc3n: String,
    cell_id: String,
}

pub fn run_pipeline(config: &PipelineConfig) -> Result<PipelineSummary> {
    fs::create_dir_all(&config.out_dir)?;
    let bc_fixed_3p = reverse_complement(&config.bc_fixed_3p);
    let umi_fixed_3p = reverse_complement(&config.umi_fixed_3p);
    let bc_fixed_5p = reverse_complement(&config.bc_fixed_5p);
    let umi_fixed_5p = reverse_complement(&config.umi_fixed_5p);

    println!("[Flora] Step 1/7: streaming FASTQ and extracting putative 3p/5p barcode tables");
    let step_t0 = Instant::now();
    let putative = extract_putative_rows(
        &config.fastq_fns,
        config.batch_size,
        &bc_fixed_3p,
        &umi_fixed_3p,
        &bc_fixed_5p,
        &umi_fixed_5p,
    )?;
    if config.save_intermediate {
        write_csv(config.out_dir.join("putative_bc.csv"), &putative)?;
    }
    log_step_elapsed(1, 7, step_t0);

    println!("[Flora] Step 2/7: reading whitelists and counting high-quality barcodes");
    let step_t0 = Instant::now();
    let full_wl3 = read_whitelist(&config.full_bc_whitelist_3p, config.revcomp_whitelist)?;
    let full_wl5 = read_whitelist(&config.full_bc_whitelist_5p, config.revcomp_whitelist)?;
    let raw3 = count_high_quality(&putative, true, config.min_q);
    let raw5 = count_high_quality(&putative, false, config.min_q);
    write_counts(config.out_dir.join("barcode_counts_3p.tsv"), &raw3)?;
    write_counts(config.out_dir.join("barcode_counts_5p.tsv"), &raw5)?;
    let observed_wl3 = select_observed_whitelist(&raw3, &full_wl3, config.exp_cells);
    let observed_wl5 = select_observed_whitelist(&raw5, &full_wl5, config.exp_cells);
    write_lines(config.out_dir.join("whitelist_3p.csv"), &observed_wl3)?;
    write_lines(config.out_dir.join("whitelist_5p.csv"), &observed_wl5)?;
    write_lines(config.out_dir.join("empty_bc_list_3p.csv"), &[])?;
    write_lines(config.out_dir.join("empty_bc_list_5p.csv"), &[])?;
    log_step_elapsed(2, 7, step_t0);

    println!("[Flora] Step 3/7: correcting reads with unique-barcode caches");
    let step_t0 = Instant::now();
    let wl3: HashSet<String> = observed_wl3.into_iter().collect();
    let wl5: HashSet<String> = observed_wl5.into_iter().collect();
    let wl3_index = BarcodeIndex::new(wl3.into_iter());
    let wl5_index = BarcodeIndex::new(wl5.into_iter());
    let correction_cache_3p =
        build_correction_cache(&putative, true, &wl3_index, config.max_ed, config.min_q);
    let correction_cache_5p =
        build_correction_cache(&putative, false, &wl5_index, config.max_ed, config.min_q);
    let corrected: Vec<CorrectedRead> = putative
        .par_iter()
        .map(|row| {
            correct_from_caches(
                row,
                &correction_cache_3p,
                &correction_cache_5p,
                config.min_q,
            )
        })
        .collect();
    if config.save_intermediate {
        write_corrected(config.out_dir.join("BC_corrected.csv"), &corrected)?;
    }
    write_correction_map(
        config.out_dir.join("correction_map_3p.tsv"),
        "3p",
        &putative,
        &corrected,
    )?;
    write_correction_map(
        config.out_dir.join("correction_map_5p.tsv"),
        "5p",
        &putative,
        &corrected,
    )?;
    let barcode_validity_stats = summarize_barcode_validity(&corrected);
    log_step_elapsed(3, 7, step_t0);

    println!("[Flora] Step 4/7: filtering corrected reads and stripping fixed barcode sequence");
    let step_t0 = Instant::now();
    let corrected_filtered = filter_corrected_rows(
        &corrected,
        config.pair_min,
        config.auto_pair_min_floor,
        config.auto_pair_min_quantile,
    );
    let filtered_rows = corrected_filtered.len();
    let clean_pre_pair = build_clean_reads(&corrected_filtered);
    let trimmed_barcode_uniques = summarize_trimmed_barcode_uniques(&clean_pre_pair);
    let (clean, pair_counts, pair_stats, _prune_stats) = filter_pairs_three_stage(
        &clean_pre_pair,
        config.pair_min,
        config.auto_pair_min_floor,
        config.auto_pair_min_quantile,
        config.top1_alpha,
        config.top1_alpha_umi,
        config.require_pass_both_ends,
        config.drop_umi_a_ratio_gt,
    );
    log_step_elapsed(4, 7, step_t0);

    println!("[Flora] Step 5/7: writing clean read and pair tables");
    let step_t0 = Instant::now();
    write_csv(config.out_dir.join("ReadIDs_UMI_BC_clean.csv"), &clean)?;
    write_csv(config.out_dir.join("pair_counts_kept.csv"), &pair_counts)?;
    log_step_elapsed(5, 7, step_t0);

    println!("[Flora] Step 6/7: assigning reads to graph components");
    let step_t0 = Instant::now();
    let barcode_to_cell = assign_cells(
        &pair_counts,
        config.include_other_components,
        config.max_other_component_barcodes,
    );
    let top1_dominance = compute_top1_dominance(&pair_counts);
    let assigned_all = assign_reads(
        &clean,
        &barcode_to_cell,
        &top1_dominance,
        config.dominance_min,
        config.absorb_unassigned_paired,
    );
    let n_cells_before_min_reads = assigned_all
        .iter()
        .map(|r| r.cell_id.as_str())
        .collect::<HashSet<_>>()
        .len();
    let assigned = filter_min_reads(assigned_all, config.min_reads_per_cell);
    write_csv(config.out_dir.join("read_assigned_cell.csv"), &assigned)?;
    let cell_read_stats = collect_cell_stats(&assigned);
    write_cell_stats(config.out_dir.join("cell_read_stats.csv"), &cell_read_stats)?;
    let kept_cells: HashSet<&str> = assigned.iter().map(|r| r.cell_id.as_str()).collect();
    write_barcode_to_cell(
        config.out_dir.join("barcode_to_cell.csv"),
        &barcode_to_cell,
        &kept_cells,
    )?;
    let core_barcodes = barcode_to_cell
        .iter()
        .filter(|(_, cell)| kept_cells.contains(cell.as_str()))
        .count();
    let assign_stats = AssignStats {
        n_total_reads: clean.len(),
        n_assigned_reads: assigned.len(),
        n_cells_before_min_reads,
        n_cells_after_min_reads: kept_cells.len(),
        min_reads_per_cell: config.min_reads_per_cell,
        n_core_barcodes: core_barcodes,
    };
    write_assign_stats(config.out_dir.join("assign_stats.tsv"), &assign_stats)?;
    write_barcode_validity(
        config.out_dir.join("barcode_validity_summary.tsv"),
        &barcode_validity_stats,
    )?;
    log_step_elapsed(6, 7, step_t0);

    println!("[Flora] Step 7/7: optional FASTQ outputs");
    let step_t0 = Instant::now();
    if !config.skip_matched_fastq {
        let assigned_ids: HashSet<String> = assigned.iter().map(|r| r.read_id.clone()).collect();
        write_filtered_fastq_gz(
            &config.fastq_fns,
            config.batch_size,
            config.out_dir.join("matched_reads.fastq.gz"),
            &assigned_ids,
            true,
        )?;
    }
    if !config.skip_unmatched_fastq {
        let assigned_ids: HashSet<String> = assigned.iter().map(|r| r.read_id.clone()).collect();
        write_filtered_fastq_gz(
            &config.fastq_fns,
            config.batch_size,
            config.out_dir.join("unmatched_reads.fastq.gz"),
            &assigned_ids,
            false,
        )?;
    }
    if !config.skip_cell_fastq {
        let assigned_ids: HashSet<String> = assigned.iter().map(|r| r.read_id.clone()).collect();
        write_filtered_fastq_gz(
            &config.fastq_fns,
            config.batch_size,
            config.out_dir.join("cell_reads.fastq.gz"),
            &assigned_ids,
            true,
        )?;
    }
    log_step_elapsed(7, 7, step_t0);

    Ok(PipelineSummary {
        fastq_files: config.fastq_fns.len(),
        reads_total: putative.len(),
        reads_demultiplexed: corrected
            .iter()
            .filter(|r| !r.bc3_corrected.is_empty() || !r.bc5_corrected.is_empty())
            .count(),
        barcode_validity_stats,
        putative_rows_3p: putative.len(),
        putative_rows_5p: putative.len(),
        merged_rows: putative.len(),
        filtered_rows,
        trimmed_barcode_uniques,
        clean_reads_rows: clean.len(),
        pair_counts_rows: pair_counts.len(),
        assigned_rows: assigned.len(),
        core_cells: kept_cells.len(),
        core_barcodes,
        pair_stats,
        assign_stats,
        cell_read_stats_head: cell_read_stats.into_iter().take(5).collect(),
    })
}

fn log_step_elapsed(step: usize, total_steps: usize, started: Instant) {
    println!(
        "[timing] Step {step}/{total_steps} elapsed: {:.2}s",
        started.elapsed().as_secs_f64()
    );
}

fn extract_putative_rows(
    paths: &[PathBuf],
    batch_size: usize,
    bc_fixed_3p: &str,
    umi_fixed_3p: &str,
    bc_fixed_5p: &str,
    umi_fixed_5p: &str,
) -> Result<Vec<PutativeRow>> {
    let mut putative = Vec::new();
    for_each_fastq_batch(paths, batch_size, |batch| {
        let mut rows: Vec<PutativeRow> = batch
            .par_iter()
            .map(|rec| extract_putative(rec, bc_fixed_3p, umi_fixed_3p, bc_fixed_5p, umi_fixed_5p))
            .collect();
        putative.append(&mut rows);
        Ok(())
    })?;
    Ok(putative)
}

fn write_csv<T: Serialize>(path: PathBuf, rows: &[T]) -> Result<()> {
    let mut writer =
        csv::Writer::from_path(&path).with_context(|| format!("write {}", path.display()))?;
    for row in rows {
        writer.serialize(row)?;
    }
    writer.flush()?;
    Ok(())
}

fn write_filtered_fastq_gz(
    paths: &[PathBuf],
    batch_size: usize,
    path: PathBuf,
    target_ids: &HashSet<String>,
    keep_target: bool,
) -> Result<()> {
    let file = File::create(path)?;
    let mut encoder = GzEncoder::new(file, Compression::default());
    for_each_fastq_batch(paths, batch_size, |batch| {
        for rec in &batch {
            let is_target = target_ids.contains(&rec.id);
            if is_target == keep_target {
                write_fastq_record(&mut encoder, rec)?;
            }
        }
        Ok(())
    })?;
    encoder.finish()?;
    Ok(())
}

fn write_corrected(path: PathBuf, rows: &[CorrectedRead]) -> Result<()> {
    let mut writer = csv::Writer::from_path(path)?;
    writer.write_record([
        "read_id",
        "putative_umi",
        "strand",
        "BC3_corrected",
        "BC5_corrected",
        "putative_umi_5p",
    ])?;
    for r in rows {
        writer.write_record([
            &r.read_id,
            &r.putative_umi,
            "+",
            &r.bc3_corrected,
            &r.bc5_corrected,
            &r.putative_umi_5p,
        ])?;
    }
    writer.flush()?;
    Ok(())
}

fn correction_status(putative: &str, corrected: &str) -> String {
    if corrected.trim().is_empty() {
        "failed".to_string()
    } else if putative == corrected {
        "exact".to_string()
    } else {
        "corrected".to_string()
    }
}

fn final_corrected_20bp_for_side(side: &str, corrected: &str) -> String {
    if corrected.trim().is_empty() {
        return String::new();
    }
    match side {
        "3p" => revcomp_upper(&strip_fixed_3p(corrected)),
        "5p" => strip_fixed_5p(corrected),
        _ => String::new(),
    }
}

fn write_correction_map(
    path: PathBuf,
    side: &str,
    putative: &[PutativeRow],
    corrected: &[CorrectedRead],
) -> Result<()> {
    let rows = putative
        .iter()
        .zip(corrected.iter())
        .map(|(p, c)| {
            let (putative_barcode, corrected_barcode) = if side == "3p" {
                (p.putative_bc.clone(), c.bc3_corrected.clone())
            } else {
                (p.putative_bc_5p.clone(), c.bc5_corrected.clone())
            };
            CorrectionMapRow {
                read_id: p.read_id.clone(),
                side: side.to_string(),
                putative_barcode: putative_barcode.clone(),
                corrected_barcode: corrected_barcode.clone(),
                status: correction_status(&putative_barcode, &corrected_barcode),
                final_corrected_20bp: final_corrected_20bp_for_side(side, &corrected_barcode),
            }
        })
        .collect::<Vec<_>>();
    write_csv(path, &rows)
}

fn read_whitelist(path: &Path, revcomp: bool) -> Result<HashSet<String>> {
    let mut values = Vec::new();
    if path.extension().and_then(|x| x.to_str()) == Some("zip") {
        let file = File::open(path)?;
        let mut archive = ZipArchive::new(file)?;
        let mut member = archive.by_index(0)?;
        let mut text = String::new();
        member.read_to_string(&mut text)?;
        values.extend(text.lines().map(str::to_string));
    } else {
        let reader: Box<dyn Read> = if path.extension().and_then(|x| x.to_str()) == Some("gz") {
            Box::new(MultiGzDecoder::new(File::open(path)?))
        } else {
            Box::new(File::open(path)?)
        };
        for line in BufReader::new(reader).lines() {
            values.push(line?);
        }
    }
    Ok(values
        .into_iter()
        .map(|x| x.split('-').next().unwrap_or("").trim().to_string())
        .filter(|x| !x.is_empty())
        .map(|x| if revcomp { revcomp_upper(&x) } else { x })
        .collect())
}

fn count_high_quality(rows: &[PutativeRow], is_3p: bool, min_q: i32) -> HashMap<String, usize> {
    let mut counts = HashMap::default();
    for row in rows {
        let (bc, q) = if is_3p {
            (&row.putative_bc, row.putative_bc_min_qs)
        } else {
            (&row.putative_bc_5p, row.putative_bc_min_qs_5p)
        };
        if q.is_some_and(|v| v >= min_q) && !bc.trim().is_empty() {
            *counts.entry(bc.clone()).or_insert(0) += 1;
        }
    }
    counts
}

fn putative_bc_for_side(row: &PutativeRow, is_3p: bool) -> (&str, Option<i32>) {
    if is_3p {
        (&row.putative_bc, row.putative_bc_min_qs)
    } else {
        (&row.putative_bc_5p, row.putative_bc_min_qs_5p)
    }
}

fn build_correction_cache(
    rows: &[PutativeRow],
    is_3p: bool,
    whitelist: &BarcodeIndex,
    max_ed: usize,
    min_q: i32,
) -> HashMap<String, String> {
    let unique: HashSet<String> = rows
        .iter()
        .filter_map(|row| {
            let (bc, q) = putative_bc_for_side(row, is_3p);
            if q.is_some_and(|v| v < min_q) || bc.trim().is_empty() {
                None
            } else {
                Some(bc.to_string())
            }
        })
        .collect();

    unique
        .par_iter()
        .map(|bc| {
            let corrected = if whitelist.contains(bc) {
                bc.clone()
            } else {
                correct_one_side_indexed(bc, whitelist, max_ed)
            };
            (bc.clone(), corrected)
        })
        .collect()
}

fn correct_from_caches(
    row: &PutativeRow,
    cache3: &HashMap<String, String>,
    cache5: &HashMap<String, String>,
    min_q: i32,
) -> CorrectedRead {
    let c3 = if row.putative_bc_min_qs.is_some_and(|q| q < min_q) {
        String::new()
    } else {
        cache3.get(&row.putative_bc).cloned().unwrap_or_default()
    };
    let c5 = if row.putative_bc_min_qs_5p.is_some_and(|q| q < min_q) {
        String::new()
    } else {
        cache5.get(&row.putative_bc_5p).cloned().unwrap_or_default()
    };
    CorrectedRead {
        read_id: row.read_id.clone(),
        putative_umi: if c3.is_empty() {
            String::new()
        } else {
            row.putative_umi.clone()
        },
        putative_umi_5p: if c5.is_empty() {
            String::new()
        } else {
            row.putative_umi_5p.clone()
        },
        bc3_corrected: c3,
        bc5_corrected: c5,
    }
}

fn numpy_linear_quantile(values_sorted_asc: &[usize], q: f64) -> f64 {
    if values_sorted_asc.is_empty() {
        return 0.0;
    }
    if values_sorted_asc.len() == 1 {
        return values_sorted_asc[0] as f64;
    }
    let q = q.clamp(0.0, 1.0);
    let pos = (values_sorted_asc.len() - 1) as f64 * q;
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    if lo == hi {
        return values_sorted_asc[lo] as f64;
    }
    let frac = pos - lo as f64;
    values_sorted_asc[lo] as f64 * (1.0 - frac) + values_sorted_asc[hi] as f64 * frac
}

fn select_observed_whitelist(
    counts: &HashMap<String, usize>,
    full: &HashSet<String>,
    exp_cells: usize,
) -> Vec<String> {
    let filtered: Vec<(&String, usize)> = counts
        .iter()
        .filter(|(bc, _)| full.contains(*bc))
        .map(|(bc, n)| (bc, *n))
        .collect();
    if filtered.is_empty() {
        return Vec::new();
    }

    let mut vals: Vec<usize> = filtered.iter().map(|(_, n)| *n).collect();
    vals.sort_unstable();
    let top_start = vals.len().saturating_sub(exp_cells.max(1));
    let threshold = numpy_linear_quantile(&vals[top_start..], 0.95) / 20.0;
    let mut selected: Vec<_> = counts
        .iter()
        .filter(|(bc, n)| full.contains(*bc) && **n as f64 > threshold)
        .map(|(bc, _)| bc.clone())
        .collect();
    selected.sort();
    selected
}

fn write_counts(path: PathBuf, counts: &HashMap<String, usize>) -> Result<()> {
    let mut rows: Vec<_> = counts.iter().collect();
    rows.sort_by(|a, b| b.1.cmp(a.1));
    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record(["barcode", "count"])?;
    for (bc, count) in rows {
        writer.write_record([bc, &count.to_string()])?;
    }
    writer.flush()?;
    Ok(())
}

fn write_lines(path: PathBuf, lines: &[String]) -> Result<()> {
    let mut file = File::create(path)?;
    for line in lines {
        writeln!(file, "{line}")?;
    }
    Ok(())
}

fn filter_corrected_rows(
    rows: &[CorrectedRead],
    pair_min: Option<usize>,
    auto_pair_min_floor: usize,
    auto_pair_min_quantile: f64,
) -> Vec<CorrectedRead> {
    let paired_keys = rows
        .iter()
        .filter(|r| !r.bc3_corrected.is_empty() && !r.bc5_corrected.is_empty())
        .map(|r| {
            ordered_pair(
                &strip_fixed_5p(&r.bc5_corrected),
                &revcomp_upper(&strip_fixed_3p(&r.bc3_corrected)),
            )
        })
        .collect::<Vec<_>>();
    let mut pair_counts: HashMap<(String, String), usize> = HashMap::default();
    for key in paired_keys {
        *pair_counts.entry(key).or_insert(0) += 1;
    }
    let (resolved_pair_min, _) = resolve_pair_min(
        pair_counts.values().copied().collect(),
        pair_min,
        auto_pair_min_floor,
        auto_pair_min_quantile,
    );
    let mut bad_read_ids = HashSet::default();
    for r in rows {
        if r.bc3_corrected.is_empty() || r.bc5_corrected.is_empty() {
            continue;
        }
        let key = ordered_pair(
            &strip_fixed_5p(&r.bc5_corrected),
            &revcomp_upper(&strip_fixed_3p(&r.bc3_corrected)),
        );
        if pair_counts.get(&key).copied().unwrap_or(0) < resolved_pair_min {
            bad_read_ids.insert(r.read_id.clone());
        }
    }
    rows.iter()
        .filter(|r| !bad_read_ids.contains(&r.read_id))
        .filter(|r| !r.bc3_corrected.is_empty() || !r.bc5_corrected.is_empty())
        .filter(|r| !(r.putative_umi.trim().is_empty() && r.putative_umi_5p.trim().is_empty()))
        .cloned()
        .collect()
}

fn build_clean_reads(rows: &[CorrectedRead]) -> Vec<CleanRead> {
    rows.iter()
        .filter_map(|r| {
            let bc5 = strip_fixed_5p(&r.bc5_corrected);
            let bc3 = revcomp_upper(&strip_fixed_3p(&r.bc3_corrected));
            Some(CleanRead {
                read_id: r.read_id.clone(),
                putative_umi: r.putative_umi.clone(),
                putative_umi_5p: r.putative_umi_5p.clone(),
                bc5n: bc5,
                bc3n: bc3,
            })
        })
        .collect()
}

fn ordered_pair(a: &str, b: &str) -> (String, String) {
    if a <= b {
        (a.to_string(), b.to_string())
    } else {
        (b.to_string(), a.to_string())
    }
}

fn canonical_edge_umi_key(r: &CleanRead) -> String {
    if r.bc5n == r.bc3n {
        if r.putative_umi <= r.putative_umi_5p {
            format!("{}|{}", r.putative_umi, r.putative_umi_5p)
        } else {
            format!("{}|{}", r.putative_umi_5p, r.putative_umi)
        }
    } else if r.bc5n <= r.bc3n {
        format!("{}|{}", r.putative_umi_5p, r.putative_umi)
    } else {
        format!("{}|{}", r.putative_umi, r.putative_umi_5p)
    }
}

fn resolve_pair_min(
    mut counts: Vec<usize>,
    pair_min: Option<usize>,
    floor: usize,
    q: f64,
) -> (usize, String) {
    if let Some(v) = pair_min {
        return (v, "manual".to_string());
    }
    if counts.is_empty() {
        return (floor, "auto_empty".to_string());
    }
    counts.sort_unstable();
    let quantile = numpy_linear_quantile(&counts, q);
    (floor.max(quantile as usize), "auto_quantile".to_string())
}

fn build_pair_counts(
    reads: &[CleanRead],
    pair_min: Option<usize>,
    floor: usize,
    q: f64,
    top1_alpha: f64,
    top1_alpha_umi: f64,
    require_pass_both_ends: bool,
) -> (Vec<PairCount>, PairStats, PairPruneStats) {
    let mut reads_by_pair: HashMap<(String, String), usize> = HashMap::default();
    let mut umis_by_pair: HashMap<(String, String), HashSet<String>> = HashMap::default();
    for r in reads {
        if r.bc5n.is_empty() || r.bc3n.is_empty() {
            continue;
        }
        let key = ordered_pair(&r.bc5n, &r.bc3n);
        *reads_by_pair.entry(key.clone()).or_insert(0) += 1;
        umis_by_pair
            .entry(key)
            .or_default()
            .insert(canonical_edge_umi_key(r));
    }
    let all_rows: Vec<_> = reads_by_pair
        .into_iter()
        .map(|((a, b), support_reads)| PairCount {
            bc5n: a.clone(),
            bc3n: b.clone(),
            support_reads,
            support_umis: umis_by_pair.get(&(a, b)).map_or(0, HashSet::len),
        })
        .collect();
    let (min_support, pair_min_mode) = resolve_pair_min(
        all_rows.iter().map(|x| x.support_reads).collect(),
        pair_min,
        floor,
        q,
    );
    let min_kept: Vec<_> = all_rows
        .into_iter()
        .filter(|row| row.support_reads >= min_support)
        .collect();
    let top1_reads_map = top1_map_from_metric(&min_kept, true);
    let top1_umis_map = top1_map_from_metric(&min_kept, false);
    let (mut rows, prune_stats) = prune_edges_structure_aware(
        &min_kept,
        &top1_reads_map,
        &top1_umis_map,
        top1_alpha,
        top1_alpha_umi,
        require_pass_both_ends,
    );
    rows.sort_by(|a, b| b.support_reads.cmp(&a.support_reads));
    (
        rows,
        PairStats {
            pair_min: min_support,
            pair_min_mode,
            top1_alpha,
            top1_alpha_umi,
        },
        prune_stats,
    )
}

fn filter_pairs_three_stage(
    reads: &[CleanRead],
    pair_min: Option<usize>,
    floor: usize,
    q: f64,
    top1_alpha: f64,
    top1_alpha_umi: f64,
    require_pass_both_ends: bool,
    drop_umi_a_ratio_gt: f64,
) -> (Vec<CleanRead>, Vec<PairCount>, PairStats, PairPruneStats) {
    let cleaned = reads
        .iter()
        .filter(|r| !r.bc5n.is_empty() || !r.bc3n.is_empty())
        .filter(|r| umi_a_ratio(&r.putative_umi) <= drop_umi_a_ratio_gt)
        .cloned()
        .collect::<Vec<_>>();
    let (paired, single) = cleaned
        .into_iter()
        .partition::<Vec<_>, _>(|r| !r.bc5n.is_empty() && !r.bc3n.is_empty());
    let (pair_counts, pair_stats, prune_stats) = build_pair_counts(
        &paired,
        pair_min,
        floor,
        q,
        top1_alpha,
        top1_alpha_umi,
        require_pass_both_ends,
    );
    let kept_pairs: HashSet<(String, String)> = pair_counts
        .iter()
        .map(|p| ordered_pair(&p.bc5n, &p.bc3n))
        .collect();
    let paired_final = paired
        .into_iter()
        .filter(|r| kept_pairs.contains(&ordered_pair(&r.bc5n, &r.bc3n)))
        .collect::<Vec<_>>();
    let mut final_reads = single;
    final_reads.extend(paired_final);
    (final_reads, pair_counts, pair_stats, prune_stats)
}

fn assign_cells(
    pairs: &[PairCount],
    include_other_components: bool,
    max_other_component_barcodes: usize,
) -> HashMap<String, String> {
    let components = graph_components(pairs);
    let mut core_components = components
        .into_iter()
        .filter(|component| {
            matches!(
                component_category(component),
                "self_only" | "pair_only" | "triangle_only" | "clique4_only" | "clique5_only"
            ) || (include_other_components
                && component_category(component) == "other"
                && component.nodes.len() <= max_other_component_barcodes)
        })
        .collect::<Vec<_>>();
    core_components.sort_by(|a, b| a.nodes.cmp(&b.nodes));
    let mut out = HashMap::default();
    for (idx, component) in core_components.into_iter().enumerate() {
        let cell = format!("cell_{:06}", idx + 1);
        for bc in component.nodes {
            out.insert(bc, cell.clone());
        }
    }
    out
}

fn edge_key(a: &str, b: &str) -> (String, String) {
    if a <= b {
        (a.to_string(), b.to_string())
    } else {
        (b.to_string(), a.to_string())
    }
}

fn graph_components(rows: &[PairCount]) -> Vec<GraphComponent> {
    let mut nodes = HashSet::default();
    let mut adjacency: HashMap<String, Vec<String>> = HashMap::default();
    let mut edges_by_key: HashMap<(String, String), EdgeRecord> = HashMap::default();
    for row in rows {
        let (u, v) = edge_key(&row.bc5n, &row.bc3n);
        nodes.insert(u.clone());
        nodes.insert(v.clone());
        adjacency.entry(u.clone()).or_default();
        adjacency.entry(v.clone()).or_default();
        if u != v {
            adjacency.entry(u.clone()).or_default().push(v.clone());
            adjacency.entry(v.clone()).or_default().push(u.clone());
        }
        edges_by_key.insert(
            (u.clone(), v.clone()),
            EdgeRecord {
                u,
                v,
                support_reads: row.support_reads,
                support_umis: row.support_umis,
            },
        );
    }

    let mut visited = HashSet::default();
    let mut ordered_nodes = nodes.into_iter().collect::<Vec<_>>();
    ordered_nodes.sort();
    let mut components = Vec::new();
    for node in ordered_nodes {
        if !visited.insert(node.clone()) {
            continue;
        }
        let mut stack = vec![node.clone()];
        let mut comp_nodes = vec![node];
        while let Some(curr) = stack.pop() {
            if let Some(neighbors) = adjacency.get(&curr) {
                for next in neighbors {
                    if visited.insert(next.clone()) {
                        stack.push(next.clone());
                        comp_nodes.push(next.clone());
                    }
                }
            }
        }
        comp_nodes.sort();
        let comp_set: HashSet<&str> = comp_nodes.iter().map(String::as_str).collect();
        let mut comp_edges = edges_by_key
            .values()
            .filter(|edge| comp_set.contains(edge.u.as_str()) && comp_set.contains(edge.v.as_str()))
            .cloned()
            .collect::<Vec<_>>();
        comp_edges.sort_by(|a, b| a.u.cmp(&b.u).then_with(|| a.v.cmp(&b.v)));
        components.push(GraphComponent {
            nodes: comp_nodes,
            edges: comp_edges,
        });
    }
    components
}

fn component_category(component: &GraphComponent) -> &'static str {
    let n = component.nodes.len();
    let m = component.edges.len();
    let m_self = component
        .edges
        .iter()
        .filter(|edge| edge.u == edge.v)
        .count();
    let m_no_self = m - m_self;
    if n == 1 && m == 1 && m_self == 1 {
        "self_only"
    } else if n == 2 && m_no_self == 1 {
        "pair_only"
    } else if n == 3 && m_no_self == 3 {
        "triangle_only"
    } else if n == 4 && m_no_self == 6 {
        "clique4_only"
    } else if n == 5 && m_no_self == 10 {
        "clique5_only"
    } else {
        "other"
    }
}

fn top1_map_from_metric(rows: &[PairCount], read_metric: bool) -> HashMap<String, usize> {
    let mut out = HashMap::default();
    for row in rows {
        let value = if read_metric {
            row.support_reads
        } else {
            row.support_umis
        };
        *out.entry(row.bc5n.clone()).or_insert(0) =
            out.get(&row.bc5n).copied().unwrap_or(0).max(value);
        *out.entry(row.bc3n.clone()).or_insert(0) =
            out.get(&row.bc3n).copied().unwrap_or(0).max(value);
    }
    out
}

fn endpoint_pass_with_umi(
    edge_reads: usize,
    edge_umis: usize,
    top_reads: usize,
    top_umis: usize,
    alpha_reads: f64,
    alpha_umis: f64,
) -> bool {
    let read_rel = if top_reads == 0 {
        0.0
    } else {
        edge_reads as f64 / top_reads as f64
    };
    let umi_rel = if top_umis == 0 {
        0.0
    } else {
        edge_umis as f64 / top_umis as f64
    };
    read_rel >= alpha_reads
        || (umi_rel >= alpha_umis && read_rel >= 0.5 * alpha_reads)
        || (read_rel >= 0.75 * alpha_reads && umi_rel >= 0.5 * alpha_umis)
}

fn median_usize(values: &[usize]) -> f64 {
    let mid = values.len() / 2;
    if values.len() % 2 == 0 {
        (values[mid - 1] + values[mid]) as f64 / 2.0
    } else {
        values[mid] as f64
    }
}

fn motif_component_passes(component: &GraphComponent) -> bool {
    let category = component_category(component);
    if category == "self_only" {
        return true;
    }
    let self_edges = component
        .edges
        .iter()
        .filter(|e| e.u == e.v)
        .collect::<Vec<_>>();
    let pair_edges = component
        .edges
        .iter()
        .filter(|e| e.u != e.v)
        .collect::<Vec<_>>();
    let mut self_map: HashMap<&str, &EdgeRecord> = HashMap::default();
    for edge in &self_edges {
        self_map.insert(edge.u.as_str(), *edge);
    }
    let pairwise_fraction = |read_metric: bool| -> f64 {
        let pair_sum: usize = pair_edges
            .iter()
            .map(|e| {
                if read_metric {
                    e.support_reads
                } else {
                    e.support_umis
                }
            })
            .sum();
        let self_sum: usize = self_edges
            .iter()
            .map(|e| {
                if read_metric {
                    e.support_reads
                } else {
                    e.support_umis
                }
            })
            .sum();
        let total = pair_sum + self_sum;
        if total == 0 {
            0.0
        } else {
            pair_sum as f64 / total as f64
        }
    };
    let weakest_to_median = |read_metric: bool| -> f64 {
        let mut vals = pair_edges
            .iter()
            .map(|e| {
                if read_metric {
                    e.support_reads
                } else {
                    e.support_umis
                }
            })
            .filter(|v| *v > 0)
            .collect::<Vec<_>>();
        if vals.is_empty() {
            return 0.0;
        }
        vals.sort_unstable();
        let med = median_usize(&vals);
        if med <= 0.0 {
            0.0
        } else {
            vals[0] as f64 / med
        }
    };
    let node_pair_fraction = |node: &str, read_metric: bool| -> f64 {
        let pair_sum: usize = pair_edges
            .iter()
            .filter(|e| e.u == node || e.v == node)
            .map(|e| {
                if read_metric {
                    e.support_reads
                } else {
                    e.support_umis
                }
            })
            .sum();
        let self_val = self_map
            .get(node)
            .map(|e| {
                if read_metric {
                    e.support_reads
                } else {
                    e.support_umis
                }
            })
            .unwrap_or(0);
        let denom = pair_sum + self_val;
        if denom == 0 {
            0.0
        } else {
            pair_sum as f64 / denom as f64
        }
    };
    let nth_incident_pair_to_self = |node: &str, read_metric: bool, nth_largest: usize| -> f64 {
        let mut vals = pair_edges
            .iter()
            .filter(|e| e.u == node || e.v == node)
            .map(|e| {
                if read_metric {
                    e.support_reads
                } else {
                    e.support_umis
                }
            })
            .collect::<Vec<_>>();
        vals.sort_unstable_by(|a, b| b.cmp(a));
        if vals.len() <= nth_largest {
            return 0.0;
        }
        let self_val = self_map
            .get(node)
            .map(|e| {
                if read_metric {
                    e.support_reads
                } else {
                    e.support_umis
                }
            })
            .unwrap_or(0)
            .max(1);
        vals[nth_largest] as f64 / self_val as f64
    };

    match category {
        "pair_only" => {
            if pair_edges.len() != 1 {
                return false;
            }
            let edge = pair_edges[0];
            let per_node = [&edge.u, &edge.v]
                .iter()
                .map(|node| {
                    let self_r = self_map
                        .get(node.as_str())
                        .map(|e| e.support_reads)
                        .unwrap_or(0);
                    let self_u = self_map
                        .get(node.as_str())
                        .map(|e| e.support_umis)
                        .unwrap_or(0);
                    let read_frac = if edge.support_reads + self_r == 0 {
                        0.0
                    } else {
                        edge.support_reads as f64 / (edge.support_reads + self_r) as f64
                    };
                    let umi_frac = if edge.support_umis + self_u == 0 {
                        0.0
                    } else {
                        edge.support_umis as f64 / (edge.support_umis + self_u) as f64
                    };
                    read_frac.max(umi_frac)
                })
                .collect::<Vec<_>>();
            per_node.into_iter().fold(f64::INFINITY, f64::min) >= 0.35 && edge.support_umis >= 2
        }
        "triangle_only" => {
            pair_edges.len() == 3
                && pairwise_fraction(true).max(pairwise_fraction(false)) >= 0.55
                && weakest_to_median(true).max(weakest_to_median(false)) >= 0.25
                && component.nodes.iter().all(|node| {
                    node_pair_fraction(node, true).max(node_pair_fraction(node, false)) >= 0.45
                })
        }
        "clique4_only" => {
            pair_edges.len() == 6
                && pairwise_fraction(true).max(pairwise_fraction(false)) >= 0.65
                && weakest_to_median(true).max(weakest_to_median(false)) >= 0.18
                && component.nodes.iter().all(|node| {
                    nth_incident_pair_to_self(node, true, 1)
                        .max(nth_incident_pair_to_self(node, false, 1))
                        >= 0.25
                })
        }
        "clique5_only" => {
            pair_edges.len() == 10
                && pairwise_fraction(true).max(pairwise_fraction(false)) >= 0.72
                && weakest_to_median(true).max(weakest_to_median(false)) >= 0.15
                && component.nodes.iter().all(|node| {
                    nth_incident_pair_to_self(node, true, 2)
                        .max(nth_incident_pair_to_self(node, false, 2))
                        >= 0.18
                })
        }
        _ => false,
    }
}

fn prune_edges_structure_aware(
    rows: &[PairCount],
    top1_reads_map: &HashMap<String, usize>,
    top1_umis_map: &HashMap<String, usize>,
    top1_alpha: f64,
    top1_alpha_umi: f64,
    require_pass_both_ends: bool,
) -> (Vec<PairCount>, PairPruneStats) {
    let components = graph_components(rows);
    let mut keep = HashSet::default();
    let mut stats = PairPruneStats {
        small_components_seen: 0,
        small_components_motif_kept: 0,
        edges_kept_by_motif: 0,
        edges_kept_by_fallback: 0,
    };
    for component in components {
        let category = component_category(&component);
        let use_motif = matches!(
            category,
            "self_only" | "pair_only" | "triangle_only" | "clique4_only" | "clique5_only"
        );
        if use_motif {
            stats.small_components_seen += 1;
        }
        if use_motif && motif_component_passes(&component) {
            stats.small_components_motif_kept += 1;
            for edge in &component.edges {
                keep.insert(edge_key(&edge.u, &edge.v));
                stats.edges_kept_by_motif += 1;
            }
            continue;
        }
        for edge in &component.edges {
            if edge.u == edge.v {
                keep.insert(edge_key(&edge.u, &edge.v));
                stats.edges_kept_by_fallback += 1;
                continue;
            }
            let pass_u = endpoint_pass_with_umi(
                edge.support_reads,
                edge.support_umis,
                *top1_reads_map.get(&edge.u).unwrap_or(&0),
                *top1_umis_map.get(&edge.u).unwrap_or(&0),
                top1_alpha,
                top1_alpha_umi,
            );
            let pass_v = endpoint_pass_with_umi(
                edge.support_reads,
                edge.support_umis,
                *top1_reads_map.get(&edge.v).unwrap_or(&0),
                *top1_umis_map.get(&edge.v).unwrap_or(&0),
                top1_alpha,
                top1_alpha_umi,
            );
            if if require_pass_both_ends {
                pass_u && pass_v
            } else {
                pass_u || pass_v
            } {
                keep.insert(edge_key(&edge.u, &edge.v));
                stats.edges_kept_by_fallback += 1;
            }
        }
    }
    (
        rows.iter()
            .filter(|row| keep.contains(&edge_key(&row.bc5n, &row.bc3n)))
            .cloned()
            .collect(),
        stats,
    )
}

fn compute_top1_dominance(rows: &[PairCount]) -> HashMap<String, Top1Dominance> {
    let mut by_barcode: HashMap<String, Vec<(String, usize)>> = HashMap::default();
    for row in rows {
        by_barcode
            .entry(row.bc5n.clone())
            .or_default()
            .push((row.bc3n.clone(), row.support_reads));
        by_barcode
            .entry(row.bc3n.clone())
            .or_default()
            .push((row.bc5n.clone(), row.support_reads));
    }
    let mut out = HashMap::default();
    for (barcode, mut vals) in by_barcode {
        vals.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
        let sum_all: usize = vals.iter().map(|x| x.1).sum();
        if let Some((top1_partner, top1_w)) = vals.first() {
            out.insert(
                barcode,
                Top1Dominance {
                    top1_partner: top1_partner.clone(),
                    dominance: if sum_all == 0 {
                        0.0
                    } else {
                        *top1_w as f64 / sum_all as f64
                    },
                },
            );
        }
    }
    out
}

fn assign_reads(
    reads: &[CleanRead],
    barcode_to_cell: &HashMap<String, String>,
    dom: &HashMap<String, Top1Dominance>,
    dominance_min: f64,
    absorb_unassigned_paired: bool,
) -> Vec<AssignedRead> {
    reads
        .iter()
        .filter_map(|r| {
            let cell_a_5 = barcode_to_cell.get(&r.bc5n).cloned();
            let cell_a_3 = barcode_to_cell.get(&r.bc3n).cloned();
            let cell_a = cell_a_5.clone().or(cell_a_3.clone());
            let has5 = !r.bc5n.is_empty();
            let has3 = !r.bc3n.is_empty();
            let try_absorb = |bc: &str| -> Option<String> {
                if bc.is_empty() {
                    return None;
                }
                let row = dom.get(bc)?;
                if row.dominance < dominance_min {
                    return None;
                }
                barcode_to_cell.get(&row.top1_partner).cloned()
            };
            let cell = if let Some(cell) = cell_a {
                Some(cell)
            } else if has5 ^ has3 {
                if has5 {
                    try_absorb(&r.bc5n)
                } else {
                    try_absorb(&r.bc3n)
                }
            } else if absorb_unassigned_paired && has5 && has3 {
                let cell_b5 = try_absorb(&r.bc5n);
                let cell_b3 = try_absorb(&r.bc3n);
                match (cell_b5, cell_b3) {
                    (Some(a), Some(b)) if a == b => Some(a),
                    (Some(a), None) => Some(a),
                    (None, Some(b)) => Some(b),
                    _ => None,
                }
            } else {
                None
            }?;
            Some(AssignedRead {
                read_id: r.read_id.clone(),
                putative_umi: r.putative_umi.clone(),
                putative_umi_5p: r.putative_umi_5p.clone(),
                bc5n: r.bc5n.clone(),
                bc3n: r.bc3n.clone(),
                cell_id: cell.clone(),
            })
        })
        .collect()
}

fn filter_min_reads(reads: Vec<AssignedRead>, min_reads: usize) -> Vec<AssignedRead> {
    let mut counts = HashMap::default();
    for r in &reads {
        *counts.entry(r.cell_id.clone()).or_insert(0usize) += 1;
    }
    reads
        .into_iter()
        .filter(|r| counts.get(&r.cell_id).copied().unwrap_or(0) >= min_reads)
        .collect()
}

fn collect_cell_stats(reads: &[AssignedRead]) -> Vec<CellReadStat> {
    let mut counts: HashMap<String, usize> = HashMap::default();
    for r in reads {
        *counts.entry(r.cell_id.clone()).or_insert(0) += 1;
    }
    let mut rows: Vec<_> = counts
        .into_iter()
        .map(|(cell_id, n_reads)| CellReadStat { cell_id, n_reads })
        .collect();
    rows.sort_by(|a, b| b.n_reads.cmp(&a.n_reads));
    rows
}

fn write_cell_stats(path: PathBuf, rows: &[CellReadStat]) -> Result<()> {
    let mut writer = csv::Writer::from_path(path)?;
    writer.write_record(["cell_id", "n_reads"])?;
    for row in rows {
        writer.write_record([row.cell_id.as_str(), &row.n_reads.to_string()])?;
    }
    writer.flush()?;
    Ok(())
}

fn write_barcode_to_cell(
    path: PathBuf,
    barcode_to_cell: &HashMap<String, String>,
    kept_cells: &HashSet<&str>,
) -> Result<()> {
    let mut cell_to_barcodes: HashMap<String, Vec<String>> = HashMap::default();
    for (barcode, cell) in barcode_to_cell {
        if kept_cells.contains(cell.as_str()) {
            cell_to_barcodes
                .entry(cell.clone())
                .or_default()
                .push(barcode.clone());
        }
    }

    let mut rows = cell_to_barcodes.into_iter().collect::<Vec<_>>();
    rows.sort_by(|a, b| a.0.cmp(&b.0));

    let mut writer = csv::Writer::from_path(path)?;
    writer.write_record(["cell", "barcode", "is_cell_barcode"])?;
    for (idx, (cell, mut barcodes)) in rows.into_iter().enumerate() {
        barcodes.sort();
        barcodes.dedup();
        let cell_number = cell_number_from_id(&cell).unwrap_or(idx + 1);
        let display_cell = format!("CELL{}_N{}", cell_number, barcodes.len());
        let barcode_joined = barcodes.join(";");
        writer.write_record([display_cell.as_str(), barcode_joined.as_str(), "1"])?;
    }
    writer.flush()?;
    Ok(())
}

fn cell_number_from_id(cell: &str) -> Option<usize> {
    cell.strip_prefix("cell_")?.parse::<usize>().ok()
}

fn write_assign_stats(path: PathBuf, stats: &AssignStats) -> Result<()> {
    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.serialize(stats)?;
    writer.flush()?;
    Ok(())
}

fn summarize_barcode_validity(rows: &[CorrectedRead]) -> BarcodeValidityStats {
    let total = rows.len();
    let both = rows
        .iter()
        .filter(|r| !r.bc3_corrected.is_empty() && !r.bc5_corrected.is_empty())
        .count();
    let only3 = rows
        .iter()
        .filter(|r| !r.bc3_corrected.is_empty() && r.bc5_corrected.is_empty())
        .count();
    let only5 = rows
        .iter()
        .filter(|r| r.bc3_corrected.is_empty() && !r.bc5_corrected.is_empty())
        .count();
    let neither = total - both - only3 - only5;
    let ratio = |n: usize| {
        if total == 0 {
            0.0
        } else {
            n as f64 / total as f64
        }
    };
    BarcodeValidityStats {
        total_rows: total,
        valid_any_n: total - neither,
        valid_any_ratio: ratio(total - neither),
        only_3p_n: only3,
        only_3p_ratio: ratio(only3),
        only_5p_n: only5,
        only_5p_ratio: ratio(only5),
        both_n: both,
        both_ratio: ratio(both),
        neither_n: neither,
        neither_ratio: ratio(neither),
    }
}

fn write_barcode_validity(path: PathBuf, stats: &BarcodeValidityStats) -> Result<()> {
    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.serialize(stats)?;
    writer.flush()?;
    Ok(())
}

fn summarize_trimmed_barcode_uniques(reads: &[CleanRead]) -> TrimmedBarcodeUniques {
    let bc3: HashSet<&str> = reads
        .iter()
        .filter_map(|r| (!r.bc3n.is_empty()).then_some(r.bc3n.as_str()))
        .collect();
    let bc5: HashSet<&str> = reads
        .iter()
        .filter_map(|r| (!r.bc5n.is_empty()).then_some(r.bc5n.as_str()))
        .collect();
    let mut union: HashSet<&str> = bc3.clone();
    union.extend(bc5.iter().copied());
    TrimmedBarcodeUniques {
        unique_bc3_20bp_rc: bc3.len(),
        unique_bc5_20bp: bc5.len(),
        unique_union_bc3_bc5: union.len(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pair_min_uses_floor_for_empty_or_low_counts() {
        assert_eq!(resolve_pair_min(vec![], None, 10, 0.1).0, 10);
        assert_eq!(resolve_pair_min(vec![1, 2, 3], None, 10, 0.1).0, 10);
        assert_eq!(resolve_pair_min(vec![1, 2, 3], Some(2), 10, 0.1).0, 2);
    }

    #[test]
    fn connected_pairs_assign_same_cell() {
        let pairs = vec![
            PairCount {
                bc5n: "A".into(),
                bc3n: "B".into(),
                support_reads: 10,
                support_umis: 1,
            },
            PairCount {
                bc5n: "B".into(),
                bc3n: "C".into(),
                support_reads: 10,
                support_umis: 1,
            },
        ];
        let map = assign_cells(&pairs, true, 8);
        assert_eq!(map["A"], map["C"]);
    }

    #[test]
    fn correction_cache_deduplicates_repeated_putative_barcodes() {
        let rows = vec![
            PutativeRow {
                putative_bc: "AAAA".into(),
                putative_bc_min_qs: Some(40),
                ..PutativeRow::default()
            },
            PutativeRow {
                putative_bc: "AAAA".into(),
                putative_bc_min_qs: Some(40),
                ..PutativeRow::default()
            },
        ];
        let whitelist = BarcodeIndex::new(["AAAA".to_string()]);
        let cache = build_correction_cache(&rows, true, &whitelist, 1, 2);
        assert_eq!(cache.len(), 1);
        assert_eq!(cache["AAAA"], "AAAA");
    }
}
