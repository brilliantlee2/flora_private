use std::cmp::Ordering;
use std::collections::HashMap;
use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::{SystemTime, UNIX_EPOCH};

use anyhow::{bail, Context, Result};
use clap::Parser;

#[derive(Clone, Debug, Parser)]
#[command(
    name = "flora glycine",
    version,
    about = "Parallel Glycine processing for one or more FASTQ files"
)]
struct BatchArgs {
    #[arg(
        short = 'f',
        long = "fastq",
        num_args = 1..,
        conflicts_with = "fastq_dir",
        required_unless_present = "fastq_dir"
    )]
    fastq: Vec<PathBuf>,
    #[arg(
        long = "fastq-dir",
        conflicts_with = "fastq",
        required_unless_present = "fastq"
    )]
    fastq_dir: Option<PathBuf>,
    #[arg(short = '5', long = "tso_seq", required = true)]
    tso_seq: String,
    #[arg(short = '3', long = "rtp_seq", required = true)]
    rtp_seq: String,
    #[arg(short = 'o', long = "outdir", required = true)]
    outdir: PathBuf,
    #[arg(short = 'n', long = "sample", required = true)]
    sample: String,
    #[arg(short = 'e', long = "err", default_value = "0.25,0.25")]
    err: String,
    #[arg(short = 's', long = "shift", default_value = "100,100")]
    shift: String,
    #[arg(short = 'L', long = "min_len", default_value_t = 100)]
    min_len: usize,
    #[arg(short = 'Q', long = "min_qual", default_value_t = 7.0)]
    min_qual: f64,
    #[arg(short = 'u', long = "trim_len", default_value_t = 0)]
    trim_len: usize,
    #[arg(short = 'l', long = "tail_len", default_value_t = 10)]
    tail_len: usize,
    #[arg(short = 'q', long = "umi_len", default_value_t = 0)]
    umi_len: usize,
    #[arg(long = "jobs", default_value_t = 10)]
    jobs: usize,
    #[arg(
        short = 't',
        long = "total-threads",
        alias = "thread",
        default_value_t = 64
    )]
    total_threads: usize,
    #[arg(long = "keep-all-outputs")]
    keep_all_outputs: bool,
}

#[derive(Clone)]
struct WorkItem {
    index: usize,
    input: PathBuf,
    outdir: PathBuf,
    sample: String,
}

pub fn run_from<I, T>(args: I) -> Result<()>
where
    I: IntoIterator<Item = T>,
    T: Into<OsString> + Clone,
{
    let args = BatchArgs::parse_from(args);
    let inputs = resolve_inputs(&args.fastq, args.fastq_dir.as_deref())?;
    if args.jobs == 0 || args.total_threads == 0 {
        bail!("--jobs and --total-threads must both be greater than zero");
    }
    fs::create_dir_all(&args.outdir)?;
    let active_jobs = inputs.len().min(args.jobs).min(args.total_threads);
    let base_threads = (args.total_threads / active_jobs).max(1);
    let extra_threads = args.total_threads % active_jobs;
    eprintln!(
        "[Glycine batch] inputs={} concurrent_jobs={} total_threads={} threads_per_job={}..{}",
        inputs.len(),
        active_jobs,
        args.total_threads,
        base_threads,
        base_threads + usize::from(extra_threads > 0)
    );

    let stamp = SystemTime::now().duration_since(UNIX_EPOCH)?.as_secs();
    let chunk_root = args.outdir.join(format!(
        ".flora-glycine-chunks-{}-{stamp}",
        std::process::id()
    ));
    fs::create_dir_all(&chunk_root)?;
    let queue = Arc::new(Mutex::new(
        inputs
            .iter()
            .enumerate()
            .map(|(index, input)| WorkItem {
                index,
                input: input.clone(),
                outdir: chunk_root.join(format!("chunk_{index:06}")),
                sample: format!("{}__chunk_{index:06}", args.sample),
            })
            .collect::<Vec<_>>()
            .into_iter(),
    ));

    let mut handles = Vec::with_capacity(active_jobs);
    for worker_id in 0..active_jobs {
        let queue = Arc::clone(&queue);
        let common = args.clone();
        let worker_threads = base_threads + usize::from(worker_id < extra_threads);
        handles.push(thread::spawn(move || -> Result<Vec<WorkItem>> {
            let mut completed = Vec::new();
            loop {
                let item = queue.lock().unwrap().next();
                let Some(item) = item else { break };
                fs::create_dir_all(&item.outdir)?;
                eprintln!(
                    "[Glycine batch] worker={} input={} threads={}",
                    worker_id + 1,
                    item.input.display(),
                    worker_threads
                );
                let mut argv = vec![
                    OsString::from("flora glycine worker"),
                    OsString::from("--fastq"),
                    item.input.as_os_str().to_owned(),
                    OsString::from("--tso_seq"),
                    OsString::from(&common.tso_seq),
                    OsString::from("--rtp_seq"),
                    OsString::from(&common.rtp_seq),
                    OsString::from("--outdir"),
                    item.outdir.as_os_str().to_owned(),
                    OsString::from("--sample"),
                    OsString::from(&item.sample),
                    OsString::from("--err"),
                    OsString::from(&common.err),
                    OsString::from("--shift"),
                    OsString::from(&common.shift),
                    OsString::from("--min_len"),
                    OsString::from(common.min_len.to_string()),
                    OsString::from("--min_qual"),
                    OsString::from(common.min_qual.to_string()),
                    OsString::from("--trim_len"),
                    OsString::from(common.trim_len.to_string()),
                    OsString::from("--tail_len"),
                    OsString::from(common.tail_len.to_string()),
                    OsString::from("--umi_len"),
                    OsString::from(common.umi_len.to_string()),
                    OsString::from("--thread"),
                    OsString::from(worker_threads.to_string()),
                ];
                if common.keep_all_outputs {
                    argv.push(OsString::from("--keep-all-outputs"));
                }
                super::run_from(argv)
                    .with_context(|| format!("Glycine failed for {}", item.input.display()))?;
                completed.push(item);
            }
            Ok(completed)
        }));
    }

    let mut completed = Vec::with_capacity(inputs.len());
    for handle in handles {
        completed.extend(
            handle
                .join()
                .map_err(|_| anyhow::anyhow!("Glycine worker panicked"))??,
        );
    }
    completed.sort_by_key(|item| item.index);
    merge_outputs(
        &completed,
        &args.outdir,
        &args.sample,
        args.keep_all_outputs,
    )?;
    fs::remove_dir_all(&chunk_root)?;
    Ok(())
}

fn resolve_inputs(explicit: &[PathBuf], directory: Option<&Path>) -> Result<Vec<PathBuf>> {
    let mut inputs = if let Some(directory) = directory {
        if !directory.is_dir() {
            bail!("FASTQ directory does not exist: {}", directory.display());
        }
        fs::read_dir(directory)?
            .filter_map(|entry| entry.ok().map(|entry| entry.path()))
            .filter(|path| path.is_file() && is_fastq_gz(path))
            .collect::<Vec<_>>()
    } else {
        explicit.to_vec()
    };
    if inputs.is_empty() {
        bail!("no FASTQ inputs found; use --fastq or --fastq-dir");
    }
    for path in &inputs {
        if !path.is_file() {
            bail!("FASTQ input does not exist: {}", path.display());
        }
    }
    inputs = inputs
        .into_iter()
        .map(|path| fs::canonicalize(&path).with_context(|| format!("resolve {}", path.display())))
        .collect::<Result<Vec<_>>>()?;
    inputs.sort_by(|left, right| natural_cmp(&left.to_string_lossy(), &right.to_string_lossy()));
    if let Some(duplicate) = inputs.windows(2).find(|pair| pair[0] == pair[1]) {
        bail!("duplicate FASTQ input: {}", duplicate[0].display());
    }
    Ok(inputs)
}

fn is_fastq_gz(path: &Path) -> bool {
    let name = path
        .file_name()
        .and_then(|v| v.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    name.ends_with(".fastq.gz") || name.ends_with(".fq.gz")
}

fn natural_cmp(left: &str, right: &str) -> Ordering {
    let mut a = left.chars().peekable();
    let mut b = right.chars().peekable();
    loop {
        match (a.peek().copied(), b.peek().copied()) {
            (Some(ac), Some(bc)) if ac.is_ascii_digit() && bc.is_ascii_digit() => {
                let mut an = String::new();
                let mut bn = String::new();
                while a.peek().is_some_and(|c| c.is_ascii_digit()) {
                    an.push(a.next().unwrap());
                }
                while b.peek().is_some_and(|c| c.is_ascii_digit()) {
                    bn.push(b.next().unwrap());
                }
                let av = an.trim_start_matches('0');
                let bv = bn.trim_start_matches('0');
                let order = av
                    .len()
                    .cmp(&bv.len())
                    .then_with(|| av.cmp(bv))
                    .then_with(|| an.len().cmp(&bn.len()));
                if order != Ordering::Equal {
                    return order;
                }
            }
            (Some(ac), Some(bc)) => {
                a.next();
                b.next();
                let order = ac.cmp(&bc);
                if order != Ordering::Equal {
                    return order;
                }
            }
            (None, None) => return Ordering::Equal,
            (None, Some(_)) => return Ordering::Less,
            (Some(_), None) => return Ordering::Greater,
        }
    }
}

fn merge_outputs(items: &[WorkItem], outdir: &Path, sample: &str, keep_all: bool) -> Result<()> {
    let suffixes: &[&str] = if keep_all {
        &[
            "full-length.fq.gz",
            "non-full-length.fq.gz",
            "discarded.fq.gz",
            "failed-filter.fq.gz",
            "rescued.fq.gz",
            "full-length-plus-rescued.fq.gz",
        ]
    } else {
        &["full-length-plus-rescued.fq.gz"]
    };
    for suffix in suffixes {
        let final_path = outdir.join(format!("{sample}.{suffix}"));
        let temporary = outdir.join(format!(".{sample}.{suffix}.tmp"));
        let mut writer = File::create(&temporary)?;
        for item in items {
            let path = item.outdir.join(format!("{}.{suffix}", item.sample));
            io::copy(&mut File::open(&path)?, &mut writer)?;
        }
        writer.flush()?;
        fs::rename(temporary, final_path)?;
    }
    let stats = items
        .iter()
        .map(|item| {
            item.outdir
                .join(format!("{}.identifying_statistic.txt", item.sample))
        })
        .collect::<Vec<_>>();
    let merged = merge_statistics(&stats)?;
    let final_stats = outdir.join(format!("{sample}.identifying_statistic.txt"));
    let temporary_stats = outdir.join(format!(".{sample}.identifying_statistic.txt.tmp"));
    fs::write(&temporary_stats, merged)?;
    fs::rename(temporary_stats, final_stats)?;
    Ok(())
}

#[derive(Default)]
struct Stats {
    total_bases: u64,
    valid_bases: u64,
    summary: HashMap<String, u64>,
    non_chimeric: HashMap<String, u64>,
    chimeric_reads: u64,
    rescued_reads: u64,
    chimeric: HashMap<String, u64>,
}

fn merge_statistics(paths: &[PathBuf]) -> Result<String> {
    let mut total = Stats::default();
    for path in paths {
        let text = fs::read_to_string(path)?;
        add_statistics(&mut total, &text).with_context(|| format!("parse {}", path.display()))?;
    }
    Ok(render_statistics(&total))
}

fn add_statistics(total: &mut Stats, text: &str) -> Result<()> {
    #[derive(Clone, Copy)]
    enum Section {
        Summary,
        NonChimeric,
        Chimeric,
    }
    let lines = text.lines().collect::<Vec<_>>();
    if lines.len() < 3 {
        bail!("incomplete identifying statistics");
    }
    let bases = lines[2].split('\t').collect::<Vec<_>>();
    total.total_bases += bases
        .first()
        .context("missing total base count")?
        .parse::<u64>()?;
    total.valid_bases += bases
        .get(1)
        .context("missing valid base count")?
        .parse::<u64>()?;
    let mut section = Section::Summary;
    let mut i = 3;
    while i < lines.len() {
        let line = lines[i].trim();
        if line.is_empty() || line.starts_with("Type\t") {
            i += 1;
            continue;
        }
        if line == "Non-chimeric" {
            section = Section::NonChimeric;
            i += 1;
            continue;
        }
        if line == "Chimeric" {
            section = Section::Chimeric;
            i += 2;
            let values = lines
                .get(i)
                .context("missing chimeric summary")?
                .split('\t')
                .collect::<Vec<_>>();
            total.chimeric_reads += values
                .first()
                .context("missing chimeric count")?
                .parse::<u64>()?;
            total.rescued_reads += values
                .get(1)
                .context("missing rescued count")?
                .parse::<u64>()?;
            i += 1;
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() >= 2 {
            let value = fields[1].parse::<u64>()?;
            let map = match section {
                Section::Summary => &mut total.summary,
                Section::NonChimeric => &mut total.non_chimeric,
                Section::Chimeric => &mut total.chimeric,
            };
            *map.entry(fields[0].to_string()).or_default() += value;
        }
        i += 1;
    }
    Ok(())
}

fn pct(value: u64, denominator: u64) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        value as f64 * 100.0 / denominator as f64
    }
}

fn render_statistics(stats: &Stats) -> String {
    const SUMMARY: &[&str] = &[
        "Total",
        "Length-filtered",
        "QC-filtered",
        "Full-length+rescued",
        "Chimeric",
        "Non-chimeric",
        "Full-length",
    ];
    const NON_CHIMERIC: &[&str] = &[
        "Total",
        "Full-length",
        "High-confidence strand:+",
        "High-confidence strand:-",
        "Low-confidence strand:+",
        "Low-confidence strand:-",
        "Non-full-length",
        "Double-rtp-double-polya/t",
        "Non-tso-polya-rtp",
        "Rtp-polyt-non-tso",
        "Tso-rtp-non-polya/t",
        "Rtp-tso-non-polya/t",
        "Single-5'-tso",
        "Single-5'-rtp",
        "Single-3'-tso",
        "Single-3'-rtp",
        "Discarded",
        "Double-tso",
        "Double-rtp",
        "Non-primer",
    ];
    const CHIMERIC: &[&str] = &[
        "Total",
        "Rescued strand:+",
        "Rescued strand:-",
        "Irrescuable",
    ];
    let total_reads = *stats.summary.get("Total").unwrap_or(&0);
    let non_total = *stats.non_chimeric.get("Total").unwrap_or(&0);
    let chim_total = *stats.chimeric.get("Total").unwrap_or(&0);
    let mut out = format!("Summary\nTotal_base_count\tValid_base_count\tValid_base_proportion(%)\n{}\t{}\t{:.2}\nType\tRead_count\tRead_proportion(%)", stats.total_bases, stats.valid_bases, pct(stats.valid_bases, stats.total_bases));
    for label in SUMMARY {
        let n = *stats.summary.get(*label).unwrap_or(&0);
        out.push_str(&format!("\n{label}\t{n}\t{:.2}", pct(n, total_reads)));
    }
    out.push_str("\n\nNon-chimeric\nType\tRead_count\tRead_proportion(%)");
    for label in NON_CHIMERIC {
        let n = *stats.non_chimeric.get(*label).unwrap_or(&0);
        out.push_str(&format!("\n{label}\t{n}\t{:.2}", pct(n, non_total)));
    }
    out.push_str(&format!("\n\nChimeric\nChimeric_read_count\tRescued_read_count\tRescued_read_proportion(%)\n{}\t{}\t{:.2}\nType\tRead_count\tRead_proportion(%)", stats.chimeric_reads, stats.rescued_reads, pct(stats.rescued_reads, stats.chimeric_reads)));
    for label in CHIMERIC {
        let n = *stats.chimeric.get(*label).unwrap_or(&0);
        out.push_str(&format!("\n{label}\t{n}\t{:.2}", pct(n, chim_total)));
    }
    out.push('\n');
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn natural_order_places_ten_after_two() {
        let mut names = vec![
            "sample_10.fastq.gz",
            "sample_2.fastq.gz",
            "sample_1.fastq.gz",
        ];
        names.sort_by(|a, b| natural_cmp(a, b));
        assert_eq!(
            names,
            vec![
                "sample_1.fastq.gz",
                "sample_2.fastq.gz",
                "sample_10.fastq.gz"
            ]
        );
    }

    #[test]
    fn merged_statistics_sum_counts_and_recompute_percentages() {
        let fixture = "Summary\nTotal_base_count\tValid_base_count\tValid_base_proportion(%)\n100\t50\t50.00\nType\tRead_count\tRead_proportion(%)\nTotal\t10\t100.00\nLength-filtered\t1\t10.00\nQC-filtered\t1\t10.00\nFull-length+rescued\t4\t40.00\nChimeric\t2\t20.00\nNon-chimeric\t6\t60.00\nFull-length\t3\t30.00\n\nNon-chimeric\nType\tRead_count\tRead_proportion(%)\nTotal\t6\t100.00\nFull-length\t3\t50.00\n\nChimeric\nChimeric_read_count\tRescued_read_count\tRescued_read_proportion(%)\n2\t1\t50.00\nType\tRead_count\tRead_proportion(%)\nTotal\t2\t100.00\nRescued strand:+\t1\t50.00\nRescued strand:-\t0\t0.00\nIrrescuable\t1\t50.00\n";
        let mut stats = Stats::default();
        add_statistics(&mut stats, fixture).unwrap();
        add_statistics(&mut stats, fixture).unwrap();
        let rendered = render_statistics(&stats);
        assert!(rendered.contains("200\t100\t50.00"));
        assert!(rendered.contains("Total\t20\t100.00"));
        assert!(rendered.contains("Full-length+rescued\t8\t40.00"));
        assert!(rendered.contains("4\t2\t50.00"));
    }

    #[test]
    fn directory_discovery_is_non_recursive_and_naturally_sorted() {
        let root = tempfile::tempdir().unwrap();
        fs::write(root.path().join("sample_10.fastq.gz"), b"").unwrap();
        fs::write(root.path().join("sample_2.fq.gz"), b"").unwrap();
        fs::write(root.path().join("notes.txt"), b"").unwrap();
        fs::create_dir(root.path().join("nested")).unwrap();
        fs::write(root.path().join("nested/ignored.fastq.gz"), b"").unwrap();

        let inputs = resolve_inputs(&[], Some(root.path())).unwrap();
        let names = inputs
            .iter()
            .map(|path| path.file_name().unwrap().to_string_lossy().into_owned())
            .collect::<Vec<_>>();
        assert_eq!(names, vec!["sample_2.fq.gz", "sample_10.fastq.gz"]);
    }

    #[test]
    fn duplicate_inputs_are_rejected_after_path_resolution() {
        let root = tempfile::tempdir().unwrap();
        let input = root.path().join("sample.fastq.gz");
        fs::write(&input, b"").unwrap();
        let error = resolve_inputs(&[input.clone(), input], None).unwrap_err();
        assert!(error.to_string().contains("duplicate FASTQ input"));
    }
}
