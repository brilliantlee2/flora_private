use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Parser;
use csv::{ReaderBuilder, StringRecord, WriterBuilder};
use flora::fastq::for_each_fastq_batch;
use rust_htslib::bam::{self, Read};
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

#[derive(Debug, Parser)]
#[command(version, about = "Compute RNA QC metrics and per-cell summaries")]
struct Cli {
    #[arg(long = "cell-umi-gene-tsv", default_value = "cell_umi_gene.tsv")]
    cell_umi_gene_tsv: PathBuf,

    #[arg(long, default_value = "filtered.sorted.bam")]
    bam: PathBuf,

    #[arg(long = "raw-fastq")]
    raw_fastq: PathBuf,

    #[arg(long = "full-length-fastq")]
    full_length_fastq: Option<PathBuf>,

    #[arg(long = "read-tags")]
    read_tags: Option<PathBuf>,

    #[arg(long = "barcode-validity-tsv")]
    barcode_validity_tsv: Option<PathBuf>,

    #[arg(long = "transcript-assigns")]
    transcript_assigns: Option<PathBuf>,

    #[arg(long = "glycine-log")]
    glycine_log: Option<PathBuf>,

    #[arg(long = "glycine-stats")]
    glycine_stats: Option<PathBuf>,

    #[arg(long = "mixed-species", default_value_t = false)]
    mixed_species: bool,

    #[arg(long = "fastq-count-file")]
    fastq_count_file: Option<PathBuf>,

    #[arg(long, default_value_t = 4)]
    threads: usize,
}

#[derive(Default)]
struct CellStats {
    read_ids: HashSet<u32>,
    known_umis: HashSet<u32>,
    known_genes: HashSet<u32>,
    mito_rows: usize,
}

struct CellUmiGeneRow {
    read_id: String,
    gene: String,
    barcode: String,
    umi: String,
}

#[derive(Default)]
struct BamSummary {
    bam_records: usize,
    mapped_primary_alignments: usize,
    unmapped_reads: usize,
    supplementary_alignments: usize,
    mapped_unique_reads: usize,
    aligned_genome_reads_in_final_cells: usize,
}

pub fn main() -> Result<()> {
    let cli = Cli::parse();
    let summarized_full_length_reads = cli
        .fastq_count_file
        .as_ref()
        .map(load_fastq_count)
        .transpose()?;
    let raw_is_full_length = cli
        .full_length_fastq
        .as_ref()
        .map(|path| same_file(&cli.raw_fastq, path))
        .transpose()?
        .unwrap_or(false);
    let glycine_raw_reads = cli
        .glycine_stats
        .as_ref()
        .map(|path| load_glycine_summary_read_count(path, "Total"))
        .transpose()?
        .flatten();

    // In multi-FASTQ Glycine runs, raw_fastq points at the merged full-length
    // output; the Glycine summary is the authoritative raw input count.
    let raw_fastq_reads = if let Some(count) = glycine_raw_reads {
        count
    } else if raw_is_full_length {
        match summarized_full_length_reads {
            Some(count) => count,
            None => count_fastq_reads(&cli.raw_fastq)?,
        }
    } else {
        count_fastq_reads(&cli.raw_fastq)?
    };
    let mut full_length_reads = if let Some(count) = summarized_full_length_reads {
        count
    } else if let Some(path) = &cli.full_length_fastq {
        if raw_is_full_length {
            raw_fastq_reads
        } else {
            count_fastq_reads(path)?
        }
    } else {
        raw_fastq_reads
    };
    if cli.mixed_species {
        if let Some(path) = &cli.glycine_log {
            if let Some(value) = load_glycine_total_full_length_reads(path)? {
                full_length_reads = value;
            }
        }
    }

    let cell_umi_gene = load_cell_umi_gene(&cli.cell_umi_gene_tsv, cli.mixed_species)?;
    let estimated_cells = cell_umi_gene.per_cell.len();
    let cell_associated_reads = cell_umi_gene.final_cell_reads.len();
    let assigned_cell_reads = match &cli.read_tags {
        Some(path) => count_read_tags(path)?,
        None => cell_associated_reads,
    };
    let barcode_valid_reads = match &cli.barcode_validity_tsv {
        Some(path) => load_barcode_valid_reads(path)?.unwrap_or(assigned_cell_reads),
        None => assigned_cell_reads,
    };
    let bam_summary =
        summarize_bam_alignments(&cli.bam, &cell_umi_gene.final_cell_reads, cli.threads)?;

    let (known_transcript_reads, unique_isoforms) = match &cli.transcript_assigns {
        Some(path) => load_transcript_assignment_summary(path, &cell_umi_gene.final_cell_reads)?,
        None => (None, None),
    };

    write_metrics_files(
        raw_fastq_reads,
        full_length_reads,
        estimated_cells,
        assigned_cell_reads,
        barcode_valid_reads,
        &cell_umi_gene,
        &bam_summary,
        known_transcript_reads,
        unique_isoforms,
    )?;
    write_per_cell_qc(&cell_umi_gene.per_cell)?;
    Ok(())
}

struct CellUmiGeneSummary {
    per_cell: HashMap<String, CellStats>,
    final_cell_reads: HashMap<String, u32>,
    known_gene_reads: usize,
    unique_genes: usize,
}

fn load_cell_umi_gene(path: &PathBuf, mixed_species: bool) -> Result<CellUmiGeneSummary> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("open {}", path.display()))?;
    let headers = reader.headers()?.clone();
    let row_idx = CellUmiGeneColumns::from_headers(&headers)?;

    let mut per_cell: HashMap<String, CellStats> = HashMap::default();
    let mut final_cell_reads = HashMap::default();
    let mut genes = HashMap::default();
    let mut umis = HashMap::default();
    let mut known_gene_reads = HashSet::default();
    let mut unique_genes = HashSet::default();

    for record in reader.records() {
        let record = record?;
        let row = row_idx.parse(&record)?;
        let read_id = intern_string(&mut final_cell_reads, row.read_id);
        let known_gene = is_known_gene(&row.gene);
        let mito_gene = is_mito_gene(&row.gene, mixed_species);
        let gene_id = known_gene.then(|| intern_string(&mut genes, row.gene));
        let umi_id = known_gene.then(|| intern_string(&mut umis, row.umi));
        if known_gene {
            known_gene_reads.insert(read_id);
            unique_genes.insert(gene_id.expect("known gene ID"));
        }
        let cell = per_cell.entry(row.barcode).or_default();
        cell.read_ids.insert(read_id);
        if known_gene {
            cell.known_umis.insert(umi_id.expect("known UMI ID"));
            cell.known_genes.insert(gene_id.expect("known gene ID"));
        }
        if mito_gene {
            cell.mito_rows += 1;
        }
    }

    Ok(CellUmiGeneSummary {
        per_cell,
        final_cell_reads,
        known_gene_reads: known_gene_reads.len(),
        unique_genes: unique_genes.len(),
    })
}

fn intern_string(values: &mut HashMap<String, u32>, value: String) -> u32 {
    if let Some(id) = values.get(value.as_str()) {
        return *id;
    }
    let id = values.len() as u32;
    values.insert(value, id);
    id
}

struct CellUmiGeneColumns {
    read_id: usize,
    gene: usize,
    barcode: usize,
    umi: usize,
}

impl CellUmiGeneColumns {
    fn from_headers(headers: &StringRecord) -> Result<Self> {
        let mut cols = HashMap::default();
        for (idx, name) in headers.iter().enumerate() {
            cols.insert(name, idx);
        }
        Ok(Self {
            read_id: *cols.get("read_id").context("missing read_id column")?,
            gene: *cols.get("gene").context("missing gene column")?,
            barcode: *cols.get("barcode").context("missing barcode column")?,
            umi: *cols.get("umi").context("missing umi column")?,
        })
    }

    fn parse(&self, record: &StringRecord) -> Result<CellUmiGeneRow> {
        Ok(CellUmiGeneRow {
            read_id: record.get(self.read_id).unwrap_or("").trim().to_string(),
            gene: record.get(self.gene).unwrap_or("").trim().to_string(),
            barcode: record.get(self.barcode).unwrap_or("").trim().to_string(),
            umi: record.get(self.umi).unwrap_or("").trim().to_string(),
        })
    }
}

fn summarize_bam_alignments(
    path: &PathBuf,
    final_cell_reads: &HashMap<String, u32>,
    threads: usize,
) -> Result<BamSummary> {
    let mut bam =
        bam::Reader::from_path(path).with_context(|| format!("open {}", path.display()))?;
    bam.set_threads(threads.max(1))?;
    let mut summary = BamSummary::default();
    let mut mapped_unique_reads = HashSet::default();
    let mut aligned_genome_reads = HashSet::default();

    for record in bam.records() {
        let rec = record?;
        summary.bam_records += 1;
        if rec.is_unmapped() {
            summary.unmapped_reads += 1;
            continue;
        }
        let read_id = String::from_utf8_lossy(rec.qname());
        if let Some(id) = final_cell_reads.get(read_id.as_ref()) {
            aligned_genome_reads.insert(*id);
        }
        mapped_unique_reads.insert(read_id.into_owned());
        if rec.is_supplementary() {
            summary.supplementary_alignments += 1;
        } else if !rec.is_secondary() {
            summary.mapped_primary_alignments += 1;
        }
    }

    summary.mapped_unique_reads = mapped_unique_reads.len();
    summary.aligned_genome_reads_in_final_cells = aligned_genome_reads.len();
    Ok(summary)
}

fn count_fastq_reads(path: &PathBuf) -> Result<usize> {
    for_each_fastq_batch(&[path], 100_000, |_batch| Ok(()))
}

fn same_file(left: &PathBuf, right: &PathBuf) -> Result<bool> {
    Ok(
        std::fs::canonicalize(left).with_context(|| format!("resolve {}", left.display()))?
            == std::fs::canonicalize(right)
                .with_context(|| format!("resolve {}", right.display()))?,
    )
}

fn load_fastq_count(path: &PathBuf) -> Result<usize> {
    let mut value = String::new();
    BufReader::new(File::open(path).with_context(|| format!("open {}", path.display()))?)
        .read_line(&mut value)?;
    value
        .trim()
        .parse()
        .with_context(|| format!("invalid FASTQ record count in {}", path.display()))
}

fn count_read_tags(path: &PathBuf) -> Result<usize> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("open {}", path.display()))?;
    let headers = reader.headers()?.clone();
    let read_idx = headers
        .iter()
        .position(|name| name == "read_id")
        .context("missing read_id column in read_tags")?;
    let mut seen = HashSet::default();
    for record in reader.records() {
        let record = record?;
        let read_id = record.get(read_idx).unwrap_or("").trim();
        if !read_id.is_empty() {
            seen.insert(read_id.to_string());
        }
    }
    Ok(seen.len())
}

fn load_barcode_valid_reads(path: &PathBuf) -> Result<Option<usize>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("open {}", path.display()))?;
    let headers = reader.headers()?.clone();
    let value_idx = headers.iter().position(|name| name == "valid_any_n");
    let Some(value_idx) = value_idx else {
        return Ok(None);
    };
    let mut records = reader.records();
    if let Some(record) = records.next() {
        let record = record?;
        let value = record.get(value_idx).unwrap_or("").trim().parse().ok();
        return Ok(value);
    }
    Ok(None)
}

fn load_transcript_assignment_summary(
    path: &PathBuf,
    final_cell_reads: &HashMap<String, u32>,
) -> Result<(Option<usize>, Option<usize>)> {
    let reader =
        BufReader::new(File::open(path).with_context(|| format!("open {}", path.display()))?);
    let mut known_reads = HashSet::default();
    let mut isoforms = HashSet::default();
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 || fields[1] != "Assigned" || fields[4].trim() == "NA" {
            continue;
        }
        let read_id = fields[0].trim();
        let Some(read_id) = final_cell_reads.get(read_id) else {
            continue;
        };
        known_reads.insert(*read_id);
        isoforms.insert(fields[4].trim().to_string());
    }
    Ok((Some(known_reads.len()), Some(isoforms.len())))
}

fn load_glycine_total_full_length_reads(path: &PathBuf) -> Result<Option<usize>> {
    load_glycine_summary_read_count(path, "Full-length")
}

fn load_glycine_summary_read_count(path: &PathBuf, row_name: &str) -> Result<Option<usize>> {
    let reader =
        BufReader::new(File::open(path).with_context(|| format!("open {}", path.display()))?);
    let mut in_type_block = false;
    for line in reader.lines() {
        let line = line?;
        let stripped = line.trim();
        if stripped.is_empty() {
            continue;
        }
        if stripped == "Type\tRead_count\tRead_proportion(%)" {
            in_type_block = true;
            continue;
        }
        if stripped == "Non-chimeric" {
            break;
        }
        if in_type_block {
            let parts: Vec<&str> = stripped.split('\t').collect();
            if parts.len() >= 2 && parts[0] == row_name {
                return Ok(parts[1].parse().ok());
            }
        }
    }
    Ok(None)
}

fn is_unannotated_region_label(gene_name: &str) -> bool {
    let bytes = gene_name.as_bytes();
    let mut i = 0usize;
    while i < bytes.len() {
        if !bytes[i].is_ascii_alphanumeric() {
            i += 1;
            continue;
        }
        let mut j = i;
        while j < bytes.len() && bytes[j].is_ascii_alphanumeric() {
            j += 1;
        }
        if j >= bytes.len() || bytes[j] != b'_' {
            i = j.saturating_add(1);
            continue;
        }
        let mut k = j + 1;
        if k >= bytes.len() || !bytes[k].is_ascii_digit() {
            i = j + 1;
            continue;
        }
        while k < bytes.len() && bytes[k].is_ascii_digit() {
            k += 1;
        }
        if k >= bytes.len() || bytes[k] != b'_' {
            i = k.saturating_add(1);
            continue;
        }
        let mut m = k + 1;
        if m >= bytes.len() || !bytes[m].is_ascii_digit() {
            i = k + 1;
            continue;
        }
        while m < bytes.len() && bytes[m].is_ascii_digit() {
            m += 1;
        }
        return true;
    }
    false
}

fn is_known_gene(gene_name: &str) -> bool {
    let gene_name = gene_name.trim();
    if matches!(gene_name, "" | "NA" | "nan" | "None") {
        return false;
    }
    !is_unannotated_region_label(gene_name)
}

fn is_mito_gene(gene_name: &str, mixed_species: bool) -> bool {
    fn starts_with_mito(s: &str) -> bool {
        s.starts_with("MT-") || s.starts_with("mt-") || s.starts_with("Mt-")
    }
    if starts_with_mito(gene_name) {
        return true;
    }
    mixed_species
        && gene_name
            .split_once('_')
            .map(|(_, suffix)| starts_with_mito(suffix))
            .unwrap_or(false)
}

fn write_metrics_files(
    raw_fastq_reads: usize,
    full_length_reads: usize,
    estimated_cells: usize,
    assigned_cell_reads: usize,
    barcode_valid_reads: usize,
    cell_umi_gene: &CellUmiGeneSummary,
    bam: &BamSummary,
    known_transcript_reads: Option<usize>,
    unique_isoforms: Option<usize>,
) -> Result<()> {
    let cell_associated_reads = cell_umi_gene.final_cell_reads.len();
    let aligned_genome_reads = bam.aligned_genome_reads_in_final_cells;
    let mean_reads_per_cell = if estimated_cells > 0 {
        cell_associated_reads as f64 / estimated_cells as f64
    } else {
        0.0
    };

    let mut reads_per_cell = Vec::with_capacity(cell_umi_gene.per_cell.len());
    let mut umis_per_cell = Vec::with_capacity(cell_umi_gene.per_cell.len());
    let mut genes_per_cell = Vec::with_capacity(cell_umi_gene.per_cell.len());
    for stats in cell_umi_gene.per_cell.values() {
        reads_per_cell.push(stats.read_ids.len() as f64);
        umis_per_cell.push(stats.known_umis.len() as f64);
        genes_per_cell.push(stats.known_genes.len() as f64);
    }

    let metrics = vec![
        ("Input reads".to_string(), MetricValue::Int(raw_fastq_reads)),
        (
            "Full length reads".to_string(),
            MetricValue::Int(full_length_reads),
        ),
        (
            "Estimated number of cells".to_string(),
            MetricValue::Int(estimated_cells),
        ),
        (
            "Raw FASTQ reads".to_string(),
            MetricValue::Int(raw_fastq_reads),
        ),
        (
            "Aligned BAM reads".to_string(),
            MetricValue::Int(bam.mapped_unique_reads),
        ),
        (
            "Aligned BAM records".to_string(),
            MetricValue::Int(bam.bam_records),
        ),
        (
            "Pass reads".to_string(),
            MetricValue::Int(full_length_reads),
        ),
        (
            "Mapped".to_string(),
            MetricValue::Int(bam.mapped_primary_alignments),
        ),
        ("Unmapped".to_string(), MetricValue::Int(bam.unmapped_reads)),
        (
            "Supplementary".to_string(),
            MetricValue::Int(bam.supplementary_alignments),
        ),
        (
            "Barcode-valid reads".to_string(),
            MetricValue::Int(barcode_valid_reads),
        ),
        (
            "Reads assigned to final cells".to_string(),
            MetricValue::Int(assigned_cell_reads),
        ),
        (
            "Gene assigned reads".to_string(),
            MetricValue::Int(cell_umi_gene.known_gene_reads),
        ),
        (
            "Reads in final cells".to_string(),
            MetricValue::Int(cell_associated_reads),
        ),
        (
            "Reads aligned to reference genome in final cells".to_string(),
            MetricValue::Int(aligned_genome_reads),
        ),
        (
            "Reads per cell (mean)".to_string(),
            MetricValue::Float(mean_reads_per_cell),
        ),
        (
            "Mean reads per cell".to_string(),
            MetricValue::Float(mean_reads_per_cell),
        ),
        (
            "Mean cell-associated reads per cell".to_string(),
            MetricValue::Float(mean(&reads_per_cell)),
        ),
        (
            "Mean UMI counts per cell".to_string(),
            MetricValue::Float(mean(&umis_per_cell)),
        ),
        (
            "Median UMI counts per cell".to_string(),
            MetricValue::Float(median_from_f64(&umis_per_cell)),
        ),
        (
            "Mean Genes per cell".to_string(),
            MetricValue::Float(mean(&genes_per_cell)),
        ),
        (
            "Median Genes per cell".to_string(),
            MetricValue::Float(median_from_f64(&genes_per_cell)),
        ),
        (
            "Genes per cell (median)".to_string(),
            MetricValue::Float(median_from_f64(&genes_per_cell)),
        ),
        (
            "Unique genes".to_string(),
            MetricValue::Int(cell_umi_gene.unique_genes),
        ),
        (
            "Total genes detected".to_string(),
            MetricValue::Int(cell_umi_gene.unique_genes),
        ),
        (
            "Fraction reads in cells".to_string(),
            MetricValue::Float(ratio(cell_associated_reads, full_length_reads)),
        ),
        (
            "Percent full length reads".to_string(),
            MetricValue::Float(ratio(full_length_reads, raw_fastq_reads)),
        ),
        (
            "Percent barcode-valid reads of full length".to_string(),
            MetricValue::Float(ratio(barcode_valid_reads, full_length_reads)),
        ),
        (
            "Percent reads assigned to final cells of full length".to_string(),
            MetricValue::Float(ratio(assigned_cell_reads, full_length_reads)),
        ),
        (
            "Percent gene assigned reads of full length".to_string(),
            MetricValue::Float(ratio(cell_umi_gene.known_gene_reads, full_length_reads)),
        ),
        (
            "Fraction reads aligned to reference genome in final cells".to_string(),
            MetricValue::Float(ratio(aligned_genome_reads, full_length_reads)),
        ),
    ];

    let mut metrics = metrics;
    if let Some(value) = known_transcript_reads {
        metrics.push((
            "Transcript assigned reads".to_string(),
            MetricValue::Int(value),
        ));
        metrics.push((
            "Percent transcript assigned reads of full length".to_string(),
            MetricValue::Float(ratio(value, full_length_reads)),
        ));
        metrics.push((
            "High-confidence known-transcript reads".to_string(),
            MetricValue::Int(value),
        ));
        metrics.push((
            "Fraction high-confidence known-transcript reads in final cells".to_string(),
            MetricValue::Float(ratio(value, full_length_reads)),
        ));
    }
    if let Some(value) = unique_isoforms {
        metrics.push(("Unique isoforms".to_string(), MetricValue::Int(value)));
    }

    let fraction_metrics: HashSet<&str> = [
        "Fraction reads in cells",
        "Percent full length reads",
        "Percent barcode-valid reads of full length",
        "Percent reads assigned to final cells of full length",
        "Percent gene assigned reads of full length",
        "Percent transcript assigned reads of full length",
        "Fraction reads aligned to reference genome in final cells",
        "Fraction high-confidence known-transcript reads in final cells",
    ]
    .into_iter()
    .collect();

    let mut writer = WriterBuilder::new()
        .delimiter(b'\t')
        .from_path("rna_qc_metrics.tsv")?;
    writer.write_record(["Metric", "Value"])?;
    for (metric, value) in &metrics {
        writer.write_record([
            metric,
            &format_metric_value(metric, value, &fraction_metrics),
        ])?;
    }
    writer.flush()?;

    let report_records = vec![
        (
            "Experiment summary",
            "Input reads",
            metric_lookup(&metrics, "Input reads"),
        ),
        (
            "Experiment summary",
            "Estimated cells",
            metric_lookup(&metrics, "Estimated number of cells"),
        ),
        (
            "Experiment summary",
            "Reads per cell (mean)",
            metric_lookup(&metrics, "Reads per cell (mean)"),
        ),
        (
            "Experiment summary",
            "UMIs per cell (median)",
            metric_lookup(&metrics, "Median UMI counts per cell"),
        ),
        (
            "Experiment summary",
            "Genes per cell (median)",
            metric_lookup(&metrics, "Genes per cell (median)"),
        ),
        (
            "Alignment / feature summary",
            "Pass reads",
            metric_lookup(&metrics, "Pass reads"),
        ),
        (
            "Alignment / feature summary",
            "Mapped",
            metric_lookup(&metrics, "Mapped"),
        ),
        (
            "Alignment / feature summary",
            "Unmapped",
            metric_lookup(&metrics, "Unmapped"),
        ),
        (
            "Alignment / feature summary",
            "Supplementary",
            metric_lookup(&metrics, "Supplementary"),
        ),
        (
            "Alignment / feature summary",
            "Unique genes",
            metric_lookup(&metrics, "Unique genes"),
        ),
        (
            "Alignment / feature summary",
            "Unique isoforms",
            metric_lookup(&metrics, "Unique isoforms"),
        ),
        (
            "Read assignment summary",
            "Full length",
            metric_lookup(&metrics, "Full length reads"),
        ),
        (
            "Read assignment summary",
            "Barcode-valid",
            metric_lookup(&metrics, "Barcode-valid reads"),
        ),
        (
            "Read assignment summary",
            "Cell-assigned",
            metric_lookup(&metrics, "Reads assigned to final cells"),
        ),
        (
            "Read assignment summary",
            "Gene assigned",
            metric_lookup(&metrics, "Gene assigned reads"),
        ),
        (
            "Read assignment summary",
            "Transcript assigned",
            metric_lookup(&metrics, "Transcript assigned reads"),
        ),
        (
            "Read assignment percentage",
            "% full length reads",
            metric_lookup(&metrics, "Percent full length reads"),
        ),
        (
            "Read assignment percentage",
            "% barcode-valid reads",
            metric_lookup(&metrics, "Percent barcode-valid reads of full length"),
        ),
        (
            "Read assignment percentage",
            "% cell-assigned reads",
            metric_lookup(
                &metrics,
                "Percent reads assigned to final cells of full length",
            ),
        ),
        (
            "Read assignment percentage",
            "% gene assigned reads",
            metric_lookup(&metrics, "Percent gene assigned reads of full length"),
        ),
        (
            "Read assignment percentage",
            "% transcript assigned reads",
            metric_lookup(&metrics, "Percent transcript assigned reads of full length"),
        ),
    ];

    let report_fraction_metrics: HashSet<&str> = [
        "% full length reads",
        "% barcode-valid reads",
        "% cell-assigned reads",
        "% gene assigned reads",
        "% transcript assigned reads",
    ]
    .into_iter()
    .collect();

    let mut report_writer = WriterBuilder::new()
        .delimiter(b'\t')
        .from_path("single_cell_report_metrics.tsv")?;
    report_writer.write_record(["Section", "Metric", "Value", "Formatted_value"])?;
    for (section, metric, value) in report_records {
        let raw = value.map(|x| x.as_f64_string()).unwrap_or_default();
        let formatted = match value {
            None => String::new(),
            Some(value) if report_fraction_metrics.contains(metric) => {
                format!("{:.2}%", value.as_f64() * 100.0)
            }
            Some(value) => value.formatted_numeric(),
        };
        report_writer.write_record([section, metric, &raw, &formatted])?;
    }
    report_writer.flush()?;

    let stdout = metrics
        .iter()
        .map(|(metric, value)| {
            format!(
                "{metric}\t{}",
                format_metric_value(metric, value, &fraction_metrics)
            )
        })
        .collect::<Vec<_>>()
        .join("\n");
    println!("\nRNA QC Metrics\nMetric\tValue\n{stdout}");
    Ok(())
}

fn write_per_cell_qc(per_cell: &HashMap<String, CellStats>) -> Result<()> {
    let mut writer = WriterBuilder::new()
        .delimiter(b'\t')
        .from_path("per_cell_qc.tsv")?;
    writer.write_record(["barcode", "reads", "umis", "genes", "mito_percent"])?;
    let mut barcodes = per_cell.keys().cloned().collect::<Vec<_>>();
    barcodes.sort();
    for barcode in barcodes {
        let stats = &per_cell[&barcode];
        let reads = stats.read_ids.len();
        let umis = stats.known_umis.len();
        let genes = stats.known_genes.len();
        let mito_percent = if reads > 0 {
            stats.mito_rows as f64 / reads as f64 * 100.0
        } else {
            0.0
        };
        writer.write_record([
            barcode,
            reads.to_string(),
            umis.to_string(),
            genes.to_string(),
            format!("{mito_percent:.6}"),
        ])?;
    }
    writer.flush()?;
    Ok(())
}

#[derive(Clone, Copy)]
enum MetricValue {
    Int(usize),
    Float(f64),
}

impl MetricValue {
    fn as_f64(self) -> f64 {
        match self {
            Self::Int(x) => x as f64,
            Self::Float(x) => x,
        }
    }

    fn as_f64_string(self) -> String {
        match self {
            Self::Int(x) => x.to_string(),
            Self::Float(x) => format!("{x}"),
        }
    }

    fn formatted_numeric(self) -> String {
        match self {
            Self::Int(x) => format_int_with_commas(x),
            Self::Float(x) => format_float_with_commas(x),
        }
    }
}

fn metric_lookup<'a>(metrics: &'a [(String, MetricValue)], key: &str) -> Option<MetricValue> {
    metrics
        .iter()
        .find(|(metric, _)| metric == key)
        .map(|(_, value)| *value)
}

fn ratio(numerator: usize, denominator: usize) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}

fn mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        0.0
    } else {
        values.iter().sum::<f64>() / values.len() as f64
    }
}

fn median_from_f64(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.total_cmp(b));
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] + sorted[mid]) / 2.0
    } else {
        sorted[mid]
    }
}

fn format_metric_value(
    metric: &str,
    value: &MetricValue,
    fraction_metrics: &HashSet<&str>,
) -> String {
    if fraction_metrics.contains(metric) {
        format!("{:.2}%", value.as_f64() * 100.0)
    } else {
        match value {
            MetricValue::Int(x) => format_int_with_commas(*x),
            MetricValue::Float(x) => format_float_with_commas(*x),
        }
    }
}

fn format_int_with_commas(value: usize) -> String {
    let s = value.to_string();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    for (i, ch) in s.chars().rev().enumerate() {
        if i != 0 && i % 3 == 0 {
            out.push(',');
        }
        out.push(ch);
    }
    out.chars().rev().collect()
}

fn format_float_with_commas(value: f64) -> String {
    let sign = if value.is_sign_negative() { "-" } else { "" };
    let abs = value.abs();
    let s = format!("{abs:.2}");
    let mut parts = s.split('.');
    let int_part = parts.next().unwrap_or("0").parse::<usize>().unwrap_or(0);
    let frac_part = parts.next().unwrap_or("00");
    format!("{sign}{}.{}", format_int_with_commas(int_part), frac_part)
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use super::*;

    #[test]
    fn mito_detection_matches_mixed_species_labels() {
        assert!(is_mito_gene("MT-CO1", false));
        assert!(!is_mito_gene("human_MT-CO1", false));
        assert!(is_mito_gene("human_MT-CO1", true));
        assert!(!is_mito_gene("ACTB", true));
    }

    #[test]
    fn unannotated_region_detection_matches_python_shape() {
        assert!(is_unannotated_region_label("chr1_1000_2000"));
        assert!(is_unannotated_region_label("NC_000001.11_1000_2000"));
        assert!(is_unannotated_region_label("human_chr1_1000_2000"));
        assert!(!is_unannotated_region_label("MALAT1"));
    }

    #[test]
    fn metric_formatting_matches_python_report_style() {
        let fraction_metrics: HashSet<&str> = ["Percent full length reads"].into_iter().collect();
        assert_eq!(
            format_metric_value(
                "Percent full length reads",
                &MetricValue::Float(0.125),
                &fraction_metrics
            ),
            "12.50%"
        );
        assert_eq!(
            format_metric_value("Input reads", &MetricValue::Int(12345), &HashSet::default()),
            "12,345"
        );
        assert_eq!(format_float_with_commas(1234.5), "1,234.50");
    }

    #[test]
    fn glycine_summary_total_and_full_length_counts_are_parsed() {
        let path =
            std::env::temp_dir().join(format!("flora-glycine-summary-{}.txt", std::process::id()));
        std::fs::write(
            &path,
            "Summary\nTotal_base_count\tValid_base_count\tValid_base_proportion(%)\n10\t9\t90.0\nType\tRead_count\tRead_proportion(%)\nTotal\t12345\t100.00\nLength-filtered\t10\t0.08\nQC-filtered\t5\t0.04\nFull-length\t10000\t81.00\n\nNon-chimeric\n",
        )
        .unwrap();

        assert_eq!(
            load_glycine_summary_read_count(&path, "Total").unwrap(),
            Some(12345)
        );
        assert_eq!(
            load_glycine_total_full_length_reads(&path).unwrap(),
            Some(10000)
        );
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn cell_summary_interns_ids_without_changing_deduplication() {
        let mut input = tempfile::NamedTempFile::new().unwrap();
        writeln!(input, "read_id\tgene\tbarcode\tumi").unwrap();
        writeln!(input, "r1\tACTB\tC1\tU1").unwrap();
        writeln!(input, "r1\tACTB\tC1\tU1").unwrap();
        writeln!(input, "r2\tchr1_100_200\tC2\tU2").unwrap();
        writeln!(input, "r3\tMT-CO1\tC1\tU3").unwrap();

        let summary = load_cell_umi_gene(&input.path().to_path_buf(), false).unwrap();
        assert_eq!(summary.final_cell_reads.len(), 3);
        assert_eq!(summary.known_gene_reads, 2);
        assert_eq!(summary.unique_genes, 2);
        assert_eq!(summary.per_cell["C1"].read_ids.len(), 2);
        assert_eq!(summary.per_cell["C1"].known_umis.len(), 2);
        assert_eq!(summary.per_cell["C1"].known_genes.len(), 2);
        assert_eq!(summary.per_cell["C1"].mito_rows, 1);
        assert_eq!(summary.per_cell["C2"].known_umis.len(), 0);
    }
}
