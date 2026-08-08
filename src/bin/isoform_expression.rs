use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Parser;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use rust_htslib::bam::{self, Read};
use flora::matrices::{add_unique_umi, matrix_axes};

#[derive(Debug, Parser)]
#[command(version, about = "Build isoform expression matrix from tagged BAM and transcript assignments")]
struct Cli {
    bam: PathBuf,
    transcripts: PathBuf,

    #[arg(long = "output", default_value = "isoform_expression.tsv")]
    output: PathBuf,

    #[arg(long = "verbosity", default_value_t = 2)]
    _verbosity: u8,

    #[arg(long = "assignment-chunk-size", default_value_t = 1_000_000)]
    _assignment_chunk_size: usize,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let read_to_tx = load_transcript_assignments(&cli.transcripts)?;
    let mut bam = bam::Reader::from_path(&cli.bam).with_context(|| format!("open {}", cli.bam.display()))?;
    let mut matrix = HashMap::default();
    for rec in bam.records() {
        let rec = rec?;
        let read_id = String::from_utf8_lossy(rec.qname()).to_string();
        let Some(tx) = read_to_tx.get(&read_id) else { continue };
        let cell = require_string_tag(&rec, b"CB", "CB")?;
        let umi = require_string_tag(&rec, b"UB", "UB")?;
        add_unique_umi(&mut matrix, tx, &cell, &umi);
    }
    write_matrix(&cli.output, &matrix)
}

fn load_transcript_assignments(path: &PathBuf) -> Result<HashMap<String, String>> {
    let reader = BufReader::new(File::open(path).with_context(|| format!("open {}", path.display()))?);
    let mut map = HashMap::default();
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 || fields[1] != "Assigned" || fields[4] == "NA" {
            continue;
        }
        map.entry(fields[0].to_string())
            .or_insert_with(|| fields[4].to_string());
    }
    Ok(map)
}

fn write_matrix(
    output: &PathBuf,
    matrix: &HashMap<(String, String), HashSet<String>>,
) -> Result<()> {
    let (rows, cols) = matrix_axes(matrix);
    let mut writer = BufWriter::new(File::create(output).with_context(|| format!("create {}", output.display()))?);
    write!(writer, "transcript_id")?;
    for col in &cols {
        write!(writer, "\t{col}")?;
    }
    writeln!(writer)?;
    for row in &rows {
        write!(writer, "{row}")?;
        for col in &cols {
            let count = matrix
                .get(&(row.clone(), col.clone()))
                .map(|umis| umis.len())
                .unwrap_or(0);
            write!(writer, "\t{count}")?;
        }
        writeln!(writer)?;
    }
    Ok(())
}

fn require_string_tag(rec: &bam::Record, tag: &[u8; 2], label: &str) -> Result<String> {
    match rec.aux(tag).with_context(|| format!("missing {label} tag on read {}", String::from_utf8_lossy(rec.qname())))? {
        bam::record::Aux::String(v) => Ok(v.to_string()),
        _ => anyhow::bail!("non-string {label} tag on read {}", String::from_utf8_lossy(rec.qname())),
    }
}
