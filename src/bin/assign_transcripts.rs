use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter};
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Parser;
use flora::annotation::{
    assign_transcript, bedtools_bamtobed_split, load_exon_gtf, parse_bed6_line, write_transcript_assignment,
    BedRecord,
};

#[derive(Debug, Parser)]
#[command(version, about = "Assign reads in BAM to transcripts using split BED blocks and GTF exons")]
struct Cli {
    bam: PathBuf,
    gtf: PathBuf,

    #[arg(long = "output", default_value = "read_transcript_assigns.tsv")]
    output: PathBuf,

    #[arg(short = 'q', long = "mapq", default_value_t = 60)]
    mapq: i32,

    #[arg(short = 'c', long = "chunk_size", default_value_t = 200_000)]
    _chunk_size: usize,

    #[arg(long = "verbosity", default_value_t = 2)]
    _verbosity: u8,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    run_cli(cli)
}

fn run_cli(cli: Cli) -> Result<()> {
    let exons = load_exon_gtf(&cli.gtf)?;
    if exons.is_empty() {
        File::create(&cli.output).with_context(|| format!("create {}", cli.output.display()))?;
        return Ok(());
    }
    let split_bed = temp_split_bed_path();
    bedtools_bamtobed_split(&cli.bam, &split_bed)?;
    let result = run(&split_bed, &cli.output, &exons, cli.mapq);
    let _ = fs::remove_file(&split_bed);
    result
}

fn run(
    split_bed: &PathBuf,
    output: &PathBuf,
    exons: &flora::annotation::ExonIndex,
    mapq: i32,
) -> Result<()> {
    let reader = BufReader::new(File::open(split_bed).with_context(|| format!("open {}", split_bed.display()))?);
    let mut writer =
        BufWriter::new(File::create(output).with_context(|| format!("create {}", output.display()))?);
    let mut current_name = String::new();
    let mut current_blocks: Vec<BedRecord> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let Some(record) = parse_bed6_line(&line) else { continue };
        if current_name.is_empty() {
            current_name = record.name.clone();
        }
        if record.name != current_name {
            let row = assign_transcript(&current_blocks, exons, mapq);
            write_transcript_assignment(&mut writer, &row)?;
            current_blocks.clear();
            current_name = record.name.clone();
        }
        current_blocks.push(record);
    }

    if !current_blocks.is_empty() {
        let row = assign_transcript(&current_blocks, exons, mapq);
        write_transcript_assignment(&mut writer, &row)?;
    }
    Ok(())
}

fn temp_split_bed_path() -> PathBuf {
    let mut path = std::env::temp_dir();
    path.push(format!("strint_assign_transcripts_{}.split.bed", std::process::id()));
    path
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_path(name: &str) -> PathBuf {
        let mut path = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        path.push(format!("strint_assign_transcripts_{name}_{nanos}"));
        path
    }

    #[test]
    fn empty_exon_gtf_matches_python_empty_output_behavior() {
        let gtf = temp_path("tx.gtf");
        let output = temp_path("out.tsv");
        fs::write(&gtf, "").unwrap();

        run_cli(Cli {
            bam: PathBuf::from("/nonexistent/test.bam"),
            gtf: gtf.clone(),
            output: output.clone(),
            mapq: 60,
            _chunk_size: 200_000,
            _verbosity: 2,
        })
        .unwrap();

        assert_eq!(fs::read_to_string(&output).unwrap(), "");

        let _ = fs::remove_file(gtf);
        let _ = fs::remove_file(output);
    }
}
