use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Parser;
use flora::annotation::{assign_gene, load_gene_gtf, parse_bed6_line, write_gene_assignment};

#[derive(Debug, Parser)]
#[command(
    version,
    about = "Assign reads in BED6 to genes using GTF gene intervals"
)]
struct Cli {
    bed: PathBuf,
    gtf: PathBuf,

    #[arg(short = 'q', long = "mapq", default_value_t = 60)]
    mapq: i32,

    #[arg(long = "output", default_value = "./read_annotations.tsv")]
    output: PathBuf,

    #[arg(short = 'c', long = "chunk_size", default_value_t = 200_000)]
    _chunk_size: usize,

    #[arg(long = "verbosity", default_value_t = 2)]
    _verbosity: u8,
}

pub fn main() -> Result<()> {
    let cli = Cli::parse();
    run(cli)
}

fn run(cli: Cli) -> Result<()> {
    let genes = load_gene_gtf(&cli.gtf)?;
    if genes.is_empty() {
        File::create(&cli.output).with_context(|| format!("create {}", cli.output.display()))?;
        return Ok(());
    }
    let reader = BufReader::new(
        File::open(&cli.bed).with_context(|| format!("open {}", cli.bed.display()))?,
    );
    let mut writer = BufWriter::new(
        File::create(&cli.output).with_context(|| format!("create {}", cli.output.display()))?,
    );

    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        if let Some(record) = parse_bed6_line(&line) {
            let row = assign_gene(&record, &genes, cli.mapq);
            write_gene_assignment(&mut writer, &row)?;
        }
    }
    Ok(())
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
        path.push(format!("strint_assign_genes_{name}_{nanos}"));
        path
    }

    #[test]
    fn empty_gtf_matches_python_empty_output_behavior() {
        let bed = temp_path("reads.bed");
        let gtf = temp_path("genes.gtf");
        let output = temp_path("out.tsv");
        fs::write(&bed, "chr1\t0\t10\tread1\t60\t+\n").unwrap();
        fs::write(&gtf, "").unwrap();

        run(Cli {
            bed: bed.clone(),
            gtf: gtf.clone(),
            mapq: 60,
            output: output.clone(),
            _chunk_size: 200_000,
            _verbosity: 2,
        })
        .unwrap();

        assert_eq!(fs::read_to_string(&output).unwrap(), "");

        let _ = fs::remove_file(bed);
        let _ = fs::remove_file(gtf);
        let _ = fs::remove_file(output);
    }
}
