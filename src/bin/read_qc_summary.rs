use std::path::PathBuf;

use anyhow::Result;
use clap::Parser;
use flora::fastq::for_each_fastq_batch;
use flora::read_qc::ReadQcAccumulator;

#[derive(Debug, Parser)]
#[command(
    version,
    about = "Summarize FASTQ read length and quality distributions"
)]
struct Cli {
    #[arg(long)]
    fastq: PathBuf,

    #[arg(long = "output-json")]
    output_json: PathBuf,

    #[arg(long = "output-fastq-count")]
    output_fastq_count: Option<PathBuf>,

    #[arg(long = "curve-points", default_value_t = 300)]
    curve_points: usize,

    #[arg(long = "batch-size", default_value_t = 100_000)]
    batch_size: usize,
}

pub fn main() -> Result<()> {
    let cli = Cli::parse();
    let mut qc = ReadQcAccumulator::new();
    for_each_fastq_batch(&[cli.fastq], cli.batch_size, |batch| {
        for rec in batch {
            qc.observe(rec.seq.as_bytes(), rec.qual.as_bytes());
        }
        Ok(())
    })?;

    if let Some(count_path) = &cli.output_fastq_count {
        qc.write_outputs(&cli.output_json, count_path, cli.curve_points)?;
    } else {
        qc.write_json(&cli.output_json, cli.curve_points)?;
    }
    println!("{}", cli.output_json.display());
    Ok(())
}
