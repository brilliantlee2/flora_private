use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Parser;
use rust_htslib::bam::{self, Read};

#[derive(Debug, Parser)]
#[command(
    version,
    about = "Extract read_id, gene, barcode, and UMI table from tagged BAM"
)]
struct Cli {
    bam: PathBuf,

    #[arg(long = "output", default_value = "cell_umi_gene.tsv")]
    output: PathBuf,

    #[arg(long = "verbosity", default_value_t = 2)]
    _verbosity: u8,
}

pub fn main() -> Result<()> {
    let cli = Cli::parse();
    let mut bam =
        bam::Reader::from_path(&cli.bam).with_context(|| format!("open {}", cli.bam.display()))?;
    let mut writer = BufWriter::new(
        File::create(&cli.output).with_context(|| format!("create {}", cli.output.display()))?,
    );
    writeln!(writer, "read_id\tgene\tbarcode\tumi")?;
    for rec in bam.records() {
        let rec = rec?;
        let read_id = String::from_utf8_lossy(rec.qname()).to_string();
        let gene = require_string_tag(&rec, b"GN", "GN")?;
        let barcode = require_string_tag(&rec, b"CB", "CB")?;
        let umi = require_string_tag(&rec, b"UB", "UB")?;
        writeln!(writer, "{read_id}\t{gene}\t{barcode}\t{umi}")?;
    }
    Ok(())
}

fn require_string_tag(rec: &bam::Record, tag: &[u8; 2], label: &str) -> Result<String> {
    match rec.aux(tag).with_context(|| {
        format!(
            "missing {label} tag on read {}",
            String::from_utf8_lossy(rec.qname())
        )
    })? {
        bam::record::Aux::String(v) => Ok(v.to_string()),
        _ => anyhow::bail!(
            "non-string {label} tag on read {}",
            String::from_utf8_lossy(rec.qname())
        ),
    }
}
