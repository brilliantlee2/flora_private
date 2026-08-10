use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Parser;
use rust_htslib::bam::{self, Read};

#[derive(Debug, Parser)]
#[command(version, about = "Add GN tags to BAM using read/gene assignment TSV")]
struct Cli {
    bam: PathBuf,
    gene_assigns: PathBuf,

    #[arg(long = "output", default_value = "gene.sorted.bam")]
    output: PathBuf,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let mut bam = bam::Reader::from_path(&cli.bam).with_context(|| format!("open {}", cli.bam.display()))?;
    let header = bam::Header::from_template(bam.header());
    let mut out = bam::Writer::from_path(&cli.output, &header, bam::Format::Bam)
        .with_context(|| format!("create {}", cli.output.display()))?;
    let mut lines = BufReader::new(File::open(&cli.gene_assigns).with_context(|| format!("open {}", cli.gene_assigns.display()))?).lines();

    for rec in bam.records() {
        let mut rec = rec?;
        let Some(line) = lines.next() else { break };
        let line = line?;
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() >= 4 {
            let rid = String::from_utf8_lossy(rec.qname()).to_string();
            if rid != fields[0] {
                anyhow::bail!("BAM and gene assignment reads not ordered: {} != {}", rid, fields[0]);
            }
            rec.update_aux(b"GN", bam::record::Aux::String(fields[3]))?;
        }
        out.write(&rec)?;
    }
    Ok(())
}
