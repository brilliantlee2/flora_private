use std::collections::BTreeSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Parser;
use flora::matrices::{add_unique_umi, is_genomic_placeholder, matrix_axes};
use rust_htslib::bam::{self, Read};
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

#[derive(Debug, Parser)]
#[command(
    version,
    about = "Extract read_id, gene, barcode, and UMI table from tagged BAM"
)]
struct Cli {
    bam: PathBuf,

    #[arg(long = "output", default_value = "cell_umi_gene.tsv")]
    output: PathBuf,

    #[arg(long = "gene-expression-output")]
    gene_expression_output: Option<PathBuf>,

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
    let mut matrix = HashMap::default();
    let mut cells = BTreeSet::new();
    for rec in bam.records() {
        let rec = rec?;
        let read_id = String::from_utf8_lossy(rec.qname()).to_string();
        let gene = require_string_tag(&rec, b"GN", "GN")?;
        let barcode = require_string_tag(&rec, b"CB", "CB")?;
        let umi = require_string_tag(&rec, b"UB", "UB")?;
        writeln!(writer, "{read_id}\t{gene}\t{barcode}\t{umi}")?;
        if cli.gene_expression_output.is_some() {
            cells.insert(barcode.clone());
            if !is_genomic_placeholder(&gene) {
                add_unique_umi(&mut matrix, &gene, &barcode, &umi);
            }
        }
    }
    writer.flush()?;
    if let Some(path) = &cli.gene_expression_output {
        write_gene_matrix(path, &matrix, &cells)?;
    }
    Ok(())
}

fn write_gene_matrix(
    output: &PathBuf,
    matrix: &HashMap<(String, String), HashSet<String>>,
    cells: &BTreeSet<String>,
) -> Result<()> {
    let (genes, _) = matrix_axes(matrix);
    let mut writer = BufWriter::new(
        File::create(output).with_context(|| format!("create {}", output.display()))?,
    );
    write!(writer, "gene")?;
    for cell in cells {
        write!(writer, "\t{cell}")?;
    }
    writeln!(writer)?;
    for gene in genes {
        write!(writer, "{gene}")?;
        for cell in cells {
            let count = matrix
                .get(&(gene.clone(), cell.clone()))
                .map(HashSet::len)
                .unwrap_or(0);
            write!(writer, "\t{count}")?;
        }
        writeln!(writer)?;
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn combined_gene_matrix_preserves_placeholder_only_cells() {
        let mut matrix = HashMap::default();
        let cells = BTreeSet::from(["CELL_A".to_string(), "CELL_B".to_string()]);
        add_unique_umi(&mut matrix, "ACTB", "CELL_A", "UMI_A");
        let tmp = tempfile::NamedTempFile::new().unwrap();
        write_gene_matrix(&tmp.path().to_path_buf(), &matrix, &cells).unwrap();
        assert_eq!(
            std::fs::read_to_string(tmp.path()).unwrap(),
            "gene\tCELL_A\tCELL_B\nACTB\t1\t0\n"
        );
    }
}
