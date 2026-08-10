use std::collections::BTreeSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Parser;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use rust_htslib::bam::{self, Read};
use flora::matrices::{add_unique_umi, is_genomic_placeholder, matrix_axes};

#[derive(Debug, Parser)]
#[command(version, about = "Build gene expression matrix from tagged BAM")]
struct Cli {
    bam: PathBuf,

    #[arg(long = "output", default_value = "gene_expression.tsv")]
    output: PathBuf,

    #[arg(long = "verbosity", default_value_t = 2)]
    _verbosity: u8,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let mut bam = bam::Reader::from_path(&cli.bam).with_context(|| format!("open {}", cli.bam.display()))?;
    let mut matrix = HashMap::default();
    let mut cells = BTreeSet::new();
    for rec in bam.records() {
        let rec = rec?;
        let gene = require_string_tag(&rec, b"GN", "GN")?;
        let cell = require_string_tag(&rec, b"CB", "CB")?;
        let umi = require_string_tag(&rec, b"UB", "UB")?;
        collect_gene_observation(&gene, &cell, &umi, &mut matrix, &mut cells);
    }
    write_matrix(&cli.output, "gene", &matrix, &cells)
}

fn collect_gene_observation(
    gene: &str,
    cell: &str,
    umi: &str,
    matrix: &mut HashMap<(String, String), HashSet<String>>,
    cells: &mut BTreeSet<String>,
) {
    cells.insert(cell.to_string());
    if !is_genomic_placeholder(gene) {
        add_unique_umi(matrix, gene, cell, umi);
    }
}

fn write_matrix(
    output: &PathBuf,
    first_col: &str,
    matrix: &HashMap<(String, String), HashSet<String>>,
    cells: &BTreeSet<String>,
) -> Result<()> {
    let (rows, _) = matrix_axes(matrix);
    let mut writer = BufWriter::new(File::create(output).with_context(|| format!("create {}", output.display()))?);
    write!(writer, "{first_col}")?;
    for col in cells {
        write!(writer, "\t{col}")?;
    }
    writeln!(writer)?;
    for row in &rows {
        write!(writer, "{row}")?;
        for col in cells {
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

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use super::*;

    #[test]
    fn placeholder_only_cells_are_preserved_in_matrix_columns() {
        let mut matrix = HashMap::default();
        let mut cells = BTreeSet::new();
        collect_gene_observation("ACTB", "CELL_A", "UMI_A", &mut matrix, &mut cells);
        collect_gene_observation(
            "chr1_1000_2000",
            "CELL_B",
            "UMI_B",
            &mut matrix,
            &mut cells,
        );

        let tmp = tempfile::NamedTempFile::new().unwrap();
        write_matrix(
            &tmp.path().to_path_buf(),
            "gene",
            &matrix,
            &cells,
        )
        .unwrap();
        let text = std::fs::read_to_string(tmp.path()).unwrap();
        let lines: Vec<_> = text.lines().collect();

        assert_eq!(lines[0], "gene\tCELL_A\tCELL_B");
        assert_eq!(lines[1], "ACTB\t1\t0");
    }
}
