use std::collections::HashSet;
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Parser;

#[derive(Debug, Parser)]
#[command(version, about = "Prepare Sockeye-style read tag table from Flora assigned reads")]
struct Cli {
    #[arg(long)]
    input: PathBuf,

    #[arg(long)]
    output: PathBuf,

    #[arg(long = "read-id-col", default_value = "read_id")]
    read_id_col: String,

    #[arg(long = "cell-col", default_value = "cell_id")]
    cell_col: String,

    #[arg(long = "umi-primary-col", default_value = "putative_umi_5p")]
    umi_primary_col: String,

    #[arg(long = "umi-backup-col", default_value = "putative_umi")]
    umi_backup_col: String,

    #[arg(long = "barcode-5p-col", default_value = "BC5n")]
    barcode_5p_col: String,

    #[arg(long = "barcode-3p-col", default_value = "BC3n")]
    barcode_3p_col: String,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    run(cli)
}

fn run(cli: Cli) -> Result<()> {
    let delimiter = if cli.input.extension().and_then(|x| x.to_str()) == Some("csv") {
        b','
    } else {
        b'\t'
    };
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(delimiter)
        .from_path(&cli.input)
        .with_context(|| format!("read {}", cli.input.display()))?;
    let headers = reader.headers()?.clone();

    let read_idx = required_idx(&headers, &cli.read_id_col)?;
    let cell_idx = required_idx(&headers, &cli.cell_col)?;
    let umi_primary_idx = required_idx(&headers, &cli.umi_primary_col)?;
    let umi_backup_idx = required_idx(&headers, &cli.umi_backup_col)?;
    let barcode_5p_idx = required_idx(&headers, &cli.barcode_5p_col)?;
    let barcode_3p_idx = required_idx(&headers, &cli.barcode_3p_col)?;

    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(&cli.output)
        .with_context(|| format!("write {}", cli.output.display()))?;
    writer.write_record([
        "read_id",
        "cell_id",
        "barcode_5p",
        "barcode_3p",
        "barcode_dual",
        "umi_primary",
        "umi_backup",
        "umi_for_clustering",
    ])?;

    let mut seen = HashSet::new();
    for record in reader.records() {
        let record = record?;
        let read_id = clean(record.get(read_idx));
        let cell_id = clean(record.get(cell_idx));
        let barcode_5p = clean(record.get(barcode_5p_idx));
        let barcode_3p = clean(record.get(barcode_3p_idx));
        let umi_primary = clean(record.get(umi_primary_idx));
        let umi_backup = clean(record.get(umi_backup_idx));
        let umi_for_clustering = if umi_primary.is_empty() {
            umi_backup.clone()
        } else {
            umi_primary.clone()
        };
        if read_id.is_empty() || cell_id.is_empty() || umi_for_clustering.is_empty() {
            continue;
        }
        if !seen.insert(read_id.clone()) {
            continue;
        }
        let barcode_dual = format!("{barcode_5p}+{barcode_3p}");
        writer.write_record([
            read_id,
            cell_id,
            barcode_5p,
            barcode_3p,
            barcode_dual,
            umi_primary,
            umi_backup,
            umi_for_clustering,
        ])?;
    }
    writer.flush()?;
    Ok(())
}

fn required_idx(headers: &csv::StringRecord, name: &str) -> Result<usize> {
    headers
        .iter()
        .position(|header| header == name)
        .with_context(|| format!("missing required column: {name}"))
}

fn clean(value: Option<&str>) -> String {
    value.unwrap_or("").trim().to_string()
}
