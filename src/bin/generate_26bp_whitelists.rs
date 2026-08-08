use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use clap::Parser;
use flora::barcode::reverse_complement;
use zip::ZipArchive;

#[derive(Debug, Parser)]
#[command(version, about = "Expand 10bp barcodes to Flora 26bp 3p/5p whitelists")]
struct Cli {
    #[arg(long = "barcode-list-10bp")]
    barcode_list_10bp: PathBuf,

    #[arg(long = "out-3p")]
    out_3p: PathBuf,

    #[arg(long = "out-5p")]
    out_5p: PathBuf,

    #[arg(long = "middle-3p", default_value = "GGTAGC")]
    middle_3p: String,

    #[arg(long = "middle-5p", default_value = "GGAAGG")]
    middle_5p: String,
}

pub fn main() -> Result<()> {
    let cli = Cli::parse();
    let barcodes = read_10bp_barcodes(&cli.barcode_list_10bp)?;
    if let Some(parent) = cli.out_3p.parent() {
        fs::create_dir_all(parent)?;
    }
    if let Some(parent) = cli.out_5p.parent() {
        fs::create_dir_all(parent)?;
    }
    let (count_3p, count_5p) = write_whitelists(
        &barcodes,
        &cli.middle_3p,
        &cli.middle_5p,
        &cli.out_3p,
        &cli.out_5p,
    )?;
    println!("Input 10bp barcode count: {}", barcodes.len());
    println!(
        "3p whitelist: {} ({count_3p} records)",
        cli.out_3p.display()
    );
    println!(
        "5p whitelist: {} ({count_5p} records)",
        cli.out_5p.display()
    );
    println!(
        "First 3 raw barcodes: {:?}",
        &barcodes[..barcodes.len().min(3)]
    );
    let first_rc: Vec<_> = barcodes
        .iter()
        .take(3)
        .map(|bc| reverse_complement(bc))
        .collect();
    println!("First 3 5p RC barcodes: {first_rc:?}");
    Ok(())
}

fn read_10bp_barcodes(path: &Path) -> Result<Vec<String>> {
    let barcodes = if path
        .extension()
        .and_then(|x| x.to_str())
        .is_some_and(|x| x.eq_ignore_ascii_case("xlsx"))
    {
        read_xlsx_barcodes(path)?
    } else {
        read_plain_barcodes(path)?
    };
    validate_barcodes(barcodes)
}

fn read_plain_barcodes(path: &Path) -> Result<Vec<String>> {
    let text = fs::read_to_string(path).with_context(|| format!("read {}", path.display()))?;
    let first_line = text
        .lines()
        .find(|line| !line.trim().is_empty())
        .unwrap_or("");
    let delimiter = if first_line.contains(',') {
        Some(',')
    } else if first_line.contains('\t') {
        Some('\t')
    } else {
        None
    };
    if let Some(delimiter) = delimiter {
        let mut lines = text.lines();
        let header_line = lines.next().unwrap_or("");
        let headers: Vec<String> = header_line
            .split(delimiter)
            .map(|x| x.trim().trim_matches('\u{feff}').to_ascii_lowercase())
            .collect();
        let seq_idx = headers.iter().position(|h| {
            matches!(
                h.as_str(),
                "sequence" | "seq" | "barcode" | "cellbarcode" | "序列"
            )
        });
        if let Some(idx) = seq_idx {
            return Ok(lines
                .filter_map(|line| line.split(delimiter).nth(idx))
                .map(|x| x.trim().to_ascii_uppercase())
                .filter(|x| !x.is_empty())
                .collect());
        }
    }
    Ok(text
        .lines()
        .filter_map(|line| line.split_whitespace().next())
        .map(|x| x.trim().trim_matches('\u{feff}').to_ascii_uppercase())
        .filter(|x| !x.is_empty() && !matches!(x.as_str(), "SEQUENCE" | "SEQ" | "BARCODE" | "序列"))
        .collect())
}

fn read_xlsx_barcodes(path: &Path) -> Result<Vec<String>> {
    let file = File::open(path)?;
    let mut zip = ZipArchive::new(file)?;
    let shared_strings = read_shared_strings(&mut zip)?;
    let mut sheet = String::new();
    zip.by_name("xl/worksheets/sheet1.xml")?
        .read_to_string(&mut sheet)?;
    let rows = parse_sheet_rows(&sheet, &shared_strings);
    if rows.is_empty() {
        return Ok(Vec::new());
    }
    let header = &rows[0];
    let seq_col = header
        .iter()
        .find_map(|(col, value)| (value.trim() == "序列").then_some(col.clone()))
        .unwrap_or_else(|| "B".to_string());
    Ok(rows
        .iter()
        .skip(1)
        .filter_map(|row| row.get(&seq_col))
        .map(|x| x.trim().to_ascii_uppercase())
        .filter(|x| !x.is_empty())
        .collect())
}

fn read_shared_strings(zip: &mut ZipArchive<File>) -> Result<Vec<String>> {
    let Ok(mut member) = zip.by_name("xl/sharedStrings.xml") else {
        return Ok(Vec::new());
    };
    let mut xml = String::new();
    member.read_to_string(&mut xml)?;
    Ok(xml
        .split("<si>")
        .skip(1)
        .map(|item| {
            item.split("<t")
                .skip(1)
                .filter_map(|part| part.split_once('>')?.1.split_once("</t>").map(|x| x.0))
                .collect::<String>()
        })
        .collect())
}

fn parse_sheet_rows(xml: &str, shared_strings: &[String]) -> Vec<HashMap<String, String>> {
    let mut rows = Vec::new();
    for row_xml in xml.split("<row").skip(1) {
        let Some(row_body) = row_xml.split_once("</row>").map(|x| x.0) else {
            continue;
        };
        let mut row = HashMap::new();
        for cell_xml in row_body.split("<c ").skip(1) {
            let cell_head = cell_xml.split_once('>').map(|x| x.0).unwrap_or("");
            let col = attr_value(cell_head, "r")
                .unwrap_or_default()
                .chars()
                .filter(|c| c.is_ascii_alphabetic())
                .collect::<String>();
            if col.is_empty() {
                continue;
            }
            let raw = cell_xml
                .split_once("<v>")
                .and_then(|x| x.1.split_once("</v>"))
                .map(|x| x.0)
                .unwrap_or("");
            let value = if cell_head.contains("t=\"s\"") {
                raw.parse::<usize>()
                    .ok()
                    .and_then(|idx| shared_strings.get(idx))
                    .cloned()
                    .unwrap_or_default()
            } else {
                raw.to_string()
            };
            row.insert(col, value);
        }
        rows.push(row);
    }
    rows
}

fn attr_value(text: &str, name: &str) -> Option<String> {
    let needle = format!("{name}=\"");
    let rest = text.split_once(&needle)?.1;
    Some(rest.split_once('"')?.0.to_string())
}

fn validate_barcodes(barcodes: Vec<String>) -> Result<Vec<String>> {
    let mut seen = HashSet::new();
    for bc in &barcodes {
        if bc.len() != 10 {
            bail!("Found non-10bp barcode example: {bc}");
        }
        if !bc.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
            bail!("Found barcode with non-ACGT bases example: {bc}");
        }
        if !seen.insert(bc.clone()) {
            bail!("Found duplicate barcode example: {bc}");
        }
    }
    Ok(barcodes)
}

fn write_whitelists(
    barcodes: &[String],
    middle_3p: &str,
    middle_5p: &str,
    out_3p: &Path,
    out_5p: &Path,
) -> Result<(usize, usize)> {
    let barcodes_5p: Vec<String> = barcodes.iter().map(|bc| reverse_complement(bc)).collect();
    if barcodes_5p.iter().collect::<HashSet<_>>().len() != barcodes_5p.len() {
        bail!("5p barcode reverse-complement list has duplicates.");
    }
    let mut writer_3p = File::create(out_3p)?;
    let mut count_3p = 0usize;
    for left in barcodes {
        for right in barcodes {
            writeln!(writer_3p, "{left}{middle_3p}{right}")?;
            count_3p += 1;
        }
    }

    let mut writer_5p = File::create(out_5p)?;
    let mut count_5p = 0usize;
    for left in &barcodes_5p {
        for right in &barcodes_5p {
            writeln!(writer_5p, "{left}{middle_5p}{right}")?;
            count_5p += 1;
        }
    }
    Ok((count_3p, count_5p))
}
