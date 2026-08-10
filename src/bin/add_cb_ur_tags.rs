use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::{ArgAction, Parser};
use rust_htslib::bam::{self, Read};

#[derive(Debug, Parser)]
#[command(
    version,
    about = "Add CB/UR and related Flora tags to BAM from read_tags.tsv"
)]
struct Cli {
    #[arg(long)]
    bam: PathBuf,

    #[arg(long)]
    tags: PathBuf,

    #[arg(long)]
    output: PathBuf,

    #[arg(long = "read-id-col", default_value = "read_id")]
    read_id_col: String,
    #[arg(long = "cb-col", default_value = "cell_id")]
    cb_col: String,
    #[arg(long = "ur-col", default_value = "umi_for_clustering")]
    ur_col: String,
    #[arg(long = "cr-col", default_value = "barcode_dual")]
    cr_col: String,
    #[arg(long = "bc5-col", default_value = "barcode_5p")]
    bc5_col: String,
    #[arg(long = "bc3-col", default_value = "barcode_3p")]
    bc3_col: String,
    #[arg(long = "umi5-col", default_value = "umi_primary")]
    umi5_col: String,
    #[arg(long = "umi3-col", default_value = "umi_backup")]
    umi3_col: String,
    #[arg(long = "keep-untagged", action = ArgAction::SetTrue)]
    keep_untagged: bool,
}

pub fn main() -> Result<()> {
    let cli = Cli::parse();
    let tag_map = load_tag_map(&cli)?;
    let mut bam =
        bam::Reader::from_path(&cli.bam).with_context(|| format!("open {}", cli.bam.display()))?;
    let header = bam::Header::from_template(bam.header());
    let mut out = bam::Writer::from_path(&cli.output, &header, bam::Format::Bam)
        .with_context(|| format!("create {}", cli.output.display()))?;

    let mut total = 0usize;
    let mut tagged = 0usize;
    let mut total_read_ids = HashSet::new();
    let mut tagged_read_ids = HashSet::new();
    for rec in bam.records() {
        let mut rec = rec?;
        total += 1;
        let rid = String::from_utf8_lossy(rec.qname()).to_string();
        total_read_ids.insert(rid.clone());
        if let Some(tags) = tag_map.get(&rid) {
            push_tag(&mut rec, b"CB", &tags.cb)?;
            push_tag(&mut rec, b"UR", &tags.ur)?;
            maybe_push_tag(&mut rec, b"CR", &tags.cr)?;
            maybe_push_tag(&mut rec, b"C5", &tags.c5)?;
            maybe_push_tag(&mut rec, b"C3", &tags.c3)?;
            maybe_push_tag(&mut rec, b"U5", &tags.u5)?;
            maybe_push_tag(&mut rec, b"U3", &tags.u3)?;
            tagged += 1;
            tagged_read_ids.insert(rid);
            out.write(&rec)?;
        } else if cli.keep_untagged {
            out.write(&rec)?;
        }
    }
    drop(out);
    drop(bam);
    bam::index::build(&cli.output, None, bam::index::Type::Bai, 1)?;
    println!("total_alignments\t{}", total_read_ids.len());
    println!("tagged_alignments\t{}", tagged_read_ids.len());
    println!(
        "unmatched_alignments\t{}",
        total_read_ids.len().saturating_sub(tagged_read_ids.len())
    );
    println!("total_alignment_records\t{total}");
    println!("tagged_alignment_records\t{tagged}");
    println!(
        "unmatched_alignment_records\t{}",
        total.saturating_sub(tagged)
    );
    Ok(())
}

#[derive(Clone)]
struct TagRow {
    cb: String,
    ur: String,
    cr: String,
    c5: String,
    c3: String,
    u5: String,
    u3: String,
}

fn load_tag_map(cli: &Cli) -> Result<HashMap<String, TagRow>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(&cli.tags)
        .with_context(|| format!("open {}", cli.tags.display()))?;
    let headers = reader.headers()?.clone();
    let idx = |name: &str| {
        headers
            .iter()
            .position(|h| h == name)
            .context(format!("missing column: {name}"))
    };
    let read_id_idx = idx(&cli.read_id_col)?;
    let cb_idx = idx(&cli.cb_col)?;
    let ur_idx = idx(&cli.ur_col)?;
    let cr_idx = idx(&cli.cr_col)?;
    let bc5_idx = idx(&cli.bc5_col)?;
    let bc3_idx = idx(&cli.bc3_col)?;
    let umi5_idx = idx(&cli.umi5_col)?;
    let umi3_idx = idx(&cli.umi3_col)?;

    let mut tag_map = HashMap::new();
    for row in reader.records() {
        let row = row?;
        let rid = clean(row.get(read_id_idx));
        if rid.is_empty() || tag_map.contains_key(&rid) {
            continue;
        }
        let cb = clean(row.get(cb_idx));
        let ur = clean(row.get(ur_idx));
        if cb.is_empty() || ur.is_empty() {
            continue;
        }
        tag_map.insert(
            rid,
            TagRow {
                cb,
                ur,
                cr: clean(row.get(cr_idx)),
                c5: clean(row.get(bc5_idx)),
                c3: clean(row.get(bc3_idx)),
                u5: clean(row.get(umi5_idx)),
                u3: clean(row.get(umi3_idx)),
            },
        );
    }
    Ok(tag_map)
}

fn clean(value: Option<&str>) -> String {
    value.unwrap_or("").trim().to_string()
}

fn push_tag(rec: &mut bam::Record, tag: &[u8; 2], value: &str) -> Result<()> {
    rec.update_aux(tag, bam::record::Aux::String(value))?;
    Ok(())
}

fn maybe_push_tag(rec: &mut bam::Record, tag: &[u8; 2], value: &str) -> Result<()> {
    if !value.is_empty() {
        push_tag(rec, tag, value)?;
    }
    Ok(())
}
