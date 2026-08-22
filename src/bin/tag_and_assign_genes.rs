use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Context, Result};
use clap::Parser;
use flora::annotation::{assign_gene, bed_record_from_bam, load_gene_gtf, write_gene_assignment};
use flora::bam_runtime::bounded_hts_threads;
use rust_htslib::bam::{self, Read};
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

#[derive(Debug, Parser)]
#[command(version, about = "Add Flora barcode/UMI and gene tags in one BAM pass")]
struct Cli {
    #[arg(long)]
    bam: PathBuf,

    #[arg(long)]
    tags: PathBuf,

    #[arg(long)]
    gtf: PathBuf,

    #[arg(long = "gene-assigns")]
    gene_assigns: PathBuf,

    #[arg(long)]
    output: PathBuf,

    #[arg(short = 'q', long = "mapq", default_value_t = 60)]
    mapq: i32,

    #[arg(short = 't', long = "threads", default_value_t = 4)]
    threads: usize,
}

pub fn main() -> Result<()> {
    run(Cli::parse())
}

fn run(cli: Cli) -> Result<()> {
    let total_phase = Instant::now();
    let phase = Instant::now();
    let tag_map = load_tag_map(&cli.tags)?;
    eprintln!(
        "[timing] tag_and_assign.load_tags: {:.2}s",
        phase.elapsed().as_secs_f64()
    );
    let phase = Instant::now();
    let genes = load_gene_gtf(&cli.gtf)?;
    eprintln!(
        "[timing] tag_and_assign.load_gtf: {:.2}s",
        phase.elapsed().as_secs_f64()
    );

    let threads = bounded_hts_threads(cli.threads);
    let mut input =
        bam::Reader::from_path(&cli.bam).with_context(|| format!("open {}", cli.bam.display()))?;
    input.set_threads(threads)?;
    let header_view = input.header().clone();
    let header = bam::Header::from_template(&header_view);
    let mut output = bam::Writer::from_path(&cli.output, &header, bam::Format::Bam)
        .with_context(|| format!("create {}", cli.output.display()))?;
    output.set_threads(threads)?;
    let mut assignments = BufWriter::new(
        File::create(&cli.gene_assigns)
            .with_context(|| format!("create {}", cli.gene_assigns.display()))?,
    );

    let phase = Instant::now();
    let mut total_records = 0usize;
    let mut tagged_records = 0usize;
    let mut written_records = 0usize;
    let mut total_read_ids = HashSet::default();
    let mut tagged_read_ids = HashSet::default();
    for record in input.records() {
        let mut record = record?;
        total_records += 1;
        let read_id = String::from_utf8_lossy(record.qname()).to_string();
        total_read_ids.insert(read_id.clone());
        let Some(tags) = tag_map.get(&read_id) else {
            continue;
        };
        add_read_tags(&mut record, tags)?;
        tagged_records += 1;
        tagged_read_ids.insert(read_id);

        // Default bedtools bamtobed skips unmapped records. The legacy
        // add_gene_tags stage therefore stops before writing them as well.
        let Some(bed_record) = bed_record_from_bam(&record, &header_view) else {
            continue;
        };
        let assignment = assign_gene(&bed_record, &genes, cli.mapq);
        write_gene_assignment(&mut assignments, &assignment)?;
        record.update_aux(b"GN", bam::record::Aux::String(&assignment.gene))?;
        output.write(&record)?;
        written_records += 1;
    }
    drop(assignments);
    drop(output);
    drop(input);
    eprintln!(
        "[timing] tag_and_assign.read_tag_assign_write: {:.2}s",
        phase.elapsed().as_secs_f64()
    );

    let phase = Instant::now();
    bam::index::build(&cli.output, None, bam::index::Type::Bai, threads as u32)?;
    eprintln!(
        "[timing] tag_and_assign.index: {:.2}s",
        phase.elapsed().as_secs_f64()
    );
    eprintln!(
        "[timing] tag_and_assign.total: {:.2}s",
        total_phase.elapsed().as_secs_f64()
    );

    println!("total_alignments\t{}", total_read_ids.len());
    println!("tagged_alignments\t{}", tagged_read_ids.len());
    println!(
        "unmatched_alignments\t{}",
        total_read_ids.len().saturating_sub(tagged_read_ids.len())
    );
    println!("total_alignment_records\t{total_records}");
    println!("tagged_alignment_records\t{tagged_records}");
    println!(
        "unmatched_alignment_records\t{}",
        total_records.saturating_sub(tagged_records)
    );
    println!("gene_tagged_alignment_records\t{written_records}");
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

fn load_tag_map(path: &Path) -> Result<HashMap<String, TagRow>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("open {}", path.display()))?;
    let headers = reader.headers()?.clone();
    let idx = |name: &str| {
        headers
            .iter()
            .position(|header| header == name)
            .with_context(|| format!("missing column: {name}"))
    };
    let read_id_idx = idx("read_id")?;
    let cb_idx = idx("cell_id")?;
    let ur_idx = idx("umi_for_clustering")?;
    let cr_idx = idx("barcode_dual")?;
    let c5_idx = idx("barcode_5p")?;
    let c3_idx = idx("barcode_3p")?;
    let u5_idx = idx("umi_primary")?;
    let u3_idx = idx("umi_backup")?;

    let mut tag_map = HashMap::default();
    for row in reader.records() {
        let row = row?;
        let read_id = clean(row.get(read_id_idx));
        if read_id.is_empty() || tag_map.contains_key(&read_id) {
            continue;
        }
        let cb = clean(row.get(cb_idx));
        let ur = clean(row.get(ur_idx));
        if cb.is_empty() || ur.is_empty() {
            continue;
        }
        tag_map.insert(
            read_id,
            TagRow {
                cb,
                ur,
                cr: clean(row.get(cr_idx)),
                c5: clean(row.get(c5_idx)),
                c3: clean(row.get(c3_idx)),
                u5: clean(row.get(u5_idx)),
                u3: clean(row.get(u3_idx)),
            },
        );
    }
    Ok(tag_map)
}

fn add_read_tags(record: &mut bam::Record, tags: &TagRow) -> Result<()> {
    push_tag(record, b"CB", &tags.cb)?;
    push_tag(record, b"UR", &tags.ur)?;
    maybe_push_tag(record, b"CR", &tags.cr)?;
    maybe_push_tag(record, b"C5", &tags.c5)?;
    maybe_push_tag(record, b"C3", &tags.c3)?;
    maybe_push_tag(record, b"U5", &tags.u5)?;
    maybe_push_tag(record, b"U3", &tags.u3)?;
    Ok(())
}

fn clean(value: Option<&str>) -> String {
    value.unwrap_or("").trim().to_string()
}

fn push_tag(record: &mut bam::Record, tag: &[u8; 2], value: &str) -> Result<()> {
    record.update_aux(tag, bam::record::Aux::String(value))?;
    Ok(())
}

fn maybe_push_tag(record: &mut bam::Record, tag: &[u8; 2], value: &str) -> Result<()> {
    if !value.is_empty() {
        push_tag(record, tag, value)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::header::HeaderRecord;
    use rust_htslib::bam::record::{Aux, Cigar, CigarString};
    use tempfile::tempdir;

    #[test]
    fn fused_stage_writes_the_same_tags_and_assignment_shape() {
        let temp = tempdir().unwrap();
        let input_path = temp.path().join("input.bam");
        let output_path = temp.path().join("output.bam");
        let tags_path = temp.path().join("tags.tsv");
        let gtf_path = temp.path().join("genes.gtf");
        let assigns_path = temp.path().join("assigns.tsv");

        let mut header = bam::Header::new();
        header.push_record(HeaderRecord::new(b"HD").push_tag(b"SO", "coordinate"));
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 1000),
        );
        let mut writer = bam::Writer::from_path(&input_path, &header, bam::Format::Bam).unwrap();
        let mut record = bam::Record::new();
        record.set(
            b"read1",
            Some(&CigarString(vec![Cigar::Match(20)])),
            &[b'A'; 20],
            &[30; 20],
        );
        record.set_tid(0);
        record.set_pos(100);
        record.set_mapq(60);
        record.unset_unmapped();
        writer.write(&record).unwrap();
        drop(writer);

        std::fs::write(
            &tags_path,
            "read_id\tcell_id\tumi_for_clustering\tbarcode_dual\tbarcode_5p\tbarcode_3p\tumi_primary\tumi_backup\nread1\tCELL1\tUMI1\tDUAL\tBC5\tBC3\tU5\tU3\n",
        )
        .unwrap();
        std::fs::write(
            &gtf_path,
            "chr1\ttest\tgene\t91\t130\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\";\n",
        )
        .unwrap();

        run(Cli {
            bam: input_path,
            tags: tags_path,
            gtf: gtf_path,
            gene_assigns: assigns_path.clone(),
            output: output_path.clone(),
            mapq: 60,
            threads: 2,
        })
        .unwrap();

        assert_eq!(
            std::fs::read_to_string(assigns_path).unwrap(),
            "read1\tAssigned\t60\tGENE1\n"
        );
        let mut reader = bam::Reader::from_path(output_path).unwrap();
        let records = reader.records().collect::<Result<Vec<_>, _>>().unwrap();
        assert_eq!(records.len(), 1);
        let record = &records[0];
        assert!(matches!(record.aux(b"CB").unwrap(), Aux::String("CELL1")));
        assert!(matches!(record.aux(b"UR").unwrap(), Aux::String("UMI1")));
        assert!(matches!(record.aux(b"GN").unwrap(), Aux::String("GENE1")));
    }
}
