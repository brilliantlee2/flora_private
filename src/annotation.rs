use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::process::{Command, Stdio};

use anyhow::{bail, Context, Result};
use rustc_hash::FxHashMap as HashMap;

#[derive(Clone, Debug)]
pub struct BedRecord {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
    pub name: String,
    pub score: i32,
    pub strand: String,
}

#[derive(Clone, Debug)]
pub struct GeneInterval {
    pub start: u32,
    pub end: u32,
    pub strand: String,
    pub label: String,
}

#[derive(Clone, Debug)]
pub struct ExonInterval {
    pub start: u32,
    pub end: u32,
    pub strand: String,
    pub transcript_id: String,
    pub gene_label: String,
}

#[derive(Clone, Debug)]
pub struct ReadAssignment {
    pub read_id: String,
    pub status: String,
    pub score: i32,
    pub gene: String,
}

#[derive(Clone, Debug)]
pub struct TranscriptAssignment {
    pub read_id: String,
    pub status: String,
    pub score: i32,
    pub gene: String,
    pub transcript_id: String,
}

pub type GeneIndex = HashMap<String, Vec<GeneInterval>>;
pub type ExonIndex = HashMap<String, Vec<ExonInterval>>;

pub fn load_gene_gtf(path: &Path) -> Result<GeneIndex> {
    let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut index: GeneIndex = HashMap::default();
    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let Some(fields) = parse_gtf_line(&line) else {
            continue;
        };
        if fields.feature != "gene" {
            continue;
        }
        let Some(start) = parse_gtf_start(fields.start) else {
            continue;
        };
        let Some(end) = parse_u32(fields.end) else {
            continue;
        };
        let Some(label) = extract_attr(fields.attrs, "gene_name")
            .or_else(|| extract_attr(fields.attrs, "gene_id"))
        else {
            continue;
        };
        index
            .entry(fields.chrom.to_string())
            .or_default()
            .push(GeneInterval {
                start,
                end,
                strand: fields.strand.to_string(),
                label,
            });
    }
    sort_gene_index(&mut index);
    Ok(index)
}

pub fn load_exon_gtf(path: &Path) -> Result<ExonIndex> {
    let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut index: ExonIndex = HashMap::default();
    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let Some(fields) = parse_gtf_line(&line) else {
            continue;
        };
        if fields.feature != "exon" {
            continue;
        }
        let Some(start) = parse_gtf_start(fields.start) else {
            continue;
        };
        let Some(end) = parse_u32(fields.end) else {
            continue;
        };
        let Some(transcript_id) = extract_attr(fields.attrs, "transcript_id") else {
            continue;
        };
        let gene_label = extract_attr(fields.attrs, "gene_name")
            .or_else(|| extract_attr(fields.attrs, "gene_id"))
            .unwrap_or_else(|| "NA".to_string());
        index
            .entry(fields.chrom.to_string())
            .or_default()
            .push(ExonInterval {
                start,
                end,
                strand: fields.strand.to_string(),
                transcript_id,
                gene_label,
            });
    }
    sort_exon_index(&mut index);
    Ok(index)
}

pub fn parse_bed6_line(line: &str) -> Option<BedRecord> {
    let mut fields = line.split('\t');
    Some(BedRecord {
        chrom: fields.next()?.to_string(),
        start: parse_u32(fields.next()?)?,
        end: parse_u32(fields.next()?)?,
        name: fields.next()?.to_string(),
        score: fields.next()?.parse().ok().unwrap_or(0),
        strand: fields.next()?.trim().to_string(),
    })
}

pub fn assign_gene(record: &BedRecord, genes: &GeneIndex, mapq: i32) -> ReadAssignment {
    if record.score < mapq {
        return ReadAssignment {
            read_id: record.name.clone(),
            status: "Unassigned_mapq".to_string(),
            score: record.score,
            gene: "NA".to_string(),
        };
    }
    let Some(intervals) = genes.get(&record.chrom) else {
        return no_gene_assignment(record);
    };
    let mut best_overlap = 0u32;
    let mut best_gene: Option<&str> = None;
    let mut ambiguous = false;
    for gene in candidate_genes(intervals, record.start, record.end) {
        if strand_mismatch(&record.strand, &gene.strand) {
            continue;
        }
        let overlap = overlap_bp(record.start, record.end, gene.start, gene.end);
        if overlap == 0 {
            continue;
        }
        if overlap > best_overlap {
            best_overlap = overlap;
            best_gene = Some(gene.label.as_str());
            ambiguous = false;
        } else if overlap == best_overlap {
            ambiguous = true;
        }
    }
    if best_overlap == 0 {
        return no_gene_assignment(record);
    }
    if ambiguous {
        return ReadAssignment {
            read_id: record.name.clone(),
            status: "Unassigned_ambiguous".to_string(),
            score: record.score,
            gene: "NA".to_string(),
        };
    }
    ReadAssignment {
        read_id: record.name.clone(),
        status: "Assigned".to_string(),
        score: record.score,
        gene: best_gene.unwrap_or("NA").to_string(),
    }
}

pub fn assign_transcript(
    read_blocks: &[BedRecord],
    exons: &ExonIndex,
    mapq: i32,
) -> TranscriptAssignment {
    let read_id = read_blocks
        .first()
        .map(|x| x.name.clone())
        .unwrap_or_default();
    let score = read_blocks.iter().map(|x| x.score).max().unwrap_or(0);
    if score < mapq {
        return TranscriptAssignment {
            read_id,
            status: "Unassigned_mapq".to_string(),
            score,
            gene: "NA".to_string(),
            transcript_id: "NA".to_string(),
        };
    }
    let mut support: HashMap<(&str, &str), u32> = HashMap::default();
    for block in read_blocks {
        let Some(intervals) = exons.get(&block.chrom) else {
            continue;
        };
        for exon in candidate_exons(intervals, block.start, block.end) {
            if strand_mismatch(&block.strand, &exon.strand) {
                continue;
            }
            let overlap = overlap_bp(block.start, block.end, exon.start, exon.end);
            if overlap == 0 {
                continue;
            }
            *support
                .entry((exon.transcript_id.as_str(), exon.gene_label.as_str()))
                .or_insert(0) += overlap;
        }
    }
    if support.is_empty() {
        return TranscriptAssignment {
            read_id,
            status: "Unassigned_no_features".to_string(),
            score,
            gene: "NA".to_string(),
            transcript_id: "NA".to_string(),
        };
    }
    let mut best_key: Option<(&str, &str)> = None;
    let mut best_overlap = 0u32;
    let mut ambiguous = false;
    for (key, overlap) in support {
        if overlap > best_overlap {
            best_overlap = overlap;
            best_key = Some(key);
            ambiguous = false;
        } else if overlap == best_overlap {
            ambiguous = true;
        }
    }
    if ambiguous {
        return TranscriptAssignment {
            read_id,
            status: "Unassigned_ambiguous".to_string(),
            score,
            gene: "NA".to_string(),
            transcript_id: "NA".to_string(),
        };
    }
    let (transcript_id, gene) = best_key.unwrap_or(("NA", "NA"));
    TranscriptAssignment {
        read_id,
        status: "Assigned".to_string(),
        score,
        gene: gene.to_string(),
        transcript_id: transcript_id.to_string(),
    }
}

pub fn bedtools_bamtobed_split(bam: &Path, output: &Path) -> Result<()> {
    let file = File::create(output)?;
    let status = Command::new("bedtools")
        .args(["bamtobed", "-split", "-i"])
        .arg(bam)
        .stdout(Stdio::from(file))
        .status()
        .with_context(|| "run bedtools bamtobed -split")?;
    if !status.success() {
        bail!("bedtools bamtobed -split failed with status {status}");
    }
    Ok(())
}

pub fn write_gene_assignment(mut writer: impl Write, row: &ReadAssignment) -> Result<()> {
    writeln!(
        writer,
        "{}\t{}\t{}\t{}",
        row.read_id, row.status, row.score, row.gene
    )?;
    Ok(())
}

pub fn write_transcript_assignment(
    mut writer: impl Write,
    row: &TranscriptAssignment,
) -> Result<()> {
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}",
        row.read_id, row.status, row.score, row.gene, row.transcript_id
    )?;
    Ok(())
}

fn parse_u32(value: &str) -> Option<u32> {
    value.parse().ok()
}

fn parse_gtf_start(value: &str) -> Option<u32> {
    let pos: u32 = value.parse().ok()?;
    pos.checked_sub(1)
}

fn extract_attr(attr: &str, key: &str) -> Option<String> {
    let needle = format!("{key} \"");
    let rest = attr.split_once(&needle)?.1;
    Some(rest.split_once('"')?.0.to_string())
}

struct GtfFields<'a> {
    chrom: &'a str,
    feature: &'a str,
    start: &'a str,
    end: &'a str,
    strand: &'a str,
    attrs: &'a str,
}

fn parse_gtf_line(line: &str) -> Option<GtfFields<'_>> {
    let mut fields = line.splitn(9, '\t');
    let chrom = fields.next()?;
    fields.next()?;
    let feature = fields.next()?;
    let start = fields.next()?;
    let end = fields.next()?;
    fields.next()?;
    let strand = fields.next()?;
    fields.next()?;
    let attrs = fields.next()?;
    Some(GtfFields {
        chrom,
        feature,
        start,
        end,
        strand,
        attrs,
    })
}

fn overlap_bp(a_start: u32, a_end: u32, b_start: u32, b_end: u32) -> u32 {
    let start = a_start.max(b_start);
    let end = a_end.min(b_end);
    end.saturating_sub(start)
}

fn strand_mismatch(read_strand: &str, feature_strand: &str) -> bool {
    matches!(read_strand, "+" | "-")
        && matches!(feature_strand, "+" | "-")
        && read_strand != feature_strand
}

fn no_gene_assignment(record: &BedRecord) -> ReadAssignment {
    ReadAssignment {
        read_id: record.name.clone(),
        status: "Unassigned_no_features".to_string(),
        score: record.score,
        gene: "NA".to_string(),
    }
}

fn sort_gene_index(index: &mut GeneIndex) {
    for intervals in index.values_mut() {
        intervals.sort_by_key(|x| x.start);
    }
}

fn sort_exon_index(index: &mut ExonIndex) {
    for intervals in index.values_mut() {
        intervals.sort_by_key(|x| x.start);
    }
}

fn candidate_genes(
    intervals: &[GeneInterval],
    start: u32,
    end: u32,
) -> impl Iterator<Item = &GeneInterval> {
    let start_idx = intervals.partition_point(|x| x.start < end);
    intervals[..start_idx]
        .iter()
        .rev()
        .take_while(move |x| x.end > start)
}

fn candidate_exons(
    intervals: &[ExonInterval],
    start: u32,
    end: u32,
) -> impl Iterator<Item = &ExonInterval> {
    let start_idx = intervals.partition_point(|x| x.start < end);
    intervals[..start_idx]
        .iter()
        .rev()
        .take_while(move |x| x.end > start)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn assign_gene_prefers_largest_same_strand_overlap() {
        let mut index = GeneIndex::default();
        index.insert(
            "chr1".to_string(),
            vec![
                GeneInterval {
                    start: 10,
                    end: 30,
                    strand: "+".to_string(),
                    label: "gene_a".to_string(),
                },
                GeneInterval {
                    start: 18,
                    end: 40,
                    strand: "-".to_string(),
                    label: "gene_b".to_string(),
                },
                GeneInterval {
                    start: 15,
                    end: 50,
                    strand: "+".to_string(),
                    label: "gene_c".to_string(),
                },
            ],
        );
        sort_gene_index(&mut index);
        let row = assign_gene(
            &BedRecord {
                chrom: "chr1".to_string(),
                start: 20,
                end: 35,
                name: "read1".to_string(),
                score: 60,
                strand: "+".to_string(),
            },
            &index,
            60,
        );
        assert_eq!(row.status, "Assigned");
        assert_eq!(row.gene, "gene_c");
    }

    #[test]
    fn assign_gene_marks_equal_overlap_ambiguous() {
        let mut index = GeneIndex::default();
        index.insert(
            "chr1".to_string(),
            vec![
                GeneInterval {
                    start: 10,
                    end: 20,
                    strand: "+".to_string(),
                    label: "gene_a".to_string(),
                },
                GeneInterval {
                    start: 20,
                    end: 30,
                    strand: "+".to_string(),
                    label: "gene_b".to_string(),
                },
            ],
        );
        sort_gene_index(&mut index);
        let row = assign_gene(
            &BedRecord {
                chrom: "chr1".to_string(),
                start: 15,
                end: 25,
                name: "read_tie".to_string(),
                score: 60,
                strand: "+".to_string(),
            },
            &index,
            60,
        );
        assert_eq!(row.status, "Unassigned_ambiguous");
        assert_eq!(row.gene, "NA");
    }

    #[test]
    fn assign_gene_low_mapq_takes_priority() {
        let mut index = GeneIndex::default();
        index.insert(
            "chr1".to_string(),
            vec![GeneInterval {
                start: 10,
                end: 50,
                strand: "+".to_string(),
                label: "gene_a".to_string(),
            }],
        );
        sort_gene_index(&mut index);
        let row = assign_gene(
            &BedRecord {
                chrom: "chr1".to_string(),
                start: 12,
                end: 18,
                name: "low_mapq".to_string(),
                score: 10,
                strand: "+".to_string(),
            },
            &index,
            60,
        );
        assert_eq!(row.status, "Unassigned_mapq");
        assert_eq!(row.gene, "NA");
    }

    #[test]
    fn assign_transcript_marks_tied_support_ambiguous() {
        let mut index = ExonIndex::default();
        index.insert(
            "chr1".to_string(),
            vec![
                ExonInterval {
                    start: 10,
                    end: 20,
                    strand: "+".to_string(),
                    transcript_id: "tx1".to_string(),
                    gene_label: "g1".to_string(),
                },
                ExonInterval {
                    start: 10,
                    end: 20,
                    strand: "+".to_string(),
                    transcript_id: "tx2".to_string(),
                    gene_label: "g1".to_string(),
                },
            ],
        );
        sort_exon_index(&mut index);
        let row = assign_transcript(
            &[BedRecord {
                chrom: "chr1".to_string(),
                start: 12,
                end: 18,
                name: "read1".to_string(),
                score: 60,
                strand: "+".to_string(),
            }],
            &index,
            60,
        );
        assert_eq!(row.status, "Unassigned_ambiguous");
        assert_eq!(row.transcript_id, "NA");
    }

    #[test]
    fn assign_transcript_sums_split_block_support_like_python() {
        let mut index = ExonIndex::default();
        index.insert(
            "chr1".to_string(),
            vec![
                ExonInterval {
                    start: 10,
                    end: 20,
                    strand: "+".to_string(),
                    transcript_id: "tx1".to_string(),
                    gene_label: "g1".to_string(),
                },
                ExonInterval {
                    start: 30,
                    end: 40,
                    strand: "+".to_string(),
                    transcript_id: "tx1".to_string(),
                    gene_label: "g1".to_string(),
                },
                ExonInterval {
                    start: 10,
                    end: 20,
                    strand: "+".to_string(),
                    transcript_id: "tx2".to_string(),
                    gene_label: "g1".to_string(),
                },
                ExonInterval {
                    start: 30,
                    end: 35,
                    strand: "+".to_string(),
                    transcript_id: "tx2".to_string(),
                    gene_label: "g1".to_string(),
                },
            ],
        );
        sort_exon_index(&mut index);
        let row = assign_transcript(
            &[
                BedRecord {
                    chrom: "chr1".to_string(),
                    start: 12,
                    end: 18,
                    name: "read_blocks".to_string(),
                    score: 60,
                    strand: "+".to_string(),
                },
                BedRecord {
                    chrom: "chr1".to_string(),
                    start: 31,
                    end: 39,
                    name: "read_blocks".to_string(),
                    score: 60,
                    strand: "+".to_string(),
                },
            ],
            &index,
            60,
        );
        assert_eq!(row.status, "Assigned");
        assert_eq!(row.gene, "g1");
        assert_eq!(row.transcript_id, "tx1");
    }

    #[test]
    fn assign_transcript_low_mapq_takes_priority() {
        let mut index = ExonIndex::default();
        index.insert(
            "chr1".to_string(),
            vec![ExonInterval {
                start: 10,
                end: 20,
                strand: "+".to_string(),
                transcript_id: "tx1".to_string(),
                gene_label: "g1".to_string(),
            }],
        );
        sort_exon_index(&mut index);
        let row = assign_transcript(
            &[BedRecord {
                chrom: "chr1".to_string(),
                start: 12,
                end: 18,
                name: "low_mapq_tx".to_string(),
                score: 10,
                strand: "+".to_string(),
            }],
            &index,
            60,
        );
        assert_eq!(row.status, "Unassigned_mapq");
        assert_eq!(row.transcript_id, "NA");
    }
}
