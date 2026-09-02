use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::process::{Command, Stdio};

use anyhow::{bail, Context, Result};
use rust_htslib::bam;
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

trait IntervalCoordinates {
    fn start(&self) -> u32;
    fn end(&self) -> u32;
}

impl IntervalCoordinates for GeneInterval {
    fn start(&self) -> u32 {
        self.start
    }

    fn end(&self) -> u32 {
        self.end
    }
}

impl IntervalCoordinates for ExonInterval {
    fn start(&self) -> u32 {
        self.start
    }

    fn end(&self) -> u32 {
        self.end
    }
}

#[derive(Clone, Debug)]
pub struct IntervalIndex<T> {
    intervals: Vec<T>,
    prefix_max_end: Vec<u32>,
}

impl<T> IntervalIndex<T> {
    fn new(mut intervals: Vec<T>) -> Self
    where
        T: IntervalCoordinates,
    {
        intervals.sort_by_key(IntervalCoordinates::start);
        let mut max_end = 0;
        let prefix_max_end = intervals
            .iter()
            .map(|interval| {
                max_end = max_end.max(interval.end());
                max_end
            })
            .collect();
        Self {
            intervals,
            prefix_max_end,
        }
    }

    fn candidates(&self, start: u32, end: u32) -> impl Iterator<Item = &T>
    where
        T: IntervalCoordinates,
    {
        let right = self.intervals.partition_point(|x| x.start() < end);
        let left = self.prefix_max_end[..right].partition_point(|x| *x <= start);
        self.intervals[left..right]
            .iter()
            .filter(move |x| x.end() > start)
    }
}

impl From<Vec<GeneInterval>> for IntervalIndex<GeneInterval> {
    fn from(intervals: Vec<GeneInterval>) -> Self {
        Self::new(intervals)
    }
}

impl From<Vec<ExonInterval>> for IntervalIndex<ExonInterval> {
    fn from(intervals: Vec<ExonInterval>) -> Self {
        Self::new(intervals)
    }
}

pub type GeneIndex = HashMap<String, IntervalIndex<GeneInterval>>;
pub type ExonIndex = HashMap<String, IntervalIndex<ExonInterval>>;

pub fn load_gene_gtf(path: &Path) -> Result<GeneIndex> {
    let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut intervals_by_chrom: HashMap<String, Vec<GeneInterval>> = HashMap::default();
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
        intervals_by_chrom
            .entry(fields.chrom.to_string())
            .or_default()
            .push(GeneInterval {
                start,
                end,
                strand: fields.strand.to_string(),
                label,
            });
    }
    Ok(intervals_by_chrom
        .into_iter()
        .map(|(chrom, intervals)| (chrom, intervals.into()))
        .collect())
}

pub fn load_exon_gtf(path: &Path) -> Result<ExonIndex> {
    let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut intervals_by_chrom: HashMap<String, Vec<ExonInterval>> = HashMap::default();
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
        intervals_by_chrom
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
    Ok(intervals_by_chrom
        .into_iter()
        .map(|(chrom, intervals)| (chrom, intervals.into()))
        .collect())
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

/// Reproduce the default, unsplit `bedtools bamtobed` BED6 projection.
pub fn bed_record_from_bam(record: &bam::Record, header: &bam::HeaderView) -> Option<BedRecord> {
    if record.is_unmapped() || record.tid() < 0 || record.pos() < 0 {
        return None;
    }
    let tid = record.tid() as u32;
    let chrom = std::str::from_utf8(header.tid2name(tid)).ok()?.to_string();
    let mut name = String::from_utf8_lossy(record.qname()).to_string();
    if record.is_first_in_template() {
        name.push_str("/1");
    }
    if record.is_last_in_template() {
        name.push_str("/2");
    }
    let start = u32::try_from(record.pos()).ok()?;
    let end = u32::try_from(record.cigar().end_pos()).ok()?;
    Some(BedRecord {
        chrom,
        start,
        end,
        name,
        score: i32::from(record.mapq()),
        strand: if record.is_reverse() { "-" } else { "+" }.to_string(),
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

fn candidate_genes(
    intervals: &IntervalIndex<GeneInterval>,
    start: u32,
    end: u32,
) -> impl Iterator<Item = &GeneInterval> {
    intervals.candidates(start, end)
}

fn candidate_exons(
    intervals: &IntervalIndex<ExonInterval>,
    start: u32,
    end: u32,
) -> impl Iterator<Item = &ExonInterval> {
    intervals.candidates(start, end)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam;
    use rust_htslib::bam::header::HeaderRecord;
    use rust_htslib::bam::record::{Cigar, CigarString};

    fn next_random(state: &mut u64) -> u32 {
        *state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1);
        (*state >> 32) as u32
    }

    fn bam_header() -> bam::HeaderView {
        let mut header = bam::Header::new();
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 1_000_000),
        );
        bam::HeaderView::from_header(&header)
    }

    fn bam_record(cigar: Vec<Cigar>, pos: i64) -> bam::Record {
        let cigar = CigarString(cigar);
        let query_len = cigar
            .iter()
            .map(|op| match op {
                Cigar::Match(n)
                | Cigar::Ins(n)
                | Cigar::SoftClip(n)
                | Cigar::Equal(n)
                | Cigar::Diff(n) => *n as usize,
                _ => 0,
            })
            .sum::<usize>();
        let mut record = bam::Record::new();
        record.set(
            b"read1",
            Some(&cigar),
            &vec![b'A'; query_len],
            &vec![30; query_len],
        );
        record.set_tid(0);
        record.set_pos(pos);
        record.set_mapq(60);
        record.unset_unmapped();
        record
    }

    #[test]
    fn bam_record_matches_unsplit_bamtobed_coordinates() {
        let header = bam_header();
        let record = bam_record(
            vec![
                Cigar::SoftClip(5),
                Cigar::Match(10),
                Cigar::Ins(3),
                Cigar::Match(5),
                Cigar::Del(4),
                Cigar::RefSkip(100),
                Cigar::Equal(6),
                Cigar::Diff(2),
            ],
            100,
        );
        let bed = bed_record_from_bam(&record, &header).unwrap();
        assert_eq!(bed.chrom, "chr1");
        assert_eq!(bed.start, 100);
        assert_eq!(bed.end, 227);
        assert_eq!(bed.name, "read1");
        assert_eq!(bed.score, 60);
        assert_eq!(bed.strand, "+");
    }

    #[test]
    fn bam_record_matches_bamtobed_name_strand_and_unmapped_rules() {
        let header = bam_header();
        let mut reverse_first = bam_record(vec![Cigar::Match(20)], 500);
        reverse_first.set_reverse();
        reverse_first.set_first_in_template();
        let bed = bed_record_from_bam(&reverse_first, &header).unwrap();
        assert_eq!(bed.name, "read1/1");
        assert_eq!(bed.strand, "-");

        reverse_first.set_unmapped();
        assert!(bed_record_from_bam(&reverse_first, &header).is_none());
    }

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
            ]
            .into(),
        );
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
            ]
            .into(),
        );
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
            }]
            .into(),
        );
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
    fn assign_gene_finds_outer_gene_past_nonoverlapping_inner_gene() {
        let mut index = GeneIndex::default();
        index.insert(
            "chr1".to_string(),
            vec![
                GeneInterval {
                    start: 10,
                    end: 1_000,
                    strand: "+".to_string(),
                    label: "outer_gene".to_string(),
                },
                GeneInterval {
                    start: 100,
                    end: 200,
                    strand: "+".to_string(),
                    label: "inner_gene".to_string(),
                },
            ]
            .into(),
        );

        let row = assign_gene(
            &BedRecord {
                chrom: "chr1".to_string(),
                start: 500,
                end: 600,
                name: "nested_gene_read".to_string(),
                score: 60,
                strand: "+".to_string(),
            },
            &index,
            60,
        );

        assert_eq!(row.status, "Assigned");
        assert_eq!(row.gene, "outer_gene");
    }

    #[test]
    fn candidate_genes_match_brute_force_for_random_intervals() {
        let mut state = 0x5eed_1234_u64;
        let intervals = IntervalIndex::new(
            (0..256)
                .map(|idx| {
                    let start = next_random(&mut state) % 10_000;
                    let length = 1 + next_random(&mut state) % 2_000;
                    GeneInterval {
                        start,
                        end: start + length,
                        strand: "+".to_string(),
                        label: format!("gene_{idx}"),
                    }
                })
                .collect::<Vec<_>>(),
        );

        for _ in 0..512 {
            let start = next_random(&mut state) % 10_000;
            let end = start + 1 + next_random(&mut state) % 500;
            let mut actual = candidate_genes(&intervals, start, end)
                .map(|x| x.label.as_str())
                .collect::<Vec<_>>();
            let mut expected = intervals
                .intervals
                .iter()
                .filter(|x| x.start < end && x.end > start)
                .map(|x| x.label.as_str())
                .collect::<Vec<_>>();
            actual.sort_unstable();
            expected.sort_unstable();
            assert_eq!(actual, expected, "query interval {start}..{end}");
        }
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
            ]
            .into(),
        );
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
            ]
            .into(),
        );
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
            }]
            .into(),
        );
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

    #[test]
    fn assign_transcript_finds_outer_exon_past_nonoverlapping_inner_exon() {
        let mut index = ExonIndex::default();
        index.insert(
            "chr1".to_string(),
            vec![
                ExonInterval {
                    start: 10,
                    end: 1_000,
                    strand: "+".to_string(),
                    transcript_id: "outer_tx".to_string(),
                    gene_label: "outer_gene".to_string(),
                },
                ExonInterval {
                    start: 100,
                    end: 200,
                    strand: "+".to_string(),
                    transcript_id: "inner_tx".to_string(),
                    gene_label: "inner_gene".to_string(),
                },
            ]
            .into(),
        );

        let row = assign_transcript(
            &[BedRecord {
                chrom: "chr1".to_string(),
                start: 500,
                end: 600,
                name: "nested_exon_read".to_string(),
                score: 60,
                strand: "+".to_string(),
            }],
            &index,
            60,
        );

        assert_eq!(row.status, "Assigned");
        assert_eq!(row.gene, "outer_gene");
        assert_eq!(row.transcript_id, "outer_tx");
    }

    #[test]
    fn candidate_exons_match_brute_force_for_random_intervals() {
        let mut state = 0x5eed_5678_u64;
        let intervals = IntervalIndex::new(
            (0..256)
                .map(|idx| {
                    let start = next_random(&mut state) % 10_000;
                    let length = 1 + next_random(&mut state) % 2_000;
                    ExonInterval {
                        start,
                        end: start + length,
                        strand: "+".to_string(),
                        transcript_id: format!("tx_{idx}"),
                        gene_label: format!("gene_{idx}"),
                    }
                })
                .collect::<Vec<_>>(),
        );

        for _ in 0..512 {
            let start = next_random(&mut state) % 10_000;
            let end = start + 1 + next_random(&mut state) % 500;
            let mut actual = candidate_exons(&intervals, start, end)
                .map(|x| x.transcript_id.as_str())
                .collect::<Vec<_>>();
            let mut expected = intervals
                .intervals
                .iter()
                .filter(|x| x.start < end && x.end > start)
                .map(|x| x.transcript_id.as_str())
                .collect::<Vec<_>>();
            actual.sort_unstable();
            expected.sort_unstable();
            assert_eq!(actual, expected, "query interval {start}..{end}");
        }
    }
}
