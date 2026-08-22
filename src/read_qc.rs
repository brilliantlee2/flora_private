use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

use anyhow::{bail, Context, Result};
use rustc_hash::FxHashMap as HashMap;
use serde::Serialize;

const READ_QUALITY_TRIM_LENGTH: usize = 10;
const READ_LENGTH_HIST_UPPER_QUANTILE: f64 = 0.995;
const ACCUMULATOR_MAGIC: &[u8; 8] = b"FLRQCA01";

static QUALITY_ERROR_PROBABILITIES: OnceLock<[f64; 256]> = OnceLock::new();

fn quality_error_probabilities() -> &'static [f64; 256] {
    QUALITY_ERROR_PROBABILITIES.get_or_init(|| {
        std::array::from_fn(|quality| {
            10_f64.powf((u8::saturating_sub(quality as u8, 33)) as f64 * -0.1)
        })
    })
}

#[derive(Default)]
pub struct ReadQcAccumulator {
    length_counts: HashMap<u32, usize>,
    mean_qualities: Vec<f32>,
    total_bases: usize,
    n_reads: usize,
}

#[derive(Debug, Serialize)]
pub struct ReadQcPayload {
    n_reads: usize,
    quality: QualityPayload,
    length: LengthPayload,
    yield_above_length: YieldPayload,
}

#[derive(Debug, Serialize)]
struct QualityPayload {
    mean: f64,
    median: f64,
    bins: Vec<f64>,
    counts: Vec<usize>,
}

#[derive(Debug, Serialize)]
struct LengthPayload {
    mean: f64,
    median: f64,
    min: usize,
    max: usize,
    bins_kb: Vec<f64>,
    counts: Vec<usize>,
}

#[derive(Debug, Serialize)]
struct YieldPayload {
    x_kb: Vec<f64>,
    y_gb: Vec<f64>,
    total_gb: f64,
    n50_bp: usize,
}

impl ReadQcAccumulator {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn observe(&mut self, seq: &[u8], qual: &[u8]) {
        self.n_reads += 1;
        if seq.is_empty() {
            return;
        }
        *self.length_counts.entry(seq.len() as u32).or_default() += 1;
        self.total_bases += seq.len();
        self.mean_qualities.push(mean_quality(qual) as f32);
    }

    pub fn merge(&mut self, mut other: Self) {
        self.n_reads += other.n_reads;
        self.total_bases += other.total_bases;
        for (length, count) in other.length_counts.drain() {
            *self.length_counts.entry(length).or_default() += count;
        }
        self.mean_qualities.append(&mut other.mean_qualities);
    }

    pub fn n_reads(&self) -> usize {
        self.n_reads
    }

    pub fn write_accumulator(&self, path: &Path) -> Result<()> {
        let mut writer = BufWriter::new(
            File::create(path).with_context(|| format!("create {}", path.display()))?,
        );
        writer.write_all(ACCUMULATOR_MAGIC)?;
        write_u64(&mut writer, self.n_reads as u64)?;
        write_u64(&mut writer, self.total_bases as u64)?;
        let mut lengths: Vec<_> = self.length_counts.iter().collect();
        lengths.sort_unstable_by_key(|(length, _)| **length);
        write_u64(&mut writer, lengths.len() as u64)?;
        for (&length, &count) in lengths {
            writer.write_all(&length.to_le_bytes())?;
            write_u64(&mut writer, count as u64)?;
        }
        write_u64(&mut writer, self.mean_qualities.len() as u64)?;
        for quality in &self.mean_qualities {
            writer.write_all(&quality.to_bits().to_le_bytes())?;
        }
        writer.flush()?;
        Ok(())
    }

    pub fn read_accumulator(path: &Path) -> Result<Self> {
        let mut reader =
            BufReader::new(File::open(path).with_context(|| format!("open {}", path.display()))?);
        let mut magic = [0u8; 8];
        reader.read_exact(&mut magic)?;
        if &magic != ACCUMULATOR_MAGIC {
            bail!("invalid Read QC accumulator: {}", path.display());
        }
        let n_reads = usize_from_u64(read_u64(&mut reader)?)?;
        let total_bases = usize_from_u64(read_u64(&mut reader)?)?;
        let n_lengths = usize_from_u64(read_u64(&mut reader)?)?;
        let mut length_counts = HashMap::default();
        for _ in 0..n_lengths {
            let mut bytes = [0u8; 4];
            reader.read_exact(&mut bytes)?;
            let length = u32::from_le_bytes(bytes);
            let count = usize_from_u64(read_u64(&mut reader)?)?;
            length_counts.insert(length, count);
        }
        let n_qualities = usize_from_u64(read_u64(&mut reader)?)?;
        let mut mean_qualities = Vec::with_capacity(n_qualities);
        for _ in 0..n_qualities {
            let mut bytes = [0u8; 4];
            reader.read_exact(&mut bytes)?;
            mean_qualities.push(f32::from_bits(u32::from_le_bytes(bytes)));
        }
        if mean_qualities.len() != length_counts.values().sum::<usize>() {
            bail!("inconsistent Read QC accumulator: {}", path.display());
        }
        Ok(Self {
            length_counts,
            mean_qualities,
            total_bases,
            n_reads,
        })
    }

    pub fn write_json(self, path: &Path, curve_points: usize) -> Result<()> {
        atomic_write(
            path,
            serde_json::to_string(&self.into_payload(curve_points))?.as_bytes(),
        )
    }

    pub fn write_outputs(
        self,
        json_path: &Path,
        count_path: &Path,
        curve_points: usize,
    ) -> Result<()> {
        let n_reads = self.n_reads();
        atomic_write(
            json_path,
            serde_json::to_string(&self.into_payload(curve_points))?.as_bytes(),
        )?;
        atomic_write(count_path, format!("{n_reads}\n").as_bytes())
    }

    pub fn into_payload(mut self, curve_points: usize) -> ReadQcPayload {
        if self.length_counts.is_empty() {
            return empty_payload(self.n_reads);
        }
        let mut lengths: Vec<(u32, usize)> = self.length_counts.into_iter().collect();
        lengths.sort_unstable_by_key(|(length, _)| *length);
        self.mean_qualities.sort_by(|a, b| a.total_cmp(b));
        let q_min = (self.mean_qualities[0].floor() as f64).max(0.0) - 1.0;
        let q_max = self.mean_qualities[self.mean_qualities.len() - 1].ceil() as f64 + 1.0;
        let (q_bins, q_counts) = histogram_f32(&self.mean_qualities, q_min, q_max, 60);
        let len_hist_max_kb =
            quantile_length_counts(&lengths, READ_LENGTH_HIST_UPPER_QUANTILE).max(0.5);
        let (len_bins, len_counts) = histogram_length_counts(&lengths, len_hist_max_kb, 80);
        let (x_kb, y_gb, n50_bp) =
            yield_curve_length_counts(&lengths, self.total_bases, curve_points);
        let observed_reads = lengths.iter().map(|(_, count)| *count).sum::<usize>();

        ReadQcPayload {
            n_reads: self.n_reads,
            quality: QualityPayload {
                mean: mean_sorted_f32(&self.mean_qualities),
                median: median_sorted_f32(&self.mean_qualities),
                bins: q_bins,
                counts: q_counts,
            },
            length: LengthPayload {
                mean: self.total_bases as f64 / observed_reads as f64,
                median: median_length_counts(&lengths),
                min: lengths[0].0 as usize,
                max: lengths[lengths.len() - 1].0 as usize,
                bins_kb: len_bins,
                counts: len_counts,
            },
            yield_above_length: YieldPayload {
                x_kb,
                y_gb,
                total_gb: self.total_bases as f64 / 1e9,
                n50_bp,
            },
        }
    }
}

fn atomic_write(path: &Path, content: &[u8]) -> Result<()> {
    let temporary = temporary_path(path);
    std::fs::write(&temporary, content)
        .with_context(|| format!("write {}", temporary.display()))?;
    std::fs::rename(&temporary, path).with_context(|| format!("publish {}", path.display()))?;
    Ok(())
}

fn temporary_path(path: &Path) -> PathBuf {
    let name = path
        .file_name()
        .and_then(|v| v.to_str())
        .unwrap_or("read-qc");
    path.with_file_name(format!(".{name}.{}.tmp", std::process::id()))
}

fn write_u64<W: Write>(writer: &mut W, value: u64) -> Result<()> {
    writer.write_all(&value.to_le_bytes())?;
    Ok(())
}

fn read_u64<R: Read>(reader: &mut R) -> Result<u64> {
    let mut bytes = [0u8; 8];
    reader.read_exact(&mut bytes)?;
    Ok(u64::from_le_bytes(bytes))
}

fn usize_from_u64(value: u64) -> Result<usize> {
    usize::try_from(value).context("Read QC accumulator is too large for this platform")
}

fn empty_payload(n_reads: usize) -> ReadQcPayload {
    ReadQcPayload {
        n_reads,
        quality: QualityPayload {
            mean: 0.0,
            median: 0.0,
            bins: vec![],
            counts: vec![],
        },
        length: LengthPayload {
            mean: 0.0,
            median: 0.0,
            min: 0,
            max: 0,
            bins_kb: vec![],
            counts: vec![],
        },
        yield_above_length: YieldPayload {
            x_kb: vec![],
            y_gb: vec![],
            total_gb: 0.0,
            n50_bp: 0,
        },
    }
}

fn mean_quality(quality: &[u8]) -> f64 {
    let (start, end) = if quality.len() >= 100 && quality.len() > 2 * READ_QUALITY_TRIM_LENGTH {
        (
            READ_QUALITY_TRIM_LENGTH,
            quality.len() - READ_QUALITY_TRIM_LENGTH,
        )
    } else {
        (0, quality.len())
    };
    if start >= end {
        return 0.0;
    }
    let lookup = quality_error_probabilities();
    let error_sum = quality[start..end]
        .iter()
        .map(|q| lookup[*q as usize])
        .sum::<f64>();
    let mean_error_rate = error_sum / (end - start) as f64;
    if mean_error_rate > 0.0 {
        -10.0 * mean_error_rate.log10()
    } else {
        0.0
    }
}

fn histogram_f32(
    values: &[f32],
    min_value: f64,
    max_value: f64,
    bins: usize,
) -> (Vec<f64>, Vec<usize>) {
    let width = ((max_value - min_value) / bins as f64).max(f64::EPSILON);
    let mut counts = vec![0usize; bins];
    for value in values {
        let idx = (((*value as f64 - min_value) / width).floor() as usize).min(bins - 1);
        counts[idx] += 1;
    }
    (
        (0..bins)
            .map(|idx| min_value + (idx as f64 + 0.5) * width)
            .collect(),
        counts,
    )
}

fn histogram_length_counts(
    values: &[(u32, usize)],
    max_kb: f64,
    bins: usize,
) -> (Vec<f64>, Vec<usize>) {
    let width = (max_kb / bins as f64).max(f64::EPSILON);
    let mut counts = vec![0usize; bins];
    for &(value, count) in values {
        let idx = (((value as f64 / 1000.0).min(max_kb) / width).floor() as usize).min(bins - 1);
        counts[idx] += count;
    }
    (
        (0..bins).map(|idx| (idx as f64 + 0.5) * width).collect(),
        counts,
    )
}

fn yield_curve_length_counts(
    lengths: &[(u32, usize)],
    total_bases: usize,
    curve_points: usize,
) -> (Vec<f64>, Vec<f64>, usize) {
    let n = lengths.iter().map(|(_, count)| *count).sum::<usize>();
    let points = curve_points.max(1).min(n);
    let mut cumulative = 0usize;
    let mut n50_bp = lengths[0].0 as usize;
    let half = total_bases.div_ceil(2);
    for &(length, count) in lengths.iter().rev() {
        let bases = length as usize * count;
        if cumulative < half && cumulative + bases >= half {
            n50_bp = length as usize;
        }
        cumulative += bases;
    }
    let mut x = Vec::new();
    let mut y = Vec::new();
    let mut last = None;
    for point in 0..points {
        let idx = if points == 1 {
            0
        } else {
            point * (n - 1) / (points - 1)
        };
        if last == Some(idx) {
            continue;
        }
        last = Some(idx);
        let (length, bases_before) = descending_length_at(lengths, idx);
        x.push(length as f64 / 1000.0);
        y.push((total_bases - bases_before) as f64 / 1e9);
    }
    (x, y, n50_bp)
}

fn mean_sorted_f32(values: &[f32]) -> f64 {
    values.iter().map(|x| *x as f64).sum::<f64>() / values.len() as f64
}

fn median_sorted_f32(values: &[f32]) -> f64 {
    let m = values.len() / 2;
    if values.len() % 2 == 0 {
        (values[m - 1] as f64 + values[m] as f64) / 2.0
    } else {
        values[m] as f64
    }
}

fn median_length_counts(values: &[(u32, usize)]) -> f64 {
    let count = values.iter().map(|(_, count)| *count).sum::<usize>();
    let left = value_at_rank(values, (count - 1) / 2);
    let right = value_at_rank(values, count / 2);
    (left + right) as f64 / 2.0
}

fn quantile_length_counts(values: &[(u32, usize)], q: f64) -> f64 {
    let count = values.iter().map(|(_, count)| *count).sum::<usize>();
    let position = (count - 1) as f64 * q.clamp(0.0, 1.0);
    let lo = position.floor() as usize;
    let hi = position.ceil() as usize;
    let weight = position - lo as f64;
    let lower = value_at_rank(values, lo) as f64 / 1000.0;
    let upper = value_at_rank(values, hi) as f64 / 1000.0;
    lower + (upper - lower) * weight
}

fn value_at_rank(values: &[(u32, usize)], rank: usize) -> u32 {
    let mut seen = 0usize;
    for &(value, count) in values {
        if rank < seen + count {
            return value;
        }
        seen += count;
    }
    values.last().map(|(value, _)| *value).unwrap_or(0)
}

fn descending_length_at(values: &[(u32, usize)], rank: usize) -> (u32, usize) {
    let mut seen = 0usize;
    let mut bases_before = 0usize;
    for &(value, count) in values.iter().rev() {
        if rank < seen + count {
            bases_before += (rank - seen) * value as usize;
            return (value, bases_before);
        }
        seen += count;
        bases_before += count * value as usize;
    }
    (values[0].0, bases_before)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn merged_payload_matches_single_pass() {
        let records: &[(&[u8], &[u8])] = &[
            (b"AAAA", b"IIII"),
            (b"AAAAAA", b"######"),
            (b"", b""),
            (&[b'A'; 120], &[b'I'; 120]),
        ];
        let mut single = ReadQcAccumulator::new();
        for (seq, qual) in records {
            single.observe(seq, qual);
        }
        let mut left = ReadQcAccumulator::new();
        let mut right = ReadQcAccumulator::new();
        for (seq, qual) in &records[..2] {
            left.observe(seq, qual);
        }
        for (seq, qual) in &records[2..] {
            right.observe(seq, qual);
        }
        left.merge(right);
        assert_eq!(
            serde_json::to_value(single.into_payload(300)).unwrap(),
            serde_json::to_value(left.into_payload(300)).unwrap()
        );
    }

    #[test]
    fn binary_accumulator_round_trip_is_lossless() {
        let mut qc = ReadQcAccumulator::new();
        qc.observe(b"AAAA", b"IIII");
        qc.observe(b"AAAAAA", b"######");
        let expected = serde_json::to_value(qc.into_payload(300)).unwrap();

        let mut qc = ReadQcAccumulator::new();
        qc.observe(b"AAAA", b"IIII");
        qc.observe(b"AAAAAA", b"######");
        let tmp = tempfile::NamedTempFile::new().unwrap();
        qc.write_accumulator(tmp.path()).unwrap();
        let restored = ReadQcAccumulator::read_accumulator(tmp.path()).unwrap();
        assert_eq!(
            expected,
            serde_json::to_value(restored.into_payload(300)).unwrap()
        );
    }

    #[test]
    fn quality_lookup_matches_original_summary_formula() {
        let lookup = quality_error_probabilities();
        for quality in 0_u16..=255 {
            let score = (quality as u8).saturating_sub(33) as f64;
            assert_eq!(lookup[quality as usize], 10_f64.powf(score * -0.1));
        }
    }
}
