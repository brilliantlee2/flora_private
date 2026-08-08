use std::path::PathBuf;

use anyhow::Result;
use clap::Parser;
use serde::Serialize;
use flora::fastq::for_each_fastq_batch;

const READ_QUALITY_TRIM_LENGTH: usize = 10;
const READ_LENGTH_HIST_UPPER_QUANTILE: f64 = 0.995;

#[derive(Debug, Parser)]
#[command(version, about = "Summarize FASTQ read length and quality distributions")]
struct Cli {
    #[arg(long)]
    fastq: PathBuf,

    #[arg(long = "output-json")]
    output_json: PathBuf,

    #[arg(long = "curve-points", default_value_t = 300)]
    curve_points: usize,

    #[arg(long = "batch-size", default_value_t = 100_000)]
    batch_size: usize,
}

#[derive(Serialize)]
struct Payload {
    n_reads: usize,
    quality: QualityPayload,
    length: LengthPayload,
    yield_above_length: YieldPayload,
}

#[derive(Serialize)]
struct QualityPayload {
    mean: f64,
    median: f64,
    bins: Vec<f64>,
    counts: Vec<usize>,
}

#[derive(Serialize)]
struct LengthPayload {
    mean: f64,
    median: f64,
    min: usize,
    max: usize,
    bins_kb: Vec<f64>,
    counts: Vec<usize>,
}

#[derive(Serialize)]
struct YieldPayload {
    x_kb: Vec<f64>,
    y_gb: Vec<f64>,
    total_gb: f64,
    n50_bp: usize,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let mut lengths: Vec<u32> = Vec::new();
    let mut mean_qualities: Vec<f32> = Vec::new();
    let mut total_bases = 0usize;

    for_each_fastq_batch(&[cli.fastq.clone()], cli.batch_size, |batch| {
        for rec in batch {
            let length = rec.seq.len();
            if length == 0 {
                continue;
            }
            lengths.push(length as u32);
            total_bases += length;
            mean_qualities.push(mean_quality(&rec.qual) as f32);
        }
        Ok(())
    })?;

    let payload = build_payload(lengths, mean_qualities, total_bases, cli.curve_points);
    std::fs::write(&cli.output_json, serde_json::to_string(&payload)?)?;
    println!("{}", cli.output_json.display());
    Ok(())
}

fn mean_quality(quality: &str) -> f64 {
    let bytes = quality.as_bytes();
    let read_length = bytes.len();
    let (start, end) = if read_length >= 100 && read_length > 2 * READ_QUALITY_TRIM_LENGTH {
        (READ_QUALITY_TRIM_LENGTH, read_length - READ_QUALITY_TRIM_LENGTH)
    } else {
        (0, read_length)
    };
    if start >= end {
        return 0.0;
    }
    let trimmed_len = end - start;
    let mut error_sum = 0.0f64;
    for &q in &bytes[start..end] {
        let score = (q.saturating_sub(33)) as f64;
        error_sum += 10_f64.powf(score * -0.1);
    }
    let mean_error_rate = error_sum / trimmed_len as f64;
    if mean_error_rate > 0.0 {
        -10.0 * mean_error_rate.log10()
    } else {
        0.0
    }
}

fn build_payload(
    mut lengths: Vec<u32>,
    mut qualities: Vec<f32>,
    total_bases: usize,
    curve_points: usize,
) -> Payload {
    if lengths.is_empty() {
        return Payload {
            n_reads: 0,
            quality: QualityPayload { mean: 0.0, median: 0.0, bins: vec![], counts: vec![] },
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
        };
    }

    lengths.sort_unstable();
    qualities.sort_by(|a, b| a.total_cmp(b));
    let q_min = (qualities[0].floor() as f64).max(0.0) - 1.0;
    let q_max = qualities[qualities.len() - 1].ceil() as f64 + 1.0;
    let (q_bins, q_counts) = histogram_f32(&qualities, q_min, q_max, 60);

    let len_hist_max_kb =
        quantile_sorted_linear_u32(&lengths, READ_LENGTH_HIST_UPPER_QUANTILE).max(0.5);
    let (len_bins, len_counts) = histogram_u32_clipped_kb(&lengths, len_hist_max_kb, 80);

    let (x_kb, y_gb, n50_bp) = yield_curve_sorted_asc(&lengths, total_bases, curve_points);

    Payload {
        n_reads: lengths.len(),
        quality: QualityPayload {
            mean: mean_sorted_f32(&qualities),
            median: median_sorted_f32(&qualities),
            bins: q_bins,
            counts: q_counts,
        },
        length: LengthPayload {
            mean: lengths.iter().map(|x| *x as f64).sum::<f64>() / lengths.len() as f64,
            median: median_u32_sorted(&lengths),
            min: lengths[0] as usize,
            max: lengths[lengths.len() - 1] as usize,
            bins_kb: len_bins,
            counts: len_counts,
        },
        yield_above_length: YieldPayload {
            x_kb,
            y_gb,
            total_gb: total_bases as f64 / 1e9,
            n50_bp,
        },
    }
}

fn histogram_f32(values: &[f32], min_value: f64, max_value: f64, bins: usize) -> (Vec<f64>, Vec<usize>) {
    let width = ((max_value - min_value) / bins as f64).max(f64::EPSILON);
    let mut counts = vec![0usize; bins];
    for value in values {
        let mut idx = (((*value as f64) - min_value) / width).floor() as usize;
        if idx >= bins {
            idx = bins - 1;
        }
        counts[idx] += 1;
    }
    let centers = (0..bins)
        .map(|idx| min_value + (idx as f64 + 0.5) * width)
        .collect();
    (centers, counts)
}

fn histogram_u32_clipped_kb(values: &[u32], max_kb: f64, bins: usize) -> (Vec<f64>, Vec<usize>) {
    let width = (max_kb / bins as f64).max(f64::EPSILON);
    let mut counts = vec![0usize; bins];
    for &value in values {
        let clipped = (value as f64 / 1000.0).min(max_kb);
        let mut idx = (clipped / width).floor() as usize;
        if idx >= bins {
            idx = bins - 1;
        }
        counts[idx] += 1;
    }
    let centers = (0..bins)
        .map(|idx| (idx as f64 + 0.5) * width)
        .collect();
    (centers, counts)
}

fn yield_curve_sorted_asc(lengths_asc: &[u32], total_bases: usize, curve_points: usize) -> (Vec<f64>, Vec<f64>, usize) {
    let n = lengths_asc.len();
    let points = curve_points.max(1).min(n);
    let mut cumulative_before_desc = vec![0usize; n];
    let mut x_kb = Vec::new();
    let mut y_gb = Vec::new();
    let mut cumulative = 0usize;
    let mut n50_bp = lengths_asc[0] as usize;
    let half = total_bases.div_ceil(2);
    for (desc_idx, &length) in lengths_asc.iter().rev().enumerate() {
        cumulative_before_desc[desc_idx] = cumulative;
        let length = length as usize;
        if cumulative < half && cumulative + length >= half {
            n50_bp = length;
        }
        cumulative += length;
    }
    let mut last_idx = None;
    for point in 0..points {
        let idx = if points == 1 {
            0
        } else {
            point * (n - 1) / (points - 1)
        };
        if last_idx == Some(idx) {
            continue;
        }
        last_idx = Some(idx);
        let asc_idx = n - 1 - idx;
        x_kb.push(lengths_asc[asc_idx] as f64 / 1000.0);
        y_gb.push((total_bases - cumulative_before_desc[idx]) as f64 / 1e9);
    }
    (x_kb, y_gb, n50_bp)
}

fn mean_sorted_f32(values: &[f32]) -> f64 {
    values.iter().map(|x| *x as f64).sum::<f64>() / values.len() as f64
}

fn median_sorted_f32(values: &[f32]) -> f64 {
    let mid = values.len() / 2;
    if values.len() % 2 == 0 {
        (values[mid - 1] as f64 + values[mid] as f64) / 2.0
    } else {
        values[mid] as f64
    }
}

fn median_u32_sorted(values: &[u32]) -> f64 {
    let mid = values.len() / 2;
    if values.len() % 2 == 0 {
        (values[mid - 1] + values[mid]) as f64 / 2.0
    } else {
        values[mid] as f64
    }
}

fn quantile_sorted_linear_u32(values: &[u32], q: f64) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let position = (values.len() - 1) as f64 * q.clamp(0.0, 1.0);
    let lower = position.floor() as usize;
    let upper = position.ceil() as usize;
    if lower == upper {
        return values[lower] as f64 / 1000.0;
    }
    let weight = position - lower as f64;
    let lower_value = values[lower] as f64 / 1000.0;
    let upper_value = values[upper] as f64 / 1000.0;
    lower_value + (upper_value - lower_value) * weight
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linear_quantile_matches_numpy_style_interpolation() {
        let values = vec![1000u32, 2000, 4000, 8000];
        let q = quantile_sorted_linear_u32(&values, 0.5);
        assert!((q - 3.0).abs() < 1e-9);
    }

    #[test]
    fn yield_curve_matches_descending_cumulative_intent() {
        let lengths = vec![100u32, 200, 300, 400];
        let (x_kb, y_gb, n50_bp) = yield_curve_sorted_asc(&lengths, 1000, 4);
        assert_eq!(n50_bp, 300);
        assert_eq!(x_kb.len(), y_gb.len());
        assert_eq!(x_kb.first().copied(), Some(0.4));
        assert!((y_gb.first().copied().unwrap_or_default() - 0.000001).abs() < 1e-12);
    }
}
