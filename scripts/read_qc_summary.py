import argparse
import gzip
import json
from array import array

import numpy as np

READ_QUALITY_TRIM_LENGTH = 10
READ_LENGTH_HIST_UPPER_QUANTILE = 0.995


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq", required=True, help="FASTQ or FASTQ.GZ to summarize")
    parser.add_argument("--output-json", required=True, help="Output JSON path")
    parser.add_argument(
        "--output-fastq-count",
        help="Optional output path for the FASTQ record count",
    )
    parser.add_argument(
        "--curve-points",
        type=int,
        default=300,
        help="Maximum points to keep for the yield-above-length curve [300]",
    )
    return parser.parse_args()


def fastq_iter(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        while True:
            title = next(handle, None)
            if title is None:
                break
            seq = next(handle).strip()
            next(handle)
            qual = next(handle).strip()
            yield seq, qual


def sampled_curve_from_lengths(lengths_sorted_desc, curve_points):
    cumulative = np.cumsum(lengths_sorted_desc, dtype=np.int64)
    total = int(cumulative[-1])
    if len(lengths_sorted_desc) <= curve_points:
        idx = np.arange(len(lengths_sorted_desc), dtype=int)
    else:
        idx = np.linspace(0, len(lengths_sorted_desc) - 1, num=curve_points, dtype=int)
        idx = np.unique(idx)
    x_bp = lengths_sorted_desc[idx]
    y_gb = (total - np.where(idx > 0, cumulative[idx - 1], 0)) / 1e9
    return x_bp, y_gb


def get_mean_quality_score(quality_scores, trim_quality_length=READ_QUALITY_TRIM_LENGTH):
    read_length = len(quality_scores)
    if read_length >= 100 and read_length > 2 * trim_quality_length:
        trimmed_scores = quality_scores[trim_quality_length:-trim_quality_length]
    else:
        trimmed_scores = quality_scores

    if len(trimmed_scores) == 0:
        return 0.0

    error_rates = 10 ** (trimmed_scores * -0.1)
    mean_error_rate = np.mean(error_rates)
    return -10 * np.log10(mean_error_rate) if mean_error_rate > 0 else 0.0


def main():
    args = parse_args()

    lengths = array("I")
    mean_qualities = array("f")
    total_bases = 0

    for seq, qual in fastq_iter(args.fastq):
        length = len(seq)
        if length == 0:
            continue
        lengths.append(length)
        total_bases += length
        q_array = np.fromiter((ord(ch) - 33 for ch in qual), dtype=np.float32, count=length)
        q_mean = get_mean_quality_score(q_array)
        mean_qualities.append(q_mean)

    if len(lengths) == 0:
        payload = {
            "n_reads": 0,
            "quality": {"mean": 0.0, "median": 0.0, "bins": [], "counts": []},
            "length": {
                "mean": 0.0,
                "median": 0.0,
                "min": 0,
                "max": 0,
                "bins_kb": [],
                "counts": [],
            },
            "yield_above_length": {
                "x_kb": [],
                "y_gb": [],
                "total_gb": 0.0,
                "n50_bp": 0,
            },
        }
    else:
        lengths_np = np.array(lengths, dtype=np.int32)
        qualities_np = np.array(mean_qualities, dtype=np.float32)

        q_min = max(0.0, float(np.floor(qualities_np.min())) - 1.0)
        q_max = float(np.ceil(qualities_np.max())) + 1.0
        q_bins = np.linspace(q_min, q_max, num=61)
        q_counts, q_edges = np.histogram(qualities_np, bins=q_bins)
        q_centers = ((q_edges[:-1] + q_edges[1:]) / 2.0).tolist()

        len_kb = lengths_np / 1000.0
        len_hist_max_kb = float(
            max(np.quantile(len_kb, READ_LENGTH_HIST_UPPER_QUANTILE), 0.5)
        )
        len_bins = np.linspace(0, len_hist_max_kb, num=81)
        len_counts, len_edges = np.histogram(np.clip(len_kb, None, len_hist_max_kb), bins=len_bins)
        len_centers = ((len_edges[:-1] + len_edges[1:]) / 2.0).tolist()

        sorted_lengths_desc = np.sort(lengths_np)[::-1]
        x_bp, y_gb = sampled_curve_from_lengths(sorted_lengths_desc, args.curve_points)
        cumulative = np.cumsum(sorted_lengths_desc, dtype=np.int64)
        n50_idx = int(np.searchsorted(cumulative, total_bases / 2.0, side="left"))
        n50_bp = int(sorted_lengths_desc[min(n50_idx, len(sorted_lengths_desc) - 1)])

        payload = {
            "n_reads": int(len(lengths_np)),
            "quality": {
                "mean": float(qualities_np.mean()),
                "median": float(np.median(qualities_np)),
                "bins": q_centers,
                "counts": q_counts.astype(int).tolist(),
            },
            "length": {
                "mean": float(lengths_np.mean()),
                "median": float(np.median(lengths_np)),
                "min": int(lengths_np.min()),
                "max": int(lengths_np.max()),
                "bins_kb": len_centers,
                "counts": len_counts.astype(int).tolist(),
            },
            "yield_above_length": {
                "x_kb": (x_bp / 1000.0).tolist(),
                "y_gb": y_gb.tolist(),
                "total_gb": float(total_bases / 1e9),
                "n50_bp": n50_bp,
            },
        }

    with open(args.output_json, "w", encoding="utf-8") as handle:
        json.dump(payload, handle)

    if args.output_fastq_count:
        with open(args.output_fastq_count, "w", encoding="utf-8") as handle:
            handle.write(f"{payload['n_reads']}\n")

    print(args.output_json)


if __name__ == "__main__":
    main()
