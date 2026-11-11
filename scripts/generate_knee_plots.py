#!/usr/bin/env python3
import argparse
import csv
import gzip
import io
import json
import math
import zipfile
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


_RC = defaultdict(lambda: ord("N"), str.maketrans("ACGTNacgtn", "TGCANTGCAN"))


def revcomp(seq):
    return seq.strip().translate(_RC)[::-1]


def read_whitelist(path, reverse_complement=True):
    values = []
    if path.endswith(".zip"):
        with zipfile.ZipFile(path) as zf:
            if not zf.namelist():
                raise ValueError(f"Empty whitelist archive: {path}")
            with io.TextIOWrapper(zf.open(zf.namelist()[0]), encoding="utf-8") as handle:
                values = [line.strip() for line in handle if line.strip()]
    elif path.endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            values = [line.strip() for line in handle if line.strip()]
    else:
        with open(path, "r", encoding="utf-8") as handle:
            values = [line.strip() for line in handle if line.strip()]
    values = [value.split("-", 1)[0].strip() for value in values]
    values = [value for value in values if value]
    if reverse_complement:
        values = [revcomp(value) for value in values]
    return set(values)


def read_counts(path):
    counts = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            barcode = (row.get("barcode") or "").strip()
            if not barcode:
                continue
            counts[barcode] = int(row.get("count") or 0)
    return counts


def threshold(counts, exp_cells):
    if not counts:
        return 0.0
    top = np.sort(np.array(list(counts), dtype=float))[::-1][: max(int(exp_cells), 1)]
    return float(np.quantile(top, 0.95) / 20.0)


def load_threshold_report(path):
    with open(path, encoding="utf-8") as handle:
        report = json.load(handle)
    if report.get("schema_version") != 1:
        raise ValueError("unsupported barcode threshold report schema")
    sides = report.get("sides", {})
    for side in ("3p", "5p"):
        item = sides.get(side, {})
        value = item.get("applied_threshold")
        if (item.get("status") != "ready" or type(value) not in (int, float)
                or not math.isfinite(value) or value < 0 or item.get("comparison") not in (">", ">=")):
            raise ValueError(f"{side}: no valid applied whitelist threshold in {path}")
    return sides


def knee_plot(counts, threshold_value, out_fn, side_report=None):
    sorted_counts = sorted(counts, reverse=True)
    plt.figure(figsize=(8, 8))
    plt.title("Barcode rank plot (from high-quality putative BC)")
    plt.loglog(range(1, len(sorted_counts) + 1), sorted_counts, marker="o", linestyle="", alpha=1, markersize=6)
    plt.xlabel("Barcodes")
    plt.ylabel("Read counts")
    if side_report is not None:
        legacy = side_report.get("legacy_threshold")
        alternative = side_report.get("steepest", {}).get("threshold")
        for value, label, color in ((legacy, "legacy threshold", "#777777"),
                                    (alternative, "steepest candidate", "#e69f00")):
            if value is not None and value > 0 and value != threshold_value:
                plt.axhline(y=value, color=color, linestyle=":", label=label)
    plt.axhline(y=threshold_value, color="r", linestyle="--", label="applied whitelist threshold")
    plt.legend()
    plt.savefig(out_fn)
    plt.close()


def plot_one(counts_path, whitelist_path, out_path, exp_cells, reverse_complement, side_report=None):
    counts = read_counts(counts_path)
    whitelist = read_whitelist(whitelist_path, reverse_complement=reverse_complement)
    filtered_counts = [count for barcode, count in counts.items() if barcode in whitelist]
    if side_report is None:
        t = threshold(filtered_counts, exp_cells)
    else:
        t = side_report["applied_threshold"]
        selected = sum(count >= t if side_report["comparison"] == ">=" else count > t for count in filtered_counts)
        if (side_report.get("observed_theoretical_barcodes") != len(filtered_counts)
                or side_report.get("selected") != selected):
            raise ValueError("barcode counts/whitelist disagree with the authoritative threshold report")
    knee_plot(filtered_counts, t, out_path, side_report)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate Flora barcode knee plots with the validated cell-calling threshold logic."
    )
    parser.add_argument("--counts-3p", required=True)
    parser.add_argument("--counts-5p", required=True)
    parser.add_argument("--full-whitelist-3p", required=True)
    parser.add_argument("--full-whitelist-5p", required=True)
    parser.add_argument("--out-3p", required=True)
    parser.add_argument("--out-5p", required=True)
    parser.add_argument("--exp-cells", type=int, default=5000)
    parser.add_argument("--threshold-report", help="authoritative barcode_thresholds.json from this run or preview")
    parser.add_argument("--no-revcomp-whitelist", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    reverse_complement = not args.no_revcomp_whitelist
    sides = load_threshold_report(args.threshold_report) if args.threshold_report else {}
    plot_one(args.counts_3p, args.full_whitelist_3p, args.out_3p, args.exp_cells, reverse_complement, sides.get("3p"))
    plot_one(args.counts_5p, args.full_whitelist_5p, args.out_5p, args.exp_cells, reverse_complement, sides.get("5p"))


if __name__ == "__main__":
    main()
