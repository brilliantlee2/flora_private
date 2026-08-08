#!/usr/bin/env python3
import argparse
import csv
import gzip
import io
import zipfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.strip().translate(_RC)[::-1]


def read_whitelist(path, reverse_complement=True):
    values = []
    if path.endswith(".zip"):
        with zipfile.ZipFile(path) as zf:
            if len(zf.namelist()) != 1:
                raise ValueError(f"Expected one file inside {path}")
            with io.TextIOWrapper(zf.open(zf.namelist()[0]), encoding="utf-8") as handle:
                values = [line.strip() for line in handle if line.strip()]
    elif path.endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            values = [line.strip() for line in handle if line.strip()]
    else:
        with open(path, "r", encoding="utf-8") as handle:
            values = [line.strip() for line in handle if line.strip()]
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


def knee_plot(counts, threshold_value, out_fn):
    sorted_counts = sorted(counts, reverse=True)
    plt.figure(figsize=(8, 8))
    plt.title("Barcode rank plot (from high-quality putative BC)")
    plt.loglog(sorted_counts, marker="o", linestyle="", alpha=1, markersize=6)
    plt.xlabel("Barcodes")
    plt.ylabel("Read counts")
    plt.axhline(y=threshold_value, color="r", linestyle="--", label="cell calling threshold")
    plt.legend()
    plt.savefig(out_fn)
    plt.close()


def plot_one(counts_path, whitelist_path, out_path, exp_cells, reverse_complement):
    counts = read_counts(counts_path)
    whitelist = read_whitelist(whitelist_path, reverse_complement=reverse_complement)
    filtered_counts = [count for barcode, count in counts.items() if barcode in whitelist]
    t = threshold(filtered_counts, exp_cells)
    knee_plot(filtered_counts, t, out_path)


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
    parser.add_argument("--no-revcomp-whitelist", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    reverse_complement = not args.no_revcomp_whitelist
    plot_one(args.counts_3p, args.full_whitelist_3p, args.out_3p, args.exp_cells, reverse_complement)
    plot_one(args.counts_5p, args.full_whitelist_5p, args.out_5p, args.exp_cells, reverse_complement)


if __name__ == "__main__":
    main()
