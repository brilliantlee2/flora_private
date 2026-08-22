#!/usr/bin/env python3
import argparse
import re

import pandas as pd

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        default="cell_umi_gene.tsv",
        help="Input cell/gene/UMI read table TSV [cell_umi_gene.tsv]",
    )
    parser.add_argument(
        "--output-tsv",
        default="saturation.tsv",
        help="Output saturation table TSV [saturation.tsv]",
    )
    parser.add_argument(
        "--output-png",
        default="saturation_curves.png",
        help="Output saturation plot PNG [saturation_curves.png]",
    )
    parser.add_argument(
        "--plot-existing-tsv",
        default=None,
        help="Plot an existing saturation TSV without loading read-level input.",
    )
    return parser.parse_args()


def is_unannotated_region_label(value):
    return bool(re.search(r"[A-Za-z0-9]+_\d+_\d+", str(value)))


def is_known_gene(value):
    gene = str(value).strip()
    if gene in {"", "NA", "nan", "None"}:
        return False
    return not is_unannotated_region_label(gene)


def median_known_genes_per_cell(df_sub):
    all_cells = pd.Index(df_sub["barcode"].dropna().astype(str).unique())
    if len(all_cells) == 0:
        return 0.0
    gene_rows = df_sub[df_sub["gene"].map(is_known_gene)]
    counts = gene_rows.groupby("barcode")["gene"].nunique()
    return float(counts.reindex(all_cells, fill_value=0).median())


def main():
    args = parse_args()
    if args.plot_existing_tsv:
        res = pd.read_csv(args.plot_existing_tsv, sep="\t")
    else:
        df = pd.read_csv(args.input, sep="\t")
        df = df.copy()
        df["barcode"] = df["barcode"].astype(str)
        df["gene"] = df["gene"].astype(str)
        df["umi"] = df["umi"].astype(str)
        df["read_id"] = df["read_id"].astype(str)

        fractions = [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        records = []

        for frac in fractions:
            df_sub = df.sample(frac=frac, random_state=42) if frac < 1.0 else df
            n_reads = df_sub.shape[0]
            if n_reads == 0:
                records.append((frac, 0, 0.0, 0.0, 0.0, 0.0))
                continue

            gene_bc_umi = df_sub["gene"] + "_" + df_sub["barcode"] + "_" + df_sub["umi"]
            n_deduped_reads = gene_bc_umi.nunique()
            saturation = 1 - (n_deduped_reads / n_reads)

            genes_per_cell = median_known_genes_per_cell(df_sub)
            umis_per_cell = float(df_sub.groupby("barcode")["umi"].nunique().median())
            reads_per_cell = float(df_sub.groupby("barcode")["read_id"].nunique().median())

            records.append(
                (frac, n_reads, reads_per_cell, genes_per_cell, umis_per_cell, saturation)
            )

        res = pd.DataFrame(
            records,
            columns=[
                "fraction",
                "reads",
                "reads_per_cell",
                "genes_per_cell",
                "umis_per_cell",
                "saturation",
            ],
        )
        res.to_csv(args.output_tsv, sep="\t", index=False)

    if plt is None:
        print("Warning: matplotlib is not installed; skipped saturation PNG output.")
        print(res.to_csv(sep="\t", index=False).strip())
        return

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    axes[0].plot(res["reads_per_cell"], res["genes_per_cell"], marker="o")
    axes[0].set_xlabel("Median reads per cell")
    axes[0].set_ylabel("Median genes per cell")
    axes[0].set_title("Genes per cell")

    axes[1].plot(res["reads_per_cell"], res["umis_per_cell"], marker="o")
    axes[1].set_xlabel("Median reads per cell")
    axes[1].set_ylabel("Median UMIs per cell")
    axes[1].set_title("UMIs per cell")

    axes[2].plot(res["reads_per_cell"], res["saturation"], marker="o")
    axes[2].set_xlabel("Median reads per cell")
    axes[2].set_ylabel("Sequencing saturation")
    axes[2].set_title("Saturation")

    plt.tight_layout()
    plt.savefig(args.output_png, dpi=300)
    plt.close(fig)

    print(res.to_csv(sep="\t", index=False).strip())


if __name__ == "__main__":
    main()
