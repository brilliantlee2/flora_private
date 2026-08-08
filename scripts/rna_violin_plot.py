import argparse

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="per_cell_qc.tsv input")
    parser.add_argument("--output", required=True, help="Output PNG path")
    return parser.parse_args()


def main():
    args = parse_args()
    per_cell_qc = pd.read_csv(args.input, sep="\t")
    plot_df = per_cell_qc.melt(
        id_vars="barcode",
        value_vars=["genes", "umis", "mito_percent"],
        var_name="metric",
        value_name="value",
    )
    plot_df["metric"] = plot_df["metric"].map(
        {
            "genes": "Genes",
            "umis": "UMIs",
            "mito_percent": "Mito percent",
        }
    )

    sns.set(style="whitegrid", font_scale=1.2)
    fig, axes = plt.subplots(1, 3, figsize=(16, 6))
    for ax, metric in zip(axes, ["Genes", "UMIs", "Mito percent"]):
        sub = plot_df[plot_df["metric"] == metric]
        sns.violinplot(
            data=sub,
            x="metric",
            y="value",
            inner="box",
            color="#4C9BD5",
            cut=0,
            ax=ax,
        )
        ax.set_xlabel("")
        ax.set_ylabel(metric)
        ax.set_title(metric)

    fig.suptitle("RNA Violin Plot", fontsize=18)
    plt.tight_layout()
    plt.savefig(args.output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(args.output)


if __name__ == "__main__":
    main()
