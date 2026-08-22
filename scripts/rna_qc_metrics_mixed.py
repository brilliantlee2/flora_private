import argparse
import os
import gzip
import re

import matplotlib.pyplot as plt
import pandas as pd
import pysam
import seaborn as sns


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cell-umi-gene-tsv",
        default="cell_umi_gene.tsv",
        help="Cell/gene/UMI read table TSV [cell_umi_gene.tsv]",
    )
    parser.add_argument(
        "--bam",
        default="filtered.sorted.bam",
        help="Aligned BAM used for QC-derived per-cell summaries [filtered.sorted.bam]",
    )
    parser.add_argument(
        "--raw-fastq",
        required=True,
        help="Input FASTQ used as the denominator for raw-read fractions.",
    )
    parser.add_argument(
        "--full-length-fastq",
        default=None,
        help="Full-length FASTQ after Glycine. Used for Sockeye-style full-length/pass read summaries.",
    )
    parser.add_argument(
        "--read-tags",
        default=None,
        help="Stint2.1 read_tags.tsv. Used to count reads assigned to final cells with usable UMIs.",
    )
    parser.add_argument(
        "--barcode-validity-tsv",
        default=None,
        help="Upstream barcode validity summary TSV. Used to count true barcode-valid reads before cell assignment.",
    )
    parser.add_argument(
        "--transcript-assigns",
        default=None,
        help="Read transcript assignments TSV. If provided, report high-confidence known-transcript read counts.",
    )
    parser.add_argument(
        "--glycine-log",
        default=None,
        help="Optional Glycine stdout/stderr log. When present, the top-level Glycine Full-length count is used for full-length/pass-read metrics.",
    )
    parser.add_argument(
        "--glycine-stats",
        default=None,
        help="Merged Glycine identifying_statistic.txt used for total raw and full-length counts.",
    )
    parser.add_argument(
        "--mixed-species",
        action="store_true",
        help="Compatibility flag supplied by the mixed-species workflow.",
    )
    return parser.parse_args()


def count_bam_reads(bam_path):
    bam = pysam.AlignmentFile(bam_path, "rb")
    stats = bam.get_index_statistics()
    n_reads = int(sum(contig.mapped for contig in stats))
    bam.close()
    return n_reads


def count_bam_aligned_unique_reads(bam_path):
    bam = pysam.AlignmentFile(bam_path, "rb")
    read_ids = set()
    for align in bam.fetch(until_eof=True):
        if align.is_unmapped:
            continue
        read_ids.add(align.query_name)
    bam.close()
    return len(read_ids)


def summarize_bam_alignments(bam_path):
    bam = pysam.AlignmentFile(bam_path, "rb")
    summary = {
        "bam_records": 0,
        "mapped_primary_alignments": 0,
        "unmapped_reads": 0,
        "supplementary_alignments": 0,
        "mapped_unique_reads": set(),
    }
    for align in bam.fetch(until_eof=True):
        summary["bam_records"] += 1
        if align.is_unmapped:
            summary["unmapped_reads"] += 1
            continue
        summary["mapped_unique_reads"].add(align.query_name)
        if align.is_supplementary:
            summary["supplementary_alignments"] += 1
        elif not align.is_secondary:
            summary["mapped_primary_alignments"] += 1
    bam.close()
    return summary


def count_fastq_reads(fastq_path):
    opener = gzip.open if str(fastq_path).endswith(".gz") else open
    n_lines = 0
    with opener(fastq_path, "rt") as handle:
        for _ in handle:
            n_lines += 1
    if n_lines % 4 != 0:
        raise ValueError(f"FASTQ line count is not divisible by 4: {fastq_path}")
    return n_lines // 4


def is_mito_gene(gene_name):
    gene_name = str(gene_name)
    if gene_name.startswith("MT-") or gene_name.startswith("mt-") or gene_name.startswith("Mt-"):
        return True
    if "_" in gene_name:
        suffix = gene_name.split("_", 1)[1]
        return suffix.startswith("MT-") or suffix.startswith("mt-") or suffix.startswith("Mt-")
    return False


def is_unannotated_region_label(gene_name):
    gene_name = str(gene_name)
    return re.search(r"[a-zA-Z0-9]+_\d+_\d+", gene_name) is not None


def is_known_gene(gene_name):
    gene_name = str(gene_name).strip()
    if gene_name in ("", "NA", "nan", "None"):
        return False
    return not is_unannotated_region_label(gene_name)


def load_transcript_assignments(path):
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["read_id", "status", "score", "gene", "transcript_id"],
        dtype=str,
    )
    df["read_id"] = df["read_id"].astype(str).str.strip()
    return df


def count_read_tags(path):
    tag_df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if "read_id" not in tag_df.columns:
        raise ValueError(f"Missing read_id column in {path}")
    tag_df["read_id"] = tag_df["read_id"].astype(str).str.strip()
    tag_df = tag_df[tag_df["read_id"] != ""]
    return tag_df["read_id"].nunique()


def load_barcode_validity_summary(path):
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return {}
    return df.iloc[0].to_dict()


def load_glycine_total_full_length_reads(path):
    if path is None:
        return None
    if not os.path.exists(path):
        return None
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        lines = [line.rstrip("\n") for line in handle]

    in_top_level_type_block = False
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped == "Type\tRead_count\tRead_proportion(%)":
            in_top_level_type_block = True
            continue
        if stripped == "Non-chimeric":
            break
        if in_top_level_type_block:
            parts = stripped.split("\t")
            if len(parts) >= 3 and parts[0] == "Full-length":
                try:
                    return int(parts[1])
                except ValueError:
                    return None
    return None


def load_glycine_summary_counts(path):
    if path is None or not os.path.exists(path):
        return {}
    counts = {}
    in_summary = False
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped == "Type\tRead_count\tRead_proportion(%)":
                in_summary = True
                continue
            if stripped == "Non-chimeric":
                break
            if in_summary and stripped:
                parts = stripped.split("\t")
                if len(parts) >= 2:
                    counts[parts[0]] = int(parts[1])
    return counts


def main():
    args = parse_args()

    df = pd.read_csv(args.cell_umi_gene_tsv, sep="\t")

    required_cols = {"read_id", "gene", "barcode", "umi"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in {args.cell_umi_gene_tsv}: {missing}")

    estimated_cells = df["barcode"].nunique()
    total_bam_reads = count_bam_reads(args.bam)
    bam_alignment_summary = summarize_bam_alignments(args.bam)
    mapped_read_ids = bam_alignment_summary["mapped_unique_reads"]
    aligned_bam_unique_reads = len(mapped_read_ids)
    glycine_counts = load_glycine_summary_counts(args.glycine_stats)
    raw_fastq_reads = (
        glycine_counts["Total"]
        if "Total" in glycine_counts
        else count_fastq_reads(args.raw_fastq)
    )
    full_length_reads = count_fastq_reads(args.full_length_fastq) if args.full_length_fastq else raw_fastq_reads
    glycine_total_full_length_reads = load_glycine_total_full_length_reads(args.glycine_log)
    if glycine_total_full_length_reads is not None:
        full_length_reads = glycine_total_full_length_reads
    elif "Full-length+rescued" in glycine_counts:
        full_length_reads = glycine_counts["Full-length+rescued"]
    final_cell_read_ids = set(df["read_id"].astype(str).str.strip())
    cell_associated_reads = len(final_cell_read_ids)
    aligned_genome_reads = len(final_cell_read_ids & mapped_read_ids)
    assigned_cell_reads = count_read_tags(args.read_tags) if args.read_tags else cell_associated_reads
    barcode_validity = load_barcode_validity_summary(args.barcode_validity_tsv) if args.barcode_validity_tsv else {}
    barcode_valid_reads = int(barcode_validity.get("valid_any_n", assigned_cell_reads))

    df_gene = df[df["gene"].apply(is_known_gene)].copy()
    known_gene_reads = df_gene["read_id"].nunique()
    unique_genes = df_gene["gene"].nunique()

    known_transcript_reads = None
    unique_isoforms = None
    if args.transcript_assigns:
        tx_df = load_transcript_assignments(args.transcript_assigns)
        tx_df = tx_df[tx_df["status"] == "Assigned"].copy()
        tx_df = tx_df[tx_df["transcript_id"].astype(str).str.strip() != "NA"].copy()
        tx_df = tx_df[tx_df["read_id"].isin(df["read_id"])].copy()
        known_transcript_reads = tx_df["read_id"].nunique()
        unique_isoforms = tx_df["transcript_id"].nunique()

    reads_per_cell = df.groupby("barcode")["read_id"].nunique()
    umis_per_cell = df_gene.groupby("barcode")["umi"].nunique().reindex(reads_per_cell.index, fill_value=0)
    genes_per_cell = df_gene.groupby("barcode")["gene"].nunique().reindex(reads_per_cell.index, fill_value=0)

    df["is_mito"] = df["gene"].apply(is_mito_gene)
    mito_reads_per_cell = df.groupby("barcode")["is_mito"].sum()
    mito_percent_per_cell = (mito_reads_per_cell / reads_per_cell) * 100

    metrics = {
        "Input reads": raw_fastq_reads,
        "Full length reads": full_length_reads,
        "Estimated number of cells": estimated_cells,
        "Raw FASTQ reads": raw_fastq_reads,
        "Aligned BAM reads": aligned_bam_unique_reads,
        "Aligned BAM records": bam_alignment_summary["bam_records"],
        "Pass reads": full_length_reads,
        "Mapped": bam_alignment_summary["mapped_primary_alignments"],
        "Unmapped": bam_alignment_summary["unmapped_reads"],
        "Supplementary": bam_alignment_summary["supplementary_alignments"],
        "Barcode-valid reads": barcode_valid_reads,
        "Reads assigned to final cells": assigned_cell_reads,
        "Gene assigned reads": known_gene_reads,
        "Reads in final cells": cell_associated_reads,
        "Reads aligned to reference genome in final cells": aligned_genome_reads,
        "Reads per cell (mean)": cell_associated_reads / estimated_cells if estimated_cells > 0 else 0,
        "Mean reads per cell": cell_associated_reads / estimated_cells if estimated_cells > 0 else 0,
        "Mean cell-associated reads per cell": reads_per_cell.mean(),
        "Mean UMI counts per cell": umis_per_cell.mean(),
        "Median UMI counts per cell": umis_per_cell.median(),
        "Mean Genes per cell": genes_per_cell.mean(),
        "Median Genes per cell": genes_per_cell.median(),
        "Genes per cell (median)": genes_per_cell.median(),
        "Unique genes": unique_genes,
        "Total genes detected": unique_genes,
        "Fraction reads in cells": cell_associated_reads / full_length_reads if full_length_reads > 0 else 0,
        "Percent full length reads": full_length_reads / raw_fastq_reads if raw_fastq_reads > 0 else 0,
        "Percent barcode-valid reads of full length": (
            barcode_valid_reads / full_length_reads if full_length_reads > 0 else 0
        ),
        "Percent reads assigned to final cells of full length": (
            assigned_cell_reads / full_length_reads if full_length_reads > 0 else 0
        ),
        "Percent gene assigned reads of full length": (
            known_gene_reads / full_length_reads if full_length_reads > 0 else 0
        ),
        "Fraction reads aligned to reference genome in final cells": (
            aligned_genome_reads / full_length_reads if full_length_reads > 0 else 0
        ),
    }

    if known_transcript_reads is not None:
        metrics["Transcript assigned reads"] = known_transcript_reads
        metrics["Percent transcript assigned reads of full length"] = (
            known_transcript_reads / full_length_reads if full_length_reads > 0 else 0
        )
        metrics["High-confidence known-transcript reads"] = known_transcript_reads
        metrics["Fraction high-confidence known-transcript reads in final cells"] = (
            known_transcript_reads / full_length_reads if full_length_reads > 0 else 0
        )
    if unique_isoforms is not None:
        metrics["Unique isoforms"] = unique_isoforms

    metrics_df = pd.DataFrame({"Metric": list(metrics.keys()), "Value": list(metrics.values())})

    fraction_metrics = {
        "Fraction reads in cells",
        "Percent full length reads",
        "Percent barcode-valid reads of full length",
        "Percent reads assigned to final cells of full length",
        "Percent gene assigned reads of full length",
        "Percent transcript assigned reads of full length",
        "Fraction reads aligned to reference genome in final cells",
        "Fraction high-confidence known-transcript reads in final cells",
    }
    formatted_metrics_df = metrics_df.copy()
    formatted_metrics_df["Value"] = formatted_metrics_df.apply(
        lambda row: f"{row['Value']:.2%}"
        if row["Metric"] in fraction_metrics
        else f"{row['Value']:,.2f}"
        if isinstance(row["Value"], float)
        else f"{row['Value']:,}",
        axis=1,
    )

    formatted_metrics_df.to_csv("rna_qc_metrics.tsv", sep="\t", index=False)

    def metric_value(name):
        return metrics.get(name, "")

    report_records = [
        ("Experiment summary", "Input reads", metric_value("Input reads")),
        ("Experiment summary", "Estimated cells", metric_value("Estimated number of cells")),
        ("Experiment summary", "Reads per cell (mean)", metric_value("Reads per cell (mean)")),
        ("Experiment summary", "UMIs per cell (median)", metric_value("Median UMI counts per cell")),
        ("Experiment summary", "Genes per cell (median)", metric_value("Genes per cell (median)")),
        ("Alignment / feature summary", "Pass reads", metric_value("Pass reads")),
        ("Alignment / feature summary", "Mapped", metric_value("Mapped")),
        ("Alignment / feature summary", "Unmapped", metric_value("Unmapped")),
        ("Alignment / feature summary", "Supplementary", metric_value("Supplementary")),
        ("Alignment / feature summary", "Unique genes", metric_value("Unique genes")),
        ("Alignment / feature summary", "Unique isoforms", metric_value("Unique isoforms")),
        ("Read assignment summary", "Full length", metric_value("Full length reads")),
        ("Read assignment summary", "Barcode-valid", metric_value("Barcode-valid reads")),
        ("Read assignment summary", "Cell-assigned", metric_value("Reads assigned to final cells")),
        ("Read assignment summary", "Gene assigned", metric_value("Gene assigned reads")),
        ("Read assignment summary", "Transcript assigned", metric_value("Transcript assigned reads")),
        ("Read assignment percentage", "% full length reads", metric_value("Percent full length reads")),
        (
            "Read assignment percentage",
            "% barcode-valid reads",
            metric_value("Percent barcode-valid reads of full length"),
        ),
        (
            "Read assignment percentage",
            "% cell-assigned reads",
            metric_value("Percent reads assigned to final cells of full length"),
        ),
        (
            "Read assignment percentage",
            "% gene assigned reads",
            metric_value("Percent gene assigned reads of full length"),
        ),
        (
            "Read assignment percentage",
            "% transcript assigned reads",
            metric_value("Percent transcript assigned reads of full length"),
        ),
    ]
    report_df = pd.DataFrame(report_records, columns=["Section", "Metric", "Value"])
    report_fraction_metrics = {
        "% full length reads",
        "% barcode-valid reads",
        "% cell-assigned reads",
        "% gene assigned reads",
        "% transcript assigned reads",
    }
    report_df["Formatted_value"] = report_df.apply(
        lambda row: ""
        if row["Value"] == ""
        else f"{row['Value']:.2%}"
        if row["Metric"] in report_fraction_metrics
        else f"{row['Value']:,.2f}"
        if isinstance(row["Value"], float)
        else f"{row['Value']:,}",
        axis=1,
    )
    report_df.to_csv("single_cell_report_metrics.tsv", sep="\t", index=False)

    print("\nRNA QC Metrics")
    print(formatted_metrics_df.to_string(index=False))

    per_cell_qc = pd.DataFrame(
        {
            "barcode": reads_per_cell.index,
            "reads": reads_per_cell.values,
            "umis": umis_per_cell.reindex(reads_per_cell.index).values,
            "genes": genes_per_cell.reindex(reads_per_cell.index).values,
            "mito_percent": mito_percent_per_cell.reindex(reads_per_cell.index).values,
        }
    )
    per_cell_qc.to_csv("per_cell_qc.tsv", sep="\t", index=False)

    plot_df = per_cell_qc.melt(
        id_vars="barcode",
        value_vars=["genes", "umis", "mito_percent"],
        var_name="metric",
        value_name="value",
    )

    label_map = {
        "genes": "Genes",
        "umis": "UMIs",
        "mito_percent": "Mito percent",
    }
    plot_df["metric"] = plot_df["metric"].map(label_map)

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
    plt.savefig("rna_violin_plot.png", dpi=300, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
