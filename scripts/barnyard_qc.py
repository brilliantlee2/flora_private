import argparse
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    plt = None


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Generate mixed-species barnyard QC summaries from Strint2 cell_umi_gene.tsv."
        )
    )
    parser.add_argument("--input", required=True, help="Input cell_umi_gene.tsv.")
    parser.add_argument("--out-dir", default="barnyard_qc")
    parser.add_argument("--human-prefix", default="hs_")
    parser.add_argument("--mouse-prefix", default="mm_")
    parser.add_argument("--genome0-name", default="human")
    parser.add_argument("--genome1-name", default="mouse")
    parser.add_argument("--singlet-threshold", type=float, default=0.9)
    parser.add_argument(
        "--ambient-umi-threshold",
        type=int,
        default=1000,
        help=(
            "Remove cells with both genome0/human and genome1/mouse UMI counts below "
            "this threshold before plotting and purity calculation [1000]."
        ),
    )
    parser.add_argument(
        "--disable-ambient-filter",
        action="store_true",
        help="Do not remove low-low barnyard points.",
    )
    parser.add_argument(
        "--core-cells-debug",
        default=None,
        help="Optional core_cells_debug.csv to annotate core_type and n_barcodes.",
    )
    parser.add_argument(
        "--filter-n-barcodes",
        default=None,
        help="Optional integer or comma-separated barcode-count filter, e.g. 1 or 1,2.",
    )
    return parser.parse_args()


def classify_gene_species(gene_series, human_prefix, mouse_prefix):
    genes = gene_series.astype(str).fillna("")
    species = np.where(
        genes.str.startswith(human_prefix),
        "human",
        np.where(genes.str.startswith(mouse_prefix), "mouse", "other"),
    )
    return pd.Series(species, index=gene_series.index, name="species")


def safe_fraction(num, den):
    num = np.asarray(num, dtype=float)
    den = np.asarray(den, dtype=float)
    out = np.zeros_like(num, dtype=float)
    mask = den > 0
    out[mask] = num[mask] / den[mask]
    return out


def parse_int_set(raw_value):
    if raw_value is None:
        return None
    values = []
    for part in str(raw_value).split(","):
        part = part.strip()
        if not part:
            continue
        try:
            value = int(part)
        except ValueError as exc:
            raise ValueError(
                f"Invalid --filter-n-barcodes value {part!r}; expected integer(s)."
            ) from exc
        if value < 1:
            raise ValueError("--filter-n-barcodes values must be >= 1")
        values.append(value)
    return sorted(set(values)) if values else None


def assignment_from_umi_fractions(human_frac, mouse_frac, singlet_threshold):
    assignments = np.full(len(human_frac), "unclassified", dtype=object)
    assignments[human_frac >= singlet_threshold] = "human_singlet"
    assignments[mouse_frac >= singlet_threshold] = "mouse_singlet"
    mixed_mask = (human_frac > 0) & (mouse_frac > 0)
    assignments[
        mixed_mask
        & (human_frac < singlet_threshold)
        & (mouse_frac < singlet_threshold)
    ] = "mixed"
    return assignments


def refresh_fraction_columns(per_cell, singlet_threshold):
    per_cell = per_cell.copy()
    per_cell["read_total"] = per_cell[["read_human", "read_mouse", "read_other"]].sum(axis=1)
    per_cell["umi_total"] = per_cell[["umi_human", "umi_mouse", "umi_other"]].sum(axis=1)
    per_cell["gene_total"] = per_cell[["gene_human", "gene_mouse", "gene_other"]].sum(axis=1)
    per_cell["human_read_fraction"] = safe_fraction(per_cell["read_human"], per_cell["read_total"])
    per_cell["mouse_read_fraction"] = safe_fraction(per_cell["read_mouse"], per_cell["read_total"])
    per_cell["human_umi_fraction"] = safe_fraction(per_cell["umi_human"], per_cell["umi_total"])
    per_cell["mouse_umi_fraction"] = safe_fraction(per_cell["umi_mouse"], per_cell["umi_total"])
    per_cell["minority_umi_fraction"] = np.minimum(
        per_cell["human_umi_fraction"], per_cell["mouse_umi_fraction"]
    )
    per_cell["assignment"] = assignment_from_umi_fractions(
        per_cell["human_umi_fraction"].to_numpy(),
        per_cell["mouse_umi_fraction"].to_numpy(),
        singlet_threshold,
    )
    return per_cell


def build_per_cell_table(df, human_prefix, mouse_prefix, singlet_threshold):
    df = df.copy()
    required = {"read_id", "gene", "barcode", "umi"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in input TSV: {sorted(missing)}")

    df["read_id"] = df["read_id"].astype(str).str.strip()
    df["gene"] = df["gene"].astype(str).str.strip()
    df["barcode"] = df["barcode"].astype(str).str.strip()
    df["umi"] = df["umi"].astype(str).str.strip()
    df = df[(df["barcode"] != "") & (df["gene"] != "") & (df["umi"] != "")]
    df["species"] = classify_gene_species(df["gene"], human_prefix, mouse_prefix)

    read_counts = (
        df.groupby(["barcode", "species"])["read_id"]
        .nunique()
        .unstack(fill_value=0)
        .rename(columns=lambda c: f"read_{c}")
    )

    umi_df = df.drop_duplicates(subset=["barcode", "gene", "umi"]).copy()
    umi_counts = (
        umi_df.groupby(["barcode", "species"])["umi"]
        .size()
        .unstack(fill_value=0)
        .rename(columns=lambda c: f"umi_{c}")
    )

    gene_counts = (
        umi_df.groupby(["barcode", "species"])["gene"]
        .nunique()
        .unstack(fill_value=0)
        .rename(columns=lambda c: f"gene_{c}")
    )

    per_cell = read_counts.join(umi_counts, how="outer").join(gene_counts, how="outer")
    per_cell = per_cell.fillna(0).astype(int).reset_index().rename(columns={"barcode": "cell_id"})

    for prefix in ("read", "umi", "gene"):
        for suffix in ("human", "mouse", "other"):
            col = f"{prefix}_{suffix}"
            if col not in per_cell.columns:
                per_cell[col] = 0

    per_cell = refresh_fraction_columns(per_cell, singlet_threshold)
    sort_cols = [
        "cell_id",
        "assignment",
        "read_total",
        "read_human",
        "read_mouse",
        "read_other",
        "umi_total",
        "umi_human",
        "umi_mouse",
        "umi_other",
        "gene_total",
        "gene_human",
        "gene_mouse",
        "gene_other",
        "human_read_fraction",
        "mouse_read_fraction",
        "human_umi_fraction",
        "mouse_umi_fraction",
        "minority_umi_fraction",
    ]
    return per_cell[sort_cols].sort_values(["assignment", "umi_total"], ascending=[True, False])


def apply_low_low_ambient_filter(per_cell, threshold=1000):
    threshold = int(threshold)
    if threshold <= 0:
        stats = {
            "ambient_filter_enabled": False,
            "ambient_umi_threshold_each_species": threshold,
            "cells_before_ambient_filter": int(len(per_cell)),
            "cells_after_ambient_filter": int(len(per_cell)),
            "cells_removed_by_ambient_filter": 0,
        }
        return per_cell.copy(), per_cell.iloc[0:0].copy(), stats

    remove_mask = (per_cell["umi_human"] < threshold) & (per_cell["umi_mouse"] < threshold)
    removed = per_cell[remove_mask].copy()
    filtered = per_cell[~remove_mask].copy()
    stats = {
        "ambient_filter_enabled": True,
        "ambient_umi_threshold_each_species": threshold,
        "cells_before_ambient_filter": int(len(per_cell)),
        "cells_after_ambient_filter": int(len(filtered)),
        "cells_removed_by_ambient_filter": int(len(removed)),
    }
    return filtered, removed, stats


def apply_min_total_umi_filter(per_cell, min_total_umi):
    """Compatibility helper for older local tests."""
    min_total_umi = int(min_total_umi)
    filtered = per_cell[per_cell["umi_total"] >= min_total_umi].copy()
    stats = {
        "min_total_umi_filter": min_total_umi,
        "cells_before_min_total_umi_filter": int(len(per_cell)),
        "cells_after_min_total_umi_filter": int(len(filtered)),
        "cells_removed_by_min_total_umi_filter": int(len(per_cell) - len(filtered)),
    }
    return filtered, stats


def quantile_or_nan(series, q):
    if len(series) == 0:
        return np.nan
    return float(series.quantile(q))


def build_summary_table(per_cell):
    total_cells = int(len(per_cell))
    assignment_counts = per_cell["assignment"].value_counts().to_dict()
    human_singlets = per_cell[per_cell["assignment"] == "human_singlet"]
    mouse_singlets = per_cell[per_cell["assignment"] == "mouse_singlet"]
    mixed = per_cell[per_cell["assignment"] == "mixed"]
    mixed_human_fraction = mixed["human_umi_fraction"]
    mixed_near_50_50 = mixed_human_fraction.between(0.4, 0.6, inclusive="both")
    mixed_mostly_90_10 = mixed_human_fraction.between(0.05, 0.10, inclusive="both") | mixed_human_fraction.between(0.90, 0.95, inclusive="both")

    classified = (
        int(assignment_counts.get("human_singlet", 0))
        + int(assignment_counts.get("mouse_singlet", 0))
        + int(assignment_counts.get("mixed", 0))
    )
    cross_species_doublet_rate = (int(assignment_counts.get("mixed", 0)) / classified) if classified else np.nan

    summary_rows = [
        ("total_cells", total_cells),
        ("human_singlet_cells", int(assignment_counts.get("human_singlet", 0))),
        ("mouse_singlet_cells", int(assignment_counts.get("mouse_singlet", 0))),
        ("mixed_cells", int(assignment_counts.get("mixed", 0))),
        ("unclassified_cells", int(assignment_counts.get("unclassified", 0))),
        ("cross_species_doublet_rate_among_classified", cross_species_doublet_rate),
        ("human_singlet_mouse_fraction_median", quantile_or_nan(human_singlets["mouse_umi_fraction"], 0.5)),
        ("human_singlet_mouse_fraction_p95", quantile_or_nan(human_singlets["mouse_umi_fraction"], 0.95)),
        ("mouse_singlet_human_fraction_median", quantile_or_nan(mouse_singlets["human_umi_fraction"], 0.5)),
        ("mouse_singlet_human_fraction_p95", quantile_or_nan(mouse_singlets["human_umi_fraction"], 0.95)),
        ("mixed_human_fraction_median", quantile_or_nan(mixed["human_umi_fraction"], 0.5)),
        ("mixed_mouse_fraction_median", quantile_or_nan(mixed["mouse_umi_fraction"], 0.5)),
        ("mixed_near_50_50_cells", int(mixed_near_50_50.sum())),
        ("mixed_near_50_50_fraction", float(mixed_near_50_50.mean()) if len(mixed_near_50_50) else np.nan),
        ("mixed_mostly_90_10_cells", int(mixed_mostly_90_10.sum())),
        ("mixed_mostly_90_10_fraction", float(mixed_mostly_90_10.mean()) if len(mixed_mostly_90_10) else np.nan),
        ("human_singlet_umi_total_median", quantile_or_nan(human_singlets["umi_total"], 0.5)),
        ("mouse_singlet_umi_total_median", quantile_or_nan(mouse_singlets["umi_total"], 0.5)),
        ("mixed_umi_total_median", quantile_or_nan(mixed["umi_total"], 0.5)),
    ]
    return pd.DataFrame(summary_rows, columns=["metric", "value"])


def build_overall_purity_table(per_cell, genome0_name="human", genome1_name="mouse"):
    human_singlets = per_cell[per_cell["assignment"] == "human_singlet"].copy()
    mouse_singlets = per_cell[per_cell["assignment"] == "mouse_singlet"].copy()
    singlets = pd.concat([human_singlets, mouse_singlets], ignore_index=True)

    human_total = float(human_singlets["umi_total"].sum())
    mouse_total = float(mouse_singlets["umi_total"].sum())
    human_own = float(human_singlets["umi_human"].sum())
    mouse_own = float(mouse_singlets["umi_mouse"].sum())
    singlet_total = float(singlets["umi_total"].sum())
    dominant_total = human_own + mouse_own

    rows = [
        ("genome0_name", genome0_name),
        ("genome1_name", genome1_name),
        ("genome0_singlet_cells", int(len(human_singlets))),
        ("genome1_singlet_cells", int(len(mouse_singlets))),
        ("effective_singlet_cells_for_overall_purity", int(len(singlets))),
        ("genome0_singlet_own_umis", int(human_own)),
        ("genome0_singlet_total_umis", int(human_total)),
        ("genome1_singlet_own_umis", int(mouse_own)),
        ("genome1_singlet_total_umis", int(mouse_total)),
        ("all_singlet_dominant_species_umis", int(dominant_total)),
        ("all_singlet_total_umis", int(singlet_total)),
        ("human_mean_purity", human_own / human_total if human_total else np.nan),
        ("mouse_mean_purity", mouse_own / mouse_total if mouse_total else np.nan),
        ("genome0_mean_purity", human_own / human_total if human_total else np.nan),
        ("genome1_mean_purity", mouse_own / mouse_total if mouse_total else np.nan),
        ("overall_purity", dominant_total / singlet_total if singlet_total else np.nan),
    ]
    return pd.DataFrame(rows, columns=["metric", "value"])


def parse_barcodes_field(value):
    if pd.isna(value):
        return []
    raw = str(value).strip()
    if raw == "":
        return []
    return [part.strip() for part in raw.split("|") if part.strip()]


def attach_core_metadata(per_cell, core_cells_debug_path):
    core = pd.read_csv(core_cells_debug_path, dtype=str).fillna("")
    if "cell_id" not in core.columns:
        raise ValueError("core_cells_debug.csv must contain cell_id")

    keep_cols = ["cell_id"]
    if "type" in core.columns:
        keep_cols.append("type")
    if "n_barcodes" in core.columns:
        keep_cols.append("n_barcodes")
    elif "barcodes" in core.columns:
        core["n_barcodes"] = core["barcodes"].apply(lambda x: len(parse_barcodes_field(x)))
        keep_cols.append("n_barcodes")

    core_small = core[keep_cols].drop_duplicates(subset=["cell_id"]).copy()
    if "type" in core_small.columns:
        core_small.rename(columns={"type": "core_type"}, inplace=True)
    if "n_barcodes" in core_small.columns:
        core_small["n_barcodes"] = (
            pd.to_numeric(core_small["n_barcodes"], errors="coerce").fillna(0).astype(int)
        )

    merged = per_cell.merge(core_small, on="cell_id", how="left", validate="one_to_one")
    ordered_cols = ["cell_id"]
    if "core_type" in merged.columns:
        ordered_cols.append("core_type")
    if "n_barcodes" in merged.columns:
        ordered_cols.append("n_barcodes")
    ordered_cols.extend([col for col in merged.columns if col not in ordered_cols])
    return merged[ordered_cols]


def filter_per_cell_by_n_barcodes(per_cell, n_barcodes_values):
    if not n_barcodes_values:
        return per_cell
    if "n_barcodes" not in per_cell.columns:
        raise ValueError(
            "--filter-n-barcodes requires --core-cells-debug with n_barcodes information"
        )
    return per_cell[per_cell["n_barcodes"].isin(n_barcodes_values)].copy()


def write_metric_table(path, values):
    pd.DataFrame(list(values.items()), columns=["metric", "value"]).to_csv(
        path, sep="\t", index=False
    )


def plot_barnyard(per_cell, out_base):
    if plt is None:
        print("[barnyard_qc] matplotlib is not available; skipping barnyard scatter plots.")
        return
    colors = {
        "human_singlet": "#2E7D32",
        "mouse_singlet": "#1565C0",
        "mixed": "#C62828",
        "unclassified": "#757575",
    }
    order = ["human_singlet", "mouse_singlet", "mixed", "unclassified"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 6), constrained_layout=True)

    for ax, x_col, y_col, title in [
        (axes[0], "umi_human", "umi_mouse", "Barnyard by UMI"),
        (axes[1], "read_human", "read_mouse", "Barnyard by Read"),
    ]:
        for label in order:
            sub = per_cell[per_cell["assignment"] == label]
            if sub.empty:
                continue
            ax.scatter(
                sub[x_col],
                sub[y_col],
                s=10,
                alpha=0.55,
                c=colors[label],
                label=f"{label} (n={len(sub)})",
                edgecolors="none",
            )
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.set_xlabel(f"Human {x_col.split('_')[0].upper()}")
        ax.set_ylabel(f"Mouse {y_col.split('_')[0].upper()}")
        ax.set_title(title)
        ax.grid(True, alpha=0.2)

    axes[1].legend(frameon=False, loc="lower right", fontsize=8)
    fig.savefig(f"{out_base}.png", dpi=220)
    fig.savefig(f"{out_base}.pdf")
    plt.close(fig)


def plot_minority_histograms(per_cell, out_base):
    if plt is None:
        print("[barnyard_qc] matplotlib is not available; skipping minority histograms.")
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    human_singlets = per_cell[per_cell["assignment"] == "human_singlet"]
    mouse_singlets = per_cell[per_cell["assignment"] == "mouse_singlet"]
    plots = [
        (axes[0], human_singlets["mouse_umi_fraction"], "Human singlets: mouse minority fraction", "#2E7D32"),
        (axes[1], mouse_singlets["human_umi_fraction"], "Mouse singlets: human minority fraction", "#1565C0"),
    ]
    bins = np.linspace(0, 0.2, 81)
    for ax, series, title, color in plots:
        if len(series) > 0:
            clipped = np.clip(series.to_numpy(dtype=float), 0, bins[-1])
            ax.hist(clipped, bins=bins, color=color, alpha=0.85)
            med = float(np.median(series))
            p95 = float(np.quantile(series, 0.95))
            ax.axvline(med, color="black", linestyle="--", linewidth=1, label=f"median={med:.4f}")
            ax.axvline(p95, color="black", linestyle=":", linewidth=1, label=f"p95={p95:.4f}")
            ax.legend(frameon=False, fontsize=8)
        ax.set_title(title)
        ax.set_xlabel("Minority-species UMI fraction")
        ax.set_ylabel("Cells")
        ax.grid(True, alpha=0.2)
    fig.savefig(f"{out_base}.png", dpi=220)
    fig.savefig(f"{out_base}.pdf")
    plt.close(fig)


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    filter_n_barcodes = parse_int_set(args.filter_n_barcodes)

    df = pd.read_csv(args.input, sep="\t", dtype=str)
    per_cell_raw = build_per_cell_table(
        df,
        human_prefix=args.human_prefix,
        mouse_prefix=args.mouse_prefix,
        singlet_threshold=args.singlet_threshold,
    )
    if args.core_cells_debug:
        per_cell_raw = attach_core_metadata(per_cell_raw, args.core_cells_debug)
    if filter_n_barcodes:
        per_cell_raw = filter_per_cell_by_n_barcodes(per_cell_raw, filter_n_barcodes)

    if args.disable_ambient_filter:
        per_cell, removed, filter_stats = apply_low_low_ambient_filter(per_cell_raw, threshold=0)
    else:
        per_cell, removed, filter_stats = apply_low_low_ambient_filter(
            per_cell_raw, threshold=args.ambient_umi_threshold
        )

    summary = build_summary_table(per_cell)
    purity = build_overall_purity_table(
        per_cell,
        genome0_name=args.genome0_name,
        genome1_name=args.genome1_name,
    )
    if filter_n_barcodes:
        filter_stats["filter_n_barcodes"] = ",".join(str(x) for x in filter_n_barcodes)

    per_cell_raw.to_csv(out_dir / "barnyard_per_cell_raw.tsv", sep="\t", index=False)
    removed.to_csv(out_dir / "barnyard_ambient_removed_cells.tsv", sep="\t", index=False)
    per_cell.to_csv(out_dir / "barnyard_per_cell.tsv", sep="\t", index=False)
    summary.to_csv(out_dir / "barnyard_summary.tsv", sep="\t", index=False)
    purity.to_csv(out_dir / "barnyard_overall_purity.tsv", sep="\t", index=False)
    write_metric_table(out_dir / "barnyard_filter_summary.tsv", filter_stats)

    plot_barnyard(per_cell, str(out_dir / "barnyard_scatter"))
    plot_minority_histograms(per_cell, str(out_dir / "barnyard_minority_fraction_hist"))

    print(out_dir / "barnyard_summary.tsv")
    print(out_dir / "barnyard_overall_purity.tsv")
    print(out_dir / "barnyard_per_cell.tsv")
    print(out_dir / "barnyard_scatter.png")


if __name__ == "__main__":
    main()
