from collections import Counter
from pathlib import Path
import time

import pandas as pd

from args_parser import set_parser
from utils import (
    assign_read,
    assign_reads_to_cells,
    build_core_cells,
    compute_top1_dominance,
    extract_reads_and_filter_df_by_raw,
    filter_pairs_three_stage,
    get_3p_features,
    get_5p_features,
    get_bc_whitelist,
    green_msg,
    is_missing,
    norm_bc,
    read_batch_generator,
    revcomp,
    reverse_complement,
    strip_fixed_3p,
    strip_fixed_5p,
)


def write_barcode_lines(path, values):
    with open(path, "w", encoding="utf-8") as handle:
        for value in values:
            handle.write(f"{value}\n")


def write_barcode_count_table(path, counts):
    rows = sorted(counts.items(), key=lambda item: item[1], reverse=True)
    pd.DataFrame(rows, columns=["barcode", "count"]).to_csv(path, sep="\t", index=False)


def count_high_quality_barcodes(df, bc_col, q_col, min_q):
    mask = df[q_col].fillna(-1) >= min_q
    series = df.loc[mask, bc_col].fillna("").astype(str).str.strip()
    series = series[series != ""]
    return Counter(series.value_counts().to_dict())


def is_valid_bc(value):
    if pd.isna(value):
        return False
    value = str(value).strip()
    return value != "" and value.upper() != "NA"


def summarize_barcode_validity(big_df):
    total_n = int(len(big_df))
    if total_n == 0:
        return {
            "total_rows": 0,
            "valid_any_n": 0,
            "valid_any_ratio": 0.0,
            "only_3p_n": 0,
            "only_3p_ratio": 0.0,
            "only_5p_n": 0,
            "only_5p_ratio": 0.0,
            "both_n": 0,
            "both_ratio": 0.0,
            "neither_n": 0,
            "neither_ratio": 0.0,
        }

    has_3p = big_df["BC3_corrected"].apply(is_valid_bc)
    has_5p = big_df["BC5_corrected"].apply(is_valid_bc)

    only_3p_n = int((has_3p & (~has_5p)).sum())
    only_5p_n = int(((~has_3p) & has_5p).sum())
    both_n = int((has_3p & has_5p).sum())
    neither_n = int(((~has_3p) & (~has_5p)).sum())
    valid_any_n = int(total_n - neither_n)

    return {
        "total_rows": total_n,
        "valid_any_n": valid_any_n,
        "valid_any_ratio": valid_any_n / total_n,
        "only_3p_n": only_3p_n,
        "only_3p_ratio": only_3p_n / total_n,
        "only_5p_n": only_5p_n,
        "only_5p_ratio": only_5p_n / total_n,
        "both_n": both_n,
        "both_ratio": both_n / total_n,
        "neither_n": neither_n,
        "neither_ratio": neither_n / total_n,
    }


def summarize_trimmed_barcode_uniques(df):
    bc3 = df["BC3_20bp_rc"].map(norm_bc) if "BC3_20bp_rc" in df.columns else pd.Series(dtype=str)
    bc5 = df["BC5_20bp"].map(norm_bc) if "BC5_20bp" in df.columns else pd.Series(dtype=str)

    bc3_nonempty = bc3[bc3 != ""]
    bc5_nonempty = bc5[bc5 != ""]
    union_n = pd.concat([bc3_nonempty, bc5_nonempty], ignore_index=True).nunique()

    return {
        "unique_bc3_20bp_rc": int(bc3_nonempty.nunique()),
        "unique_bc5_20bp": int(bc5_nonempty.nunique()),
        "unique_union_bc3_bc5": int(union_n),
    }


def build_putative_tables(args):
    read_ids = []
    putative_bcs = []
    putative_bc_min_qs = []
    bc_fixed_locs = []
    umis = []
    umi_fixed_locs = []
    post_umi_flankings = []
    polyA_starts = []
    read_types = []

    read_ids_5p = []
    putative_bcs_5p = []
    putative_bc_min_qs_5p = []
    bc_fixed_locs_5p = []
    umis_5p = []
    umi_fixed_locs_5p = []

    bc_fixed_3p = reverse_complement(args.BC_fixed_3p)
    umi_fixed_3p = reverse_complement(args.umi_fixed_3p)
    bc_fixed_5p = reverse_complement(args.BC_fixed_5p)
    umi_fixed_5p = reverse_complement(args.umi_fixed_5p)

    for batch in read_batch_generator(args.fastq_fns, args.batch_size):
        for read_info in batch:
            get_3p_features(
                read_info=read_info,
                read_ids=read_ids,
                putative_bcs=putative_bcs,
                bc_fixed_locs=bc_fixed_locs,
                putative_bc_min_qs=putative_bc_min_qs,
                umis=umis,
                umi_fixed_locs=umi_fixed_locs,
                post_umi_flankings=post_umi_flankings,
                polyA_starts=polyA_starts,
                read_types=read_types,
                BC_fixed=bc_fixed_3p,
                umi_fixed=umi_fixed_3p,
            )
            get_5p_features(
                read_info=read_info,
                read_ids_5p=read_ids_5p,
                putative_bcs_5p=putative_bcs_5p,
                bc_fixed_locs_5p=bc_fixed_locs_5p,
                putative_bc_min_qs_5p=putative_bc_min_qs_5p,
                umis_5p=umis_5p,
                umi_fixed_locs_5p=umi_fixed_locs_5p,
                BC_fixed_5p=bc_fixed_5p,
                umi_fixed_5p=umi_fixed_5p,
            )

    rst_df_3p = pd.DataFrame(
        {
            "read_id": read_ids,
            "putative_bc": putative_bcs,
            "bc_fixed_locs": bc_fixed_locs,
            "putative_bc_min_qs": putative_bc_min_qs,
            "putative_umi": umis,
            "umi_fixed_locs": umi_fixed_locs,
            "post_umi_flankings": post_umi_flankings,
            "polyA_starts": polyA_starts,
            "read_types": read_types,
        }
    )
    rst_df_5p = pd.DataFrame(
        {
            "read_id": read_ids_5p,
            "putative_bc_5p": putative_bcs_5p,
            "bc_fixed_locs_50": bc_fixed_locs_5p,
            "putative_bc_min_qs_5p": putative_bc_min_qs_5p,
            "putative_umi_5p": umis_5p,
            "umi_fixed_locs_5p": umi_fixed_locs_5p,
        }
    )
    df_merge = rst_df_3p.merge(rst_df_5p, on="read_id", how="inner")

    if args.save_intermediate:
        rst_df_3p.to_csv(args.putative_bc_3p_out, index=False)
        rst_df_5p.to_csv(args.putative_bc_5p_out, index=False)
    df_merge.to_csv(args.putative_bc_out, index=False)
    return rst_df_3p, rst_df_5p, df_merge


def filter_corrected_reads(
    big_df,
    pair_min=None,
    auto_pair_min_floor=10,
    auto_pair_min_quantile=0.1,
):
    df = big_df.copy()
    df["BC3c"] = df["BC3_corrected"].map(norm_bc)
    df["BC5c"] = df["BC5_corrected"].map(norm_bc)

    paired = df[(df["BC3c"] != "") & (df["BC5c"] != "")].copy()
    if len(paired) > 0:
        paired["pair_u"] = paired.apply(
            lambda r: min(strip_fixed_5p(r["BC5c"]), revcomp(strip_fixed_3p(r["BC3c"]))),
            axis=1,
        )
        paired["pair_v"] = paired.apply(
            lambda r: max(strip_fixed_5p(r["BC5c"]), revcomp(strip_fixed_3p(r["BC3c"]))),
            axis=1,
        )
    else:
        paired["pair_u"] = ""
        paired["pair_v"] = ""
    pair_counts = paired.groupby(["pair_u", "pair_v"]).size().reset_index(name="pair_n_reads")
    resolved_pair_min = pair_min
    pair_min_mode = "manual"
    if pair_min is None:
        if len(pair_counts) == 0:
            resolved_pair_min = int(auto_pair_min_floor)
            pair_min_mode = "auto_empty"
        else:
            q = float(pair_counts["pair_n_reads"].quantile(auto_pair_min_quantile))
            resolved_pair_min = max(int(auto_pair_min_floor), int(q))
            pair_min_mode = "auto_quantile"
    paired2 = paired.merge(pair_counts, on=["pair_u", "pair_v"], how="left")
    bad_read_ids = set(paired2.loc[paired2["pair_n_reads"] < resolved_pair_min, "read_id"])

    df_keep = df[~df["read_id"].isin(bad_read_ids)].copy()
    df_keep = df_keep[(df_keep["BC3c"] != "") | (df_keep["BC5c"] != "")].copy()

    if ("putative_umi" not in df_keep.columns) or ("putative_umi_5p" not in df_keep.columns):
        raise KeyError("df missing putative_umi or putative_umi_5p.")

    miss3 = df_keep["putative_umi"].map(is_missing)
    miss5 = df_keep["putative_umi_5p"].map(is_missing)
    n_removed_both_missing = int((miss3 & miss5).sum())
    df_keep = df_keep[~(miss3 & miss5)].copy()
    df_keep = df_keep.drop(columns=["BC3c", "BC5c"], errors="ignore")

    drop_cols = [
        "read_types",
        "putative_bc",
        "bc_fixed_locs",
        "putative_bc_min_qs",
        "umi_fixed_locs",
        "polyA_starts",
        "post_umi_flankings",
        "umi_fixed_locs_5p",
        "bc_fixed_locs_50",
        "putative_bc_min_qs_5p",
        "putative_bc_5p",
        "strand",
    ]
    big_df_filtered = df_keep.drop(columns=drop_cols, errors="ignore")

    stats = {
        "original_rows": int(len(big_df)),
        "bad_paired_removed": int(len(bad_read_ids)),
        "both_empty_barcode_rows": int(((df["BC3c"] == "") & (df["BC5c"] == "")).sum()),
        "both_umi_missing_removed": int(n_removed_both_missing),
        "filtered_rows": int(len(big_df_filtered)),
        "PAIR_MIN": int(resolved_pair_min),
        "PAIR_MIN_mode": pair_min_mode,
        "auto_pair_min_floor": int(auto_pair_min_floor),
        "auto_pair_min_quantile": float(auto_pair_min_quantile),
    }
    return big_df_filtered, stats


def prepare_final_read_table(df):
    df = df.copy()
    df["BC5_20bp"] = df["BC5_corrected"].map(strip_fixed_5p)
    df["BC3_20bp_rc"] = df["BC3_corrected"].map(lambda x: revcomp(strip_fixed_3p(x)))
    return df.drop(columns=["BC3_corrected", "BC5_corrected", "has_3p", "has_5p"], errors="ignore")


def build_clean_exports(df_final, pc_final):
    df_out = df_final.copy()
    df_out["BC5n"] = df_out["BC5_20bp"].map(norm_bc) if "BC5_20bp" in df_out.columns else ""
    df_out["BC3n"] = df_out["BC3_20bp_rc"].map(norm_bc) if "BC3_20bp_rc" in df_out.columns else ""

    cols_read = ["read_id", "putative_umi", "putative_umi_5p", "BC5n", "BC3n"]
    df_out = df_out[[col for col in cols_read if col in df_out.columns]].copy()

    pair_cols = ["BC5n", "BC3n", "support_reads"]
    if "support_umis" in pc_final.columns:
        pair_cols.append("support_umis")
    pair_counts_kept = pc_final[pair_cols].copy()
    return df_out, pair_counts_kept


def write_single_row_tsv(path, data):
    pd.DataFrame([data]).to_csv(path, sep="\t", index=False)


def serialize_core_cells_df(core_cells_df):
    if core_cells_df is None or len(core_cells_df) == 0:
        return core_cells_df
    df = core_cells_df.copy()
    if "barcodes" in df.columns:
        df["barcodes"] = df["barcodes"].apply(
            lambda x: "|".join(map(str, x)) if isinstance(x, (list, tuple, set)) else str(x)
        )
    return df


def build_assigned_reads(
    df_out,
    pair_counts_kept,
    dominance_min,
    include_other_components=False,
    max_other_component_barcodes=8,
    absorb_unassigned_paired=True,
    min_reads_per_cell=20,
):
    _, _, core_cells_df, barcode2cell = build_core_cells(
        pair_counts_kept,
        include_other_components=include_other_components,
        max_other_component_barcodes=max_other_component_barcodes,
    )
    dom_table = compute_top1_dominance(pair_counts_kept)
    df_assigned, assign_stats = assign_reads_to_cells(
        df_out,
        barcode2cell,
        dom_table,
        dominance_min=dominance_min,
        absorb_unassigned_paired=absorb_unassigned_paired,
    )

    keep_cols = ["read_id", "putative_umi", "putative_umi_5p", "BC5n", "BC3n", "cell_id"]
    df_assigned_slim = df_assigned[[col for col in keep_cols if col in df_assigned.columns]].copy()
    df_assigned_slim = df_assigned_slim[
        df_assigned_slim["cell_id"].notna() & (df_assigned_slim["cell_id"] != "")
    ].copy()

    cell_read_stats_all = (
        df_assigned.dropna(subset=["cell_id"])
        .groupby("cell_id")["read_id"]
        .size()
        .sort_values(ascending=False)
        .reset_index(name="n_reads")
    )

    kept_cell_ids = set(
        cell_read_stats_all.loc[cell_read_stats_all["n_reads"] >= int(min_reads_per_cell), "cell_id"]
    )
    if kept_cell_ids:
        df_assigned_slim = df_assigned_slim[df_assigned_slim["cell_id"].isin(kept_cell_ids)].copy()
        core_cells_df = core_cells_df[core_cells_df["cell_id"].isin(kept_cell_ids)].copy()
        cell_read_stats = cell_read_stats_all[cell_read_stats_all["cell_id"].isin(kept_cell_ids)].copy()
    else:
        df_assigned_slim = df_assigned_slim.iloc[0:0].copy()
        core_cells_df = core_cells_df.iloc[0:0].copy()
        cell_read_stats = cell_read_stats_all.iloc[0:0].copy()

    if "type" in core_cells_df.columns:
        core_type_counts = core_cells_df["type"].value_counts()
    else:
        core_type_counts = pd.Series(dtype=int)

    filtered_core_barcode_count = 0
    if "barcodes" in core_cells_df.columns and len(core_cells_df) > 0:
        filtered_core_barcode_count = int(sum(len(x) for x in core_cells_df["barcodes"]))

    assign_stats["n_cells_before_min_reads"] = int(cell_read_stats_all["cell_id"].nunique())
    assign_stats["n_cells_after_min_reads"] = int(cell_read_stats["cell_id"].nunique())
    assign_stats["min_reads_per_cell"] = int(min_reads_per_cell)

    return (
        core_cells_df,
        df_assigned,
        df_assigned_slim,
        assign_stats,
        core_type_counts,
        cell_read_stats,
        filtered_core_barcode_count,
    )


def run_pipeline(args):
    if args.barcode_extract_mode != "fixed_seq":
        raise NotImplementedError(
            "--barcode-extract-mode probe is reserved for the future Sockeye-style "
            "dual-end probe extractor. Use --barcode-extract-mode fixed_seq for the "
            "current validated Flora dual-end barcode/UMI extraction logic."
        )

    step_t0 = time.perf_counter()
    green_msg("Step 1/7: splitting 3p and 5p putative barcode tables", printit=True)
    rst_df_3p, rst_df_5p, df_merge = build_putative_tables(args)
    print(f"[timing] Step 1/7 elapsed: {time.perf_counter() - step_t0:.2f}s")

    step_t0 = time.perf_counter()
    green_msg("Step 2/7: counting high-quality putative barcodes and building whitelists", printit=True)
    raw_bc_count_3p = count_high_quality_barcodes(rst_df_3p, "putative_bc", "putative_bc_min_qs", args.minQ)
    raw_bc_count_5p = count_high_quality_barcodes(rst_df_5p, "putative_bc_5p", "putative_bc_min_qs_5p", args.minQ)

    bc_whitelist_3p, ept_bc_3p = get_bc_whitelist(
        raw_bc_count=raw_bc_count_3p,
        full_bc_whitelist=args.full_bc_whitelist_3p,
        exp_cells=args.exp_cells,
        out_plot_fn=args.knee_plot_3p_out,
        DEFAULT_EMPTY_DROP_MIN_ED=args.DEFAULT_EMPTY_DROP_MIN_ED,
        DEFAULT_EMPTY_DROP_NUM=args.DEFAULT_EMPTY_DROP_NUM,
        reverse_complement_whitelist=args.revcomp_whitelist,
    )
    bc_whitelist_5p, ept_bc_5p = get_bc_whitelist(
        raw_bc_count=raw_bc_count_5p,
        full_bc_whitelist=args.full_bc_whitelist_5p,
        exp_cells=args.exp_cells,
        out_plot_fn=args.knee_plot_5p_out,
        DEFAULT_EMPTY_DROP_MIN_ED=args.DEFAULT_EMPTY_DROP_MIN_ED,
        DEFAULT_EMPTY_DROP_NUM=args.DEFAULT_EMPTY_DROP_NUM,
        reverse_complement_whitelist=args.revcomp_whitelist,
    )
    write_barcode_count_table(args.barcode_counts_3p_out, raw_bc_count_3p)
    write_barcode_count_table(args.barcode_counts_5p_out, raw_bc_count_5p)
    write_barcode_lines(args.whitelist_3p_out, bc_whitelist_3p.keys())
    write_barcode_lines(args.whitelist_5p_out, bc_whitelist_5p.keys())
    write_barcode_lines(args.emptydrop_3p_out, ept_bc_3p)
    write_barcode_lines(args.emptydrop_5p_out, ept_bc_5p)
    print(f"[timing] Step 2/7 elapsed: {time.perf_counter() - step_t0:.2f}s")

    step_t0 = time.perf_counter()
    green_msg("Step 3/7: correcting reads and writing matched FASTQ", printit=True)
    demul_count_tot, count_tot, big_df = assign_read(
        fastq_fns=args.fastq_fns,
        fastq_out=args.fastq_out,
        putative_bc_csv=args.putative_bc_out,
        whitelsit_3p=args.whitelist_3p_out,
        whitelsit_5p=args.whitelist_5p_out,
        max_ed=args.max_ed,
        n_process=args.threads,
        batchsize=args.assign_batchsize,
        minQ=args.minQ,
        write_fastq_out=not args.skip_matched_fastq,
        write_unmatched_fastq=not args.skip_unmatched_fastq,
    )
    if args.save_intermediate:
        big_df.to_csv(args.corrected_bc_out, index=False)
    barcode_validity_stats = summarize_barcode_validity(big_df)
    print(f"[timing] Step 3/7 elapsed: {time.perf_counter() - step_t0:.2f}s")

    step_t0 = time.perf_counter()
    green_msg("Step 4/7: filtering corrected reads", printit=True)
    big_df_filtered, filter_stats = filter_corrected_reads(
        big_df,
        pair_min=args.PAIR_MIN,
        auto_pair_min_floor=args.auto_pair_min_floor,
        auto_pair_min_quantile=args.auto_pair_min_quantile,
    )
    print(f"[timing] Step 4/7 elapsed: {time.perf_counter() - step_t0:.2f}s")

    step_t0 = time.perf_counter()
    green_msg("Step 5/7: stripping fixed sequences and filtering barcode pairs", printit=True)
    df_for_pairs = prepare_final_read_table(big_df_filtered)
    trimmed_barcode_uniques = summarize_trimmed_barcode_uniques(df_for_pairs)
    (
        df_final,
        _df_single,
        _df_paired_final,
        pc_all,
        pc_min,
        pc_final,
        pc_dropped,
        pair_stats,
    ) = filter_pairs_three_stage(
        df_for_pairs,
        bc5_col="BC5_20bp",
        bc3_col="BC3_20bp_rc",
        umi3_col="putative_umi",
        umi5_col="putative_umi_5p",
        PAIR_MIN=args.PAIR_MIN,
        auto_pair_min_floor=args.auto_pair_min_floor,
        auto_pair_min_quantile=args.auto_pair_min_quantile,
        TOP1_ALPHA=args.TOP1_ALPHA,
        TOP1_ALPHA_UMI=args.TOP1_ALPHA_UMI,
        require_pass_both_ends=args.require_pass_both_ends,
        drop_umiA_ratio_gt=args.drop_umiA_ratio_gt,
    )
    print(f"[timing] Step 5/7 elapsed: {time.perf_counter() - step_t0:.2f}s")

    step_t0 = time.perf_counter()
    green_msg("Step 6/7: exporting clean read table and pair table", printit=True)
    df_out, pair_counts_kept = build_clean_exports(df_final, pc_final)
    df_out.to_csv(args.clean_reads_out, index=False)
    pair_counts_kept.to_csv(args.pair_counts_out, index=False)
    if args.save_merge_debug:
        pc_all.to_csv(args.pair_counts_all_out, index=False)
        pc_min.to_csv(args.pair_counts_pairmin_kept_out, index=False)
        pc_dropped.to_csv(args.dropped_pairs_out, index=False)
    print(f"[timing] Step 6/7 elapsed: {time.perf_counter() - step_t0:.2f}s")

    step_t0 = time.perf_counter()
    green_msg("Step 7/7: assigning reads to cells and optionally extracting cell FASTQ", printit=True)
    (
        core_cells_df,
        df_assigned,
        df_assigned_slim,
        assign_stats,
        core_type_counts,
        cell_read_stats,
        core_barcode_count,
    ) = build_assigned_reads(
        df_out,
        pair_counts_kept,
        dominance_min=args.dominance_min,
        include_other_components=args.include_other_components,
        max_other_component_barcodes=args.max_other_component_barcodes,
        absorb_unassigned_paired=args.absorb_unassigned_paired,
        min_reads_per_cell=args.min_reads_per_cell,
    )
    cell_read_stats.to_csv(args.cell_read_stats_out, index=False)
    write_single_row_tsv(args.barcode_validity_summary_out, barcode_validity_stats)
    write_single_row_tsv(args.assign_stats_out, assign_stats)
    if args.save_merge_debug:
        serialize_core_cells_df(core_cells_df).to_csv(args.core_cells_debug_out, index=False)
        df_assigned.to_csv(args.read_assigned_debug_out, index=False)

    if args.skip_cell_fastq:
        df_kept = df_assigned_slim.copy()
        extract_stats = {
            "skipped_cell_fastq": True,
            "target_rows": int(len(df_assigned_slim)),
            "target_unique_ids": int(df_assigned_slim["read_id"].nunique()),
            "found_unique_ids": int(df_assigned_slim["read_id"].nunique()),
            "written_reads": 0,
            "scanned_reads": 0,
            "dropped_rows": 0,
            "missing_unique_ids": 0,
            "out_fastq_gz": "",
        }
    else:
        df_kept, extract_stats = extract_reads_and_filter_df_by_raw(
            df_assigned_slim,
            raw_fastq_gz=args.fastq_out,
            out_fastq_gz=args.cell_fastq_out,
            read_id_col="read_id",
            remove_found=True,
        )
    df_kept.to_csv(args.read_assigned_out, index=False, encoding="utf-8")
    print(f"[timing] Step 7/7 elapsed: {time.perf_counter() - step_t0:.2f}s")

    return {
        "fastq_files": len(args.fastq_fns),
        "reads_total": int(count_tot),
        "reads_demultiplexed": int(demul_count_tot),
        "barcode_validity_stats": barcode_validity_stats,
        "putative_rows_3p": int(len(rst_df_3p)),
        "putative_rows_5p": int(len(rst_df_5p)),
        "merged_rows": int(len(df_merge)),
        "filtered_rows": int(len(big_df_filtered)),
        "trimmed_barcode_uniques": trimmed_barcode_uniques,
        "clean_reads_rows": int(len(df_out)),
        "pair_counts_rows": int(len(pair_counts_kept)),
        "assigned_rows": int(len(df_kept)),
        "core_cells": int(len(core_cells_df)),
        "core_barcodes": int(core_barcode_count),
        "core_type_counts": core_type_counts.to_dict(),
        "cell_read_stats_head": cell_read_stats.head().to_dict(orient="records"),
        "filter_stats": filter_stats,
        "pair_stats": pair_stats,
        "assign_stats": assign_stats,
        "extract_stats": extract_stats,
    }


def print_summary(summary, out_dir):
    print("\n=== Flora Summary ===")
    print(f"Output directory: {out_dir}")
    print(f"FASTQ files: {summary['fastq_files']}")
    print(f"Reads total: {summary['reads_total']}")
    print(f"Reads demultiplexed: {summary['reads_demultiplexed']}")
    print("\nBarcode validity in corrected reads:")
    print(
        f"Any barcode valid: {summary['barcode_validity_stats']['valid_any_n']} "
        f"({summary['barcode_validity_stats']['valid_any_ratio']:.2%})"
    )
    print(
        f"Only 3' valid: {summary['barcode_validity_stats']['only_3p_n']} "
        f"({summary['barcode_validity_stats']['only_3p_ratio']:.2%})"
    )
    print(
        f"Only 5' valid: {summary['barcode_validity_stats']['only_5p_n']} "
        f"({summary['barcode_validity_stats']['only_5p_ratio']:.2%})"
    )
    print(
        f"Both ends valid: {summary['barcode_validity_stats']['both_n']} "
        f"({summary['barcode_validity_stats']['both_ratio']:.2%})"
    )
    print(
        f"Neither valid: {summary['barcode_validity_stats']['neither_n']} "
        f"({summary['barcode_validity_stats']['neither_ratio']:.2%})"
    )
    print(f"Merged putative rows: {summary['merged_rows']}")
    print(f"Filtered rows: {summary['filtered_rows']}")
    print(
        f"Unique BC3_20bp_rc: {summary['trimmed_barcode_uniques']['unique_bc3_20bp_rc']}"
    )
    print(
        f"Unique BC5_20bp: {summary['trimmed_barcode_uniques']['unique_bc5_20bp']}"
    )
    print(
        "Unique union (BC3 ∪ BC5): "
        f"{summary['trimmed_barcode_uniques']['unique_union_bc3_bc5']}"
    )
    print(f"Clean read rows: {summary['clean_reads_rows']}")
    print(f"Pair rows kept: {summary['pair_counts_rows']}")
    print(
        f"PAIR_MIN used: {summary['pair_stats']['PAIR_MIN']} "
        f"(mode={summary['pair_stats']['PAIR_MIN_mode']})"
    )
    print(
        f"TOP1 alpha (reads/UMI): {summary['pair_stats']['TOP1_ALPHA']} / "
        f"{summary['pair_stats']['TOP1_ALPHA_UMI']}"
    )
    print(f"Core cells: {summary['core_cells']} core barcodes: {summary['core_barcodes']}")
    print(f"Assigned read rows: {summary['assigned_rows']}")
    print("\nCore cell type counts:")
    for cell_type, count in summary["core_type_counts"].items():
        print(f"{cell_type:16} {count}")
    print("\nassign_stats:")
    print(summary["assign_stats"])
    print("\nTop cell_read_stats:")
    for row in summary["cell_read_stats_head"]:
        print(row)


def main():
    total_t0 = time.perf_counter()
    args = set_parser()
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    summary = run_pipeline(args)
    print_summary(summary, args.out_dir)
    print(f"\n[Flora] Total elapsed: {time.perf_counter() - total_t0:.2f}s")


if __name__ == "__main__":
    main()
