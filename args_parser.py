import argparse
import os
import sys
from pathlib import Path

from utils import err_msg


def check_files_exist(file_list):
    if isinstance(file_list, str):
        file_list = [file_list]
    exit_code = 0
    for fn in file_list:
        if not os.path.exists(fn):
            exit_code = 1
            err_msg(f"Error: can not find {fn}", printit=True)
    if exit_code == 1:
        sys.exit(1)
    return True


def get_files_by_suffix(search_dir, suffix, recursive=True):
    files = []
    if isinstance(suffix, str):
        suffix = [suffix]
    for item in suffix:
        globber = Path(search_dir).rglob if recursive else Path(search_dir).glob
        files.extend(globber(item))
    return sorted(files)


def get_files_from_dir(fastq_dir):
    check_files_exist(fastq_dir)
    if os.path.isdir(fastq_dir):
        fastq_fns = get_files_by_suffix(
            search_dir=fastq_dir,
            suffix=["*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"],
            recursive=True,
        )
        if not fastq_fns:
            err_msg(f"No FASTQ files found in directory: {fastq_dir}", printit=True)
            sys.exit(1)
    elif os.path.isfile(fastq_dir):
        fastq_fns = [fastq_dir]
    else:
        err_msg(f"File type of input {fastq_dir} is not supported.", printit=True)
        sys.exit(1)
    return list(map(str, fastq_fns))


def _to_outdir(out_dir, path):
    return path if os.path.isabs(path) else os.path.join(out_dir, os.path.basename(path))


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Flora: full-length read barcode splitting, whitelist generation, "
            "read correction, pair filtering, and cell assignment."
        )
    )
    parser.add_argument(
        "fastq_fns",
        metavar="<input full-length fastq filename/directory>",
        help="Full-length FASTQ file (.fq/.fastq/.gz) or a directory containing FASTQ files.",
        type=get_files_from_dir,
    )
    parser.add_argument("--full-bc-whitelist-3p", dest="full_bc_whitelist_3p", required=True)
    parser.add_argument("--full-bc-whitelist-5p", dest="full_bc_whitelist_5p", required=True)
    parser.add_argument(
        "--no-revcomp-whitelist",
        dest="revcomp_whitelist",
        action="store_false",
        help="Use whitelist sequences as-is instead of reverse-complementing them before matching.",
    )
    parser.set_defaults(revcomp_whitelist=True)

    parser.add_argument("--out_dir", type=str, default="Flora")
    parser.add_argument("--save-intermediate", dest="save_intermediate", action="store_true")
    parser.add_argument("--light-output", dest="light_output", action="store_true")
    parser.add_argument("--full-output", dest="light_output", action="store_false")
    parser.set_defaults(light_output=True)
    parser.add_argument(
        "--barcode_extract_mode",
        "--barcode-extract-mode",
        choices=["fixed_seq", "probe"],
        default="fixed_seq",
        help=(
            "Barcode/UMI extraction mode. fixed_seq uses the current Flora dual-end "
            "fixed-sequence extractor. probe is reserved for a future Sockeye-style "
            "local-alignment extractor."
        ),
    )
    parser.add_argument("--batch_size", type=int, default=100000)
    parser.add_argument("--assign_batchsize", type=int, default=10000)
    parser.add_argument("--threads", type=int, default=32)
    parser.add_argument("--minQ", type=int, default=2)
    parser.add_argument("--exp_cells", type=int, default=5000)
    parser.add_argument("--max_ed", type=int, default=2)
    parser.add_argument("--DEFAULT_EMPTY_DROP_MIN_ED", type=int, default=5)
    parser.add_argument("--DEFAULT_EMPTY_DROP_NUM", type=int, default=2000)
    parser.add_argument("--PAIR_MIN", type=int, default=None)
    parser.add_argument("--auto_pair_min_floor", type=int, default=10)
    parser.add_argument("--auto_pair_min_quantile", type=float, default=0.1)
    parser.add_argument("--TOP1_ALPHA", type=float, default=0.1)
    parser.add_argument("--TOP1_ALPHA_UMI", type=float, default=0.3)
    parser.add_argument("--dominance_min", type=float, default=0.80)
    parser.add_argument("--drop_umiA_ratio_gt", type=float, default=0.5)
    parser.add_argument("--require_pass_both_ends", action="store_true")

    parser.add_argument("--BC_fixed_3p", type=str, default="GGTAGC")
    parser.add_argument("--umi_fixed_3p", type=str, default="GATCT")
    parser.add_argument("--BC_fixed_5p", type=str, default="GGAAGG")
    parser.add_argument("--umi_fixed_5p", type=str, default="CAGCA")

    parser.add_argument("--putative_bc_out", type=str, default="putative_bc.csv")
    parser.add_argument("--putative_bc_3p_out", type=str, default="putative_bc_p3.csv")
    parser.add_argument("--putative_bc_5p_out", type=str, default="putative_bc_p5.csv")
    parser.add_argument("--whitelist_3p_out", type=str, default="whitelist_3p.csv")
    parser.add_argument("--whitelist_5p_out", type=str, default="whitelist_5p.csv")
    parser.add_argument("--barcode_counts_3p_out", type=str, default="barcode_counts_3p.tsv")
    parser.add_argument("--barcode_counts_5p_out", type=str, default="barcode_counts_5p.tsv")
    parser.add_argument("--emptydrop_3p_out", type=str, default="empty_bc_list_3p.csv")
    parser.add_argument("--emptydrop_5p_out", type=str, default="empty_bc_list_5p.csv")
    parser.add_argument("--knee_plot_3p_out", type=str, default="knee_plot_3p.png")
    parser.add_argument("--knee_plot_5p_out", type=str, default="knee_plot_5p.png")
    parser.add_argument("--fastq_out", type=str, default="matched_reads.fastq.gz")
    parser.add_argument("--corrected_bc_out", type=str, default="BC_corrected.csv")
    parser.add_argument("--clean_reads_out", type=str, default="ReadIDs_UMI_BC_clean.csv")
    parser.add_argument("--pair_counts_out", type=str, default="pair_counts_kept.csv")
    parser.add_argument("--read_assigned_out", type=str, default="read_assigned_cell.csv")
    parser.add_argument("--cell_fastq_out", type=str, default="cell_reads.fastq.gz")
    parser.add_argument("--cell_read_stats_out", type=str, default="cell_read_stats.csv")
    parser.add_argument("--barcode_validity_summary_out", type=str, default="barcode_validity_summary.tsv")
    parser.add_argument("--assign_stats_out", type=str, default="assign_stats.tsv")
    parser.add_argument("--save_merge_debug", action="store_true")
    parser.add_argument("--pair_counts_all_out", type=str, default="pair_counts_all.csv")
    parser.add_argument("--pair_counts_pairmin_kept_out", type=str, default="pair_counts_pairmin_kept.csv")
    parser.add_argument("--dropped_pairs_out", type=str, default="dropped_pairs.csv")
    parser.add_argument("--core_cells_debug_out", type=str, default="core_cells_debug.csv")
    parser.add_argument("--read_assigned_debug_out", type=str, default="read_assigned_debug.csv")
    parser.add_argument("--include_other_components", action="store_true", default=False)
    parser.add_argument("--exclude_other_components", dest="include_other_components", action="store_false")
    parser.add_argument("--max_other_component_barcodes", type=int, default=8)
    parser.add_argument("--absorb_unassigned_paired", action="store_true", default=True)
    parser.add_argument("--disable_absorb_unassigned_paired", dest="absorb_unassigned_paired", action="store_false")
    parser.add_argument("--min_reads_per_cell", type=int, default=20)
    return parser


def set_parser(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    check_files_exist([args.full_bc_whitelist_3p, args.full_bc_whitelist_5p])

    if args.light_output:
        args.skip_matched_fastq = True
        args.skip_unmatched_fastq = True
        args.skip_cell_fastq = True
    else:
        args.skip_matched_fastq = False
        args.skip_unmatched_fastq = False
        args.skip_cell_fastq = False

    os.makedirs(args.out_dir, exist_ok=True)
    for attr in [
        "putative_bc_out",
        "putative_bc_3p_out",
        "putative_bc_5p_out",
        "whitelist_3p_out",
        "whitelist_5p_out",
        "barcode_counts_3p_out",
        "barcode_counts_5p_out",
        "emptydrop_3p_out",
        "emptydrop_5p_out",
        "knee_plot_3p_out",
        "knee_plot_5p_out",
        "fastq_out",
        "corrected_bc_out",
        "clean_reads_out",
        "pair_counts_out",
        "read_assigned_out",
        "cell_fastq_out",
        "cell_read_stats_out",
        "barcode_validity_summary_out",
        "assign_stats_out",
        "pair_counts_all_out",
        "pair_counts_pairmin_kept_out",
        "dropped_pairs_out",
        "core_cells_debug_out",
        "read_assigned_debug_out",
    ]:
        setattr(args, attr, _to_outdir(args.out_dir, getattr(args, attr)))

    print(f"[Flora] Collected {len(args.fastq_fns)} FASTQ files.")
    print(f"[Flora] Output directory: {args.out_dir}")
    return args


if __name__ == "__main__":
    parsed_args = set_parser()
    print("\n=== Parsed arguments ===")
    for key, value in vars(parsed_args).items():
        print(f"{key:24}: {value}")
