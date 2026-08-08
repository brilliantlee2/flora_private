import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="输入表格，支持 .tsv/.csv")
    parser.add_argument("--output", required=True, help="输出 TSV")
    parser.add_argument("--read-id-col", default="read_id")
    parser.add_argument("--cell-col", default="cell_id")
    parser.add_argument("--umi-primary-col", default="putative_umi_5p")
    parser.add_argument("--umi-backup-col", default="putative_umi")
    parser.add_argument("--barcode-5p-col", default="BC5n")
    parser.add_argument("--barcode-3p-col", default="BC3n")
    return parser.parse_args()


def read_table(path):
    if path.endswith(".csv"):
        return pd.read_csv(path)
    return pd.read_csv(path, sep="\t")


def main():
    args = parse_args()
    df = read_table(args.input).copy()

    required = [
        args.read_id_col,
        args.cell_col,
        args.umi_primary_col,
        args.umi_backup_col,
        args.barcode_5p_col,
        args.barcode_3p_col,
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"缺少列: {missing}")

    df[args.umi_primary_col] = df[args.umi_primary_col].fillna("").astype(str).str.strip()
    df[args.umi_backup_col] = df[args.umi_backup_col].fillna("").astype(str).str.strip()
    df[args.barcode_5p_col] = df[args.barcode_5p_col].fillna("").astype(str).str.strip()
    df[args.barcode_3p_col] = df[args.barcode_3p_col].fillna("").astype(str).str.strip()
    df[args.read_id_col] = df[args.read_id_col].astype(str).str.strip()
    df[args.cell_col] = df[args.cell_col].astype(str).str.strip()

    # Strint2.1 keeps both barcode ends, while still emitting Sockeye-style
    # single CB/CR/UR fields for downstream BAM tag tools.
    df["umi_for_clustering"] = df[args.umi_primary_col]
    mask_empty = df["umi_for_clustering"].eq("")
    df.loc[mask_empty, "umi_for_clustering"] = df.loc[mask_empty, args.umi_backup_col]
    df["barcode_dual"] = df[args.barcode_5p_col] + "+" + df[args.barcode_3p_col]
    df.loc[df[args.barcode_5p_col].eq(""), "barcode_dual"] = "+" + df.loc[
        df[args.barcode_5p_col].eq(""), args.barcode_3p_col
    ]
    df.loc[df[args.barcode_3p_col].eq(""), "barcode_dual"] = df.loc[
        df[args.barcode_3p_col].eq(""), args.barcode_5p_col
    ] + "+"

    out = df[
        [
            args.read_id_col,
            args.cell_col,
            args.barcode_5p_col,
            args.barcode_3p_col,
            "barcode_dual",
            args.umi_primary_col,
            args.umi_backup_col,
            "umi_for_clustering",
        ]
    ].copy()
    out.columns = [
        "read_id",
        "cell_id",
        "barcode_5p",
        "barcode_3p",
        "barcode_dual",
        "umi_primary",
        "umi_backup",
        "umi_for_clustering",
    ]

    out = out[out["read_id"] != ""]
    out = out[out["cell_id"] != ""]
    out = out[out["umi_for_clustering"] != ""]

    # 如果同一个 read_id 有重复记录，保留第一条
    out = out.drop_duplicates(subset=["read_id"], keep="first")

    out.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
