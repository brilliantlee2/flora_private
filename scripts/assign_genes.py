import argparse
import logging
import os
import re

import bioframe as bf
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "bed",
        help="BED file of alignments intervals",
        type=str,
    )

    parser.add_argument(
        "gtf",
        help="GTF file of gene annotations",
        type=str,
    )

    parser.add_argument(
        "-q",
        "--mapq",
        help="Minimum mapping quality to use for feature assignment [60]",
        type=int,
        default=60,
    )

    parser.add_argument(
        "--output",
        help="Output file [./read_annotations.tsv]",
        type=str,
        default="./read_annotations.tsv",
    )

    parser.add_argument(
        "-c",
        "--chunk_size",
        help="BED alignments per chunk to process [200000]",
        type=int,
        default=200000,
    )

    parser.add_argument(
        "--verbosity",
        help="logging level: <=2 logs info, <=3 logs warnings",
        type=int,
        default=2,
    )

    return parser.parse_args()


def init_logger(args):
    logging.basicConfig(
        format="%(asctime)s -- %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    logging_level = args.verbosity * 10
    logging.root.setLevel(logging_level)
    logging.root.handlers[0].addFilter(lambda x: "NumExpr" not in x.msg)


def extract_gene_label(attr):
    if pd.isna(attr):
        return np.nan

    attr = str(attr)

    m = re.search(r'gene_name "([^"]+)"', attr)
    if m:
        return m.group(1)

    m = re.search(r'gene_id "([^"]+)"', attr)
    if m:
        return m.group(1)

    return np.nan


def load_gtf(args):
    cols = [
        "chrom",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ]

    df = pd.read_csv(
        args.gtf,
        sep="\t",
        comment="#",
        header=None,
        names=cols,
        dtype=str,
        low_memory=False,
    )

    if df.shape[0] > 0:
        # 坐标列转成数值，避免 bioframe 校验失败
        df["start"] = pd.to_numeric(df["start"], errors="coerce")
        df["end"] = pd.to_numeric(df["end"], errors="coerce")

        df = df.dropna(subset=["chrom", "start", "end", "feature", "attribute"]).copy()

        # 只保留 gene 级别注释
        df = df[df["feature"] == "gene"].copy()

        if df.shape[0] == 0:
            raise ValueError(
                "No 'gene' features found in the input GTF. "
                "Please check whether this is a gene-level GTF."
            )

        df["attribute"] = df["attribute"].apply(extract_gene_label)
        df = df.dropna(subset=["attribute"]).copy()

        if df.shape[0] == 0:
            raise ValueError(
                "Failed to extract gene labels from GTF attribute column. "
                "Expected gene_name or gene_id in the 9th column."
            )

        assert bf.is_bedframe(df), "GTF file not loading as a valid dataframe!"

    return df


def normalize_bed_chunk(df, block_offset=0):
    cols = [
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
    ]

    df = df[cols].copy()
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["chrom", "start", "end", "name"]).copy()

    if df.shape[0] > 0:
        assert bf.is_bedframe(df), "BED file not loading as a valid dataframe!"

    df = df.reset_index(drop=True)
    df["index_bed"] = np.arange(block_offset, block_offset + len(df), dtype=np.int64)
    df["score"] = pd.to_numeric(df["score"], errors="coerce").fillna(0).astype(int)
    df["aln_len"] = df["end"] - df["start"]

    return df


def iter_bed_chunks(path, chunk_size):
    cols = ["chrom", "start", "end", "name", "score", "strand"]
    block_offset = 0
    reader = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=cols,
        chunksize=max(1, int(chunk_size)),
    )
    for raw_chunk in reader:
        chunk = normalize_bed_chunk(raw_chunk, block_offset=block_offset)
        block_offset += len(chunk)
        if len(chunk) > 0:
            yield chunk


def assign_status_low_mapq(df, args):
    df.loc[df["score"] < args.mapq, "status"] = "Unassigned_mapq"
    df.loc[df["score"] < args.mapq, "gene"] = "NA"
    return df


def assign_status_ambiguous_overlap(df):
    is_ambiguous = df[df["status"] == "Unknown"].duplicated(
        subset=["index_bed", "overlap_bp"], keep=False
    )
    ambiguous_idx = is_ambiguous.index[is_ambiguous]
    df.loc[ambiguous_idx, "status"] = "Unassigned_ambiguous"
    df.loc[ambiguous_idx, "gene"] = "NA"
    df = df.drop_duplicates(subset=["index_bed", "overlap_bp", "status"])

    return df


def assign_status_no_features(df):
    df.loc[df["gene"] == 0, "status"] = "Unassigned_no_features"
    df.loc[df["gene"] == 0, "gene"] = "NA"
    return df


def find_largest_overlap(df):
    max_ovlp_idx = df.groupby(["index_bed"])["overlap_bp"].idxmax().sort_values().values
    df = df.loc[max_ovlp_idx, :]
    unknown_idx = [i for i in max_ovlp_idx if df.loc[i, "status"] == "Unknown"]
    df.loc[unknown_idx, "status"] = "Assigned"
    return df


def get_overlaps(bed, gtf):
    df = bf.overlap(
        bed,
        gtf,
        how="left",
        suffixes=("_bed", "_gtf"),
        return_overlap=True,
        return_index=True,
    )

    df = df[
        [
            "index_bed",
            "name_bed",
            "chrom_bed",
            "score_bed",
            "strand_bed",
            "strand_gtf",
            "attribute_gtf",
            "overlap_start",
            "overlap_end",
        ]
    ].fillna(0)

    df = df.rename(
        columns={
            "name_bed": "read",
            "chrom_bed": "chrom",
            "score_bed": "score",
            "strand_bed": "read_strand",
            "strand_gtf": "gene_strand",
            "attribute_gtf": "gene",
        }
    )
    df["score"] = df["score"].astype(int)
    df["overlap_bp"] = df["overlap_end"] - df["overlap_start"]

    # Reads are already oriented to a biological direction upstream, so
    # only same-strand annotation overlaps should be eligible for gene calls.
    mismatch = (
        (df["gene"] != 0)
        & df["read_strand"].isin(["+", "-"])
        & df["gene_strand"].isin(["+", "-"])
        & (df["read_strand"] != df["gene_strand"])
    )
    if mismatch.any():
        df.loc[mismatch, "gene"] = 0
        df.loc[mismatch, "overlap_bp"] = 0

    return df


def process_bed_chunk(bed_chunk, gtf, args):
    gtf_sub = gtf[gtf["chrom"].isin(bed_chunk["chrom"].unique())].copy()
    if gtf_sub.empty:
        df_chunk = bed_chunk[["name", "score", "index_bed"]].copy()
        df_chunk["status"] = np.where(df_chunk["score"] < args.mapq, "Unassigned_mapq", "Unassigned_no_features")
        df_chunk["gene"] = "NA"
        df_chunk = df_chunk.rename(columns={"name": "read"})
        return df_chunk[["read", "status", "score", "gene"]].reset_index(drop=True)

    df_chunk = get_overlaps(bed_chunk, gtf_sub)
    df_chunk["status"] = "Unknown"
    df_chunk = assign_status_low_mapq(df_chunk, args)
    df_chunk = assign_status_ambiguous_overlap(df_chunk)
    df_chunk = assign_status_no_features(df_chunk)
    df_chunk = find_largest_overlap(df_chunk)

    df_chunk = df_chunk[["read", "status", "score", "gene", "index_bed"]]
    df_chunk = df_chunk.reset_index(drop=True)
    df_chunk = df_chunk.drop(["index_bed"], axis=1)

    return df_chunk


def main(args):
    init_logger(args)

    gtf = load_gtf(args)
    open(args.output, "w").close()
    if gtf.shape[0] > 0:
        for i, bed_chunk in enumerate(iter_bed_chunks(args.bed, args.chunk_size), 1):
            df_chunk = process_bed_chunk(bed_chunk, gtf, args)
            df_chunk.to_csv(args.output, sep="\t", index=False, header=False, mode="a")
            logger.info(
                "Processed gene-assignment chunk %s: bed_rows=%s reads=%s",
                i,
                len(bed_chunk),
                len(df_chunk),
            )


if __name__ == "__main__":
    args = parse_args()
    main(args)
