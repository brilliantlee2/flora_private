import argparse
import logging
import os
import re
import subprocess
import tempfile

import bioframe as bf
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("bam", help="Input BAM file", type=str)
    parser.add_argument("gtf", help="Input GTF file", type=str)
    parser.add_argument(
        "--output",
        help="Output transcript assignment TSV",
        type=str,
        default="read_transcript_assigns.tsv",
    )
    parser.add_argument(
        "-q",
        "--mapq",
        help="Minimum MAPQ to use for transcript assignment [60]",
        type=int,
        default=60,
    )
    parser.add_argument(
        "-c",
        "--chunk_size",
        help="BED alignment blocks per chunk to process [200000]",
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
    logging.root.setLevel(args.verbosity * 10)


def bam_to_split_bed(bam_path):
    tmp = tempfile.NamedTemporaryFile(suffix=".split.bed", delete=False)
    tmp.close()
    cmd = ["bedtools", "bamtobed", "-split", "-i", bam_path]
    with open(tmp.name, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)
    return tmp.name


def extract_attr(attr, key):
    m = re.search(rf'{key} "([^"]+)"', str(attr))
    if m:
        return m.group(1)
    return np.nan


def load_gtf(gtf_path):
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
        gtf_path,
        sep="\t",
        comment="#",
        header=None,
        names=cols,
        dtype=str,
        low_memory=False,
    )
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["chrom", "start", "end", "feature", "attribute"]).copy()

    # Use exon blocks instead of transcript span so splice structure contributes
    # to transcript assignment.
    df = df[df["feature"] == "exon"].copy()
    if df.shape[0] == 0:
        raise ValueError("No exon features found in GTF")

    df["transcript_id"] = df["attribute"].apply(
        lambda x: extract_attr(x, "transcript_id")
    )
    df["gene_name"] = df["attribute"].apply(lambda x: extract_attr(x, "gene_name"))
    df["gene_id"] = df["attribute"].apply(lambda x: extract_attr(x, "gene_id"))
    df["gene_label"] = df["gene_name"].fillna(df["gene_id"])
    df = df.dropna(subset=["transcript_id"]).copy()

    assert bf.is_bedframe(df), "GTF exon table is not a valid bedframe"
    return df


def normalize_bed_chunk(df, block_offset=0):
    cols = ["chrom", "start", "end", "name", "score", "strand"]
    df = df[cols].copy()
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["chrom", "start", "end", "name"]).copy()
    if df.shape[0] > 0:
        assert bf.is_bedframe(df), "BED file not loading as a valid bedframe"
    df["score"] = pd.to_numeric(df["score"], errors="coerce").fillna(0).astype(int)
    df = df.reset_index(drop=True)
    df["block_index"] = np.arange(block_offset, block_offset + len(df), dtype=np.int64)
    return df


def iter_bed_chunks(bed_path, chunk_size):
    cols = ["chrom", "start", "end", "name", "score", "strand"]
    chunk_size = max(1, int(chunk_size))
    carry = pd.DataFrame(columns=cols)
    block_offset = 0

    reader = pd.read_csv(
        bed_path,
        sep="\t",
        header=None,
        names=cols,
        chunksize=chunk_size,
    )
    for raw_chunk in reader:
        if len(carry) > 0:
            raw_chunk = pd.concat([carry, raw_chunk], ignore_index=True)
            carry = carry.iloc[0:0].copy()

        if raw_chunk.empty:
            continue

        # bedtools emits split blocks from one alignment consecutively. Keep the
        # final read_id for the next chunk so a multi-exon read is never assigned
        # using only a subset of its blocks.
        last_read = str(raw_chunk.iloc[-1]["name"])
        tail_mask = raw_chunk["name"].astype(str) == last_read
        process = raw_chunk.loc[~tail_mask].copy()
        carry = raw_chunk.loc[tail_mask].copy()

        if process.empty:
            continue

        chunk = normalize_bed_chunk(process, block_offset=block_offset)
        block_offset += len(chunk)
        yield chunk

    if len(carry) > 0:
        chunk = normalize_bed_chunk(carry, block_offset=block_offset)
        if len(chunk) > 0:
            yield chunk


def get_overlaps(bed, exons):
    df = bf.overlap(
        bed,
        exons,
        how="left",
        suffixes=("_bed", "_gtf"),
        return_overlap=True,
        return_index=True,
    )
    df = df[
        [
            "index_bed",
            "name_bed",
            "score_bed",
            "strand_bed",
            "transcript_id_gtf",
            "gene_label_gtf",
            "strand_gtf",
            "overlap_start",
            "overlap_end",
        ]
    ].fillna(0)

    df = df.rename(
        columns={
            "name_bed": "read_id",
            "score_bed": "score",
            "strand_bed": "read_strand",
            "transcript_id_gtf": "transcript_id",
            "gene_label_gtf": "gene",
            "strand_gtf": "tx_strand",
        }
    )
    df["score"] = df["score"].astype(int)
    df["overlap_bp"] = df["overlap_end"] - df["overlap_start"]

    mismatch = (
        (df["transcript_id"] != 0)
        & df["read_strand"].isin(["+", "-"])
        & df["tx_strand"].isin(["+", "-"])
        & (df["read_strand"] != df["tx_strand"])
    )
    if mismatch.any():
        df.loc[mismatch, "transcript_id"] = 0
        df.loc[mismatch, "gene"] = 0
        df.loc[mismatch, "overlap_bp"] = 0

    return df


def assign_transcripts_for_reads(bed, exons, mapq):
    read_meta = (
        bed.groupby("name", as_index=False)
        .agg(score=("score", "max"), read_order=("block_index", "min"))
        .rename(columns={"name": "read_id"})
        .sort_values("read_order")
        .reset_index(drop=True)
    )
    if read_meta.empty:
        return pd.DataFrame(columns=["read_id", "status", "score", "gene", "transcript_id"])

    exons_sub = exons[exons["chrom"].isin(bed["chrom"].unique())].copy()
    if exons_sub.empty:
        out = read_meta[["read_id", "score"]].copy()
        out["status"] = np.where(out["score"] < mapq, "Unassigned_mapq", "Unassigned_no_features")
        out["gene"] = "NA"
        out["transcript_id"] = "NA"
        return out[["read_id", "status", "score", "gene", "transcript_id"]]

    overlaps = get_overlaps(bed, exons_sub)
    overlaps = overlaps[overlaps["transcript_id"] != 0].copy()

    if overlaps.empty:
        out = read_meta[["read_id", "score"]].copy()
        out["status"] = np.where(out["score"] < mapq, "Unassigned_mapq", "Unassigned_no_features")
        out["gene"] = "NA"
        out["transcript_id"] = "NA"
        return out[["read_id", "status", "score", "gene", "transcript_id"]]

    tx_support = (
        overlaps.groupby(["read_id", "transcript_id", "gene"], as_index=False)
        .agg(score=("score", "max"), overlap_bp=("overlap_bp", "sum"))
    )
    tx_support = tx_support.merge(read_meta, on="read_id", how="left", suffixes=("", "_meta"))
    tx_support["score"] = tx_support["score_meta"].fillna(tx_support["score"]).astype(int)
    tx_support = tx_support.drop(columns=["score_meta"])

    rows = []
    grouped = tx_support.groupby("read_id", sort=False)
    for read in read_meta.itertuples(index=False):
        read_id = read.read_id
        read_score = int(read.score)

        if read_score < mapq:
            rows.append((read_id, "Unassigned_mapq", read_score, "NA", "NA"))
            continue

        if read_id not in grouped.groups:
            rows.append((read_id, "Unassigned_no_features", read_score, "NA", "NA"))
            continue

        g = grouped.get_group(read_id).copy()
        max_overlap = g["overlap_bp"].max()
        top = g[g["overlap_bp"] == max_overlap].copy()

        if len(top) > 1:
            rows.append((read_id, "Unassigned_ambiguous", read_score, "NA", "NA"))
            continue

        best = top.iloc[0]
        rows.append(
            (
                read_id,
                "Assigned",
                read_score,
                best["gene"],
                best["transcript_id"],
            )
        )

    return pd.DataFrame(
        rows, columns=["read_id", "status", "score", "gene", "transcript_id"]
    )


def main(args):
    init_logger(args)

    exons = load_gtf(args.gtf)
    bed_path = bam_to_split_bed(args.bam)
    n_chunks = 0
    n_rows = 0

    try:
        open(args.output, "w").close()
        if exons.shape[0] > 0:
            for n_chunks, bed_chunk in enumerate(iter_bed_chunks(bed_path, args.chunk_size), 1):
                out = assign_transcripts_for_reads(bed_chunk, exons, args.mapq)
                n_rows += len(out)
                out.to_csv(args.output, sep="\t", index=False, header=False, mode="a")
                logger.info(
                    "Processed transcript-assignment chunk %s: bed_blocks=%s reads=%s",
                    n_chunks,
                    len(bed_chunk),
                    len(out),
                )
        logger.info("Finished transcript assignment: chunks=%s reads=%s", n_chunks, n_rows)
    finally:
        os.remove(bed_path)


if __name__ == "__main__":
    args = parse_args()
    main(args)
