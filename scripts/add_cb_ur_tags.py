import argparse
import pysam
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True, help="输入 BAM")
    parser.add_argument("--tags", required=True, help="read_tags.tsv")
    parser.add_argument("--output", required=True, help="输出 BAM")
    parser.add_argument("--read-id-col", default="read_id")
    parser.add_argument("--cb-col", default="cell_id")
    parser.add_argument("--ur-col", default="umi_for_clustering")
    parser.add_argument("--cr-col", default="barcode_dual")
    parser.add_argument("--bc5-col", default="barcode_5p")
    parser.add_argument("--bc3-col", default="barcode_3p")
    parser.add_argument("--umi5-col", default="umi_primary")
    parser.add_argument("--umi3-col", default="umi_backup")
    parser.add_argument(
        "--keep-untagged",
        action="store_true",
        help="保留无法写入 CB/UR 的 alignments；默认按 Sockeye 逻辑丢弃。",
    )
    return parser.parse_args()


def load_tag_map(path, args):
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    required = [args.read_id_col, args.cb_col, args.ur_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"缺少列: {missing}")

    df = df[df[args.read_id_col].str.strip() != ""].copy()
    df = df.drop_duplicates(subset=[args.read_id_col], keep="first")

    tag_map = {}
    for _, row in df.iterrows():
        rid = row[args.read_id_col].strip()
        cb = row[args.cb_col].strip()
        ur = row[args.ur_col].strip()
        if cb and ur:
            tag_map[rid] = {
                "CB": cb,
                "UR": ur,
                "CR": row.get(args.cr_col, "").strip(),
                "C5": row.get(args.bc5_col, "").strip(),
                "C3": row.get(args.bc3_col, "").strip(),
                "U5": row.get(args.umi5_col, "").strip(),
                "U3": row.get(args.umi3_col, "").strip(),
            }
    return tag_map


def main():
    args = parse_args()
    tag_map = load_tag_map(args.tags, args)

    bam = pysam.AlignmentFile(args.bam, "rb")
    out = pysam.AlignmentFile(args.output, "wb", template=bam)

    total = 0
    tagged = 0
    total_read_ids = set()
    tagged_read_ids = set()

    for aln in bam.fetch(until_eof=True):
        total += 1
        rid = aln.query_name
        total_read_ids.add(rid)
        if rid in tag_map:
            tags = tag_map[rid]
            aln.set_tag("CB", tags["CB"], value_type="Z")
            aln.set_tag("UR", tags["UR"], value_type="Z")
            for tag in ["CR", "C5", "C3", "U5", "U3"]:
                if tags.get(tag):
                    aln.set_tag(tag, tags[tag], value_type="Z")
            tagged += 1
            tagged_read_ids.add(rid)
            out.write(aln)
        elif args.keep_untagged:
            out.write(aln)

    bam.close()
    out.close()
    pysam.index(args.output)

    print(f"total_alignments\t{len(total_read_ids)}")
    print(f"tagged_alignments\t{len(tagged_read_ids)}")
    print(f"unmatched_alignments\t{len(total_read_ids - tagged_read_ids)}")
    print(f"total_alignment_records\t{total}")
    print(f"tagged_alignment_records\t{tagged}")
    print(f"unmatched_alignment_records\t{total - tagged}")


if __name__ == "__main__":
    main()
