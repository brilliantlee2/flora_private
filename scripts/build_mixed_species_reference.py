#!/usr/bin/env python3
import argparse
import os
import re
import shutil
import subprocess
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build a prefixed human/mouse mixed-species reference from two "
            "single-species reference directories."
        )
    )
    parser.add_argument("--human-ref-dir", required=True, help="Directory containing GRCh38-style reference files")
    parser.add_argument("--mouse-ref-dir", required=True, help="Directory containing GRCm39-style reference files")
    parser.add_argument("--out-dir", required=True, help="Output directory for the merged reference")
    parser.add_argument("--human-prefix", default="hs", help="Prefix for human contigs/IDs/names [hs]")
    parser.add_argument("--mouse-prefix", default="mm", help="Prefix for mouse contigs/IDs/names [mm]")
    parser.add_argument("--human-label", default="human", help="Label used in auxiliary filenames [human]")
    parser.add_argument("--mouse-label", default="mouse", help="Label used in auxiliary filenames [mouse]")
    return parser.parse_args()


def require_file(path: Path, label: str):
    if not path.is_file():
        raise FileNotFoundError(f"Missing {label}: {path}")


def prefix_value(raw: str, prefix: str) -> str:
    raw = str(raw)
    return raw if raw.startswith(f"{prefix}_") else f"{prefix}_{raw}"


def prefix_gtf_attributes(attr: str, prefix: str) -> str:
    patterns = [
        r'(gene_id ")([^"]+)(")',
        r'(transcript_id ")([^"]+)(")',
        r'(gene_name ")([^"]+)(")',
    ]
    for pattern in patterns:
        attr = re.sub(
            pattern,
            lambda m: f'{m.group(1)}{prefix_value(m.group(2), prefix)}{m.group(3)}',
            attr,
        )
    return attr


def prefix_bed_name(name: str, prefix: str) -> str:
    parts = str(name).split("|")
    if len(parts) >= 1 and parts[0]:
        parts[0] = prefix_value(parts[0], prefix)
    if len(parts) >= 3 and parts[2]:
        parts[2] = prefix_value(parts[2], prefix)
    return "|".join(parts)


def process_fasta(src: Path, dst: Path, prefix: str):
    with src.open("r") as fin, dst.open("w") as fout:
        for line in fin:
            if line.startswith(">"):
                header = line[1:].rstrip("\n")
                if not header:
                    fout.write(line)
                    continue
                parts = header.split(maxsplit=1)
                contig = prefix_value(parts[0], prefix)
                if len(parts) == 2:
                    fout.write(f">{contig} {parts[1]}\n")
                else:
                    fout.write(f">{contig}\n")
            else:
                fout.write(line)


def process_chrom_sizes(src: Path, dst: Path, prefix: str):
    with src.open("r") as fin, dst.open("w") as fout:
        for line in fin:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"Invalid chrom_sizes line in {src}: {line.rstrip()}")
            fields[0] = prefix_value(fields[0], prefix)
            fout.write("\t".join(fields[:2]) + "\n")


def process_gtf(src: Path, dst: Path, prefix: str):
    with src.open("r") as fin, dst.open("w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Invalid GTF line in {src}: {line.rstrip()}")
            fields[0] = prefix_value(fields[0], prefix)
            fields[8] = prefix_gtf_attributes(fields[8], prefix)
            fout.write("\t".join(fields) + "\n")


def process_bed(src: Path, dst: Path, prefix: str):
    with src.open("r") as fin, dst.open("w") as fout:
        for line in fin:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise ValueError(f"Invalid BED line in {src}: {line.rstrip()}")
            fields[0] = prefix_value(fields[0], prefix)
            fields[3] = prefix_bed_name(fields[3], prefix)
            fout.write("\t".join(fields) + "\n")


def concatenate_text(inputs, output):
    with output.open("w") as fout:
        for path in inputs:
            with path.open("r") as fin:
                shutil.copyfileobj(fin, fout)


def run_samtools_faidx(genome_fa: Path):
    if shutil.which("samtools") is None:
        print(f"[WARN] samtools not found; skipped faidx for {genome_fa}")
        return
    subprocess.run(["samtools", "faidx", str(genome_fa)], check=True)


def write_manifest(out_dir: Path, human_ref_dir: Path, mouse_ref_dir: Path, human_prefix: str, mouse_prefix: str):
    manifest = out_dir / "manifest.tsv"
    with manifest.open("w") as fout:
        fout.write("key\tvalue\n")
        fout.write(f"human_ref_dir\t{human_ref_dir}\n")
        fout.write(f"mouse_ref_dir\t{mouse_ref_dir}\n")
        fout.write(f"human_prefix\t{human_prefix}\n")
        fout.write(f"mouse_prefix\t{mouse_prefix}\n")
        fout.write(f"genome_fa\t{out_dir / 'genome.fa'}\n")
        fout.write(f"genome_fa_fai\t{out_dir / 'genome.fa.fai'}\n")
        fout.write(f"chrom_sizes\t{out_dir / 'chrom_sizes.tsv'}\n")
        fout.write(f"genes_gtf\t{out_dir / 'genes.gtf'}\n")
        fout.write(f"genes_bed\t{out_dir / 'genes.bed'}\n")


def main():
    args = parse_args()

    human_ref_dir = Path(args.human_ref_dir).resolve()
    mouse_ref_dir = Path(args.mouse_ref_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    expected = {
        "genome.fa": "genome FASTA",
        "chrom_sizes.tsv": "chrom sizes",
        "genes.gtf": "gene GTF",
        "genes.bed": "gene BED",
    }

    for ref_dir, label in ((human_ref_dir, args.human_label), (mouse_ref_dir, args.mouse_label)):
        for filename, desc in expected.items():
            require_file(ref_dir / filename, f"{label} {desc}")

    prefixed_dir = out_dir / "prefixed_inputs"
    prefixed_dir.mkdir(exist_ok=True)

    human_files = {
        "fasta": prefixed_dir / f"{args.human_label}.genome.fa",
        "chrom_sizes": prefixed_dir / f"{args.human_label}.chrom_sizes.tsv",
        "gtf": prefixed_dir / f"{args.human_label}.genes.gtf",
        "bed": prefixed_dir / f"{args.human_label}.genes.bed",
    }
    mouse_files = {
        "fasta": prefixed_dir / f"{args.mouse_label}.genome.fa",
        "chrom_sizes": prefixed_dir / f"{args.mouse_label}.chrom_sizes.tsv",
        "gtf": prefixed_dir / f"{args.mouse_label}.genes.gtf",
        "bed": prefixed_dir / f"{args.mouse_label}.genes.bed",
    }

    process_fasta(human_ref_dir / "genome.fa", human_files["fasta"], args.human_prefix)
    process_chrom_sizes(human_ref_dir / "chrom_sizes.tsv", human_files["chrom_sizes"], args.human_prefix)
    process_gtf(human_ref_dir / "genes.gtf", human_files["gtf"], args.human_prefix)
    process_bed(human_ref_dir / "genes.bed", human_files["bed"], args.human_prefix)

    process_fasta(mouse_ref_dir / "genome.fa", mouse_files["fasta"], args.mouse_prefix)
    process_chrom_sizes(mouse_ref_dir / "chrom_sizes.tsv", mouse_files["chrom_sizes"], args.mouse_prefix)
    process_gtf(mouse_ref_dir / "genes.gtf", mouse_files["gtf"], args.mouse_prefix)
    process_bed(mouse_ref_dir / "genes.bed", mouse_files["bed"], args.mouse_prefix)

    merged_genome = out_dir / "genome.fa"
    merged_chrom_sizes = out_dir / "chrom_sizes.tsv"
    merged_gtf = out_dir / "genes.gtf"
    merged_bed = out_dir / "genes.bed"

    concatenate_text([human_files["fasta"], mouse_files["fasta"]], merged_genome)
    concatenate_text([human_files["chrom_sizes"], mouse_files["chrom_sizes"]], merged_chrom_sizes)
    concatenate_text([human_files["gtf"], mouse_files["gtf"]], merged_gtf)
    concatenate_text([human_files["bed"], mouse_files["bed"]], merged_bed)

    run_samtools_faidx(merged_genome)
    write_manifest(out_dir, human_ref_dir, mouse_ref_dir, args.human_prefix, args.mouse_prefix)

    print("[OK] Mixed-species reference created")
    print(f"  genome.fa       : {merged_genome}")
    print(f"  genome.fa.fai   : {merged_genome}.fai")
    print(f"  chrom_sizes.tsv : {merged_chrom_sizes}")
    print(f"  genes.gtf       : {merged_gtf}")
    print(f"  genes.bed       : {merged_bed}")
    print(f"  manifest.tsv    : {out_dir / 'manifest.tsv'}")


if __name__ == "__main__":
    main()
