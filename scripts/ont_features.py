"""Flora feature assignment using the documented wf-single-cell v3.3.4 rules.

Independent implementation. Keep read-sized indexes on disk and process one
chromosome at a time; transcriptome secondary alignments are needed for ranking.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import gzip
import itertools
import json
from pathlib import Path
import re
import shlex
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import time


REFERENCE_VERSION = "wf-single-cell v3.3.4"
EXCLUDED_CLASSES = {"i", "y", "p", "s"}


@dataclass(frozen=True)
class Hit:
    gene: str
    transcript: str
    score: int
    query_coverage: float
    transcript_coverage: float
    mapq: int


def select_features(hits, genome_mapq, min_mapq=30, min_read_cov=0.4, min_tr_cov=0.4):
    if not hits:
        return "-", "-"
    ranked = sorted(hits, key=lambda h: (h.score, h.query_coverage, h.transcript_coverage), reverse=True)
    best = ranked[0]
    gene = best.gene if genome_mapq >= min_mapq else "-"
    if len(ranked) == 1:
        return (gene, best.transcript) if best.mapq > 0 else ("-", "-")
    second = ranked[1]
    score_or_query_diff = best.score != second.score or best.query_coverage != second.query_coverage
    sufficient_coverage = best.query_coverage >= min_read_cov and best.transcript_coverage >= min_tr_cov
    coverage_tiebreak = best.transcript_coverage != second.transcript_coverage and best.transcript_coverage > 0.8
    keep = sufficient_coverage if score_or_query_diff else coverage_tiebreak
    return gene, best.transcript if keep else "-"


def log(message):
    print(f"[ont_features] {message}", file=sys.stderr, flush=True)


def timed_run(command):
    started = time.monotonic()
    log(shlex.join([str(x) for x in command]))
    subprocess.run(command, check=True)
    log(f"Stage {command[0]} elapsed: {time.monotonic() - started:.2f}s")


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def attributes(text):
    return dict(re.findall(r'(\w+)\s+"([^"]*)"', text))


def split_annotation(gtf, root):
    paths = {}
    last_chrom = None
    handle = None
    try:
        with open_text(gtf) as source:
            for line in source:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.split("\t")
                if len(fields) != 9:
                    raise ValueError("Reference annotation must be a nine-column GTF")
                chrom = fields[0]
                if chrom not in paths:
                    paths[chrom] = root / f"chrom_{len(paths):05d}.gtf"
                if chrom != last_chrom:
                    if handle:
                        handle.close()
                    handle = paths[chrom].open("a")
                    last_chrom = chrom
                handle.write(line)
    finally:
        if handle:
            handle.close()
    return paths


def reference_transcripts(gtf):
    result = {}
    with open_text(gtf) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            attr = attributes(fields[8])
            tx = attr.get("transcript_id")
            if tx and fields[2].lower() in {"transcript", "mrna"}:
                result[tx] = attr.get("gene_name", attr.get("genename", "-")) or "-"
    return result


def load_candidates(tmap, gtf):
    reference = reference_transcripts(gtf)
    result = {}
    with open(tmap) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            tx = row["ref_id"]
            # ONT counts only assembled candidates linked to reference transcripts.
            if tx not in reference:
                continue
            assigned_tx = "-" if row["class_code"] in EXCLUDED_CLASSES else tx
            result[row["qry_id"]] = (reference[tx], assigned_tx)
    return result


def make_hit(record, candidates, lengths):
    pair = candidates.get(record.reference_name)
    if pair is None:
        return None
    read_length = record.infer_read_length() or 0
    if not read_length:
        return None
    return Hit(
        pair[0], pair[1], record.get_tag("AS") if record.has_tag("AS") else 0,
        record.query_alignment_length / read_length,
        (record.reference_end - record.reference_start) / lengths[record.reference_name],
        record.mapping_quality,
    )


def assign_chromosome(db, chrom, transcript_bam, candidates, min_mapq):
    import pysam

    with pysam.AlignmentFile(transcript_bam, "rb", check_sq=False) as bam:
        lengths = dict(zip(bam.references, bam.lengths))
        pending = []
        for name, records in itertools.groupby(bam, key=lambda a: a.query_name):
            hits = [make_hit(a, candidates, lengths) for a in records
                    if not a.is_unmapped and not a.is_supplementary]
            hits = [h for h in hits if h is not None]
            row = db.execute("SELECT mapq FROM genome WHERE read_id=?", (name,)).fetchone()
            if row is None:
                continue
            gene, tx = select_features(hits, row[0], min_mapq)
            reason = "candidate" if hits else "no_reference_candidate"
            if hits and gene == "-" and tx == "-":
                reason = "ambiguous_or_zero_transcript_mapq"
            pending.append((chrom, name, gene, tx, reason))
            if len(pending) >= 10000:
                db.executemany("INSERT OR REPLACE INTO features VALUES (?,?,?,?,?)", pending)
                db.commit()
                pending.clear()
        db.executemany("INSERT OR REPLACE INTO features VALUES (?,?,?,?,?)", pending)
        db.commit()


def assemble(bam, chrom, gtf, genome, directory, threads, options, label):
    started = time.monotonic()
    view = subprocess.Popen(["samtools", "view", "-h", str(bam), chrom], stdout=subprocess.PIPE)
    command = ["stringtie", "-L", *shlex.split(options), "-p", str(threads),
               "-G", str(gtf), "-l", label, "-o", str(directory / "stringtie.gff"), "-"]
    log(shlex.join(command))
    try:
        result = subprocess.run(command, stdin=view.stdout)
    finally:
        view.stdout.close()
        view_code = view.wait()
    if view_code:
        raise subprocess.CalledProcessError(view_code, view.args)
    result.check_returncode()
    log(f"Stage stringtie ({chrom}) elapsed: {time.monotonic() - started:.2f}s")
    # An empty assembly is a valid outcome for sparse chromosomes.
    with (directory / "stringtie.gff").open() as handle:
        has_exons = any(not line.startswith("#") and "\texon\t" in line for line in handle)
    if not has_exons:
        log(f"No assembled transcripts on {chrom}")
        return
    timed_run(["gffread", "-g", str(genome), "-w", str(directory / "transcriptome.fa"),
               str(directory / "stringtie.gff")])
    timed_run(["gffcompare", "-o", str(directory / "gffcompare"), "-r", str(gtf),
               str(directory / "stringtie.gff")])


def remap_chromosome(db, bam_path, chrom, fasta, output, threads):
    import pysam

    started = time.monotonic()
    mm = ["minimap2", "-ax", "map-ont", "--cap-kalloc", "100m", "--cap-sw-mem", "50m",
          "--end-bonus", "10", "-p", "0.9", "-N", "3", "-t", str(max(threads - 4, 4)), str(fasta), "-"]
    log(shlex.join(mm))
    processes = []
    try:
        aligner = subprocess.Popen(mm, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        processes.append(aligner)
        viewer = subprocess.Popen(["samtools", "view", "-h", "-@", "1", "-b", "-F", "2052", "-"],
                                  stdin=aligner.stdout, stdout=subprocess.PIPE)
        processes.append(viewer)
        aligner.stdout.close()
        sorter = subprocess.Popen(["samtools", "sort", "-n", "-@", "3", "--no-PG", "-o", str(output), "-"],
                                  stdin=viewer.stdout)
        processes.append(sorter)
        viewer.stdout.close()
        pending = []
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for record in bam.fetch(chrom):
                if record.is_unmapped or record.is_secondary or record.is_supplementary:
                    continue
                name = record.query_name
                sequence = record.get_forward_sequence()
                quals = record.get_forward_qualities()
                if sequence is None:
                    raise ValueError(f"Missing sequence for {name}")
                quality = "".join(chr(q + 33) for q in quals) if quals is not None else "B" * len(sequence)
                aligner.stdin.write(f"@{name}\n{sequence}\n+\n{quality}\n".encode())
                pending.append((name, record.mapping_quality))
                if len(pending) >= 10000:
                    db.executemany("INSERT OR REPLACE INTO genome VALUES (?,?)", pending)
                    pending.clear()
        db.executemany("INSERT OR REPLACE INTO genome VALUES (?,?)", pending)
        db.commit()
        aligner.stdin.close()
        codes = [p.wait() for p in processes]
        for p, code in zip(processes, codes):
            if code:
                raise subprocess.CalledProcessError(code, p.args)
    except BaseException:
        for p in processes:
            if p.poll() is None:
                p.terminate()
        for p in processes:
            p.wait()
        raise
    log(f"Stage transcript_remapping ({chrom}) elapsed: {time.monotonic() - started:.2f}s")


def write_tagged_assignments(db, tagged_bam, output_bam, gene_tsv, transcript_tsv, min_mapq, threads):
    import pysam
    from collections import Counter

    stats = Counter()
    db.execute("CREATE TABLE IF NOT EXISTS counted_reads(read_id TEXT PRIMARY KEY) WITHOUT ROWID")
    with pysam.AlignmentFile(tagged_bam, "rb") as source, \
            pysam.AlignmentFile(output_bam, "wb", template=source, threads=max(1, min(threads, 8))) as dest, \
            open(gene_tsv, "w") as genes, open(transcript_tsv, "w") as transcripts:
        for record in source:
            # One primary genomic record per read avoids cross-chromosome SA duplication.
            if record.is_secondary or record.is_supplementary or record.is_unmapped:
                continue
            try:
                db.execute("INSERT INTO counted_reads VALUES (?)", (record.query_name,))
            except sqlite3.IntegrityError as error:
                raise ValueError(f"Duplicate primary read name: {record.query_name}; input FASTQ read IDs must be unique") from error
            row = db.execute("SELECT gene,transcript,reason FROM features WHERE chrom=? AND read_id=?",
                             (record.reference_name, record.query_name)).fetchone()
            gene, tx, reason = row if row else ("-", "-", "no_reference_candidate")
            status = "Assigned" if gene != "-" else (
                "Unassigned_mapq" if record.mapping_quality < min_mapq else "Unassigned_no_features")
            if gene == "-" and reason == "ambiguous_or_zero_transcript_mapq" and record.mapping_quality >= min_mapq:
                status = "Unassigned_ambiguous"
            tx_status = "Assigned" if tx != "-" else (
                "Unassigned_ambiguous" if reason != "no_reference_candidate" else "Unassigned_no_features")
            rid = record.query_name
            gene_out, tx_out = (gene if gene != "-" else "NA"), (tx if tx != "-" else "NA")
            genes.write(f"{rid}\t{status}\t{record.mapping_quality}\t{gene_out}\n")
            transcripts.write(f"{rid}\t{tx_status}\t{record.mapping_quality}\t{gene_out}\t{tx_out}\n")
            record.set_tag("GN", gene_out, value_type="Z")
            record.set_tag("TX", tx_out, value_type="Z")
            dest.write(record)
            stats[f"gene:{status}"] += 1
            stats[f"transcript:{tx_status}"] += 1
            stats[f"reason:{reason}"] += 1
    pysam.index("-@", str(max(1, min(threads, 8))), str(output_bam))
    return dict(stats)


def run(args):
    import pysam

    for name in ("minimap2", "samtools", "stringtie", "gffread", "gffcompare"):
        if shutil.which(name) is None:
            raise RuntimeError(f"Missing required tool: {name}; update the Flora conda environment")
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    versions = {"behavior_reference": REFERENCE_VERSION}
    for name, flag in [("minimap2", "--version"), ("samtools", "--version"),
                       ("stringtie", "--version"), ("gffread", "--version"), ("gffcompare", "--version")]:
        p = subprocess.run([name, flag], capture_output=True, text=True, errors="replace")
        versions[name] = ((p.stdout + p.stderr).strip().splitlines() or ["unknown"])[0]
    (out / "tool_versions.json").write_text(json.dumps(versions, indent=2) + "\n")
    log(f"Behavior reference: {REFERENCE_VERSION}; {json.dumps(versions)}")
    with tempfile.TemporaryDirectory(prefix="ont-work-", dir=out) as temp:
        root = Path(temp)
        genome = args.genome.resolve()
        if not Path(str(genome) + ".fai").exists():
            # Each invocation gets its own link/index, including reruns with a new reference.
            linked_genome = root / "reference.fa"
            linked_genome.symlink_to(genome)
            timed_run(["samtools", "faidx", str(linked_genome)])
            genome = linked_genome
        saved_run = out / "intermediates" / root.name
        if args.save_intermediate:
            saved_run.mkdir(parents=True)
        db = sqlite3.connect(root / "assignments.sqlite")
        db.execute("PRAGMA cache_size=-32768")
        db.execute("CREATE TABLE features(chrom TEXT,read_id TEXT,gene TEXT,transcript TEXT,reason TEXT,PRIMARY KEY(chrom,read_id)) WITHOUT ROWID")
        db.execute("CREATE TABLE genome(read_id TEXT PRIMARY KEY,mapq INTEGER) WITHOUT ROWID")
        try:
            annotations = split_annotation(args.gtf, root)
            if not reference_transcripts(args.gtf):
                raise ValueError("ONT feature assignment requires reference transcript/mRNA records with transcript_id and gene_name in the GTF")
            with pysam.AlignmentFile(args.bam, "rb") as bam:
                mapped = {s.contig for s in bam.get_index_statistics() if s.mapped}
            transcriptome_fa = out / f"{args.sample}.transcriptome.fa.gz"
            transcriptome_gtf = out / f"{args.sample}.transcriptome.gtf.gz"
            with gzip.open(transcriptome_fa, "wt") as fa_out, gzip.open(transcriptome_gtf, "wt") as gtf_out:
                for number, (chrom, gtf) in enumerate(annotations.items()):
                    if chrom not in mapped:
                        continue
                    directory = root / f"chrom_{number:05d}"
                    directory.mkdir()
                    assemble(args.bam, chrom, gtf, genome, directory, args.threads, args.stringtie_opts, f"flora_chr{number}")
                    fasta = directory / "transcriptome.fa"
                    annotated = directory / "gffcompare.annotated.gtf"
                    if fasta.exists():
                        with fasta.open() as handle:
                            shutil.copyfileobj(handle, fa_out)
                    if annotated.exists():
                        with annotated.open() as handle:
                            shutil.copyfileobj(handle, gtf_out)
                    if fasta.exists() and fasta.stat().st_size:
                        db.execute("DELETE FROM genome")
                        remap_chromosome(db, args.bam, chrom, fasta, directory / "transcript.bam", args.threads)
                        tmap = directory / "gffcompare.stringtie.gff.tmap"
                        candidates = load_candidates(tmap, gtf)
                        assign_chromosome(db, chrom, directory / "transcript.bam", candidates, args.mapq)
                    if args.save_intermediate:
                        shutil.copyfile(gtf, directory / "reference.gtf")
                        shutil.move(str(directory), str(saved_run / directory.name))
                    else:
                        shutil.rmtree(directory)
            stats = write_tagged_assignments(db, args.tagged_bam, args.output_bam, args.gene_assigns,
                                             args.transcript_assigns, args.mapq, args.threads)
            (out / "assignment_summary.json").write_text(json.dumps(stats, indent=2) + "\n")
        finally:
            db.close()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("bam", "tagged-bam", "genome", "gtf", "output-dir", "output-bam", "gene-assigns", "transcript-assigns"):
        p.add_argument(f"--{name}", required=True, type=Path)
    p.add_argument("--sample", required=True)
    p.add_argument("--threads", type=int, default=16)
    p.add_argument("--mapq", type=int, default=30)
    p.add_argument("--stringtie-opts", default="-c 2")
    p.add_argument("--save-intermediate", action="store_true")
    args = p.parse_args()
    if args.threads < 1 or not 0 <= args.mapq <= 255:
        p.error("threads must be positive and MAPQ must be between 0 and 255")
    run(args)


if __name__ == "__main__":
    main()
