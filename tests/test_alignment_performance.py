import gzip
import importlib.util
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


def shell_function(path, name):
    text = path.read_text()
    start = text.index(name + "() {")
    end = text.index("\n}\n", start) + 3
    return text[start:end]


class AlignmentPerformanceTests(unittest.TestCase):
    def pipeline(self, script, compressed=False, invalid=False):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fastq = root / ("reads.fq.gz" if compressed else "reads.fq")
            payload = b"@r1\nACGT\n+\nIIII\n"
            fastq.write_bytes(b"not gzip" if invalid else gzip.compress(payload) if compressed else payload)
            output, capture = root / "out.bam", root / "args"
            # Protocol fixtures validate plumbing; they are not alignment algorithms.
            stubs = '''
minimap2() { printf '%s\\n' "$@" > "$CAPTURE"; local last="${!#}"; if [[ "$last" == - ]]; then cat; else cat "$last"; fi; }
samtools() { local sub="$1"; shift; if [[ "$sub" == sort ]]; then while [[ $# -gt 0 ]]; do if [[ "$1" == -o ]]; then cat > "$2"; return; fi; shift; done; else cat; fi; }
log() { :; }
'''
            command = "\n".join(["set -euo pipefail", stubs,
                shell_function(ROOT / script, "alignment_fastq_stream"),
                shell_function(ROOT / script, "run_alignment_pipeline"),
                'THREADS=4; MINIMAP2_BATCH=500M; MINIMAP2_IO_THREADS=2; ALIGNMENT_GZIP_INPUT=stream',
                'JUNCTION_BED=genes.bed; GENE_FASTA=genome.fa',
                'ALIGNED_BAM="$1"; CAPTURE="$2"; run_alignment_pipeline "$3"'])
            result = subprocess.run(["bash", "-c", command, "test", str(output), str(capture), str(fastq)], capture_output=True)
            return result.returncode, output.read_bytes() if output.exists() else None, capture.read_text().splitlines() if capture.exists() else []

    def test_plain_and_gzip_are_identical_uncompressed_streams(self):
        for script in ("run_all.sh", "run_all_mixed_species.sh"):
            for compressed in (False, True):
                with self.subTest(script=script, compressed=compressed):
                    code, data, argv = self.pipeline(script, compressed)
                    self.assertEqual(code, 0)
                    self.assertEqual(data, b"@r1\nACGT\n+\nIIII\n")
                    self.assertEqual(argv[argv.index("-K") + 1], "500M")
                    self.assertIn("-2", argv)
                    self.assertEqual(argv[-1], "-")

    def test_bad_gzip_fails_the_entire_pipeline(self):
        for script in ("run_all.sh", "run_all_mixed_species.sh"):
            code, _, _ = self.pipeline(script, compressed=True, invalid=True)
            self.assertNotEqual(code, 0)

    def test_both_workflows_expose_reversible_knobs(self):
        for script in ("run_all.sh", "run_all_mixed_species.sh"):
            text = (ROOT / script).read_text()
            for option in ("--minimap2-batch", "--minimap2-io-threads", "--alignment-gzip-input", "--stringtie-bin"):
                self.assertIn(option, text)
            self.assertIn('MINIMAP2_BATCH="500M"', text)
            self.assertIn('ALIGNMENT_GZIP_INPUT="stream"', text)

    def test_explicit_stringtie_path_flows_to_assembly(self):
        spec = importlib.util.spec_from_file_location("ont_perf_tools", ROOT / "scripts/ont_features.py")
        mod = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = mod
        spec.loader.exec_module(mod)
        command = mod.assembly_command("/isolated/bin/stringtie", "genes.gtf", "out.gtf", 4, "-c 2 -v", "flora_chr1")
        self.assertEqual(command[0], "/isolated/bin/stringtie")
        self.assertEqual(mod.stringtie_executable(Path("./stringtie")), str((Path.cwd() / "stringtie").resolve()))
        self.assertEqual(mod.stringtie_executable(None), "stringtie")
        self.assertEqual(command[1:5], ["-L", "-c", "2", "-v"])
        self.assertEqual(command[-1], "-")
        self.assertNotIn("-j", command)

    def test_invalid_controls_fail_before_any_work(self):
        for script in ("run_all.sh", "run_all_mixed_species.sh"):
            for option, value in (("--minimap2-batch", "0"), ("--minimap2-batch", "1M;false"),
                                  ("--minimap2-io-threads", "4"), ("--alignment-gzip-input", "unknown")):
                done = subprocess.run(["bash", str(ROOT / script), option, value], capture_output=True, text=True)
                self.assertNotEqual(done.returncode, 0)
                self.assertIn(option, done.stderr)

    def test_stringtie_preflight_checks_selected_binary_not_path_stock(self):
        for script in ("run_all.sh", "run_all_mixed_species.sh"):
            with tempfile.TemporaryDirectory() as tmp:
                selected = Path(tmp) / "stringtie"
                selected.write_text("#!/bin/sh\nprintf '2.2.2\\n'\n")
                selected.chmod(0o755)
                code = "\n".join(["set -euo pipefail", "stringtie() { echo 2.2.3; }",
                    'require_cmd() { command -v "$1" >/dev/null; }',
                    shell_function(ROOT / script, "require_stringtie"),
                    'STRINGTIE_BIN="$1"; require_stringtie'])
                done = subprocess.run(["bash", "-c", code, "check", str(selected)], capture_output=True)
                self.assertEqual(done.returncode, 0)
                done = subprocess.run(["bash", "-c", code, "check", ""], capture_output=True)
                self.assertNotEqual(done.returncode, 0)

    @unittest.skipUnless(shutil.which("minimap2") and shutil.which("samtools"), "requires actual minimap2 and samtools")
    def test_real_minimap2_spliced_alignment_matches_previous_pipeline(self):
        import pysam
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rng = random.Random(42)
            dna = "".join(rng.choices("ACGT", k=10000))
            dna = dna[:1300] + "GT" + dna[1302:1798] + "AG" + dna[1800:]
            genome, bed, fastq = root / "genome.fa", root / "genes.bed", root / "reads.fq.gz"
            genome.write_text(">chr1\n" + dna + "\n")
            bed.write_text("chr1\t1000\t2200\tT\t0\t+\t1000\t2200\t0\t2\t300,400\t0,800\n")
            with gzip.open(fastq, "wt") as handle:
                for i in range(60):
                    seq = dna[1000 + i % 20:1300] + dna[1800:2200 - i % 10]
                    if i % 3 == 0:
                        seq = seq[:50] + ("A" if seq[50] != "A" else "T") + seq[51:]
                    handle.write(f"@read{i}\n{seq}\n+\n{'I' * len(seq)}\n")
            previous = root / "previous.bam"
            baseline = '''set -euo pipefail
minimap2 -ax splice -uf --MD -I 16G -t 2 -K 10M --junc-bed "$1" --cap-kalloc 100m "$2" "$3" |
samtools view -uh --no-PG - | samtools view -uh -F 256 - |
samtools sort --no-PG -@ 2 -o "$4" -
'''
            subprocess.run(["bash", "-c", baseline, "baseline", str(bed), str(genome), str(fastq), str(previous)], check=True, capture_output=True)
            with pysam.AlignmentFile(previous) as handle:
                expected = sorted(r.to_string() for r in handle)
            for script in ("run_all.sh", "run_all_mixed_species.sh"):
                for batch, io_threads, supplier in (("500M", 2, "stream"), ("1K", 1, "stream"), ("10M", 1, "direct")):
                    output = root / "new.bam"
                    code = "\n".join(["set -euo pipefail", "log() { :; }",
                        shell_function(ROOT / script, "alignment_fastq_stream"),
                        shell_function(ROOT / script, "run_alignment_pipeline"),
                        'THREADS=2; JUNCTION_BED="$1"; GENE_FASTA="$2"; ALIGNED_BAM="$4"',
                        f'MINIMAP2_BATCH={batch}; MINIMAP2_IO_THREADS={io_threads}; ALIGNMENT_GZIP_INPUT={supplier}',
                        'run_alignment_pipeline "$3"'])
                    subprocess.run(["bash", "-c", code, "new", str(bed), str(genome), str(fastq), str(output)], check=True, capture_output=True)
                    with pysam.AlignmentFile(output) as handle:
                        actual = sorted(r.to_string() for r in handle)
                    self.assertEqual(actual, expected)


if __name__ == "__main__":
    unittest.main()
