import os
import signal
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class Step2AlignmentTest(unittest.TestCase):
    def test_failed_streaming_alignment_drains_producer_and_keeps_status(self):
        for runner in ("run_all.sh", "run_all_mixed_species.sh"):
            source = (ROOT / runner).read_text()
            cleanup = source[source.index("cleanup_alignment_stream() {"):source.index("cleanup_large_intermediates() {")]
            start = source.index('if [[ "${UPSTREAM_ONLY}" -eq 0 && "${ALIGNMENT_STREAMING}" -eq 1 ]]; then')
            streaming = source[start:source.index("\nstep_start", start)]
            settings = r'''
UPSTREAM_ONLY=0
ALIGNMENT_STREAMING=1
GENOME_ALIGNER=bonefish
ALIGNMENT_FASTQ=input.pipe
ALIGN_LOG=alignment.log
log() { :; }
run_alignment_pipeline() {
  if [[ "$1" == input.pipe && "$FAIL_AFTER_EOF" == 1 ]]; then cat "$1" >/dev/null; fi
  echo 'synthetic mapper failure' >&2
  return 31
}
'''
            producer = r'''
dd if=/dev/zero of="$ALIGNMENT_FASTQ" bs=65536 count=4 2>/dev/null
exec 9>&-
wait "$ALIGNMENT_STREAM_PID"
'''
            for fail_after_eof in (False, True):
                with self.subTest(runner=runner, fail_after_eof=fail_after_eof), tempfile.TemporaryDirectory() as work:
                    script = (f"set -euo pipefail\nFAIL_AFTER_EOF={int(fail_after_eof)}\n" + cleanup +
                              settings + streaming + producer)
                    with subprocess.Popen(["bash", "-c", script], cwd=work, start_new_session=True,
                                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as process:
                        try:
                            output, errors = process.communicate(timeout=2)
                        except subprocess.TimeoutExpired:
                            os.killpg(process.pid, signal.SIGKILL)
                            process.communicate()
                            self.fail("alignment stream blocked after mapper failure")
                    self.assertEqual(process.returncode, 31, errors)
                    self.assertIn("synthetic mapper failure", output)
                    self.assertFalse((Path(work) / "input.pipe").exists())

    def test_bonefish_pin_only_applies_when_genome_alignment_uses_it(self):
        for runner in ("run_all.sh", "run_all_mixed_species.sh"):
            source = (ROOT / runner).read_text()
            end = source.index('require_file "${WHITELIST_3P}"')
            start = source.rindex("require_cmd python3\n", 0, end)
            preflight = source[start:end]
            settings = r'''
UPSTREAM_ONLY=$1
GENOME_ALIGNER=$2
FEATURE_BACKEND=legacy
GENE_TAGGING_BACKEND=rust
SKIP_ISOFORM=1
GENOME_MMI=''
BONEFISH_BIN=candidate
BONEFISH_VERSION=''
require_cmd() { command -v "$1" >/dev/null; }
samtools() { :; }
minimap2() { :; }
candidate() { printf '%s\n' 'bonefish 0.4.0'; }
log() { :; }
'''
            for upstream, mapper, expected in (("0", "bonefish", 1), ("1", "bonefish", 0), ("0", "minimap2", 0)):
                with self.subTest(runner=runner, upstream=upstream, mapper=mapper):
                    result = subprocess.run(["bash", "-c", "set -eu\n" + settings + preflight,
                                             "test", upstream, mapper], capture_output=True, text=True, timeout=5)
                    self.assertEqual(result.returncode, expected, result.stderr)
                    if expected:
                        self.assertIn("Bonefish 0.5.0", result.stderr)
                        self.assertIn("--genome-aligner minimap2", result.stderr)

    def test_both_workflows_use_same_filter_and_inline_bai_for_each_mapper(self):
        for runner in ("run_all.sh", "run_all_mixed_species.sh"):
            source = (ROOT / runner).read_text()
            pipeline = source[source.index("run_alignment_pipeline() {"):source.index("cleanup_alignment_stream() {")]
            for mapper, stream, index in (("bonefish", False, ""), ("bonefish", True, "splice idx.mmi"),
                                          ("minimap2", False, "splice idx.mmi"), ("minimap2", True, "")):
                with self.subTest(runner=runner, mapper=mapper, stream=stream), tempfile.TemporaryDirectory() as work:
                    result = subprocess.run(["bash", "-c", "set -euo pipefail\n" + pipeline + r'''
GENOME_ALIGNER=$1
GENOME_MMI=$2
BONEFISH_BIN=candidate
THREADS=8
GENE_FASTA='gene ref.fa'
JUNCTION_BED='genes bed'
ALIGNED_BAM='out.bam'
record_args() { printf '%s\n' "$@" > "$1.args"; }
candidate() {
  record_args bonefish "$@"
  if [[ "${@: -1}" == - ]]; then cat; else cat "${@: -1}"; fi
}
minimap2() { record_args minimap2 "$@"; cat "${@: -1}"; }
samtools() {
  case "$*" in
    'view -uh --no-PG -') record_args samtools-view "$@" ;;
    'view -uh -F 256 -') record_args samtools-filter "$@" ;;
    *) record_args samtools-sort "$@" ;;
  esac
  cat
}
if [[ "$3" == 1 ]]; then
  mkfifo input.fq
  printf 'reads\n' > input.fq &
else
  printf 'reads\n' > input.fq
fi
run_alignment_pipeline input.fq
''', "test", mapper, index, str(int(stream))], cwd=work, capture_output=True, text=True, timeout=5)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertEqual(result.stdout, "reads\n")
                    def argv(name):
                        path = Path(work) / f"{name}.args"
                        self.assertTrue(path.exists(), f"{name} was not called")
                        return path.read_text().splitlines()[1:]
                    prefix = ["-a", "-x", "cdna", "-u", "f"] if mapper == "bonefish" else ["-ax", "splice", "-uf"]
                    flags = ["--MD", "-I", "16G", "-t", "8", "-K", "10M", "--junc-bed", "genes bed"]
                    if mapper == "minimap2":
                        flags += ["--cap-kalloc", "100m"]
                    fastq_arg = "-" if mapper == "bonefish" and stream else "input.fq"
                    self.assertEqual(argv(mapper), prefix + flags + [index or "gene ref.fa", fastq_arg])
                    self.assertEqual(argv("samtools-view"), ["view", "-uh", "--no-PG", "-"])
                    self.assertEqual(argv("samtools-filter"), ["view", "-uh", "-F", "256", "-"])
                    self.assertEqual(argv("samtools-sort"), ["sort", "--no-PG", "-@", "8", "--write-index",
                                                            "-o", "out.bam##idx##out.bam.bai", "-"])
