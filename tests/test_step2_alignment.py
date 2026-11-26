import json
import os
import signal
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def workflow_preflight(source):
    end = source.index('require_file "${WHITELIST_3P}"')
    start = source.rindex("require_cmd python3\n", 0, end)
    first = source.index("resolve_workflow_tool() {")
    helper = source[first:source.index("\n}\n", first) + 3]
    return helper + "\n" + source[start:end]


def version_executable(path, version):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"#!/bin/sh\nprintf '%s\\n' '{version}'\n")
    path.chmod(0o755)


class Step2AlignmentTest(unittest.TestCase):
    def test_workflow_routes_child_tempfiles_to_selected_scratch(self):
        for name in ("run_all.sh", "run_all_mixed_species.sh"):
            with tempfile.TemporaryDirectory() as work:
                root = Path(work).resolve()
                runner = root / name
                runner.write_bytes((ROOT / name).read_bytes())
                (root / "reference").mkdir()
                (root / "barcodes.txt").write_text("ACGTACGTAC\n")
                (root / "input.fastq").touch()
                (root / "scripts").mkdir()
                child = root / "scripts/generate_26bp_whitelists.py"
                child.write_text(
                    "import json, os, tempfile\n"
                    "with tempfile.NamedTemporaryFile() as handle:\n"
                    "    print('child scratch ' + json.dumps([os.environ.get('TMPDIR'), handle.name]))\n"
                    "raise SystemExit(73)\n"
                )
                original_tmp = root / "parent scratch"
                original_tmp.mkdir()
                for override in (None, "explicit scratch"):
                    with self.subTest(runner=name, override=override):
                        output = root / ("output default" if override is None else "output explicit")
                        args = ["--fastq", "input.fastq", "--ref-dir", "reference",
                                "--barcode-list-10bp", "barcodes.txt", "--out-dir", output.name]
                        expected = output / ".flora-scratch"
                        if override is not None:
                            args += ["--ont-temp-dir", override]
                            expected = root / override
                        result = subprocess.run(
                            ["/bin/bash", str(runner), *args], cwd=root, timeout=5,
                            env={**os.environ, "TMPDIR": str(original_tmp)},
                            capture_output=True, text=True,
                        )
                        self.assertEqual(result.returncode, 73, result.stdout + result.stderr)
                        line = next(line for line in result.stdout.splitlines() if line.startswith("child scratch "))
                        env_path, file_path = json.loads(line.removeprefix("child scratch "))
                        self.assertIsNotNone(env_path)
                        child_tmp = Path(env_path).resolve()
                        self.assertIn(expected, (child_tmp, *child_tmp.parents))
                        self.assertEqual(Path(file_path).resolve().parent, child_tmp)

    def test_managed_tools_match_install_paths_and_preserve_explicit_choices(self):
        settings = r'''
UPSTREAM_ONLY=0
GENOME_ALIGNER=bonefish
FEATURE_BACKEND=ont
GENE_TAGGING_BACKEND=rust
SKIP_ISOFORM=0
GENOME_MMI=''
BONEFISH_BIN=$1
BONEFISH_VERSION=''
STRINGTIE_BIN=$2
DOWNSTREAM_DIR=unused
require_cmd() { command -v "$1" >/dev/null; }
python_asset() { :; }
python3() { "$TEST_PYTHON" "$@"; }
samtools() { :; }
minimap2() { :; }
gffread() { :; }
gffcompare() { :; }
log() { printf '%s\n' "$*"; }
'''
        for runner in ("run_all.sh", "run_all_mixed_species.sh"):
            preflight = workflow_preflight((ROOT / runner).read_text())
            with tempfile.TemporaryDirectory() as work:
                work = Path(work).resolve()
                stock = work / "stock bin"
                data_home = work / "data home"
                home = work / "home"
                override = work / "managed override"
                home_override = home / "other install" / "bin"
                for bin_dir in (stock, home / ".local/share/flora/tools/bin", data_home / "flora/tools/bin",
                                override / "bin", home_override):
                    version_executable(bin_dir / "bonefish", "bonefish 0.5.0")
                    stringtie = "stringtie" if bin_dir == stock else "stringtie-2.2.2-same-start"
                    version_executable(bin_dir / stringtie, "2.2.2")

                env = {**os.environ, "HOME": str(home), "PATH": str(stock), "TEST_PYTHON": sys.executable}
                for name in ("FLORA_TOOLS_DIR", "XDG_DATA_HOME", "STRINGTIE_INPUT_LOOKUP"):
                    env.pop(name, None)
                xdg = {"XDG_DATA_HOME": data_home.name}
                custom = {**xdg, "FLORA_TOOLS_DIR": "managed override"}
                cases = (
                    ({}, ("", ""), home / ".local/share/flora/tools/bin", "indexed"),
                    (xdg, ("", ""), data_home / "flora/tools/bin", "indexed"),
                    ({**custom, "STRINGTIE_INPUT_LOOKUP": "linear"}, ("", ""), override / "bin", "linear"),
                    ({"FLORA_TOOLS_DIR": "~/other install"}, ("", ""), home_override, "indexed"),
                    (custom, ("bonefish", "stringtie"), stock, "indexed"),
                    ({"FLORA_TOOLS_DIR": "missing"}, ("", ""), stock, "indexed"),
                )
                output = '\nprintf "tools:%s:%s:%s:%s\\n" "$BONEFISH_BIN" "$STRINGTIE_BIN" "$STRINGTIE_INPUT_LOOKUP" "$PATH"\n'
                for settings_env, requested, selected, mode in cases:
                    with self.subTest(runner=runner, environment=settings_env, requested=requested):
                        result = subprocess.run(["/bin/bash", "-c", "set -eu\n" + settings + preflight + output,
                                                 "test", *requested], cwd=work, env={**env, **settings_env},
                                                capture_output=True, text=True, timeout=5)
                        self.assertEqual(result.returncode, 0, result.stderr)
                        bin_dir = selected
                        if requested[0]:
                            bin_dir = override / "bin"
                        expected_path = str(stock) if bin_dir == stock else f"{bin_dir}:{stock}"
                        stringtie = "stringtie" if selected == stock else "stringtie-2.2.2-same-start"
                        self.assertIn(f"tools:{selected / 'bonefish'}:{selected / stringtie}:{mode}:{expected_path}",
                                      result.stdout)
                        self.assertIn(f"StringTie: {selected / stringtie}", result.stdout)

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
            preflight = workflow_preflight(source)
            settings = r'''
UPSTREAM_ONLY=$1
GENOME_ALIGNER=$2
FEATURE_BACKEND=legacy
GENE_TAGGING_BACKEND=rust
SKIP_ISOFORM=1
GENOME_MMI=''
BONEFISH_BIN=$3
BONEFISH_VERSION=''
require_cmd() { command -v "$1" >/dev/null; }
samtools() { :; }
minimap2() { :; }
log() { :; }
'''
            with tempfile.TemporaryDirectory() as work:
                candidate = Path(work) / "candidate"
                version_executable(candidate, "bonefish 0.4.0")
                for upstream, mapper, expected in (("0", "bonefish", 1), ("1", "bonefish", 0), ("0", "minimap2", 0)):
                    with self.subTest(runner=runner, upstream=upstream, mapper=mapper):
                        result = subprocess.run(["bash", "-c", "set -eu\n" + settings + preflight,
                                                 "test", upstream, mapper, str(candidate)],
                                                capture_output=True, text=True, timeout=5)
                        self.assertEqual(result.returncode, expected, result.stderr)
                        if expected:
                            self.assertIn("Bonefish 0.5.0", result.stderr)
                            self.assertIn("bonefish 0.4.0", result.stderr)
                            self.assertIn("flora setup-tools", result.stderr)
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
