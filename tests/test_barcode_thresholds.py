import csv
import importlib.util
import json
import math
import os
import random
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest


ROOT = Path(__file__).resolve().parents[1]
BINARY = Path(os.environ.get("FLORA_TEST_BINARY", ROOT / "target/debug/flora"))


def load_plotter():
    spec = importlib.util.spec_from_file_location("threshold_plotter", ROOT / "scripts/generate_knee_plots.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def signal_counts():
    return [max(1, int(100000 * math.exp(-.1 * math.log(i) -
                4 / (1 + math.exp(-(i - 3000) / 100))))) for i in range(1, 10001)]


class PlotThresholdTests(unittest.TestCase):
    def test_plot_whitelist_normalization_matches_rust(self):
        plotter = load_plotter()
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "whitelist.txt"
            path.write_text("acgaaa-1\nACGT-1\n")
            self.assertEqual(plotter.read_whitelist(str(path), True), {"TTTCGT", "ACGT"})
            self.assertEqual(plotter.read_whitelist(str(path), False), {"acgaaa", "ACGT"})

    def test_plot_uses_authoritative_applied_threshold_not_recomputed_legacy(self):
        plotter = load_plotter()
        self.assertTrue(hasattr(plotter, "load_threshold_report"))
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "thresholds.json"
            path.write_text(json.dumps({"schema_version": 1, "sides": {
                "3p": {"applied_threshold": 456, "comparison": ">=", "status": "ready"},
                "5p": {"applied_threshold": 789, "comparison": ">=", "status": "ready"}}}))
            self.assertEqual(plotter.load_threshold_report(path)["3p"]["applied_threshold"], 456)

    def test_plot_rejects_unavailable_or_nonfinite_applied_threshold(self):
        plotter = load_plotter()
        self.assertTrue(hasattr(plotter, "load_threshold_report"))
        for threshold in (None, float("nan"), -1):
            with self.subTest(threshold=threshold), tempfile.TemporaryDirectory() as temporary:
                path = Path(temporary) / "thresholds.json"
                path.write_text(json.dumps({"schema_version": 1, "sides": {
                    side: {"applied_threshold": threshold, "comparison": ">=", "status": "ready"}
                    for side in ("3p", "5p")}}))
                with self.assertRaises(ValueError):
                    plotter.load_threshold_report(path)

    def test_plot_checks_source_counts_and_uses_rank_one(self):
        from unittest.mock import patch
        plotter = load_plotter()
        with patch.object(plotter.plt, "loglog") as loglog, patch.object(plotter.plt, "savefig"):
            plotter.knee_plot([300, 200, 100], 150, "unused.png")
        self.assertEqual(list(loglog.call_args.args[0]), [1, 2, 3])

    def test_both_workflows_pass_threshold_options_and_authoritative_report(self):
        for name in ("run_all.sh", "run_all_mixed_species.sh"):
            runner = (ROOT / name).read_text()
            for fragment in ('WHITELIST_METHOD="legacy"', '--whitelist-method "${WHITELIST_METHOD}"',
                             '--threshold-report "${UPSTREAM_DIR}/barcode_thresholds.json"',
                             'printf "whitelist_method\\t%s\\n"'):
                self.assertIn(fragment, runner)
            invalidate = 'rm -f "${UPSTREAM_DIR}/barcode_thresholds.json" "${UPSTREAM_DIR}/barcode_thresholds.tsv"'
            self.assertIn(invalidate, runner)
            self.assertLess(runner.index(invalidate), runner.index('FLORA_TRIM_CACHE_DIR="${ONT_TEMP_DIR}" run_stage flora'))


@unittest.skipUnless(BINARY.is_file(), "set FLORA_TEST_BINARY to a built Flora binary")
class ThresholdPreviewTests(unittest.TestCase):
    def invoke(self, counts3, counts5=None, *options):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        for side, counts in (("3p", counts3), ("5p", counts3 if counts5 is None else counts5)):
            with (root / (side + ".tsv")).open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["barcode", "count"])
                writer.writerows(("BC%d" % i, count) for i, count in enumerate(counts))
            (root / (side + ".txt")).write_text("".join("BC%d\n" % i for i in range(len(counts))))
        result = subprocess.run([str(BINARY), "barcode-thresholds",
            "--counts-3p", str(root / "3p.tsv"), "--counts-5p", str(root / "5p.tsv"),
            "--full-whitelist-3p", str(root / "3p.txt"),
            "--full-whitelist-5p", str(root / "5p.txt"),
            "--no-revcomp-whitelist", "--out-dir", str(root / "preview"), *options],
            capture_output=True, text=True)
        return result, root / "preview"

    def test_preview_defaults_to_legacy_and_preserves_strict_greater_than(self):
        result, out = self.invoke([200, 10, 9])
        self.assertEqual(result.returncode, 0, result.stderr)
        data = json.loads((out / "barcode_thresholds.json").read_text())
        side = data["sides"]["3p"]
        self.assertEqual(data["method"], "legacy")
        self.assertEqual(side["comparison"], ">")
        self.assertEqual(side["applied_threshold"], side["legacy_threshold"])

    def test_steepest_matches_ngs_signal_and_keeps_count_ties(self):
        counts = signal_counts()
        result, out = self.invoke(counts, None, "--whitelist-method", "steepest")
        self.assertEqual(result.returncode, 0, result.stderr)
        data = json.loads((out / "barcode_thresholds.json").read_text())
        side = data["sides"]["3p"]
        rank = side["steepest"]["center_rank"]
        # Frozen from NGS call_cells.py, not an approximate expected location.
        self.assertEqual(rank, 3007)
        self.assertEqual(side["applied_threshold"], 5665)
        self.assertEqual([item["center_rank"] for item in side["steepest"]["windows"]], [3008, 3007, 3015, 3171])
        self.assertEqual(side["applied_threshold"], counts[rank - 1])
        self.assertEqual(side["comparison"], ">=")
        expected = {"BC%d" % i for i, count in enumerate(counts) if count >= side["applied_threshold"]}
        self.assertEqual(set((out / "whitelist_3p.csv").read_text().splitlines()), expected)

    def test_one_failed_side_writes_diagnostics_but_no_whitelists(self):
        result, out = self.invoke(signal_counts(), [1000] * 1000, "--whitelist-method", "steepest")
        self.assertNotEqual(result.returncode, 0)
        data = json.loads((out / "barcode_thresholds.json").read_text())
        self.assertEqual(data["sides"]["3p"]["status"], "ready")
        self.assertEqual(data["sides"]["5p"]["status"], "needs_manual_review")
        self.assertFalse((out / "whitelist_3p.csv").exists())
        self.assertFalse((out / "whitelist_5p.csv").exists())

    def test_invalid_window_fails_before_creating_outputs(self):
        result, out = self.invoke([1000] * 1000, None, "--whitelist-slope-window", "0")
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(out.exists())

    def test_invalid_spread_fails_before_outputs(self):
        for spread in ("NaN", "inf", "0", "1", "-0.1"):
            result, out = self.invoke([1000] * 1000, None, "--whitelist-slope-max-rank-spread", spread)
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(out.exists())

    def test_legacy_empty_and_asymmetric_sides_are_not_fatal(self):
        result, out = self.invoke([], [200, 10])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((out / "whitelist_3p.csv").read_text(), "")

    def test_exp_cells_does_not_change_steepest(self):
        results = []
        for exp in ("0", "1", "5000"):
            result, out = self.invoke(signal_counts(), None, "--whitelist-method", "steepest", "--exp-cells", exp)
            self.assertEqual(result.returncode, 0, result.stderr)
            results.append(json.loads((out / "barcode_thresholds.json").read_text())["sides"]["3p"]["applied_threshold"])
        self.assertEqual(results, [5665] * 3)

    def test_preview_does_not_overwrite_existing_directory(self):
        result, out = self.invoke([200, 10])
        self.assertEqual(result.returncode, 0, result.stderr)
        before = {path.name: path.read_bytes() for path in out.iterdir()}
        result = subprocess.run([str(BINARY), "barcode-thresholds", "--counts-3p", str(out.parent / "3p.tsv"),
            "--counts-5p", str(out.parent / "5p.tsv"), "--full-whitelist-3p", str(out.parent / "3p.txt"),
            "--full-whitelist-5p", str(out.parent / "5p.txt"), "--out-dir", str(out)], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(before, {path.name: path.read_bytes() for path in out.iterdir()})

    def test_orientation_and_off_whitelist_noise_use_the_same_selection_domain(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            # Reverse complement of ACGAAA is TTTCGT; the 1M-count noise is not eligible.
            for side in ("3p", "5p"):
                (root / (side + ".txt")).write_text("ACGAAA\n")
                (root / (side + ".tsv")).write_text("barcode\tcount\nTTTCGT\t200\nACGAAA\t3\nGGGGGG\t1000000\n")
            for reverse in (True, False):
                out = root / str(reverse)
                command = [str(BINARY), "barcode-thresholds", "--counts-3p", str(root / "3p.tsv"),
                    "--counts-5p", str(root / "5p.tsv"), "--full-whitelist-3p", str(root / "3p.txt"),
                    "--full-whitelist-5p", str(root / "5p.txt"), "--out-dir", str(out)]
                if not reverse:
                    command += ["--no-revcomp-whitelist"]
                result = subprocess.run(command, capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                report = json.loads((out / "barcode_thresholds.json").read_text())["sides"]["3p"]
                self.assertEqual(report["observed_theoretical_barcodes"], 1)
                self.assertEqual((out / "whitelist_3p.csv").read_text(), "TTTCGT\n" if reverse else "ACGAAA\n")
                self.assertEqual(report["legacy_threshold"], 10 if reverse else .15)


@unittest.skipUnless(BINARY.is_file(), "set FLORA_TEST_BINARY to a built Flora binary")
class ProductionThresholdTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        rng = random.Random(42)
        self.codes = []
        while len(self.codes) < 250:
            code = "".join(rng.choice("ACGT") for _ in range(20))
            if code not in self.codes and "CCTTCC" not in code and "GCTACC" not in code:
                self.codes.append(code)
        counts = [max(1, int(100 * math.exp(-.1 * math.log(i) - 4 / (1 + math.exp(-(i - 60) / 5))))) for i in range(1, 251)]
        self.bc5 = [code[:10] + "CCTTCC" + code[10:] for code in self.codes]
        rc = lambda code: code.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        self.bc3 = [rc(code)[:10] + "GCTACC" + rc(code)[10:] for code in self.codes]
        for side, barcodes in (("3p", self.bc3), ("5p", self.bc5)):
            (self.root / (side + ".txt")).write_text("\n".join(barcodes) + "\n")
        with (self.root / "input.fq").open("w") as handle:
            for i, count in enumerate(counts):
                sequence = self.bc5[i] + "CAGCA" + "C" * 10 + "CGTC" * 40 + "C" * 10 + "AGATC" + self.bc3[i]
                for j in range(count):
                    handle.write("@r%d_%d\n%s\n+\n%s\n" % (i, j, sequence, "I" * len(sequence)))

    def analyze(self, name, *options, binary=BINARY):
        return subprocess.run([str(binary), "analyze", str(self.root / "input.fq"),
            "--full-bc-whitelist-3p", str(self.root / "3p.txt"), "--full-bc-whitelist-5p", str(self.root / "5p.txt"),
            "--no-revcomp-whitelist", "--out_dir", str(self.root / name), "--threads", "1", *options],
            capture_output=True, text=True)

    def test_production_and_preview_generate_identical_whitelists_and_reports(self):
        options = ["--whitelist-method", "steepest", "--whitelist-slope-window", "5",
                   "--whitelist-slope-min-count", "0", "--whitelist-slope-max-rank-spread", ".3"]
        result = self.analyze("actual", *options)
        self.assertEqual(result.returncode, 0, result.stderr)
        actual = self.root / "actual"
        result = subprocess.run([str(BINARY), "barcode-thresholds", "--counts-3p", str(actual / "barcode_counts_3p.tsv"),
            "--counts-5p", str(actual / "barcode_counts_5p.tsv"), "--full-whitelist-3p", str(self.root / "3p.txt"),
            "--full-whitelist-5p", str(self.root / "5p.txt"), "--no-revcomp-whitelist",
            "--out-dir", str(self.root / "preview"), *options], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        for name in ("barcode_thresholds.json", "barcode_thresholds.tsv", "whitelist_3p.csv", "whitelist_5p.csv"):
            self.assertEqual((actual / name).read_bytes(), (self.root / "preview" / name).read_bytes(), name)

    def test_steepest_refuses_reused_production_outputs_without_modifying_them(self):
        result = self.analyze("actual")
        self.assertEqual(result.returncode, 0, result.stderr)
        before = {path.name: path.read_bytes() for path in (self.root / "actual").iterdir() if path.is_file()}
        result = self.analyze("actual", "--whitelist-method", "steepest")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("new/empty output directory", result.stderr)
        self.assertEqual(before, {path.name: path.read_bytes() for path in (self.root / "actual").iterdir() if path.is_file()})

    @unittest.skipUnless(hasattr(os, "mkfifo"), "requires POSIX FIFO")
    def test_steepest_claims_output_before_blocking_on_fastq(self):
        source = self.root / "input.fq"
        source.unlink()
        os.mkfifo(source)
        out = self.root / "reserved"
        command = [str(BINARY), "analyze", str(source), "--full-bc-whitelist-3p", str(self.root / "3p.txt"),
                   "--full-bc-whitelist-5p", str(self.root / "5p.txt"), "--no-revcomp-whitelist",
                   "--out_dir", str(out), "--whitelist-method", "steepest"]
        first = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        try:
            deadline = time.monotonic() + 3
            while time.monotonic() < deadline and not (out / ".barcode-threshold-reservation").exists():
                if first.poll() is not None:
                    self.fail(first.communicate()[1])
                time.sleep(.01)
            self.assertTrue((out / ".barcode-threshold-reservation").exists())
            second = subprocess.run(command, capture_output=True, text=True, timeout=3)
            self.assertNotEqual(second.returncode, 0)
            self.assertIn("output directory", second.stderr)
        finally:
            first.terminate()
            first.communicate(timeout=3)

    def test_default_legacy_matches_original_main_outputs(self):
        baseline = os.environ.get("FLORA_BASELINE_BINARY")
        if not baseline:
            self.skipTest("set FLORA_BASELINE_BINARY for original main comparison")
        result = self.analyze("baseline", binary=Path(baseline))
        self.assertEqual(result.returncode, 0, result.stderr)
        result = self.analyze("current")
        self.assertEqual(result.returncode, 0, result.stderr)
        for name in ("whitelist_3p.csv", "whitelist_5p.csv", "read_assigned_cell.csv", "barcode_to_cell.csv",
                     "pair_counts_kept.csv", "ReadIDs_UMI_BC_clean.csv", "cell_read_stats.csv", "assign_stats.tsv"):
            self.assertEqual((self.root / "baseline" / name).read_bytes(), (self.root / "current" / name).read_bytes(), name)

    def test_embedded_workflow_rejects_old_output_before_external_stages(self):
        fixture = ROOT / "tests/fixtures/workflows/single"
        existing = self.root / "existing"
        existing.mkdir()
        sentinel = existing / "old_report.txt"
        sentinel.write_text("do not modify")
        for command in ([], ["mixed"]):
            args = [str(BINARY), *command, "--python", sys.executable, "--fastq", str(fixture / "reads.fastq"),
                    "--barcode-list-10bp", str(fixture / "barcodes_10bp.txt"),
                    "--out-dir", str(existing), "--sample-id", "test", "--whitelist-method", "steepest"]
            args += ["--ref-dir", str(fixture / "ref")]
            result = subprocess.run(args, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("NEW output directory", result.stderr)
            self.assertEqual(sentinel.read_text(), "do not modify")
            self.assertEqual(list(existing.iterdir()), [sentinel])


if __name__ == "__main__":
    unittest.main()
