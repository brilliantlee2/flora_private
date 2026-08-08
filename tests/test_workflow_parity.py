import importlib.util
import json
import os
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WORKFLOW_FIXTURES = ROOT / "tests" / "fixtures" / "workflows"
COMPARATOR_PATH = ROOT / "tests" / "compare_workflow_outputs.py"
CLASSIFICATIONS = {
    "byte_exact",
    "parsed_exact",
    "numeric_tolerance",
    "canonicalized_html_log",
    "intentionally_absent",
}
REQUIRED_SCENARIOS = {
    "single": {
        "full",
        "light",
        "skip_glycine",
        "skip_isoform",
        "upstream_only",
        "stale_output",
        "malformed_input",
        "forced_failure",
    },
    "mixed": {
        "full",
        "light",
        "skip_glycine",
        "skip_isoform",
        "upstream_only",
        "stale_output",
        "malformed_input",
        "forced_failure",
    },
}
REQUIRED_INPUTS = {
    "reads.fastq",
    "barcodes_10bp.txt",
    "ref/genome.fa",
    "ref/genes.gtf",
    "ref/isoforms.gtf",
    "ref/genes.bed",
    "ref/chrom_sizes.tsv",
}


def load_comparator():
    spec = importlib.util.spec_from_file_location("workflow_comparator", COMPARATOR_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("unable to load workflow comparator")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class WorkflowOracleFixtureTests(unittest.TestCase):
    def test_fixture_inputs_are_small_complete_and_synthetic(self):
        for workflow in REQUIRED_SCENARIOS:
            fixture = WORKFLOW_FIXTURES / workflow
            with self.subTest(workflow=workflow):
                paths = {
                    path.relative_to(fixture).as_posix()
                    for path in fixture.rglob("*")
                    if path.is_file()
                }
                self.assertTrue(REQUIRED_INPUTS.issubset(paths))
                for relative in REQUIRED_INPUTS:
                    path = fixture / relative
                    self.assertLess(path.stat().st_size, 100_000)
                    self.assertNotIn("/Users/", path.read_text(encoding="utf-8"))

    def test_fixture_tree_is_modest_and_contains_no_private_paths(self):
        import gzip

        files = [
            path
            for path in WORKFLOW_FIXTURES.rglob("*")
            if path.is_file() and not path.is_symlink()
        ]
        self.assertLess(sum(path.stat().st_size for path in files), 20_000_000)
        self.assertTrue(all(path.stat().st_size < 5_000_000 for path in files))
        needles = ("/Users/", "/private/", "/tmp/")
        text_suffixes = {
            ".bed",
            ".csv",
            ".fa",
            ".fastq",
            ".gtf",
            ".html",
            ".json",
            ".log",
            ".sh",
            ".tsv",
            ".txt",
        }
        for path in files:
            if path.suffix == ".gz":
                try:
                    with gzip.open(path, "rt", encoding="utf-8") as handle:
                        text = handle.read()
                except (gzip.BadGzipFile, UnicodeDecodeError):
                    continue
            elif path.suffix in text_suffixes:
                text = path.read_text(encoding="utf-8")
            else:
                continue
            for needle in needles:
                self.assertNotIn(needle, text, str(path))
        for path in WORKFLOW_FIXTURES.rglob("*"):
            if path.is_symlink():
                target = os.readlink(path)
                self.assertFalse(Path(target).is_absolute(), str(path))

    def test_manifests_are_complete_frozen_and_classified(self):
        for workflow, required_scenarios in REQUIRED_SCENARIOS.items():
            manifest_path = ROOT / "tests" / f"artifact_manifest_{workflow}.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            with self.subTest(workflow=workflow):
                self.assertEqual(1, manifest["schema_version"])
                self.assertEqual(workflow, manifest["workflow"])
                self.assertEqual(required_scenarios, set(manifest["scenarios"]))
                self.assertEqual(
                    {"absolute", "relative"}, set(manifest["numeric_tolerance"])
                )
                self.assertGreater(len(manifest["artifacts"]), 20)
                for artifact in manifest["artifacts"]:
                    self.assertIn(artifact["classification"], CLASSIFICATIONS)
                    self.assertNotIn("..", Path(artifact["path"]).parts)
                    if artifact["classification"] != "intentionally_absent":
                        self.assertRegex(artifact["sha256"], r"^[0-9a-f]{64}$")

                present = {
                    (artifact["scenario"], artifact["path"])
                    for artifact in manifest["artifacts"]
                    if artifact["classification"] != "intentionally_absent"
                }
                for artifact in manifest["artifacts"]:
                    if not artifact["path"].endswith(".bam") or artifact[
                        "classification"
                    ] != "parsed_exact":
                        continue
                    index_key = (
                        artifact["scenario"],
                        artifact["path"] + ".bai",
                    )
                    if index_key in present:
                        self.assertTrue(
                            artifact.get("require_index"),
                            f"indexed BAM is not validated: {index_key}",
                        )

    def test_full_output_failure_is_approved_only_for_legacy_orchestration(self):
        expected_commands = {
            "single": ["flora", "--full-output"],
            "mixed": ["flora", "mixed", "--full-output"],
        }
        required_outputs = {
            "upstream/matched_reads.fastq.gz",
            "upstream/unmatched_reads.fastq.gz",
            "upstream/cell_reads.fastq.gz",
        }
        for workflow, future_command in expected_commands.items():
            manifest = json.loads(
                (ROOT / "tests" / f"artifact_manifest_{workflow}.json").read_text(
                    encoding="utf-8"
                )
            )
            deviation = manifest["known_deviations"][
                "legacy_full_output_orchestration"
            ]
            with self.subTest(workflow=workflow):
                self.assertTrue(deviation["approved"])
                self.assertEqual("full", deviation["legacy_scenario"])
                self.assertEqual(1, deviation["legacy_expected_exit"])
                self.assertEqual(
                    "Missing cell_reads.fastq.gz",
                    deviation["legacy_diagnostic_contains"],
                )
                self.assertEqual(future_command, deviation["future_command"])
                self.assertEqual(0, deviation["future_expected_exit"])
                self.assertEqual(
                    required_outputs, set(deviation["future_required_outputs"])
                )

                full_oracle = WORKFLOW_FIXTURES / workflow / "oracles" / "full"
                self.assertEqual(
                    "status\t1",
                    (full_oracle / "exit_status.tsv")
                    .read_text(encoding="utf-8")
                    .strip(),
                )
                self.assertIn(
                    deviation["legacy_diagnostic_contains"],
                    (full_oracle / "workflow.log").read_text(encoding="utf-8"),
                )

    def test_tool_versions_record_all_pinned_environment_dimensions(self):
        version_text = (WORKFLOW_FIXTURES / "tool_versions.tsv").read_text(
            encoding="utf-8"
        )
        versions = {
            line.split("\t", 1)[0]: line.split("\t", 1)[1]
            for line in version_text.splitlines()[1:]
        }
        keys = set(versions)
        self.assertTrue(
            {
                "python",
                "python_packages",
                "flora_build",
                "minimap2",
                "samtools",
                "bedtools",
                "libc",
                "architecture",
                "os",
            }.issubset(keys)
        )
        self.assertIn("current version", versions["libc"])

    def test_every_scenario_records_exact_rust_release_stage_selection(self):
        full_pipeline_stages = {
            "flora",
            "generate_26bp_whitelists",
            "prepare_read_tags",
            "add_cb_ur_tags",
            "assign_genes",
            "add_gene_tags",
            "cluster_umis_allbam",
            "cell_umi_gene_table",
            "gene_expression",
            "assign_transcripts",
            "isoform_expression",
            "rna_qc_metrics",
            "read_qc_summary",
        }
        expected_by_scenario = {
            "light": full_pipeline_stages,
            "skip_glycine": full_pipeline_stages,
            "stale_output": full_pipeline_stages,
            "skip_isoform": full_pipeline_stages
            - {"assign_transcripts", "isoform_expression"},
            "upstream_only": {"flora", "generate_26bp_whitelists"},
            "full": {"flora", "generate_26bp_whitelists"},
            "malformed_input": {
                "flora",
                "generate_26bp_whitelists",
                "prepare_read_tags",
            },
            "forced_failure": {
                "flora",
                "generate_26bp_whitelists",
                "prepare_read_tags",
                "add_cb_ur_tags",
            },
        }
        for workflow, scenarios in REQUIRED_SCENARIOS.items():
            for scenario in scenarios:
                selection_log = (
                    WORKFLOW_FIXTURES
                    / workflow
                    / "oracles"
                    / scenario
                    / "rust_stage_selection.log"
                )
                with self.subTest(workflow=workflow, scenario=scenario):
                    rows = selection_log.read_text(encoding="utf-8").splitlines()
                    self.assertEqual(
                        "stage\trelease_binary\tsha256\tinvocation_count", rows[0]
                    )
                    self.assertEqual(
                        expected_by_scenario[scenario],
                        {row.split("\t")[0] for row in rows[1:]},
                    )
                    for row in rows[1:]:
                        stage, binary, digest, invocation_count = row.split("\t")
                        self.assertTrue(binary.endswith(f"/target/release/{stage}"))
                        self.assertRegex(digest, r"^[0-9a-f]{64}$")
                        self.assertGreaterEqual(int(invocation_count), 1)
                    self.assertEqual(
                        "fallback\tcount\npython\t0",
                        (
                            selection_log.parent / "python_fallback_audit.tsv"
                        ).read_text(encoding="utf-8").strip(),
                    )

    def test_frozen_oracles_self_validate_for_every_scenario(self):
        comparator = load_comparator()
        for workflow, scenarios in REQUIRED_SCENARIOS.items():
            manifest = json.loads(
                (ROOT / "tests" / f"artifact_manifest_{workflow}.json").read_text(
                    encoding="utf-8"
                )
            )
            for scenario in sorted(scenarios):
                expected = WORKFLOW_FIXTURES / workflow / "oracles" / scenario
                with self.subTest(workflow=workflow, scenario=scenario):
                    with tempfile.TemporaryDirectory() as temp:
                        actual = Path(temp) / "actual"
                        actual.mkdir()
                        artifacts = [
                            artifact
                            for artifact in manifest["artifacts"]
                            if artifact["scenario"] == scenario
                            and artifact["classification"]
                            != "intentionally_absent"
                        ]
                        for artifact in artifacts:
                            source = expected / artifact.get(
                                "oracle_path", artifact["path"]
                            )
                            destination = actual / artifact["path"]
                            destination.parent.mkdir(parents=True, exist_ok=True)
                            if source.is_symlink():
                                destination.symlink_to(os.readlink(source))
                            elif "oracle_path" in artifact:
                                import gzip

                                with gzip.open(source, "rb") as input_handle:
                                    destination.write_bytes(input_handle.read())
                            else:
                                shutil.copyfile(source, destination)
                        comparator.compare_tree(
                            expected, actual, manifest, scenario=scenario
                        )

    def test_regenerator_is_explicit_pinned_and_guards_rust_stages(self):
        script_path = WORKFLOW_FIXTURES / "generate_legacy_oracles.sh"
        script = script_path.read_text(
            encoding="utf-8"
        )
        self.assertIn("--refresh-environment", script)
        self.assertIn(".conda-env/bin", script)
        self.assertIn("cargo build --release --locked --bins", script)
        self.assertIn("run_all.sh", script)
        self.assertIn("run_all_mixed_species.sh", script)
        self.assertIn("forbidden Python fallback", script)
        self.assertIn("FLORA_ORACLE_INVOCATION_LOG", script)
        self.assertIn("python_fallback_audit.tsv", script)
        self.assertIn('exec "${real_binary}" "$@"', script)
        self.assertIn("tool_versions.tsv", script)
        for scenario in sorted(set().union(*REQUIRED_SCENARIOS.values())):
            self.assertIn(scenario, script)
        self.assertNotRegex(script, r"/Users/(?!\$|\{)")
        syntax = subprocess.run(
            ["bash", "-n", str(script_path)], text=True, capture_output=True
        )
        self.assertEqual(0, syntax.returncode, syntax.stderr)

    def test_numbered_qc_collision_artifacts_are_rejected_and_cleaned(self):
        collision_name = re.compile(
            r"^(?:(?:cell_umi_gene|filtered\.sorted|.+\.single_cell_report) [0-9]+"
            r"(?:\.tsv|\.bam|\.html(?:\.gz)?)|filtered\.sorted\.bam [0-9]+\.bai)$"
        )
        collisions = [
            str(path.relative_to(WORKFLOW_FIXTURES))
            for path in WORKFLOW_FIXTURES.rglob("*")
            if collision_name.search(path.name)
        ]
        self.assertEqual([], collisions)

        script = (
            WORKFLOW_FIXTURES / "generate_legacy_oracles.sh"
        ).read_text(encoding="utf-8")
        self.assertIn("remove_numbered_qc_collision_artifacts", script)
        post_compression_cleanup = (
            'remove_numbered_qc_collision_artifacts "${oracle_dir}"'
        )
        self.assertLess(
            script.index("gzip -n -9"),
            script.index(post_compression_cleanup),
        )
        self.assertLess(
            script.index(post_compression_cleanup),
            script.index('import hashlib\nimport json'),
        )


class WorkflowComparatorSecurityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.comparator = load_comparator()

    def test_rejects_unclassified_extra_and_missing_paths(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            expected = root / "expected"
            actual = root / "actual"
            expected.mkdir()
            actual.mkdir()
            (expected / "value.tsv").write_text("key\tvalue\na\t1\n", encoding="utf-8")
            manifest = {
                "numeric_tolerance": {"absolute": 1e-9, "relative": 1e-7},
                "artifacts": [
                    {
                        "path": "value.tsv",
                        "classification": "parsed_exact",
                        "sha256": self.comparator.sha256_path(expected / "value.tsv"),
                    }
                ],
            }
            with self.assertRaisesRegex(self.comparator.ComparisonError, "missing"):
                self.comparator.compare_tree(expected, actual, manifest)
            (actual / "value.tsv").write_text("key\tvalue\na\t1\n", encoding="utf-8")
            (actual / "extra.txt").write_text("extra\n", encoding="utf-8")
            with self.assertRaisesRegex(self.comparator.ComparisonError, "extra"):
                self.comparator.compare_tree(expected, actual, manifest)

    def test_missing_pinned_samtools_falls_back_to_active_environment(self):
        configured = {"samtools": ".conda-env/bin/samtools"}
        previous = os.environ.get("FLORA_SAMTOOLS")
        os.environ["FLORA_SAMTOOLS"] = "/active/conda/bin/samtools"
        try:
            self.assertEqual(
                "/active/conda/bin/samtools",
                self.comparator._samtools(configured),
            )
        finally:
            if previous is None:
                os.environ.pop("FLORA_SAMTOOLS", None)
            else:
                os.environ["FLORA_SAMTOOLS"] = previous

    def test_rejects_traversal_and_escaping_symlink(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            expected = root / "expected"
            actual = root / "actual"
            expected.mkdir()
            actual.mkdir()
            manifest = {
                "numeric_tolerance": {"absolute": 1e-9, "relative": 1e-7},
                "artifacts": [
                    {
                        "path": "../escape",
                        "classification": "byte_exact",
                        "sha256": "0" * 64,
                    }
                ],
            }
            with self.assertRaisesRegex(self.comparator.ComparisonError, "traversal"):
                self.comparator.compare_tree(expected, actual, manifest)

            outside = root / "outside.txt"
            outside.write_text("outside\n", encoding="utf-8")
            (expected / "link").symlink_to(outside)
            (actual / "link").symlink_to(outside)
            manifest["artifacts"][0]["path"] = "link"
            manifest["artifacts"][0]["sha256"] = self.comparator.sha256_path(
                expected / "link"
            )
            with self.assertRaisesRegex(self.comparator.ComparisonError, "symlink"):
                self.comparator.compare_tree(expected, actual, manifest)

    def test_numeric_tolerance_is_fieldwise_and_finite(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            expected = root / "expected"
            actual = root / "actual"
            expected.mkdir()
            actual.mkdir()
            (expected / "metrics.tsv").write_text(
                "metric\tvalue\nscore\t1.0\n", encoding="utf-8"
            )
            (actual / "metrics.tsv").write_text(
                "metric\tvalue\nscore\t1.00000001\n", encoding="utf-8"
            )
            manifest = {
                "numeric_tolerance": {"absolute": 1e-9, "relative": 1e-7},
                "artifacts": [
                    {
                        "path": "metrics.tsv",
                        "classification": "numeric_tolerance",
                        "sha256": self.comparator.sha256_path(expected / "metrics.tsv"),
                    }
                ],
            }
            self.comparator.compare_tree(expected, actual, manifest)
            (actual / "metrics.tsv").write_text(
                "metric\tvalue\nscore\tnan\n", encoding="utf-8"
            )
            with self.assertRaises(self.comparator.ComparisonError):
                self.comparator.compare_tree(expected, actual, manifest)

    def test_compressed_oracle_maps_to_uncompressed_actual_path(self):
        import gzip

        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            expected = root / "expected"
            actual = root / "actual"
            expected.mkdir()
            actual.mkdir()
            html = "<html><body>report</body></html>\n"
            with gzip.GzipFile(
                expected / "report.html.gz", mode="wb", mtime=0
            ) as handle:
                handle.write(html.encode("utf-8"))
            (actual / "report.html").write_text(html, encoding="utf-8")
            manifest = {
                "numeric_tolerance": {"absolute": 1e-9, "relative": 1e-7},
                "artifacts": [
                    {
                        "path": "report.html",
                        "oracle_path": "report.html.gz",
                        "classification": "canonicalized_html_log",
                        "canonicalize": [],
                        "sha256": self.comparator.sha256_path(
                            expected / "report.html.gz"
                        ),
                    }
                ],
            }
            self.comparator.compare_tree(expected, actual, manifest)

    def test_internal_symlink_targets_compare_relative_to_each_tree(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            expected = root / "expected"
            actual = root / "actual"
            for tree in (expected, actual):
                (tree / "data").mkdir(parents=True)
                (tree / "data" / "value.txt").write_text("value\n", encoding="utf-8")
            (expected / "link").symlink_to("data/value.txt")
            (actual / "link").symlink_to(actual / "data" / "value.txt")
            manifest = {
                "numeric_tolerance": {"absolute": 1e-9, "relative": 1e-7},
                "artifacts": [
                    {
                        "path": "data/value.txt",
                        "classification": "parsed_exact",
                        "sha256": self.comparator.sha256_path(
                            expected / "data" / "value.txt"
                        ),
                    },
                    {
                        "path": "link",
                        "classification": "byte_exact",
                        "sha256": self.comparator.sha256_path(expected / "link"),
                    },
                ],
            }
            self.comparator.compare_tree(expected, actual, manifest)


if __name__ == "__main__":
    unittest.main()
