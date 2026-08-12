import importlib.util
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
BUILD_REPORT_PATH = PROJECT_ROOT / "scripts" / "build_report.py"
TEMPLATE_PATH = PROJECT_ROOT / "scripts" / "report_template.html"
PLOTLY_PATH = PROJECT_ROOT / "scripts" / "plotly-2.26.0.min.js"


def load_build_report_module():
    spec = importlib.util.spec_from_file_location("strint_build_report", BUILD_REPORT_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class ReportTemplateContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.build_report = load_build_report_module()

    def test_parse_args_accepts_barnyard_paths_with_spaces(self):
        summary_path = "/tmp/report inputs/barnyard summary.tsv"
        per_cell_path = "/tmp/report inputs/barnyard per cell.tsv"
        argv = [
            "build_report.py",
            "--sample-id",
            "sample",
            "--output-html",
            "report.html",
            "--report-metrics-tsv",
            "report_metrics.tsv",
            "--rna-qc-metrics-tsv",
            "rna_qc_metrics.tsv",
            "--saturation-tsv",
            "saturation.tsv",
            "--read-qc-json",
            "read_qc.json",
            "--parameters-tsv",
            "parameters.tsv",
            "--per-cell-qc-tsv",
            "per_cell_qc.tsv",
            "--rna-cluster-tsv",
            "rna_cluster.tsv",
            "--barnyard-summary-tsv",
            summary_path,
            "--barnyard-per-cell-tsv",
            per_cell_path,
        ]

        with patch.object(sys, "argv", argv):
            try:
                args = self.build_report.parse_args()
            except SystemExit as exc:
                self.fail(f"Barnyard report arguments were rejected: {exc}")

        self.assertEqual(args.barnyard_summary_tsv, summary_path)
        self.assertEqual(args.barnyard_per_cell_tsv, per_cell_path)

    def write_barnyard_inputs(self, directory, summary_rows=None, per_cell_rows=None):
        summary_path = Path(directory) / "barnyard_summary.tsv"
        per_cell_path = Path(directory) / "barnyard_per_cell.tsv"
        if summary_rows is None:
            summary_rows = [
                ("total_cells", "12"),
                ("human_singlet_cells", "5"),
                ("mouse_singlet_cells", "4"),
                ("mixed_cells", "2"),
                ("unclassified_cells", "1"),
                ("cross_species_doublet_rate_among_classified", "0.1818"),
            ]
        if per_cell_rows is None:
            per_cell_rows = [
                ("CELL_A", "human_singlet", "7", "1", "9"),
                ("CELL_B", "mouse_singlet", "0", "8", "10"),
            ]
        summary_path.write_text(
            "metric\tvalue\n"
            + "".join(f"{metric}\t{value}\n" for metric, value in summary_rows),
            encoding="utf-8",
        )
        per_cell_path.write_text(
            "cell_id\tassignment\tumi_human\tumi_mouse\tumi_total\n"
            + "".join("\t".join(row) + "\n" for row in per_cell_rows),
            encoding="utf-8",
        )
        return summary_path, per_cell_path

    def report_sections(self, **overrides):
        section_names = [
            "sample_title",
            "summary_cards",
            "summary_bar",
            "summary_rows",
            "summary_violin_plots",
            "read_qc_bar",
            "beads_bar",
            "rna_cluster_bar",
            "sequencing_bar",
            "read_summary_rows",
            "mapping_bar",
            "mapping_rows",
            "saturation_bar",
        ]
        sections = {name: f"SECTION_{name}" for name in section_names}
        sections.update(overrides)
        return sections

    def minimal_payload(self, **overrides):
        payload = {
            "readSummary": {"labels": [], "counts": []},
            "readQc": {},
            "perCell": {"reads": [], "umis": [], "genes": []},
            "saturation": {},
            "barcodeRank5p": {"rank": [], "count": [], "is_true": []},
            "beads": {"x": [], "y": [], "n_cells": 0},
            "rnaCluster": [],
        }
        payload.update(overrides)
        return payload

    def test_barnyard_summary_selects_and_strictly_formats_six_metrics(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary_path, _ = self.write_barnyard_inputs(
                tmp,
                summary_rows=[
                    ("total_cells", "1234"),
                    ("human_singlet_cells", "-1"),
                    ("mouse_singlet_cells", "2.5"),
                    ("mixed_cells", "7"),
                    ("unclassified_cells", "nan"),
                    ("cross_species_doublet_rate_among_classified", "1.2"),
                ],
            )
            stderr = io.StringIO()
            with redirect_stderr(stderr):
                rows = self.build_report.barnyard_summary_rows(summary_path)

        self.assertEqual(
            rows,
            [
                {"Metric": "Cells after ambient filter", "Value": "1,234"},
                {"Metric": "Human singlets", "Value": "NA"},
                {"Metric": "Mouse singlets", "Value": "NA"},
                {"Metric": "Mixed cells", "Value": "7"},
                {"Metric": "Unclassified cells", "Value": "NA"},
                {"Metric": "Cross-species doublet rate", "Value": "NA"},
            ],
        )
        self.assertIn("4 invalid Barnyard summary values", stderr.getvalue())

    def test_barnyard_summary_missing_values_are_na(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary_path, _ = self.write_barnyard_inputs(
                tmp,
                summary_rows=[("total_cells", "3")],
            )
            with redirect_stderr(io.StringIO()):
                rows = self.build_report.barnyard_summary_rows(summary_path)

        self.assertEqual(rows[0]["Value"], "3")
        self.assertTrue(all(row["Value"] == "NA" for row in rows[1:]))

    def test_barnyard_summary_formats_a_valid_rate_as_a_percentage(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary_path, _ = self.write_barnyard_inputs(tmp)
            rows = self.build_report.barnyard_summary_rows(summary_path)

        self.assertEqual(rows[-1]["Value"], "18.18%")

    def test_barnyard_summary_preserves_an_integer_larger_than_float_precision(self):
        exact_count = "9007199254740993"
        with tempfile.TemporaryDirectory() as tmp:
            summary_path, _ = self.write_barnyard_inputs(
                tmp,
                summary_rows=[
                    ("total_cells", exact_count),
                    ("human_singlet_cells", "0"),
                    ("mouse_singlet_cells", "0"),
                    ("mixed_cells", "0"),
                    ("unclassified_cells", "0"),
                    ("cross_species_doublet_rate_among_classified", "0"),
                ],
            )
            rows = self.build_report.barnyard_summary_rows(summary_path)

        self.assertEqual(rows[0]["Value"], "9,007,199,254,740,993")

    def test_barnyard_summary_treats_serialization_overflow_as_invalid(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary_path, per_cell_path = self.write_barnyard_inputs(
                tmp,
                summary_rows=[
                    ("total_cells", "1e4300"),
                    ("human_singlet_cells", "0"),
                    ("mouse_singlet_cells", "0"),
                    ("mixed_cells", "0"),
                    ("unclassified_cells", "0"),
                    ("cross_species_doublet_rate_among_classified", "0"),
                ],
            )
            stderr = io.StringIO()
            with redirect_stderr(stderr):
                payload = self.build_report.barnyard_payload(summary_path, per_cell_path)

        self.assertIsNotNone(payload)
        self.assertEqual(payload["summaryRows"][0]["Value"], "NA")
        self.assertIn("1 invalid Barnyard summary values", stderr.getvalue())
        json.dumps(payload)

    def test_barnyard_integer_parser_is_exact_and_rejects_invalid_numbers(self):
        parser = self.build_report._finite_nonnegative_integer

        self.assertEqual(parser("9007199254740993"), 9007199254740993)
        self.assertEqual(parser(" 7.0 "), 7)
        self.assertEqual(parser("1e3"), 1000)
        for invalid in ("1.5", "-1", "nan", "inf", "-inf", ""):
            with self.subTest(value=invalid):
                self.assertIsNone(parser(invalid))

    def test_barnyard_payload_rejects_duplicate_selected_summary_metric(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary_path, per_cell_path = self.write_barnyard_inputs(tmp)
            with summary_path.open("a", encoding="utf-8") as handle:
                handle.write("total_cells\t12\n")
            stderr = io.StringIO()
            with redirect_stderr(stderr):
                payload = self.build_report.barnyard_payload(summary_path, per_cell_path)

        self.assertIsNone(payload)
        self.assertIn("duplicate selected metric", stderr.getvalue())

    def test_barnyard_payload_requires_both_existing_well_formed_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary_path, per_cell_path = self.write_barnyard_inputs(tmp)
            malformed = Path(tmp) / "malformed.tsv"
            malformed.write_text("wrong\theader\n1\t2\n", encoding="utf-8")
            cases = [
                (summary_path, None),
                (None, per_cell_path),
                (Path(tmp) / "missing.tsv", per_cell_path),
                (malformed, per_cell_path),
                (summary_path, malformed),
            ]
            for summary_arg, per_cell_arg in cases:
                with self.subTest(summary=summary_arg, per_cell=per_cell_arg):
                    stderr = io.StringIO()
                    with redirect_stderr(stderr):
                        payload = self.build_report.barnyard_payload(summary_arg, per_cell_arg)
                    self.assertIsNone(payload)
                    self.assertIn("WARNING", stderr.getvalue())

    def test_barnyard_payload_validates_before_dedup_and_preserves_umi_total(self):
        rows = [
            ("CELL_A", "human_singlet", "bad", "1", "2"),
            ("CELL_A", "human_singlet", "7", "1", "11"),
            ("CELL_A", "mouse_singlet", "0", "9", "12"),
            ("", "mixed", "1", "1", "2"),
            ("CELL_UNKNOWN", "other", "1", "1", "2"),
            ("CELL_NEG", "mixed", "-1", "2", "3"),
            ("CELL_FRAC", "mixed", "1.5", "2", "4"),
            ("CELL_INF", "mixed", "inf", "2", "4"),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            summary_path, per_cell_path = self.write_barnyard_inputs(
                tmp, per_cell_rows=rows
            )
            stderr = io.StringIO()
            with redirect_stderr(stderr):
                payload = self.build_report.barnyard_payload(summary_path, per_cell_path)

        trace = payload["traces"]["human_singlet"]
        self.assertEqual(trace["cell_id"], ["CELL_A"])
        self.assertEqual(trace["umi_human"], [7])
        self.assertEqual(trace["umi_mouse"], [1])
        self.assertEqual(trace["umi_total"], [11])
        self.assertEqual(payload["classCounts"]["human_singlet"], 1)
        self.assertEqual(payload["totalValid"], 1)
        self.assertEqual(payload["displayed"], 1)
        self.assertFalse(payload["sampled"])
        warnings = stderr.getvalue()
        self.assertIn("1 duplicate cell IDs", warnings)
        self.assertIn("1 unknown assignments", warnings)
        self.assertIn("5 invalid per-cell rows", warnings)

    def test_barnyard_payload_preserves_large_per_cell_integers_exactly(self):
        exact_total = "9007199254740991"
        with tempfile.TemporaryDirectory() as tmp:
            summary_path, per_cell_path = self.write_barnyard_inputs(
                tmp,
                per_cell_rows=[
                    ("CELL_BIG", "human_singlet", exact_total, "1", exact_total),
                ],
            )
            payload = self.build_report.barnyard_payload(summary_path, per_cell_path)

        trace = payload["traces"]["human_singlet"]
        self.assertEqual(trace["umi_human"], [9007199254740991])
        self.assertEqual(trace["umi_total"], [9007199254740991])

    def test_barnyard_payload_skips_js_unsafe_umi_before_dedup(self):
        exact_total = "9007199254740991"
        with tempfile.TemporaryDirectory() as tmp:
            summary_path, per_cell_path = self.write_barnyard_inputs(
                tmp,
                per_cell_rows=[
                    ("CELL_BIG", "human_singlet", "9007199254740992", "1", "9007199254740992"),
                    ("CELL_BIG", "human_singlet", exact_total, "1", exact_total),
                ],
            )
            stderr = io.StringIO()
            with redirect_stderr(stderr):
                payload = self.build_report.barnyard_payload(summary_path, per_cell_path)

        self.assertIsNotNone(payload)
        trace = payload["traces"]["human_singlet"]
        self.assertEqual(trace["cell_id"], ["CELL_BIG"])
        self.assertEqual(trace["umi_human"], [9007199254740991])
        self.assertEqual(trace["umi_total"], [9007199254740991])
        self.assertIn("1 invalid per-cell rows", stderr.getvalue())
        self.assertNotIn("duplicate cell IDs", stderr.getvalue())
        json.dumps(payload)

    def test_barnyard_payload_rejects_js_unsafe_value_in_any_umi_column(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary_path, per_cell_path = self.write_barnyard_inputs(
                tmp,
                per_cell_rows=[
                    ("CELL_H", "human_singlet", "9007199254740992", "1", "2"),
                    ("CELL_M", "mouse_singlet", "1", "9007199254740992", "2"),
                    ("CELL_T", "mixed", "1", "1", "9007199254740992"),
                ],
            )
            stderr = io.StringIO()
            with redirect_stderr(stderr):
                payload = self.build_report.barnyard_payload(summary_path, per_cell_path)

        self.assertEqual(payload["totalValid"], 0)
        self.assertEqual(payload["displayed"], 0)
        self.assertIn("3 invalid per-cell rows", stderr.getvalue())

    def test_barnyard_header_only_input_is_valid(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary_path, per_cell_path = self.write_barnyard_inputs(
                tmp, per_cell_rows=[]
            )
            payload = self.build_report.barnyard_payload(summary_path, per_cell_path)

        self.assertEqual(payload["totalValid"], 0)
        self.assertEqual(payload["displayed"], 0)
        self.assertFalse(payload["sampled"])
        self.assertTrue(all(not trace["cell_id"] for trace in payload["traces"].values()))

    def test_barnyard_sampler_retains_all_through_exact_threshold(self):
        sampler = self.build_report.BarnyardPointSampler()
        for index in range(self.build_report.BARNYARD_MAX_POINTS):
            sampler.add((f"CELL_{index}", "human_singlet", index, 0, index))

        traces = sampler.finish()

        self.assertEqual(sampler.retained_count, 100_000)
        self.assertEqual(sampler.max_retained_count, 100_000)
        self.assertFalse(sampler.sampled)
        self.assertEqual(len(traces["human_singlet"]["cell_id"]), 100_000)

    def test_barnyard_sampler_converts_at_100001_and_is_deterministic(self):
        def sampled_cells():
            sampler = self.build_report.BarnyardPointSampler()
            for index in range(self.build_report.BARNYARD_MAX_POINTS + 1):
                assignment = self.build_report.BARNYARD_ASSIGNMENTS[index % 4]
                sampler.add((f"CELL_{index}", assignment, index, index + 1, index + 2))
            return sampler, sampler.finish()

        first_sampler, first = sampled_cells()
        second_sampler, second = sampled_cells()

        self.assertTrue(first_sampler.sampled)
        self.assertLessEqual(first_sampler.retained_count, 100_000)
        self.assertLessEqual(first_sampler.max_retained_count, 100_000)
        self.assertTrue(
            all(size <= 25_000 for size in first_sampler.heap_sizes_by_class.values())
        )
        self.assertEqual(first, second)
        self.assertEqual(first_sampler.retained_count, second_sampler.retained_count)

    def test_barnyard_sampler_caps_an_unbalanced_class_after_conversion(self):
        sampler = self.build_report.BarnyardPointSampler()
        for index in range(self.build_report.BARNYARD_MAX_POINTS + 1):
            sampler.add((f"CELL_{index}", "human_singlet", index, 0, index))

        traces = sampler.finish()

        self.assertTrue(sampler.sampled)
        self.assertEqual(sampler.class_counts["human_singlet"], 100_001)
        self.assertEqual(sampler.retained_count, 25_000)
        self.assertEqual(len(traces["human_singlet"]["cell_id"]), 25_000)
        self.assertTrue(
            all(
                not traces[assignment]["cell_id"]
                for assignment in self.build_report.BARNYARD_ASSIGNMENTS[1:]
            )
        )

    def test_barnyard_payload_uses_chunked_selected_columns_and_cleans_sqlite(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary_path, per_cell_path = self.write_barnyard_inputs(tmp)
            observed_reads = []
            observed_databases = []
            real_read_csv = self.build_report.pd.read_csv
            real_mkstemp = self.build_report.tempfile.mkstemp

            def recording_read_csv(*args, **kwargs):
                if Path(args[0]) == per_cell_path:
                    observed_reads.append(dict(kwargs))
                return real_read_csv(*args, **kwargs)

            def recording_mkstemp(*args, **kwargs):
                descriptor, path = real_mkstemp(*args, **kwargs)
                observed_databases.append(path)
                return descriptor, path

            with patch.object(self.build_report.pd, "read_csv", side_effect=recording_read_csv), patch.object(
                self.build_report.tempfile, "mkstemp", side_effect=recording_mkstemp
            ):
                payload = self.build_report.barnyard_payload(summary_path, per_cell_path)

        self.assertIsNotNone(payload)
        self.assertTrue(observed_reads)
        self.assertEqual(
            set(observed_reads[0]["usecols"]),
            {"cell_id", "assignment", "umi_human", "umi_mouse", "umi_total"},
        )
        self.assertEqual(observed_reads[0]["chunksize"], self.build_report.BARNYARD_CHUNK_ROWS)
        self.assertTrue(observed_databases)
        self.assertTrue(all(not os.path.exists(path) for path in observed_databases))
        serialized = json.dumps(payload)
        self.assertNotIn("max_retained_count", serialized)
        self.assertNotIn("heap_sizes_by_class", serialized)

    def test_barnyard_payload_cleans_sqlite_after_a_read_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary_path, per_cell_path = self.write_barnyard_inputs(tmp)
            observed_databases = []
            real_mkstemp = self.build_report.tempfile.mkstemp
            real_read_csv = self.build_report.pd.read_csv

            def recording_mkstemp(*args, **kwargs):
                descriptor, path = real_mkstemp(*args, **kwargs)
                observed_databases.append(path)
                return descriptor, path

            def failing_per_cell_read(*args, **kwargs):
                if Path(args[0]) == per_cell_path:
                    raise ValueError("broken TSV")
                return real_read_csv(*args, **kwargs)

            with patch.object(
                self.build_report.tempfile, "mkstemp", side_effect=recording_mkstemp
            ), patch.object(
                self.build_report.pd, "read_csv", side_effect=failing_per_cell_read
            ), redirect_stderr(io.StringIO()):
                payload = self.build_report.barnyard_payload(summary_path, per_cell_path)

        self.assertIsNone(payload)
        self.assertTrue(observed_databases)
        self.assertTrue(all(not os.path.exists(path) for path in observed_databases))

    def test_valid_barnyard_section_is_between_beads_and_rna_cluster(self):
        barnyard = {
            "summaryRows": [
                {"Metric": "Cells after ambient filter", "Value": "12"},
                {"Metric": "Human singlets", "Value": "5"},
                {"Metric": "Mouse singlets", "Value": "4"},
                {"Metric": "Mixed cells", "Value": "2"},
                {"Metric": "Unclassified cells", "Value": "1"},
                {"Metric": "Cross-species doublet rate", "Value": "18.18%"},
            ],
            "traces": {
                "human_singlet": {"cell_id": ["CELL_A"], "umi_human": [7], "umi_mouse": [1], "umi_total": [9]},
                "mouse_singlet": {"cell_id": [], "umi_human": [], "umi_mouse": [], "umi_total": []},
                "mixed": {"cell_id": [], "umi_human": [], "umi_mouse": [], "umi_total": []},
                "unclassified": {"cell_id": [], "umi_human": [], "umi_mouse": [], "umi_total": []},
            },
            "classCounts": {"human_singlet": 5, "mouse_singlet": 4, "mixed": 2, "unclassified": 1},
            "totalValid": 12,
            "displayed": 1,
            "sampled": True,
        }
        sections = self.report_sections(
            barnyard_section=self.build_report.barnyard_report_section(barnyard)
        )
        rendered = self.build_report.build_html(
            SimpleNamespace(sample_id="MIXED"),
            self.minimal_payload(barnyard=barnyard),
            sections,
        )
        cells_html = rendered.split('id="cells-tab"', 1)[1].split('id="library-tab"', 1)[0]

        self.assertLess(cells_html.index("SECTION_beads_bar"), cells_html.index("Barnyard QC"))
        self.assertLess(cells_html.index("Barnyard QC"), cells_html.index("SECTION_rna_cluster_bar"))
        self.assertIn('id="barnyard-umi"', cells_html)
        for label in [
            "Cells after ambient filter",
            "Human singlets",
            "Mouse singlets",
            "Mixed cells",
            "Unclassified cells",
            "Cross-species doublet rate",
        ]:
            self.assertIn(label, cells_html)
        for help_text in [
            "Ambient filtering",
            "Singlets",
            "Mixed cells",
            "Cross-species doublet rate",
            "Human UMI",
            "Mouse UMI",
            "Class colors",
            "display sampling",
        ]:
            self.assertIn(help_text, cells_html)

    def test_barnyard_renderer_has_expected_plotly_contract(self):
        barnyard = {
            "summaryRows": [],
            "traces": {
                assignment: {"cell_id": [], "umi_human": [], "umi_mouse": [], "umi_total": []}
                for assignment in self.build_report.BARNYARD_ASSIGNMENTS
            },
            "classCounts": {"human_singlet": 8, "mouse_singlet": 7, "mixed": 3, "unclassified": 2},
            "totalValid": 20,
            "displayed": 2001,
            "sampled": True,
        }
        rendered = self.build_report.build_html(
            SimpleNamespace(sample_id="MIXED"),
            self.minimal_payload(barnyard=barnyard),
            self.report_sections(
                barnyard_section=self.build_report.barnyard_report_section(barnyard)
            ),
        )

        self.assertIn('const barnyardTraceType = barnyard.displayed > 2000 ? "scattergl" : "scatter";', rendered)
        for color in ["#2E7D32", "#1565C0", "#C62828", "#757575"]:
            self.assertIn(color, rendered)
        self.assertIn('x: trace.umi_human', rendered)
        self.assertIn('y: trace.umi_mouse', rendered)
        self.assertIn('name: `${barnyardClassLabels[assignment]} (${fullCount.toLocaleString()})`', rendered)
        self.assertIn('Cell ID: %{customdata[0]}', rendered)
        self.assertIn('Assignment: %{customdata[1]}', rendered)
        self.assertIn('Human UMI: %{x:,.0f}', rendered)
        self.assertIn('Mouse UMI: %{y:,.0f}', rendered)
        self.assertIn('Total UMI: %{customdata[2]:,.0f}', rendered)
        self.assertIn('text: "Human UMI"', rendered)
        self.assertIn('text: "Mouse UMI"', rendered)
        self.assertIn('rangemode: "nonnegative"', rendered)
        self.assertIn('displayModeBar: true', rendered)
        self.assertIn('scrollZoom: true', rendered)
        self.assertIn('modeBarButtonsToRemove:', rendered)
        self.assertIn('orientation: "v"', rendered)
        self.assertIn('x: 0.98', rendered)
        self.assertIn('y: 0.98', rendered)
        self.assertIn('xanchor: "right"', rendered)
        self.assertIn('yanchor: "top"', rendered)
        self.assertIn('bgcolor: "rgba(255,255,255,0.88)"', rendered)
        self.assertIn('bordercolor: "rgba(44,62,80,0.25)"', rendered)
        self.assertNotIn('y: 1.02', self.build_report.barnyard_report_script(self.minimal_payload(barnyard=barnyard)))
        self.assertIn('${barnyard.displayed.toLocaleString()} of ${barnyard.totalValid.toLocaleString()}', rendered)
        self.assertIn('function plotIf(id, data, layout, config)', rendered)
        self.assertIn('Object.assign({}, plotConfig, config || {})', rendered)
        self.assertIn('displayModeBar: false', rendered)

    def test_header_only_barnyard_input_renders_empty_state_without_plot_call(self):
        barnyard = {
            "summaryRows": [],
            "traces": {
                assignment: {"cell_id": [], "umi_human": [], "umi_mouse": [], "umi_total": []}
                for assignment in self.build_report.BARNYARD_ASSIGNMENTS
            },
            "classCounts": {assignment: 0 for assignment in self.build_report.BARNYARD_ASSIGNMENTS},
            "totalValid": 0,
            "displayed": 0,
            "sampled": False,
        }
        rendered = self.build_report.build_html(
            SimpleNamespace(sample_id="MIXED"),
            self.minimal_payload(barnyard=barnyard),
            self.report_sections(
                barnyard_section=self.build_report.barnyard_report_section(barnyard)
            ),
        )

        self.assertIn("No valid Barnyard cells were available for plotting.", rendered)
        self.assertNotIn('plotIf("barnyard-umi"', rendered)

    def test_report_without_barnyard_has_no_barnyard_artifacts(self):
        rendered = self.build_report.build_html(
            SimpleNamespace(sample_id="SINGLE"),
            self.minimal_payload(),
            self.report_sections(),
        )

        for marker in ["Barnyard QC", 'id="barnyard-umi"', '"barnyard":', "payload.barnyard", 'plotIf("barnyard-umi"']:
            self.assertNotIn(marker, rendered)

    def test_template_is_rna_only_new_report_shell(self):
        template = TEMPLATE_PATH.read_text(encoding="utf-8")

        self.assertIn('data-template-source="template_multi"', template)
        self.assertIn("__SAMPLE__", template)
        self.assertIn("__REPORT_BODY__", template)
        self.assertIn("__PLOTLY_LOADER__", template)
        self.assertIn('class="sidebar-navigation"', template)
        self.assertIn('data-tab="summary"', template)
        self.assertIn('data-tab="cells"', template)
        self.assertIn('data-tab="library"', template)
        self.assertIn(">RNA<", template)
        self.assertRegex(
            template,
            r'\[data-library-content="gene-expression"\]\s*\{\s*display:\s*block',
        )
        self.assertNotIn("P2026042903", template)
        self.assertNotIn("VDJ-T", template)
        self.assertNotIn("VDJ-B", template)
        self.assertNotIn(">ATAC<", template)
        self.assertLess(TEMPLATE_PATH.stat().st_size, 1_500_000)

    def test_barcode_rank_plot_uses_only_the_5prime_curve(self):
        source = BUILD_REPORT_PATH.read_text(encoding="utf-8")

        self.assertNotIn("payload.barcodeRank3p", source)
        self.assertNotIn("3' barcode", source)
        self.assertNotIn("3' points", source)
        self.assertIn("payload.barcodeRank5p", source)
        self.assertIn("5' barcode", source)

    def test_barcode_rank_payload_labels_true_and_noise_from_observed_whitelist(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            counts = tmp_path / "counts.tsv"
            whitelist = tmp_path / "whitelist.csv"
            counts.write_text(
                "barcode\tcount\nBC_A\t100\nBC_B\t50\nBC_C\t10\n",
                encoding="utf-8",
            )
            whitelist.write_text("BC_A\nBC_C\n", encoding="utf-8")

            payload = self.build_report.barcode_rank_payload(counts, whitelist)

        self.assertEqual(payload["rank"], [1, 2, 3])
        self.assertEqual(payload["count"], [100.0, 50.0, 10.0])
        self.assertEqual(payload["is_true"], [True, False, True])
        self.assertNotIn("threshold", payload)

    def test_beads_payload_prefers_one_row_per_cell_barcode_map(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            barcode_to_cell = tmp_path / "barcode_to_cell.csv"
            read_assigned_cell = tmp_path / "read_assigned_cell.csv"
            barcode_to_cell.write_text(
                "cell,barcode,is_cell_barcode\n"
                "CELL1_N2,BC_A;BC_B,1\n"
                "CELL2_N1,BC_C,1\n",
                encoding="utf-8",
            )
            read_assigned_cell.write_text(
                "cell_id,BC5n,BC3n\n"
                "should_not_be_used,WRONG_A,WRONG_B\n",
                encoding="utf-8",
            )

            payload = self.build_report.beads_per_droplet_payload(
                read_assigned_cell,
                barcode_to_cell,
            )

        self.assertEqual(payload, {"x": [1, 2], "y": [1, 1], "n_cells": 2})

    def test_beads_payload_falls_back_to_read_assignments(self):
        with tempfile.TemporaryDirectory() as tmp:
            read_assigned_cell = Path(tmp) / "read_assigned_cell.csv"
            read_assigned_cell.write_text(
                "cell_id,BC5n,BC3n\n"
                "cell_1,BC_A,BC_B\n"
                "cell_1,BC_A,BC_B\n"
                "cell_2,BC_C,BC_C\n",
                encoding="utf-8",
            )

            payload = self.build_report.beads_per_droplet_payload(
                read_assigned_cell,
                None,
            )

        self.assertEqual(payload, {"x": [1, 2], "y": [1, 1], "n_cells": 2})

    def test_reference_plot_colors_and_labels_are_used(self):
        source = BUILD_REPORT_PATH.read_text(encoding="utf-8")

        self.assertIn('name: "TRUE"', source)
        self.assertIn('name: "NOISE"', source)
        self.assertIn('color: "#1358A2"', source)
        self.assertIn('color: "#DDDDDD"', source)
        self.assertIn('"rgba(15, 76, 129, 0.8)"', source)
        self.assertIn('"rgba(163, 194, 255, 0.8)"', source)
        self.assertIn('type: "violin"', source)
        self.assertIn('fillcolor: "rgba(41, 128, 185, 0.8)"', source)
        self.assertNotIn("5' threshold", source)

    def test_bundled_plotly_starts_with_a_valid_javascript_comment(self):
        self.assertTrue(PLOTLY_PATH.exists())
        self.assertTrue(PLOTLY_PATH.read_text(encoding="utf-8").startswith("/*!"))

    def test_build_html_preserves_all_existing_report_content(self):
        section_names = [
            "sample_title",
            "summary_cards",
            "summary_bar",
            "summary_rows",
            "summary_violin_plots",
            "read_qc_bar",
            "beads_bar",
            "rna_cluster_bar",
            "sequencing_bar",
            "read_summary_rows",
            "mapping_bar",
            "mapping_rows",
            "saturation_bar",
        ]
        sections = {name: f"SECTION_{name}" for name in section_names}
        payload = {
            "readSummary": {"labels": [], "counts": []},
            "readQc": {},
            "saturation": {},
            "barcodeRank3p": {"ranks": [], "counts": [], "threshold": None},
            "barcodeRank5p": {"ranks": [], "counts": [], "threshold": None},
            "beads": {"labels": [], "counts": []},
            "rnaCluster": [],
        }

        rendered = self.build_report.build_html(
            SimpleNamespace(sample_id="TEST_SAMPLE"), payload, sections
        )

        self.assertIn("TEST_SAMPLE", rendered)
        self.assertNotIn("__SAMPLE__", rendered)
        self.assertNotIn("__REPORT_BODY__", rendered)
        self.assertNotIn("P2026042903", rendered)
        self.assertNotIn("VDJ-T", rendered)
        self.assertNotIn("VDJ-B", rendered)
        for marker in sections.values():
            self.assertIn(marker, rendered)
        for plot_id in [
            "read-qc-quality",
            "read-qc-length",
            "read-qc-yield",
            "barcode-rank",
            "beads-per-droplet",
            "rna-cluster-assignment",
            "rna-umi-counts",
            "read-assignment-plot",
            "saturation-genes",
            "saturation-umis",
            "saturation-rate",
        ]:
            self.assertIn(f'id="{plot_id}"', rendered)

        violin_markup = self.build_report.summary_violin_cards([
            {"title": "Reads per cell", "plot_id": "violin-reads"},
            {"title": "UMIs per cell", "plot_id": "violin-umis"},
            {"title": "Genes per cell", "plot_id": "violin-genes"},
        ])
        for plot_id in ["violin-reads", "violin-umis", "violin-genes"]:
            self.assertIn(f'id="{plot_id}"', violin_markup)

    def test_summary_details_are_on_cells_tab_before_beads(self):
        section_names = [
            "sample_title",
            "summary_cards",
            "summary_bar",
            "summary_rows",
            "summary_violin_plots",
            "read_qc_bar",
            "beads_bar",
            "rna_cluster_bar",
            "sequencing_bar",
            "read_summary_rows",
            "mapping_bar",
            "mapping_rows",
            "saturation_bar",
        ]
        sections = {name: f"SECTION_{name}" for name in section_names}

        markup = self.build_report.new_report_markup(sections)
        summary_tab = markup.split('id="cells-tab"', 1)[0]
        cells_tab = markup.split('id="cells-tab"', 1)[1].split('id="library-tab"', 1)[0]

        self.assertNotIn("SECTION_summary_bar", summary_tab)
        self.assertNotIn("SECTION_summary_rows", summary_tab)
        self.assertNotIn("SECTION_summary_violin_plots", summary_tab)
        self.assertLess(cells_tab.index("SECTION_summary_bar"), cells_tab.index("SECTION_beads_bar"))
        self.assertLess(
            cells_tab.index("SECTION_beads_bar"),
            cells_tab.index("SECTION_rna_cluster_bar"),
        )

    def test_embedded_payload_is_safe_inside_script_element(self):
        section_names = [
            "sample_title",
            "summary_cards",
            "summary_bar",
            "summary_rows",
            "summary_violin_plots",
            "read_qc_bar",
            "beads_bar",
            "rna_cluster_bar",
            "sequencing_bar",
            "read_summary_rows",
            "mapping_bar",
            "mapping_rows",
            "saturation_bar",
        ]
        sections = {name: "" for name in section_names}
        rendered = self.build_report.build_html(
            SimpleNamespace(sample_id="TEST"),
            {
                "readSummary": {"labels": [], "counts": []},
                "readQc": {},
                "saturation": {},
                "barcodeRank3p": {},
                "barcodeRank5p": {},
                "beads": {"x": [], "y": []},
                "rnaCluster": [{"cell": "</script>&", "UMAP_1": 0, "UMAP_2": 0}],
            },
            sections,
        )

        self.assertNotIn('"</script>&"', rendered)
        self.assertIn("\\u003c/script\\u003e\\u0026", rendered)

    def test_plots_are_container_sized_without_fixed_pixel_widths(self):
        source = BUILD_REPORT_PATH.read_text(encoding="utf-8")

        self.assertNotRegex(source, r"width:\s*(360|450),")
        self.assertIn("overflow: hidden", source)
        self.assertIn("Plotly.Plots.resize", source)

    def test_sequencing_table_has_fixed_columns_and_single_line_values(self):
        table = self.build_report.pbmc_metric_rows_3col(
            [{"Metric": "Full length", "Read count": "123,456", "Percent": "98.7%"}]
        )

        self.assertIn('class="stats-table sequencing-table"', table)
        self.assertIn("<colgroup>", table)
        self.assertIn('class="metric-value read-count-value"', table)
        self.assertIn('class="metric-value percent-value"', table)
        source = BUILD_REPORT_PATH.read_text(encoding="utf-8")
        self.assertRegex(source, r"\.stats-table td\.metric-value\s*\{[^}]*display:\s*table-cell")

    def test_read_summary_uses_raw_reads_as_denominator_with_glycine(self):
        report_df = pd.DataFrame(
            [
                ("Experiment summary", "Input reads", 1000),
                ("Read assignment summary", "Full length", 600),
                ("Read assignment summary", "Barcode-valid", 500),
                ("Read assignment summary", "Cell-assigned", 400),
                ("Read assignment summary", "Gene assigned", 300),
                ("Read assignment summary", "Transcript assigned", 200),
            ],
            columns=["Section", "Metric", "Value"],
        )
        with tempfile.TemporaryDirectory() as tmp:
            glycine_stats = Path(tmp) / "identifying_statistics.txt"
            glycine_stats.write_text(
                "Summary\n"
                "Total_base_count\tValid_base_count\tValid_base_proportion(%)\n"
                "100000\t90000\t90.00\n"
                "Type\tRead_count\tRead_proportion(%)\n"
                "Total\t1000\t100.00\n"
                "Length-filtered\t100\t10.00\n"
                "QC-filtered\t50\t5.00\n"
                "Full-length+rescued\t600\t60.00\n"
                "Full-length\t550\t55.00\n\n"
                "Non-chimeric\n"
                "Type\tRead_count\tRead_proportion(%)\n"
                "Total\t700\t100.00\n"
                "Full-length\t400\t57.14\n",
                encoding="utf-8",
            )
            rows, denominator = self.build_report.build_read_summary(
                report_df, pd.DataFrame(), False, glycine_stats
            )

        self.assertEqual(denominator, 1000)
        self.assertEqual(
            rows,
            [
                {"Metric": "Raw reads", "Read count": "1,000", "Percent": "100.00%"},
                {"Metric": "Clean reads", "Read count": "850", "Percent": "85.00%"},
                {"Metric": "Full length", "Read count": "600", "Percent": "60.00%"},
                {"Metric": "Barcode-valid", "Read count": "500", "Percent": "50.00%"},
                {"Metric": "Cell-assigned", "Read count": "400", "Percent": "40.00%"},
                {"Metric": "Gene assigned", "Read count": "300", "Percent": "30.00%"},
                {"Metric": "Transcript assigned", "Read count": "200", "Percent": "20.00%"},
            ],
        )

    def test_clean_reads_uses_raw_fastq_count_when_it_differs_from_glycine_total(self):
        report_df = pd.DataFrame(
            [
                ("Experiment summary", "Input reads", 1010),
                ("Read assignment summary", "Full length", 600),
            ],
            columns=["Section", "Metric", "Value"],
        )
        with tempfile.TemporaryDirectory() as tmp:
            glycine_stats = Path(tmp) / "sample.identifying_statistic.txt"
            glycine_stats.write_text(
                "Summary\n"
                "Type\tRead_count\tRead_proportion(%)\n"
                "Total\t1000\t100.00\n"
                "Length-filtered\t100\t10.00\n"
                "QC-filtered\t50\t5.00\n"
                "Full-length\t600\t60.00\n",
                encoding="utf-8",
            )
            rows, denominator = self.build_report.build_read_summary(
                report_df, pd.DataFrame(), False, glycine_stats
            )

        self.assertEqual(denominator, 1010)
        self.assertEqual(rows[0], {"Metric": "Raw reads", "Read count": "1,010", "Percent": "100.00%"})
        self.assertEqual(rows[1], {"Metric": "Clean reads", "Read count": "860", "Percent": "85.15%"})

    def test_non_skipped_glycine_rejects_malformed_statistics(self):
        report_df = pd.DataFrame(
            [
                ("Experiment summary", "Input reads", 1000),
                ("Read assignment summary", "Full length", 600),
            ],
            columns=["Section", "Metric", "Value"],
        )
        with tempfile.TemporaryDirectory() as tmp:
            glycine_stats = Path(tmp) / "sample.identifying_statistic.txt"
            glycine_stats.write_text("Summary\nnot a valid table\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Glycine Summary"):
                self.build_report.build_read_summary(
                    report_df, pd.DataFrame(), False, glycine_stats
                )

    def test_read_summary_raw_clean_and_full_length_match_when_skipping_glycine(self):
        report_df = pd.DataFrame(
            [
                ("Experiment summary", "Input reads", 600),
                ("Read assignment summary", "Full length", 600),
                ("Read assignment summary", "Barcode-valid", 450),
                ("Read assignment summary", "Cell-assigned", 300),
                ("Read assignment summary", "Gene assigned", 240),
                ("Read assignment summary", "Transcript assigned", 180),
            ],
            columns=["Section", "Metric", "Value"],
        )

        rows, denominator = self.build_report.build_read_summary(
            report_df, pd.DataFrame(), True, None
        )

        self.assertEqual(denominator, 600)
        self.assertEqual(rows[:3], [
            {"Metric": "Raw reads", "Read count": "600", "Percent": "100.00%"},
            {"Metric": "Clean reads", "Read count": "600", "Percent": "100.00%"},
            {"Metric": "Full length", "Read count": "600", "Percent": "100.00%"},
        ])

    def test_cli_populates_report_from_pipeline_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            report_metrics = tmp_path / "report_metrics.tsv"
            qc_metrics = tmp_path / "qc_metrics.tsv"
            saturation = tmp_path / "saturation.tsv"
            read_qc = tmp_path / "read_qc.json"
            parameters = tmp_path / "parameters.tsv"
            per_cell = tmp_path / "per_cell.tsv"
            rna_cluster = tmp_path / "rna_cluster.tsv"
            barnyard_summary, barnyard_per_cell = self.write_barnyard_inputs(tmp_path)
            output = tmp_path / "report.html"

            metrics = [
                ("Summary", "Estimated cells", 123),
                ("Summary", "Input reads", 1000),
                ("Summary", "Reads per cell (mean)", 8.13),
                ("Summary", "UMIs per cell (median)", 10),
                ("Summary", "Genes per cell (median)", 7),
                ("Summary", "Unique genes", 99),
                ("Summary", "Unique isoforms", 44),
                ("Summary", "Aligned BAM reads", 800),
                ("Summary", "Unmapped", 200),
                ("Read assignment summary", "Full length", 1000),
                ("Read assignment summary", "Barcode-valid", 900),
                ("Read assignment summary", "Cell assigned", 700),
                ("Read assignment summary", "Gene assigned", 600),
                ("Read assignment summary", "Transcript assigned", 500),
            ]
            report_metrics.write_text(
                "Section\tMetric\tValue\n"
                + "".join(f"{section}\t{metric}\t{value}\n" for section, metric, value in metrics),
                encoding="utf-8",
            )
            qc_metrics.write_text("Metric\tValue\n", encoding="utf-8")
            saturation.write_text(
                "reads_per_cell\tgenes_per_cell\tumis_per_cell\tsaturation\n"
                "100\t5\t8\t0.1\n200\t7\t10\t0.2\n",
                encoding="utf-8",
            )
            read_qc.write_text(json.dumps({}), encoding="utf-8")
            parameters.write_text(
                "Parameter\tValue\nsample_id\tPIPELINE_SAMPLE\n",
                encoding="utf-8",
            )
            per_cell.write_text(
                "reads\tumis\tgenes\tmito_percent\n10\t4\t3\t1.0\n20\t7\t5\t2.0\n",
                encoding="utf-8",
            )
            rna_cluster.write_text(
                "cell\tUMAP_1\tUMAP_2\tleiden\ttotal_counts\tstatus\n"
                "CELL_A\t1.5\t-2\t0\t11\tscanpy\n"
                "CELL_ZERO\t0\t0\tunassigned\t0\tunassigned\n",
                encoding="utf-8",
            )

            command = [
                sys.executable,
                str(BUILD_REPORT_PATH),
                "--sample-id",
                "PIPELINE_SAMPLE",
                "--output-html",
                str(output),
                "--report-metrics-tsv",
                str(report_metrics),
                "--rna-qc-metrics-tsv",
                str(qc_metrics),
                "--saturation-tsv",
                str(saturation),
                "--read-qc-json",
                str(read_qc),
                "--parameters-tsv",
                str(parameters),
                "--per-cell-qc-tsv",
                str(per_cell),
                "--rna-cluster-tsv",
                str(rna_cluster),
                "--barnyard-summary-tsv",
                str(barnyard_summary),
                "--barnyard-per-cell-tsv",
                str(barnyard_per_cell),
                "--skip-glycine",
            ]
            subprocess.run(command, check=True, capture_output=True, text=True)
            rendered = output.read_text(encoding="utf-8")

            self.assertIn("PIPELINE_SAMPLE", rendered)
            self.assertIn(">123<", rendered)
            self.assertIn(">Raw reads<", rendered)
            self.assertIn('"reads_per_cell": [100, 200]', rendered)
            self.assertIn('"cell": "CELL_A"', rendered)
            self.assertIn('"status": "unassigned"', rendered)
            self.assertIn('id="rna-cluster-assignment"', rendered)
            self.assertIn('id="rna-umi-counts"', rendered)
            self.assertIn('id="barnyard-umi"', rendered)
            self.assertIn('"barnyard": {', rendered)
            self.assertIn('plotIf("barnyard-umi"', rendered)
            self.assertNotIn("P2026042903", rendered)
            self.assertNotIn("16,682", rendered)
            self.assertNotIn("VDJ-T", rendered)
            self.assertNotIn("VDJ-B", rendered)
            self.assertNotIn("Sample information", rendered)
            self.assertRegex(rendered, r"<h2>\s*Summary\s*</h2>")


if __name__ == "__main__":
    unittest.main()
