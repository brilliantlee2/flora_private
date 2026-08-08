import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

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
                "Read_count: 1000\nLength-filtered: 100\nQC-filtered: 50\n",
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
            self.assertNotIn("P2026042903", rendered)
            self.assertNotIn("16,682", rendered)
            self.assertNotIn("VDJ-T", rendered)
            self.assertNotIn("VDJ-B", rendered)
            self.assertNotIn("Sample information", rendered)
            self.assertRegex(rendered, r"<h2>\s*Summary\s*</h2>")


if __name__ == "__main__":
    unittest.main()
