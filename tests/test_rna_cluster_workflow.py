import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


class RnaClusterWorkflowTests(unittest.TestCase):
    def test_both_entrypoints_run_clustering_and_pass_result_to_report(self):
        for script_name in ["run_all.sh", "run_all_mixed_species.sh"]:
            with self.subTest(script=script_name):
                source = (PROJECT_ROOT / script_name).read_text(encoding="utf-8")
                self.assertIn("rna_cluster_analysis.py", source)
                self.assertIn(
                    '--input "${SAMPLE_ID}.gene_expression.tsv"', source
                )
                self.assertIn(
                    '--output "${SAMPLE_ID}.rna_cluster.tsv"', source
                )
                self.assertIn(
                    '--rna-cluster-tsv "${MATRIX_DIR}/${SAMPLE_ID}.rna_cluster.tsv"',
                    source,
                )
                self.assertIn(
                    '${MATRIX_DIR}/${SAMPLE_ID}.rna_cluster.tsv',
                    source,
                )

                gene_position = source.index(
                    '--gene-expression-output "${SAMPLE_ID}.gene_expression.tsv"'
                )
                self.assertNotIn(
                    'run_stage gene_expression "${DOWNSTREAM_DIR}/gene_expression.py"',
                    source,
                )
                cluster_position = source.index(
                    'python3 "$(python_asset "${DOWNSTREAM_DIR}/rna_cluster_analysis.py")"'
                )
                report_position = source.index(
                    'run_stage build_report "${DOWNSTREAM_DIR}/build_report.py"'
                )
                self.assertLess(gene_position, cluster_position)
                self.assertLess(cluster_position, report_position)

    def test_mixed_entrypoint_defines_barnyard_report_outputs_after_qc(self):
        source = (PROJECT_ROOT / "run_all_mixed_species.sh").read_text(
            encoding="utf-8"
        )

        barnyard_qc_position = source.index(
            'run_stage barnyard_qc "${DOWNSTREAM_DIR}/barnyard_qc.py"'
        )
        summary_definition = (
            'BARNYARD_SUMMARY_TSV="barnyard_qc/barnyard_summary.tsv"'
        )
        per_cell_definition = (
            'BARNYARD_PER_CELL_TSV="barnyard_qc/barnyard_per_cell.tsv"'
        )

        self.assertIn(summary_definition, source)
        self.assertIn(per_cell_definition, source)
        summary_position = source.index(summary_definition)
        per_cell_position = source.index(per_cell_definition)

        self.assertLess(barnyard_qc_position, summary_position)
        self.assertLess(summary_position, per_cell_position)

    def test_mixed_entrypoint_forwards_complete_barnyard_pair_nonfatally(self):
        source = (PROJECT_ROOT / "run_all_mixed_species.sh").read_text(
            encoding="utf-8"
        )
        expected_block = '''if [[ -f "${BARNYARD_SUMMARY_TSV}" && -f "${BARNYARD_PER_CELL_TSV}" ]]; then
  BUILD_REPORT_ARGS+=(
    --barnyard-summary-tsv "${BARNYARD_SUMMARY_TSV}"
    --barnyard-per-cell-tsv "${BARNYARD_PER_CELL_TSV}"
  )
else
  log "WARNING: Barnyard report inputs are incomplete; omitting Barnyard QC from the HTML report."
fi'''

        self.assertIn(expected_block, source)
        self.assertEqual(source.count("--barnyard-summary-tsv"), 1)
        self.assertEqual(source.count("--barnyard-per-cell-tsv"), 1)
        self.assertNotIn('require_file "${BARNYARD_SUMMARY_TSV}"', source)
        self.assertNotIn('require_file "${BARNYARD_PER_CELL_TSV}"', source)

    def test_single_species_entrypoint_has_no_barnyard_report_arguments(self):
        source = (PROJECT_ROOT / "run_all.sh").read_text(encoding="utf-8")

        self.assertNotIn("--barnyard-summary-tsv", source)
        self.assertNotIn("--barnyard-per-cell-tsv", source)


if __name__ == "__main__":
    unittest.main()
