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
                    'run_stage gene_expression "${DOWNSTREAM_DIR}/gene_expression.py"'
                )
                cluster_position = source.index(
                    'python3 "$(python_asset "${DOWNSTREAM_DIR}/rna_cluster_analysis.py")"'
                )
                report_position = source.index(
                    'run_stage build_report "${DOWNSTREAM_DIR}/build_report.py"'
                )
                self.assertLess(gene_position, cluster_position)
                self.assertLess(cluster_position, report_position)


if __name__ == "__main__":
    unittest.main()
