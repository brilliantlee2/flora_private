import importlib.util
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy import sparse


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = PROJECT_ROOT / "scripts" / "rna_cluster_analysis.py"


def load_module():
    spec = importlib.util.spec_from_file_location("strint_rna_cluster", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class MatrixFixtureMixin:
    def write_matrix(self, directory, text):
        path = Path(directory) / "gene_expression.tsv"
        path.write_text(text, encoding="utf-8")
        return path


class SparseLoaderTests(MatrixFixtureMixin, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = load_module()

    def test_loader_preserves_orientation_order_totals_and_sparse_type(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self.write_matrix(
                tmp,
                "gene\tCELL_B\tCELL_A\n"
                "G1\t3\t5\n"
                "ZERO\t0\t0\n"
                "G2\t0\t1\n",
            )
            loaded = self.module.load_expression_matrix(path, chunk_size=1)

        self.assertEqual(loaded.cells, ["CELL_B", "CELL_A"])
        self.assertEqual(loaded.genes, ["G1", "G2"])
        self.assertTrue(sparse.isspmatrix_csr(loaded.matrix))
        np.testing.assert_array_equal(loaded.matrix.toarray(), [[3, 0], [5, 1]])
        np.testing.assert_array_equal(
            loaded.total_counts, np.array([3, 6], dtype=np.uint64)
        )

    def test_uint64_totals_do_not_overflow_at_uint32_boundary(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self.write_matrix(
                tmp,
                "gene\tCELL_A\n"
                "G1\t4294967295\n"
                "G2\t1\n",
            )
            loaded = self.module.load_expression_matrix(path, chunk_size=1)

        self.assertEqual(int(loaded.total_counts[0]), 4294967296)

    def test_checked_uint64_addition_rejects_overflow(self):
        current = np.array([np.iinfo(np.uint64).max], dtype=np.uint64)
        increment = np.array([1], dtype=np.uint64)

        with self.assertRaisesRegex(ValueError, "uint64"):
            self.module.checked_add_uint64(current, increment)

    def test_invalid_matrix_values_fail_clearly(self):
        cases = {
            "wrong first column": "feature\tCELL_A\nG1\t1\n",
            "no cell columns": "gene\nG1\n",
            "blank cell": "gene\t\nG1\t1\n",
            "duplicate cell": "gene\tCELL_A\tCELL_A\nG1\t1\t2\n",
            "blank gene": "gene\tCELL_A\n\t1\n",
            "nonnumeric": "gene\tCELL_A\nG1\tx\n",
            "fractional": "gene\tCELL_A\nG1\t1.5\n",
            "negative": "gene\tCELL_A\nG1\t-1\n",
            "nan": "gene\tCELL_A\nG1\tNaN\n",
            "inf": "gene\tCELL_A\nG1\tInf\n",
            "uint32 overflow": "gene\tCELL_A\nG1\t4294967296\n",
        }
        for label, text in cases.items():
            with self.subTest(label=label), tempfile.TemporaryDirectory() as tmp:
                path = self.write_matrix(tmp, text)
                with self.assertRaises((ValueError, OverflowError)):
                    self.module.load_expression_matrix(path)


class FallbackAnalysisTests(MatrixFixtureMixin, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = load_module()

    def analyze(self, text):
        with tempfile.TemporaryDirectory() as tmp:
            loaded = self.module.load_expression_matrix(
                self.write_matrix(tmp, text), chunk_size=1
            )
            return self.module.analyze_expression(loaded)

    def test_small_inputs_preserve_every_cell_deterministically(self):
        for count in [1, 2, 3]:
            cells = [f"CELL_{i}" for i in range(count)]
            header = "gene\t" + "\t".join(cells)
            values = "\t".join(str(i + 1) for i in range(count))
            text = f"{header}\nG1\t{values}\nG2\t{values}\n"
            with self.subTest(cells=count):
                first = self.analyze(text)
                second = self.analyze(text)
                self.assertEqual(first, second)
                self.assertEqual([row.cell for row in first], cells)
                self.assertEqual({row.leiden for row in first}, {"0"})
                self.assertEqual({row.status for row in first}, {"fallback"})

    def test_zero_count_cells_are_retained_as_unassigned(self):
        rows = self.analyze(
            "gene\tCELL_A\tCELL_ZERO\tCELL_B\n"
            "G1\t2\t0\t1\n"
            "G2\t1\t0\t3\n"
        )

        self.assertEqual([row.cell for row in rows], ["CELL_A", "CELL_ZERO", "CELL_B"])
        zero = rows[1]
        self.assertEqual(zero.total_counts, 0)
        self.assertEqual(zero.leiden, "unassigned")
        self.assertEqual(zero.status, "unassigned")

    def test_all_zero_cells_are_emitted_as_unassigned(self):
        rows = self.analyze("gene\tCELL_A\tCELL_B\nZERO\t0\t0\n")

        self.assertEqual([row.cell for row in rows], ["CELL_A", "CELL_B"])
        self.assertEqual({row.status for row in rows}, {"unassigned"})


@unittest.skipUnless(importlib.util.find_spec("scanpy"), "Scanpy is not installed")
class ScanpyIntegrationTests(MatrixFixtureMixin, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = load_module()

    def test_standard_path_clusters_informative_cells_and_retains_zero_cell(self):
        rng = np.random.default_rng(7)
        cells = [f"CELL_{index}" for index in range(12)] + ["CELL_ZERO"]
        counts = rng.poisson(2, size=(40, 12))
        counts[:12, :6] += 6
        counts[12:24, 6:] += 6
        counts = np.column_stack((counts, np.zeros(40, dtype=int)))
        text = (
            "gene\t"
            + "\t".join(cells)
            + "\n"
            + "\n".join(
                f"G{index}\t" + "\t".join(str(value) for value in row)
                for index, row in enumerate(counts)
            )
            + "\n"
        )

        with tempfile.TemporaryDirectory() as tmp:
            loaded = self.module.load_expression_matrix(
                self.write_matrix(tmp, text), chunk_size=7
            )
            rows = self.module.analyze_expression(loaded)

        self.assertEqual([row.cell for row in rows], cells)
        self.assertEqual({row.status for row in rows[:-1]}, {"scanpy"})
        self.assertGreaterEqual(len({row.leiden for row in rows[:-1]}), 2)
        self.assertEqual(rows[-1].status, "unassigned")
        self.assertEqual(rows[-1].leiden, "unassigned")


class CliTests(MatrixFixtureMixin, unittest.TestCase):
    def run_cli(self, input_path, output_path):
        return subprocess.run(
            [
                sys.executable,
                str(SCRIPT_PATH),
                "--input",
                str(input_path),
                "--output",
                str(output_path),
            ],
            text=True,
            capture_output=True,
            env=os.environ.copy(),
            check=False,
        )

    def test_cli_writes_exact_columns_and_is_deterministic(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_path = self.write_matrix(
                tmp,
                "gene\tCELL_B\tCELL_A\nG1\t2\t1\nG2\t1\t3\n",
            )
            output = Path(tmp) / "cluster.tsv"
            result = self.run_cli(input_path, output)
            self.assertEqual(result.returncode, 0, result.stderr)
            first = output.read_bytes()
            result = self.run_cli(input_path, output)
            self.assertEqual(result.returncode, 0, result.stderr)

            self.assertEqual(first, output.read_bytes())
            lines = output.read_text(encoding="utf-8").splitlines()

        self.assertEqual(
            lines[0], "cell\tUMAP_1\tUMAP_2\tleiden\ttotal_counts\tstatus"
        )
        self.assertEqual(
            [line.split("\t", 1)[0] for line in lines[1:]],
            ["CELL_B", "CELL_A"],
        )

    def test_cli_returns_nonzero_without_partial_output_for_invalid_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_path = self.write_matrix(tmp, "feature\tCELL_A\nG1\t1\n")
            output = Path(tmp) / "cluster.tsv"
            output.write_text("old output\n", encoding="utf-8")
            result = self.run_cli(input_path, output)

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("gene", result.stderr)
            self.assertEqual(output.read_text(encoding="utf-8"), "old output\n")


if __name__ == "__main__":
    unittest.main()
