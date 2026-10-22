import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


class MatrixTests(unittest.TestCase):
    def test_filters_are_simultaneous_then_normalized_log10(self):
        script = ROOT / "scripts/ont_matrices.py"
        self.assertTrue(script.exists(), "ONT sparse matrix exporter is missing")
        spec = importlib.util.spec_from_file_location("ont_matrices", script)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        import numpy as np
        from scipy import sparse
        raw = sparse.csr_matrix([[1, 1, 0], [1, 0, 0], [0, 1, 5]])
        matrix, features, cells = mod.process_matrix(raw, ["G1", "G2", "MT-G"], ["C1", "C2", "C3"],
                                                     min_features=2, min_cells=2, max_mito=40, norm_count=100)
        # Initial masks retain G1,MT-G and C1,C2. C2 then fails mito filtering.
        self.assertEqual(features, ["G1", "MT-G"])
        self.assertEqual(cells, ["C1"])
        np.testing.assert_allclose(matrix.toarray(), [[np.log10(101)], [0]])


if __name__ == "__main__":
    unittest.main()
