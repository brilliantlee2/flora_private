import importlib.util
import importlib
import sys
import tempfile
import types
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
GENE_EXPRESSION_PATH = PROJECT_ROOT / "scripts" / "gene_expression.py"


def load_gene_expression_module():
    injected = []
    try:
        try:
            importlib.import_module("pysam")
        except ImportError:
            sys.modules["pysam"] = types.ModuleType("pysam")
            injected.append("pysam")
        try:
            importlib.import_module("tqdm")
        except ImportError:
            tqdm_module = types.ModuleType("tqdm")
            tqdm_module.tqdm = lambda iterable, **_kwargs: iterable
            sys.modules["tqdm"] = tqdm_module
            injected.append("tqdm")

        spec = importlib.util.spec_from_file_location(
            "strint_gene_expression", GENE_EXPRESSION_PATH
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        for name in injected:
            sys.modules.pop(name, None)


class FakeAlignment:
    def __init__(self, gene, cell, umi):
        self.tags = {"GN": gene, "CB": cell, "UB": umi}

    def get_tag(self, tag):
        return self.tags[tag]


class FakeBam:
    def __init__(self, alignments):
        self.alignments = alignments

    def fetch(self):
        return iter(self.alignments)


class GeneExpressionCellParityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = load_gene_expression_module()

    def test_python_writer_preserves_placeholder_only_cell(self):
        bam = FakeBam(
            [
                FakeAlignment("ACTB", "CELL_A", "UMI_A"),
                FakeAlignment("chr1_1000_2000", "CELL_B", "UMI_B"),
            ]
        )

        genes, cells, umi_sets = self.module.read_bam_entries(bam, 2)

        self.assertEqual(genes, {"ACTB"})
        self.assertEqual(cells, {"CELL_A", "CELL_B"})
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "gene_expression.tsv"
            self.module.write_matrix(output, genes, cells, umi_sets)
            lines = output.read_text(encoding="utf-8").splitlines()

        self.assertEqual(lines[0], "gene\tCELL_A\tCELL_B")
        self.assertEqual(lines[1], "ACTB\t1\t0")


if __name__ == "__main__":
    unittest.main()
