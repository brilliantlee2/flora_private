"""Decoded BAM/matrix parity against an exhaustive directional-UMI oracle."""

from collections import Counter, deque
import csv
import importlib.util
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


BIN_DIR = Path(os.environ.get("FLORA_TEST_BIN_DIR", "target/debug"))
AVAILABLE = importlib.util.find_spec("pysam") is not None and all(
    (BIN_DIR / name).is_file() for name in ("cluster_umis_allbam", "cell_umi_gene_table")
)


def distance(left, right):
    row = list(range(len(right) + 1))
    for i, a in enumerate(left, 1):
        following = [i]
        for j, b in enumerate(right, 1):
            following.append(min(following[-1] + 1, row[j] + 1, row[j - 1] + (a != b)))
        row = following
    return row[-1]


def exhaustive_directional(counts):
    graph = {umi: [] for umi in counts}
    for a in counts:
        for b in counts:
            if a != b and counts[a] >= max(2 * counts[b] - 1, 0) and distance(a, b) <= 3:
                graph[a].append(b)
    corrected = {}
    for root in sorted(counts, key=lambda umi: (-counts[umi], umi)):
        if root in corrected:
            continue
        corrected[root] = root
        queue = deque([root])
        while queue:
            for child in graph[queue.popleft()]:
                if child not in corrected:
                    corrected[child] = root
                    queue.append(child)
    return corrected


def dna(value, length=10):
    sequence = []
    for _ in range(length):
        sequence.append("ACGT"[value % 4])
        value //= 4
    return "".join(sequence)


@unittest.skipUnless(AVAILABLE, "requires pysam and built Rust stage binaries")
class UmiHotspotParityTests(unittest.TestCase):
    def test_constant_prefix_and_mixed_lengths_preserve_tags_and_matrices(self):
        import pysam

        groups = {
            ("CELL_A", "GENE_ONE"): {"TGCTG" + dna(i): 1 for i in range(128)},
            ("CELL_B", "GENE_ONE"): {"TGCTG" + dna(i): 2 for i in range(128)},
            ("CELL_A", "GENE_TWO"): {"TGCTG" + dna(i, 6 + i % 5): 1 + i % 7 for i in range(120)},
            ("CELL_B", "GENE_TWO"): {dna(i): 1 + i % 5 for i in range(64)},
            ("CELL_C", "GENE_THREE"): {dna(i): 1 + i % 5 for i in range(65)},
        }
        expected = {}
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source, output = root / "input.bam", root / "output.bam"
            matrix, table = root / "matrix.tsv", root / "table.tsv"
            header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 100000}]}
            with pysam.AlignmentFile(source, "wb", header=header) as handle:
                index = 0
                for (cell, gene), counts in groups.items():
                    corrected = exhaustive_directional(counts)
                    for umi, count in counts.items():
                        for _ in range(count):
                            record = pysam.AlignedSegment(handle.header)
                            record.query_name = f"read_{index:06d}"
                            record.query_sequence = "ACGT" * 25
                            record.query_qualities = pysam.qualitystring_to_array("I" * 100)
                            record.reference_id, record.reference_start = 0, index
                            record.flag, record.mapping_quality, record.cigarstring = 0, 60, "100M"
                            record.set_tags([("CB", cell), ("GN", gene), ("UR", umi)])
                            expected[record.query_name] = (record.to_string(), corrected[umi])
                            handle.write(record)
                            index += 1
            pysam.index(str(source))
            decoded = []
            for threads in (1, 4):
                subprocess.run([
                    str(BIN_DIR / "cluster_umis_allbam"), str(source), "--output", str(output),
                    "--cell_gene_max_reads", "0", "--threads", str(threads),
                ], check=True, capture_output=True, text=True)
                with pysam.AlignmentFile(output) as handle:
                    records = list(handle)
                self.assertEqual(len(records), len(expected))
                for record in records:
                    original, umi = expected[record.query_name]
                    self.assertEqual(record.get_tag("UB"), umi)
                    record.set_tag("UB", None)
                    self.assertEqual(record.to_string(), original)
                decoded.append([record.to_string() for record in records])
            self.assertEqual(decoded[0], decoded[1])
            subprocess.run([
                str(BIN_DIR / "cell_umi_gene_table"), str(output), "--output", str(table),
                "--gene-expression-output", str(matrix),
            ], check=True, capture_output=True, text=True)
            with matrix.open() as handle:
                reader = csv.reader(handle, delimiter="\t")
                cells = next(reader)[1:]
                actual = {row[0]: dict(zip(cells, map(int, row[1:]))) for row in reader}
            self.assertEqual(set(cells), {cell for cell, _ in groups})
            self.assertEqual(set(actual), {gene for _, gene in groups})
            for gene, counts_by_cell in actual.items():
                for cell, number in counts_by_cell.items():
                    counts = groups.get((cell, gene), {})
                    expected_number = len(Counter(exhaustive_directional(counts).values()))
                    self.assertEqual(number, expected_number, (cell, gene))


if __name__ == "__main__":
    unittest.main()
