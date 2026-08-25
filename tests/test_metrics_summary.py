import importlib.util
import json
import tempfile
import unittest
import zipfile
from pathlib import Path
from xml.etree import ElementTree


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = PROJECT_ROOT / "scripts" / "metrics_summary.py"


def load_module():
    spec = importlib.util.spec_from_file_location("metrics_summary", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def workbook_rows(path):
    namespace = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with zipfile.ZipFile(path) as archive:
        required = {
            "[Content_Types].xml",
            "_rels/.rels",
            "xl/workbook.xml",
            "xl/_rels/workbook.xml.rels",
            "xl/styles.xml",
            "xl/worksheets/sheet1.xml",
        }
        assert required <= set(archive.namelist())
        root = ElementTree.fromstring(archive.read("xl/worksheets/sheet1.xml"))
    rows = []
    for row in root.findall(".//x:sheetData/x:row", namespace):
        values = []
        for cell in row.findall("x:c", namespace):
            if cell.attrib.get("t") == "inlineStr":
                value = cell.findtext("x:is/x:t", default="", namespaces=namespace)
            else:
                value = cell.findtext("x:v", default="", namespaces=namespace)
            values.append(value)
        rows.append(values)
    return rows


def workbook_styles(path):
    with zipfile.ZipFile(path) as archive:
        return archive.read("xl/styles.xml").decode("utf-8")


class MetricsSummaryTests(unittest.TestCase):
    def test_writes_requested_metrics_with_glycine_clean_read_definition(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            read_qc = root / "read_qc.json"
            read_qc.write_text(
                json.dumps(
                    {
                        "quality": {"mean": 18.25, "median": 19.0},
                        "length": {"mean": 1234.5, "median": 1200.0},
                        "yield_above_length": {"n50_bp": 1500},
                    }
                ),
                encoding="utf-8",
            )
            rna_qc = root / "rna_qc.tsv"
            rna_qc.write_text(
                "Metric\tValue\n"
                "Input reads\t1,000\n"
                "Full length reads\t600\n"
                "Estimated number of cells\t2\n"
                "Mean reads per cell\t250.00\n"
                "Mean UMI counts per cell\t120.00\n"
                "Median UMI counts per cell\t110.00\n"
                "Total genes detected\t99\n"
                "Mean Genes per cell\t50.50\n"
                "Median Genes per cell\t48.00\n"
                "Fraction reads in cells\t0.500000\n"
                "Barcode-valid reads\t550\n"
                "Reads assigned to final cells\t500\n"
                "Gene assigned reads\t450\n"
                "Transcript assigned reads\t400\n"
                "Aligned BAM reads\t800\n"
                "Unmapped\t200\n"
                "Unique genes\t99\n"
                "Unique isoforms\t44\n",
                encoding="utf-8",
            )
            saturation = root / "saturation.tsv"
            saturation.write_text(
                "fraction\treads\treads_per_cell\tgenes_per_cell\tumis_per_cell\tsaturation\n"
                "0.5\t300\t150\t40\t90\t0.25\n"
                "1.0\t600\t300\t48\t110\t0.40\n",
                encoding="utf-8",
            )
            glycine = root / "identifying_statistic.txt"
            glycine.write_text(
                "Summary\n"
                "Total_base_count\tValid_base_count\tValid_base_proportion(%)\n"
                "100\t50\t50.00\n"
                "Type\tRead_count\tRead_proportion(%)\n"
                "Total\t1000\t100.00\n"
                "Length-filtered\t100\t10.00\n"
                "QC-filtered\t50\t5.00\n"
                "Full-length+rescued\t600\t60.00\n",
                encoding="utf-8",
            )
            output = root / "metrics_summary.xlsx"

            module.build_metrics_workbook(
                sample_id="sample-A",
                species="GRCH38",
                read_qc_json=read_qc,
                rna_qc_metrics_tsv=rna_qc,
                saturation_tsv=saturation,
                output_xlsx=output,
                glycine_stats=glycine,
                skip_glycine=False,
            )

            rows = workbook_rows(output)
            self.assertEqual(len(rows), 2)
            values = dict(zip(rows[0], rows[1]))
            self.assertEqual(values["SampleName"], "sample-A")
            self.assertEqual(values["species"], "GRCH38")
            self.assertEqual(values["Mean read quality"], "18.25")
            self.assertEqual(values["N50(b)"], "1500")
            self.assertEqual(values["Mean reads per cell"], "500.0")
            self.assertEqual(values["Clean reads"], "850")
            self.assertEqual(values["Sequencing saturation"], "0.4")
            self.assertEqual(values["Aligned BAM reads / total reads"], "0.8")
            self.assertEqual(values["Unmapped / total reads"], "0.2")
            self.assertNotIn('patternType="solid"', workbook_styles(output))

    def test_skip_glycine_sets_clean_reads_to_full_length(self):
        module = load_module()
        metrics = {"Input reads": 600, "Full length reads": 600}
        self.assertEqual(module.clean_reads(metrics, None, skip_glycine=True), 600)


if __name__ == "__main__":
    unittest.main()
