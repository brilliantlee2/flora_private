#!/usr/bin/env python3
"""Build Flora's compact, dependency-free Excel metrics summary."""

import argparse
import csv
import json
import math
import zipfile
from pathlib import Path
from xml.sax.saxutils import escape


SHEET_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--read-qc-json", required=True, type=Path)
    parser.add_argument("--rna-qc-metrics-tsv", required=True, type=Path)
    parser.add_argument("--saturation-tsv", required=True, type=Path)
    parser.add_argument("--output-xlsx", required=True, type=Path)
    parser.add_argument("--glycine-stats", type=Path)
    parser.add_argument("--skip-glycine", action="store_true")
    return parser.parse_args()


def parse_number(value):
    text = str(value or "").strip().replace(",", "")
    if not text:
        return None
    is_percent = text.endswith("%")
    if is_percent:
        text = text[:-1]
    try:
        number = float(text)
    except ValueError:
        return None
    return number / 100.0 if is_percent else number


def load_rna_metrics(path):
    values = {}
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            metric = (row.get("Metric") or "").strip()
            if metric:
                values[metric] = parse_number(row.get("Value"))
    return values


def parse_glycine_counts(path):
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    start = next(
        (
            index
            for index, line in enumerate(lines)
            if line.strip().split("\t")[:2] == ["Type", "Read_count"]
        ),
        None,
    )
    if start is None:
        raise ValueError(f"Glycine statistics table not found in {path}")
    values = {}
    for line in lines[start + 1 :]:
        if not line.strip():
            break
        fields = line.split("\t")
        if len(fields) < 2:
            break
        try:
            values[fields[0].strip()] = int(fields[1].replace(",", ""))
        except ValueError:
            break
    required = {"Total", "Length-filtered", "QC-filtered"}
    missing = required - values.keys()
    if missing:
        raise ValueError(f"Glycine statistics missing rows: {', '.join(sorted(missing))}")
    return values


def clean_reads(metrics, glycine_stats, skip_glycine):
    full_length = metrics.get("Full length reads")
    if skip_glycine:
        return full_length
    if glycine_stats is None:
        raise ValueError("--glycine-stats is required unless --skip-glycine is set")
    counts = parse_glycine_counts(Path(glycine_stats))
    return max(0, counts["Total"] - counts["Length-filtered"] - counts["QC-filtered"])


def full_depth_saturation(path):
    best_fraction = -1.0
    best_value = None
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            fraction = parse_number(row.get("fraction"))
            saturation = parse_number(row.get("saturation"))
            if fraction is not None and saturation is not None and fraction >= best_fraction:
                best_fraction = fraction
                best_value = saturation
    return best_value


def ratio(numerator, denominator):
    if numerator is None or denominator in (None, 0):
        return None
    return numerator / denominator


def metric_rows(sample_id, species, read_qc, metrics, saturation, clean_count):
    quality = read_qc.get("quality", {})
    length = read_qc.get("length", {})
    yield_data = read_qc.get("yield_above_length", {})
    raw_reads = metrics.get("Input reads")
    rows = [
        ("SampleName", sample_id, "text", False),
        ("species", species, "text", False),
        ("Mean read quality", quality.get("mean"), "decimal", True),
        ("Median read quality", quality.get("median"), "decimal", True),
        ("Mean read length (b)", length.get("mean"), "decimal", True),
        ("Median read length(b)", length.get("median"), "decimal", True),
        ("N50(b)", yield_data.get("n50_bp"), "integer", True),
        ("Estimated number of cell", metrics.get("Estimated number of cells"), "integer", False),
        ("Mean reads per cell", metrics.get("Mean reads per cell"), "decimal", False),
        ("Mean UMI counts per cell", metrics.get("Mean UMI counts per cell"), "decimal", False),
        ("Median UMI counts per cell", metrics.get("Median UMI counts per cell"), "decimal", False),
        ("Total genes detected", metrics.get("Total genes detected"), "integer", False),
        ("Mean genes per cell", metrics.get("Mean Genes per cell"), "decimal", False),
        ("Median genes per cell", metrics.get("Median Genes per cell"), "decimal", False),
        ("Sequencing saturation", saturation, "percent", False),
        ("Fraction Reads in cell", metrics.get("Fraction reads in cells"), "percent", False),
        ("Raw reads", raw_reads, "integer", False),
        ("Clean reads", clean_count, "integer", False),
        ("Full length", metrics.get("Full length reads"), "integer", False),
        ("Barcode-valid", metrics.get("Barcode-valid reads"), "integer", False),
        ("Cell-assigned", metrics.get("Reads assigned to final cells"), "integer", False),
        ("Gene assigned", metrics.get("Gene assigned reads"), "integer", False),
        ("Transcript assigned", metrics.get("Transcript assigned reads"), "integer", False),
        (
            "Aligned BAM reads / total reads",
            ratio(metrics.get("Aligned BAM reads"), raw_reads),
            "percent",
            False,
        ),
        ("Unmapped / total reads", ratio(metrics.get("Unmapped"), raw_reads), "percent", False),
        ("Unique genes", metrics.get("Unique genes"), "integer", False),
        ("Unique isoforms", metrics.get("Unique isoforms"), "integer", False),
    ]
    return rows


def cell(reference, value, kind, style):
    if value is None or (isinstance(value, float) and not math.isfinite(value)):
        return f'<c r="{reference}" s="{style}" t="inlineStr"><is><t></t></is></c>'
    if kind == "text":
        return (
            f'<c r="{reference}" s="{style}" t="inlineStr"><is><t>'
            f"{escape(str(value))}</t></is></c>"
        )
    number = int(value) if kind == "integer" else float(value)
    return f'<c r="{reference}" s="{style}"><v>{number}</v></c>'


def worksheet_xml(rows):
    xml_rows = [
        '<row r="1" ht="24" customHeight="1">'
        + cell("A1", "Metric", "text", 1)
        + cell("B1", "Value", "text", 1)
        + "</row>"
    ]
    styles = {"text": 2, "integer": 4, "decimal": 5, "percent": 6}
    for row_number, (label, value, kind, highlight) in enumerate(rows, start=2):
        label_style = 3 if highlight else 2
        value_style = 7 if highlight else styles[kind]
        xml_rows.append(
            f'<row r="{row_number}">'
            + cell(f"A{row_number}", label, "text", label_style)
            + cell(f"B{row_number}", value, kind, value_style)
            + "</row>"
        )
    last_row = len(rows) + 1
    return f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<worksheet xmlns="{SHEET_NS}">
  <sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/></sheetView></sheetViews>
  <cols><col min="1" max="1" width="38" customWidth="1"/><col min="2" max="2" width="22" customWidth="1"/></cols>
  <sheetData>{''.join(xml_rows)}</sheetData>
  <autoFilter ref="A1:B{last_row}"/>
</worksheet>'''


def write_xlsx(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    parts = {
        "[Content_Types].xml": '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>
  <Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>
  <Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>
</Types>''',
        "_rels/.rels": '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
</Relationships>''',
        "xl/workbook.xml": f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="{SHEET_NS}" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
  <sheets><sheet name="Metrics Summary" sheetId="1" r:id="rId1"/></sheets>
</workbook>''',
        "xl/_rels/workbook.xml.rels": '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>
</Relationships>''',
        "xl/styles.xml": f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="{SHEET_NS}">
  <fonts count="3"><font><sz val="11"/><name val="Calibri"/></font><font><b/><color rgb="FFFFFFFF"/><sz val="12"/><name val="Calibri"/></font><font><b/><sz val="11"/><name val="Calibri"/></font></fonts>
  <fills count="4"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill><fill><patternFill patternType="solid"><fgColor rgb="FF4472C4"/><bgColor indexed="64"/></patternFill></fill><fill><patternFill patternType="solid"><fgColor rgb="FFFFC000"/><bgColor indexed="64"/></patternFill></fill></fills>
  <borders count="2"><border><left/><right/><top/><bottom/><diagonal/></border><border><left style="thin"><color rgb="FFD9D9D9"/></left><right style="thin"><color rgb="FFD9D9D9"/></right><top style="thin"><color rgb="FFD9D9D9"/></top><bottom style="thin"><color rgb="FFD9D9D9"/></bottom><diagonal/></border></borders>
  <cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>
  <cellXfs count="8">
    <xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>
    <xf numFmtId="0" fontId="1" fillId="2" borderId="1" xfId="0" applyAlignment="1"><alignment horizontal="center"/></xf>
    <xf numFmtId="0" fontId="2" fillId="0" borderId="1" xfId="0"/>
    <xf numFmtId="0" fontId="2" fillId="3" borderId="1" xfId="0"/>
    <xf numFmtId="3" fontId="0" fillId="0" borderId="1" xfId="0" applyNumberFormat="1"/>
    <xf numFmtId="4" fontId="0" fillId="0" borderId="1" xfId="0" applyNumberFormat="1"/>
    <xf numFmtId="10" fontId="0" fillId="0" borderId="1" xfId="0" applyNumberFormat="1"/>
    <xf numFmtId="4" fontId="0" fillId="3" borderId="1" xfId="0" applyNumberFormat="1"/>
  </cellXfs>
  <cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>
</styleSheet>''',
        "xl/worksheets/sheet1.xml": worksheet_xml(rows),
    }
    temporary = path.with_name(f".{path.name}.tmp")
    with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, content in parts.items():
            archive.writestr(name, content)
    temporary.replace(path)


def build_metrics_workbook(
    sample_id,
    species,
    read_qc_json,
    rna_qc_metrics_tsv,
    saturation_tsv,
    output_xlsx,
    glycine_stats=None,
    skip_glycine=False,
):
    read_qc = json.loads(Path(read_qc_json).read_text(encoding="utf-8"))
    metrics = load_rna_metrics(Path(rna_qc_metrics_tsv))
    saturation = full_depth_saturation(Path(saturation_tsv))
    clean_count = clean_reads(metrics, glycine_stats, skip_glycine)
    rows = metric_rows(sample_id, species, read_qc, metrics, saturation, clean_count)
    write_xlsx(Path(output_xlsx), rows)


def main():
    args = parse_args()
    build_metrics_workbook(
        sample_id=args.sample_id,
        species=args.species,
        read_qc_json=args.read_qc_json,
        rna_qc_metrics_tsv=args.rna_qc_metrics_tsv,
        saturation_tsv=args.saturation_tsv,
        output_xlsx=args.output_xlsx,
        glycine_stats=args.glycine_stats,
        skip_glycine=args.skip_glycine,
    )
    print(args.output_xlsx)


if __name__ == "__main__":
    main()
