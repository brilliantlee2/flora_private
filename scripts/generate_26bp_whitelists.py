#!/usr/bin/env python3
import argparse
import csv
import re
import sys
import xml.etree.ElementTree as ET
from itertools import product
from pathlib import Path
from zipfile import ZipFile


VALID_BC_RE = re.compile(r"^[ACGT]{10}$")


def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def read_plain_barcodes(path):
    barcodes = []
    with open(path, "r", encoding="utf-8-sig") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        if "," in sample or "\t" in sample:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
            reader = csv.DictReader(handle, dialect=dialect)
            if reader.fieldnames:
                seq_col = next(
                    (
                        col
                        for col in reader.fieldnames
                        if col and col.strip().lower() in {"sequence", "seq", "barcode", "cellbarcode", "序列"}
                    ),
                    None,
                )
                if seq_col:
                    for row in reader:
                        value = (row.get(seq_col) or "").strip().upper()
                        if value:
                            barcodes.append(value)
                    return barcodes
            handle.seek(0)
        for line in handle:
            value = line.strip().split()[0].upper() if line.strip() else ""
            if value and value not in {"SEQUENCE", "SEQ", "BARCODE", "序列"}:
                barcodes.append(value)
    return barcodes


def _xlsx_shared_strings(zf):
    if "xl/sharedStrings.xml" not in zf.namelist():
        return []
    root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
    ns = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    values = []
    for item in root.findall("x:si", ns):
        texts = [node.text or "" for node in item.findall(".//x:t", ns)]
        values.append("".join(texts))
    return values


def read_xlsx_barcodes(path, sequence_column_name="序列"):
    with ZipFile(path) as zf:
        shared_strings = _xlsx_shared_strings(zf)
        sheet_name = "xl/worksheets/sheet1.xml"
        root = ET.fromstring(zf.read(sheet_name))
    ns = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}

    rows = []
    for row in root.findall(".//x:sheetData/x:row", ns):
        values = {}
        for cell in row.findall("x:c", ns):
            ref = cell.attrib.get("r", "")
            col = re.sub(r"\d+", "", ref)
            value_node = cell.find("x:v", ns)
            if value_node is None:
                value = ""
            elif cell.attrib.get("t") == "s":
                value = shared_strings[int(value_node.text)]
            else:
                value = value_node.text or ""
            values[col] = str(value).strip()
        rows.append(values)

    if not rows:
        return []
    header = rows[0]
    seq_col = next((col for col, name in header.items() if name.strip() == sequence_column_name), None)
    if seq_col is None:
        seq_col = "B"

    return [row.get(seq_col, "").strip().upper() for row in rows[1:] if row.get(seq_col, "").strip()]


def read_10bp_barcodes(path):
    path = Path(path)
    if path.suffix.lower() == ".xlsx":
        barcodes = read_xlsx_barcodes(path)
    else:
        barcodes = read_plain_barcodes(path)

    bad_len = [bc for bc in barcodes if len(bc) != 10]
    bad_base = [bc for bc in barcodes if len(bc) == 10 and not VALID_BC_RE.match(bc)]
    if bad_len:
        raise ValueError(f"Found non-10bp barcode examples: {bad_len[:10]}")
    if bad_base:
        raise ValueError(f"Found barcode with non-ACGT bases examples: {bad_base[:10]}")
    if len(barcodes) != len(set(barcodes)):
        duplicates = [bc for bc, n in csv_counter(barcodes).items() if n > 1]
        raise ValueError(f"Found duplicate barcodes examples: {duplicates[:10]}")
    return barcodes


def csv_counter(values):
    counts = {}
    for value in values:
        counts[value] = counts.get(value, 0) + 1
    return counts


def write_whitelists(barcodes, middle_3p, middle_5p, out_3p, out_5p):
    barcodes_5p = [reverse_complement(bc) for bc in barcodes]
    if len(barcodes_5p) != len(set(barcodes_5p)):
        raise ValueError("5p barcode reverse-complement list has duplicates.")

    count_3p = 0
    with open(out_3p, "w", encoding="utf-8") as handle:
        for left, right in product(barcodes, repeat=2):
            handle.write(f"{left}{middle_3p}{right}\n")
            count_3p += 1

    count_5p = 0
    with open(out_5p, "w", encoding="utf-8") as handle:
        for left, right in product(barcodes_5p, repeat=2):
            handle.write(f"{left}{middle_5p}{right}\n")
            count_5p += 1

    return count_3p, count_5p


def main():
    parser = argparse.ArgumentParser(description="Expand a 10bp barcode list to Strint2 26bp 3p/5p whitelists.")
    parser.add_argument("--barcode-list-10bp", required=True)
    parser.add_argument("--out-3p", required=True)
    parser.add_argument("--out-5p", required=True)
    parser.add_argument("--middle-3p", default="GGTAGC")
    parser.add_argument("--middle-5p", default="GGAAGG")
    args = parser.parse_args()

    barcodes = read_10bp_barcodes(args.barcode_list_10bp)
    Path(args.out_3p).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_5p).parent.mkdir(parents=True, exist_ok=True)
    count_3p, count_5p = write_whitelists(
        barcodes=barcodes,
        middle_3p=args.middle_3p,
        middle_5p=args.middle_5p,
        out_3p=args.out_3p,
        out_5p=args.out_5p,
    )

    print(f"Input 10bp barcode count: {len(barcodes)}")
    print(f"3p whitelist: {args.out_3p} ({count_3p} records)")
    print(f"5p whitelist: {args.out_5p} ({count_5p} records)")
    print(f"First 3 raw barcodes: {barcodes[:3]}")
    print(f"First 3 5p RC barcodes: {[reverse_complement(bc) for bc in barcodes[:3]]}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        sys.exit(1)
