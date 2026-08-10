#!/usr/bin/env python3
import argparse
import base64
import csv
import hashlib
import heapq
import html
import io
import json
import math
import os
import re
import sqlite3
import sys
import tempfile
from collections import Counter, defaultdict
from decimal import Decimal, InvalidOperation
from pathlib import Path

import pandas as pd


BARCODE_RANK_MAX_POINTS = 20000
BEAD_CHUNK_ROWS = 1_000_000
BARNYARD_ASSIGNMENTS = (
    "human_singlet",
    "mouse_singlet",
    "mixed",
    "unclassified",
)
BARNYARD_REQUIRED_COLUMNS = (
    "cell_id",
    "assignment",
    "umi_human",
    "umi_mouse",
    "umi_total",
)
BARNYARD_MAX_POINTS = 100_000
BARNYARD_MAX_POINTS_PER_CLASS = 25_000
BARNYARD_CHUNK_ROWS = 100_000
BARNYARD_JS_MAX_SAFE_INTEGER = 9_007_199_254_740_991
_RUNTIME_INTEGER_DIGIT_LIMIT = sys.get_int_max_str_digits()
BARNYARD_MAX_INTEGER_DIGITS = min(
    _RUNTIME_INTEGER_DIGIT_LIMIT or sys.int_info.default_max_str_digits,
    sys.int_info.default_max_str_digits,
)
BARNYARD_SUMMARY_METRICS = (
    ("total_cells", "Cells after ambient filter", "count"),
    ("human_singlet_cells", "Human singlets", "count"),
    ("mouse_singlet_cells", "Mouse singlets", "count"),
    ("mixed_cells", "Mixed cells", "count"),
    ("unclassified_cells", "Unclassified cells", "count"),
    (
        "cross_species_doublet_rate_among_classified",
        "Cross-species doublet rate",
        "rate",
    ),
)
SCRIPT_DIR = Path(__file__).resolve().parent
REPORT_TEMPLATE_PATH = SCRIPT_DIR / "report_template.html"
PLOTLY_JS_PATH = SCRIPT_DIR / "plotly-2.26.0.min.js"
PARAMETER_DISPLAY_ORDER = [
    "sample_id",
    "fastq",
    "full_length_fastq",
    "tso_seq",
    "rtp_seq",
    "ref_dir",
    "out_dir",
    "threads",
    "cluster_threads",
    "exp_cells",
]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output-html", required=True)
    parser.add_argument("--report-metrics-tsv", required=True)
    parser.add_argument("--rna-qc-metrics-tsv", required=True)
    parser.add_argument("--saturation-tsv", required=True)
    parser.add_argument("--read-qc-json", required=True)
    parser.add_argument("--parameters-tsv", required=True)
    parser.add_argument("--per-cell-qc-tsv", required=True)
    parser.add_argument("--rna-cluster-tsv", required=True)
    parser.add_argument("--barcode-counts-3p-tsv", default=None)
    parser.add_argument("--barcode-counts-5p-tsv", default=None)
    parser.add_argument("--whitelist-3p", default=None)
    parser.add_argument("--whitelist-5p", default=None)
    parser.add_argument("--read-assigned-cell", default=None)
    parser.add_argument("--glycine-stats", default=None)
    parser.add_argument("--skip-glycine", action="store_true")
    parser.add_argument("--knee-plot-3p", default=None)
    parser.add_argument("--knee-plot-5p", default=None)
    parser.add_argument("--saturation-png", default=None)
    parser.add_argument("--rna-violin-png", default=None)
    parser.add_argument("--barnyard-summary-tsv", default=None)
    parser.add_argument("--barnyard-per-cell-tsv", default=None)
    return parser.parse_args()


def read_tsv(path, required=True):
    if not path or not Path(path).exists():
        if required:
            raise FileNotFoundError(path)
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t").fillna("")


def value_to_float(value):
    if value is None:
        return None
    text = str(value).strip().replace(",", "")
    if text == "":
        return None
    if text.endswith("%"):
        text = text[:-1]
    try:
        return float(text)
    except ValueError:
        return None


def metric(report_df, qc_df, name, section=None):
    if not report_df.empty and {"Metric", "Value"}.issubset(report_df.columns):
        sub = report_df[report_df["Metric"].astype(str) == name]
        if section is not None and "Section" in sub.columns:
            sub = sub[sub["Section"].astype(str) == section]
        if not sub.empty:
            return value_to_float(sub.iloc[0]["Value"])
    if not qc_df.empty and {"Metric", "Value"}.issubset(qc_df.columns):
        sub = qc_df[qc_df["Metric"].astype(str) == name]
        if not sub.empty:
            return value_to_float(sub.iloc[0]["Value"])
    return None


def metric_text(report_df, qc_df, name, section=None):
    if not report_df.empty and {"Metric", "Formatted_value"}.issubset(report_df.columns):
        sub = report_df[report_df["Metric"].astype(str) == name]
        if section is not None and "Section" in sub.columns:
            sub = sub[sub["Section"].astype(str) == section]
        if not sub.empty and str(sub.iloc[0]["Formatted_value"]).strip():
            return str(sub.iloc[0]["Formatted_value"])
    value = metric(report_df, qc_df, name, section)
    return format_number(value)


def format_number(value, digits=2):
    value = value_to_float(value)
    if value is None or not math.isfinite(value):
        return "NA"
    if float(value).is_integer():
        return f"{int(value):,}"
    return f"{value:,.{digits}f}"


def format_integer(value):
    value = value_to_float(value)
    if value is None or not math.isfinite(value):
        return "NA"
    return f"{int(round(value)):,}"


def format_percent(value, digits=2):
    value = value_to_float(value)
    if value is None or not math.isfinite(value):
        return "NA"
    return f"{value * 100:.{digits}f}%"


def compact_number(value):
    value = value_to_float(value)
    if value is None or not math.isfinite(value):
        return "NA"
    abs_value = abs(value)
    if abs_value >= 1_000_000:
        return f"{value / 1_000_000:.2f}M"
    if abs_value >= 1_000:
        return f"{value / 1_000:.1f}K"
    if float(value).is_integer():
        return str(int(value))
    return f"{value:.2f}"


def read_parameters(path):
    df = read_tsv(path, required=False)
    if df.empty or not {"Parameter", "Value"}.issubset(df.columns):
        return {}
    return {
        str(row["Parameter"]): str(row["Value"])
        for _, row in df.iterrows()
        if str(row["Parameter"]).strip()
    }


def html_table(rows, headers=("Metric", "Value"), class_name="report-table"):
    header_html = "".join(f"<th>{html.escape(str(h))}</th>" for h in headers)
    body = []
    for row in rows:
        body.append(
            "<tr>"
            + "".join(f"<td>{html.escape(str(row.get(h, '')))}</td>" for h in headers)
            + "</tr>"
        )
    return f"<table class=\"{class_name}\"><thead><tr>{header_html}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def load_report_template():
    text = REPORT_TEMPLATE_PATH.read_text(encoding="utf-8", errors="ignore")

    # The PBMC-derived template includes one malformed style/script splice and a
    # large "Relocated from body" CSS blob that breaks HTML parsing in some
    # browsers, causing a blank white report page.
    text = re.sub(
        r"""
        \s*/\*\s*Relocated\ from\ body\s*\*/.*?
        </script><script>"use\ strict";\s*
        </style>
        """,
        "\n    </style>\n",
        text,
        count=1,
        flags=re.S | re.X,
    )

    # The template shell should have a real <body>; some captured PBMC exports
    # omit it and rely on browser recovery, which is fragile once we inject our
    # own report body and scripts.
    head_end = text.find("</head>")
    if head_end != -1:
        body_after_head = re.search(r"<body\b", text[head_end:], flags=re.I)
        if body_after_head is None:
            text = text.replace("</head>", "</head>\n<body>", 1)

    required_markers = ["__SAMPLE__", "__REPORT_BODY__", "__PLOTLY_LOADER__"]
    missing = [marker for marker in required_markers if marker not in text]
    if missing:
        raise ValueError("report template is missing: " + ", ".join(missing))
    return text


def embedded_plotly_loader():
    if not PLOTLY_JS_PATH.is_file():
        raise FileNotFoundError(f"missing bundled Plotly runtime: {PLOTLY_JS_PATH}")
    return f"<script>{PLOTLY_JS_PATH.read_text(encoding='utf-8')}</script>"


def pbmc_metric_rows(rows):
    body = "".join(
        f"""
        <tr>
          <td>{html.escape(str(row.get("Metric", "")))}</td>
          <td class="metric-value">{html.escape(str(row.get("Value", "")))}</td>
        </tr>
        """
        for row in rows
    )
    return f'<div class="stats-table-container"><table class="stats-table"><tbody>{body}</tbody></table></div>'


def pbmc_metric_rows_3col(rows):
    body = "".join(
        f"""
        <tr>
          <td>{html.escape(str(row.get("Metric", "")))}</td>
          <td class="metric-value read-count-value">{html.escape(str(row.get("Read count", "")))}</td>
          <td class="metric-value percent-value">{html.escape(str(row.get("Percent", "")))}</td>
        </tr>
        """
        for row in rows
    )
    return (
        '<div class="stats-table-container"><table class="stats-table sequencing-table">'
        '<colgroup><col class="metric-column"><col class="read-count-column">'
        '<col class="percent-column"></colgroup>'
        '<thead><tr><th>Metric</th><th>Read count</th><th>Percent</th></tr></thead>'
        f"<tbody>{body}</tbody></table></div>"
    )


def pbmc_summary_cards(card_rows):
    return "".join(
        f"""
        <div class="summary-card">
          <div class="summary-card-content">
            <div class="value">{html.escape(str(row.get("value", "NA")))}</div>
            <div class="label">{html.escape(str(row.get("label", "")))}</div>
          </div>
        </div>
        """
        for row in card_rows
    )


def summary_violin_cards(cards):
    blocks = []
    for card in cards:
        plot_id = str(card.get("plot_id", "")).strip()
        image_html = (
            f'<div id="{html.escape(plot_id)}" class="summary-violin-plot"></div>'
            if plot_id
            else f'<div class="empty-plot">{html.escape(card.get("title", "NA"))}</div>'
        )
        blocks.append(
            f"""
            <div class="summary-plot-col">
              <div class="violin-wrap">
                {image_html}
              </div>
            </div>
            """
        )
    return "".join(blocks)


def pbmc_title_block(title, help_id, help_html, sample_id=None):
    title_text = f"{title} : {sample_id}" if sample_id else title
    return f"""
        <div class="section-heading">
            <h2>{html.escape(title_text)}</h2>
            <button class="help-button" type="button" onclick="show('{html.escape(help_id)}')" aria-label="Show help">?</button>
            <div id="{html.escape(help_id)}" class="help-panel">
                {help_html}
            </div>
        </div>
    """


def pbmc_section_bar(title, help_id, help_html, width_px=240):
    return f"""
        <div class="section-heading">
            <h2>{html.escape(title)}</h2>
            <button class="help-button" type="button" onclick="show('{html.escape(help_id)}')" aria-label="Show help">?</button>
            <div id="{html.escape(help_id)}" class="help-panel">
                {help_html}
            </div>
        </div>
    """


def help_dl(items):
    parts = ['<dl class="help-list">']
    for title, desc in items:
        parts.append(f"<dt>{html.escape(title)}</dt>")
        parts.append(f"<dd>{html.escape(desc)}</dd>")
    parts.append("</dl>")
    return "".join(parts)


def barnyard_report_section(payload):
    if payload is None:
        return ""
    help_html = help_dl([
        ("UMI axes", "Human UMI is shown on the x-axis and Mouse UMI on the y-axis."),
        ("Class colors", "Human singlets are green, mouse singlets are blue, mixed cells are red, and unclassified cells are gray."),
        ("Ambient filtering", "Cells with both species' UMI counts below the configured ambient threshold are excluded before this summary."),
        ("Singlets", "Cells whose UMI fraction meets the configured human- or mouse-dominance threshold."),
        ("Mixed cells", "Classified cells containing both species that do not meet either singlet threshold."),
        ("Cross-species doublet rate", "Mixed cells divided by classified human singlet, mouse singlet, and mixed cells."),
        ("Display sampling", "For more than 100,000 valid cells, deterministic per-class display sampling limits browser memory while legend counts retain all valid cells."),
    ])
    summary_table = pbmc_metric_rows(payload.get("summaryRows", []))
    if payload.get("displayed", 0):
        plot_content = """
          <div id="barnyard-umi" class="dynamic-plot-wide barnyard-plot"></div>
          <div class="dynamic-note" id="barnyard-sampling-note"></div>
        """
    else:
        plot_content = """
          <div class="barnyard-empty-state" role="status">
            No valid Barnyard cells were available for plotting.
          </div>
        """
    return f"""
    <section class="section-card report-section barnyard-section" data-library-content="gene-expression">
      {pbmc_section_bar("Barnyard QC", "barnyard-detail", help_html, width_px=220)}
      <div class="report-grid two barnyard-grid">
        <div class="table-panel barnyard-summary-table">{summary_table}</div>
        <div class="plot-panel barnyard-plot-panel">
          <h4 class="plot-title">Barnyard by UMI</h4>
          {plot_content}
        </div>
      </div>
    </section>
    """


def build_sample_rows(params):
    rows = []
    for key in PARAMETER_DISPLAY_ORDER:
        if key in params and params[key] != "":
            rows.append({"Metric": key, "Value": params[key]})
    return rows


def parse_glycine_read_counts(path):
    if not path or not Path(path).exists():
        return None, None
    text = Path(path).read_text(encoding="utf-8", errors="ignore")
    values = {}
    for line in text.splitlines():
        parts = re.split(r"[\t: ]+", line.strip())
        if len(parts) < 2:
            continue
        key = parts[0].strip()
        for token in parts[1:]:
            if re.fullmatch(r"\d+(?:\.\d+)?", token):
                values.setdefault(key, float(token))
                break
    read_count = values.get("Read_count")
    length_filtered = values.get("Length-filtered", 0.0)
    qc_filtered = values.get("QC-filtered", 0.0)
    if read_count is None:
        return None, None
    return int(read_count), max(0, int(read_count - length_filtered - qc_filtered))


def parse_glycine_clean_reads(path):
    return parse_glycine_read_counts(path)[1]


def build_read_summary(report_df, qc_df, skip_glycine, glycine_stats):
    ordered = [
        ("Full length", "Full length"),
        ("Barcode-valid", "Barcode-valid"),
        ("Cell-assigned", "Cell-assigned"),
        ("Gene assigned", "Gene assigned"),
        ("Transcript assigned", "Transcript assigned"),
    ]
    full_length = metric(report_df, qc_df, "Full length", "Read assignment summary")
    if full_length is None:
        full_length = metric(report_df, qc_df, "Full length reads")
    glycine_raw_reads, glycine_clean_reads = parse_glycine_read_counts(glycine_stats)
    clean_reads = int(full_length) if skip_glycine and full_length is not None else glycine_clean_reads
    if clean_reads is None and full_length is not None:
        clean_reads = int(full_length)

    raw_reads = metric(report_df, qc_df, "Input reads")
    if raw_reads is None:
        raw_reads = full_length if skip_glycine else glycine_raw_reads
    if raw_reads is None:
        raw_reads = clean_reads or full_length

    rows = [{"Metric": "Raw reads", "Read count": format_number(raw_reads), "Percent": "100.00%"}]
    if clean_reads is not None:
        clean_ratio = clean_reads / raw_reads if raw_reads else None
        rows.append({
            "Metric": "Clean reads",
            "Read count": format_number(clean_reads),
            "Percent": format_percent(clean_ratio),
        })
    denominator = raw_reads
    for label, metric_name in ordered:
        count = metric(report_df, qc_df, metric_name, "Read assignment summary")
        percent = count / denominator if count is not None and denominator else None
        rows.append(
            {
                "Metric": label,
                "Read count": format_number(count),
                "Percent": format_percent(percent),
            }
        )
    return rows, denominator


def load_whitelist(path):
    if not path or not Path(path).exists():
        return set()
    return {
        line.strip().split(",")[0]
        for line in Path(path).read_text(encoding="utf-8", errors="ignore").splitlines()
        if line.strip()
    }


def count_data_rows(path):
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        total = sum(1 for _ in handle)
    return max(0, total - 1)


def barcode_rank_payload(path, whitelist_path):
    if not path or not Path(path).exists():
        return None
    total = count_data_rows(path)
    if total == 0:
        return {
            "rank": [],
            "count": [],
            "is_true": [],
            "displayed_points": 0,
            "original_points": 0,
            "true_points": 0,
            "noise_points": 0,
        }
    if total <= BARCODE_RANK_MAX_POINTS:
        keep_ranks = set(range(1, total + 1))
    else:
        log_start = 0.0
        log_end = math.log(total)
        keep_ranks = {
            max(1, min(total, int(round(math.exp(log_start + (log_end - log_start) * i / (BARCODE_RANK_MAX_POINTS - 1))))))
            for i in range(BARCODE_RANK_MAX_POINTS)
        }
        keep_ranks.add(1)
        keep_ranks.add(total)

    whitelist = load_whitelist(whitelist_path)
    ranks = []
    counts = []
    is_true = []
    true_points = 0
    noise_points = 0
    with open(path, "r", encoding="utf-8", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for rank, row in enumerate(reader, start=1):
            barcode = str(row.get("barcode", "")).strip()
            count = value_to_float(row.get("count", 0)) or 0.0
            called = not whitelist or barcode in whitelist
            if called:
                true_points += 1
            else:
                noise_points += 1
            if rank in keep_ranks:
                ranks.append(rank)
                counts.append(count)
                is_true.append(called)
    return {
        "rank": ranks,
        "count": counts,
        "is_true": is_true,
        "displayed_points": len(ranks),
        "original_points": total,
        "true_points": true_points,
        "noise_points": noise_points,
    }


def _barnyard_warning(message):
    print(f"WARNING: {message}", file=sys.stderr)


def _finite_nonnegative_integer(value):
    try:
        number = Decimal(str(value).strip())
    except (InvalidOperation, TypeError, ValueError):
        return None
    if not number.is_finite() or number < 0 or number != number.to_integral_value():
        return None
    integer_digits = 1 if number.is_zero() else number.adjusted() + 1
    if integer_digits > BARNYARD_MAX_INTEGER_DIGITS:
        return None
    return int(number)


def barnyard_summary_rows(path):
    path = Path(path) if path else None
    if path is None or not path.is_file():
        _barnyard_warning(f"Barnyard summary TSV is missing: {path}")
        return None
    try:
        frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except (OSError, UnicodeError, pd.errors.ParserError, ValueError) as exc:
        _barnyard_warning(f"cannot read Barnyard summary TSV {path}: {exc}")
        return None
    if not {"metric", "value"}.issubset(frame.columns):
        _barnyard_warning("Barnyard summary TSV must contain metric and value columns")
        return None

    selected_names = {name for name, _, _ in BARNYARD_SUMMARY_METRICS}
    selected = frame[frame["metric"].astype(str).isin(selected_names)]
    duplicate_names = sorted(
        name for name, count in selected["metric"].value_counts().items() if count > 1
    )
    if duplicate_names:
        _barnyard_warning(
            "Barnyard summary contains duplicate selected metric(s): "
            + ", ".join(duplicate_names)
        )
        return None

    values = {
        str(row.metric): row.value
        for row in selected[["metric", "value"]].itertuples(index=False)
    }
    invalid_count = 0
    rows = []
    for name, label, value_type in BARNYARD_SUMMARY_METRICS:
        raw_value = values.get(name)
        if value_type == "count":
            value = _finite_nonnegative_integer(raw_value)
            formatted = f"{value:,}" if value is not None else "NA"
        else:
            try:
                value = float(str(raw_value).strip())
            except (TypeError, ValueError):
                value = None
            if value is None or not math.isfinite(value) or not 0 <= value <= 1:
                formatted = "NA"
            else:
                formatted = format_percent(value)
        if formatted == "NA":
            invalid_count += 1
        rows.append({"Metric": label, "Value": formatted})

    if invalid_count:
        _barnyard_warning(f"{invalid_count} invalid Barnyard summary values displayed as NA")
    return rows


class BarnyardPointSampler:
    def __init__(self):
        self._points = []
        self._heaps = {assignment: [] for assignment in BARNYARD_ASSIGNMENTS}
        self.class_counts = {assignment: 0 for assignment in BARNYARD_ASSIGNMENTS}
        self.sampled = False
        self.max_retained_count = 0

    @staticmethod
    def _priority(cell_id):
        digest = hashlib.blake2b(cell_id.encode("utf-8"), digest_size=16).digest()
        return int.from_bytes(digest, "big")

    @property
    def retained_count(self):
        if self.sampled:
            return sum(len(heap) for heap in self._heaps.values())
        return len(self._points)

    @property
    def heap_sizes_by_class(self):
        return {assignment: len(heap) for assignment, heap in self._heaps.items()}

    def _heap_add(self, point):
        cell_id, assignment, _, _, _ = point
        priority = self._priority(cell_id)
        heap = self._heaps[assignment]
        entry = (-priority, point)
        if len(heap) < BARNYARD_MAX_POINTS_PER_CLASS:
            heapq.heappush(heap, entry)
        elif priority < -heap[0][0]:
            heapq.heapreplace(heap, entry)

    def add(self, point):
        assignment = point[1]
        self.class_counts[assignment] += 1
        if not self.sampled and len(self._points) < BARNYARD_MAX_POINTS:
            self._points.append(point)
        else:
            if not self.sampled:
                self.sampled = True
                for retained_point in self._points:
                    self._heap_add(retained_point)
                self._points = []
            self._heap_add(point)
        self.max_retained_count = max(self.max_retained_count, self.retained_count)

    def finish(self):
        grouped = {assignment: [] for assignment in BARNYARD_ASSIGNMENTS}
        if self.sampled:
            for assignment, heap in self._heaps.items():
                grouped[assignment] = [
                    entry[1]
                    for entry in sorted(heap, key=lambda entry: (-entry[0], entry[1][0]))
                ]
        else:
            for point in self._points:
                grouped[point[1]].append(point)

        traces = {}
        for assignment in BARNYARD_ASSIGNMENTS:
            points = grouped[assignment]
            traces[assignment] = {
                "cell_id": [point[0] for point in points],
                "umi_human": [point[2] for point in points],
                "umi_mouse": [point[3] for point in points],
                "umi_total": [point[4] for point in points],
            }
        return traces


def barnyard_payload(summary_path, per_cell_path):
    if not summary_path and not per_cell_path:
        return None
    if not summary_path or not per_cell_path:
        _barnyard_warning("both Barnyard summary and per-cell TSVs are required")
        return None

    summary_rows = barnyard_summary_rows(summary_path)
    if summary_rows is None:
        return None

    per_cell_path = Path(per_cell_path)
    if not per_cell_path.is_file():
        _barnyard_warning(f"Barnyard per-cell TSV is missing: {per_cell_path}")
        return None

    database_fd = None
    database_path = None
    connection = None
    sampler = BarnyardPointSampler()
    duplicate_count = 0
    unknown_assignment_count = 0
    invalid_row_count = 0
    try:
        database_fd, database_path = tempfile.mkstemp(
            prefix="flora-barnyard-seen-", suffix=".sqlite3"
        )
        os.close(database_fd)
        database_fd = None
        connection = sqlite3.connect(database_path)
        connection.execute("PRAGMA journal_mode=OFF")
        connection.execute("PRAGMA synchronous=OFF")
        connection.execute("CREATE TABLE seen (cell_id TEXT PRIMARY KEY)")

        chunks = pd.read_csv(
            per_cell_path,
            sep="\t",
            dtype=str,
            keep_default_na=False,
            usecols=list(BARNYARD_REQUIRED_COLUMNS),
            chunksize=BARNYARD_CHUNK_ROWS,
        )
        for chunk in chunks:
            for raw_point in chunk.itertuples(index=False, name=None):
                row = dict(zip(chunk.columns, raw_point))
                cell_id = str(row["cell_id"]).strip()
                assignment = str(row["assignment"]).strip()
                if not cell_id:
                    invalid_row_count += 1
                    continue
                if assignment not in BARNYARD_ASSIGNMENTS:
                    unknown_assignment_count += 1
                    continue
                umi_values = [
                    _finite_nonnegative_integer(row[column])
                    for column in ("umi_human", "umi_mouse", "umi_total")
                ]
                if any(
                    value is None or value > BARNYARD_JS_MAX_SAFE_INTEGER
                    for value in umi_values
                ):
                    invalid_row_count += 1
                    continue

                try:
                    connection.execute("INSERT INTO seen(cell_id) VALUES (?)", (cell_id,))
                except sqlite3.IntegrityError:
                    duplicate_count += 1
                    continue
                sampler.add((cell_id, assignment, *umi_values))
        connection.commit()
    except (OSError, UnicodeError, ValueError, KeyError, sqlite3.Error, pd.errors.ParserError) as exc:
        _barnyard_warning(f"cannot build Barnyard per-cell payload from {per_cell_path}: {exc}")
        return None
    finally:
        if connection is not None:
            connection.close()
        if database_fd is not None:
            os.close(database_fd)
        if database_path is not None:
            try:
                os.unlink(database_path)
            except FileNotFoundError:
                pass

    if duplicate_count:
        _barnyard_warning(f"skipped {duplicate_count} duplicate cell IDs")
    if unknown_assignment_count:
        _barnyard_warning(f"skipped {unknown_assignment_count} unknown assignments")
    if invalid_row_count:
        _barnyard_warning(f"skipped {invalid_row_count} invalid per-cell rows")

    return {
        "summaryRows": summary_rows,
        "traces": sampler.finish(),
        "classCounts": sampler.class_counts,
        "totalValid": sum(sampler.class_counts.values()),
        "displayed": sampler.retained_count,
        "sampled": sampler.sampled,
    }


def beads_per_droplet_payload(path):
    if not path or not Path(path).exists():
        return {"x": [], "y": [], "n_cells": 0}
    cell_barcodes = defaultdict(set)
    usecols = None
    header = pd.read_csv(path, nrows=0)
    candidates = [c for c in ["cell_id", "BC5n", "BC3n"] if c in header.columns]
    if {"cell_id", "BC5n", "BC3n"}.issubset(candidates):
        usecols = candidates
    for chunk in pd.read_csv(path, chunksize=BEAD_CHUNK_ROWS, usecols=usecols):
        if not {"cell_id", "BC5n", "BC3n"}.issubset(chunk.columns):
            return {"x": [], "y": [], "n_cells": 0}
        for cell_id, bc5, bc3 in zip(chunk["cell_id"], chunk["BC5n"], chunk["BC3n"]):
            cell = str(cell_id).strip()
            if not cell:
                continue
            for barcode in (bc5, bc3):
                bc = str(barcode).strip()
                if bc and bc.lower() != "nan":
                    cell_barcodes[cell].add(bc)
    hist = Counter(len(v) for v in cell_barcodes.values())
    xs = sorted(hist)
    return {"x": xs, "y": [hist[x] for x in xs], "n_cells": len(cell_barcodes)}


def dataframe_payload(df):
    return {col: df[col].tolist() for col in df.columns}


def rna_cluster_payload(path):
    required_columns = [
        "cell",
        "UMAP_1",
        "UMAP_2",
        "leiden",
        "total_counts",
        "status",
    ]
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if frame.columns.tolist() != required_columns:
        raise ValueError(
            "RNA cluster TSV must contain exactly: " + ", ".join(required_columns)
        )

    rows = []
    seen_cells = set()
    for row in frame.itertuples(index=False, name=None):
        cell, umap_1, umap_2, leiden, total_counts, status = row
        if not cell or cell in seen_cells:
            raise ValueError("RNA cluster TSV contains a blank or duplicate cell")
        seen_cells.add(cell)
        coordinates = [float(umap_1), float(umap_2)]
        if not all(math.isfinite(value) for value in coordinates):
            raise ValueError("RNA cluster TSV contains a non-finite UMAP coordinate")
        count = int(total_counts)
        if count < 0:
            raise ValueError("RNA cluster TSV contains a negative total_counts value")
        rows.append(
            {
                "cell": cell,
                "UMAP_1": coordinates[0],
                "UMAP_2": coordinates[1],
                "leiden": leiden,
                "total_counts": count,
                "status": status,
            }
        )
    return rows


def script_safe_json(value):
    return (
        json.dumps(value, ensure_ascii=False)
        .replace("&", "\\u0026")
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
        .replace("\u2028", "\\u2028")
        .replace("\u2029", "\\u2029")
    )


def safe_read_json(path):
    if not path or not Path(path).exists():
        return {}
    return json.loads(Path(path).read_text(encoding="utf-8"))


def summary_rows(report_df, qc_df):
    items = [
        ("Estimated cells", "Estimated cells"),
        ("Input reads", "Input reads"),
        ("Reads per cell (mean)", "Reads per cell (mean)"),
        ("UMIs per cell (median)", "UMIs per cell (median)"),
        ("Genes per cell (median)", "Genes per cell (median)"),
        ("Unique genes", "Unique genes"),
        ("Unique isoforms", "Unique isoforms"),
    ]
    return [
        {"Metric": label, "Value": format_integer(metric(report_df, qc_df, name))}
        for label, name in items
    ]


def mapping_rows(report_df, qc_df, denominator):
    aligned = metric(report_df, qc_df, "Aligned BAM reads")
    unmapped = metric(report_df, qc_df, "Unmapped")
    if denominator is None:
        denominator = metric(report_df, qc_df, "Input reads")
    aligned_ratio = aligned / denominator if aligned is not None and denominator else None
    unmapped_ratio = unmapped / denominator if unmapped is not None and denominator else None
    return [
        {"Metric": "Aligned BAM reads / total reads", "Value": format_percent(aligned_ratio)},
        {"Metric": "Unmapped / total reads", "Value": format_percent(unmapped_ratio)},
        {"Metric": "Unique genes", "Value": metric_text(report_df, qc_df, "Unique genes")},
        {"Metric": "Unique isoforms", "Value": metric_text(report_df, qc_df, "Unique isoforms")},
    ]


def per_cell_payload(path):
    df = read_tsv(path, required=False)
    if df.empty:
        return {"reads": [], "umis": [], "genes": [], "mito_percent": []}
    keep = min(len(df), 20000)
    if len(df) > keep:
        df = df.sample(n=keep, random_state=1)
    payload = {}
    for col in ["reads", "umis", "genes", "mito_percent"]:
        payload[col] = pd.to_numeric(df.get(col, pd.Series(dtype=float)), errors="coerce").fillna(0).tolist()
    return payload


def violin_data_uri(values, color, title, ylabel):
    series = pd.to_numeric(pd.Series(values, dtype=float), errors="coerce")
    series = series[series.notna()]
    if series.empty:
        return ""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return ""

    fig, ax = plt.subplots(figsize=(3.5, 5.6), dpi=150)
    violin = ax.violinplot(series.tolist(), showmeans=True, showmedians=True, showextrema=False, widths=0.75)
    for body in violin["bodies"]:
        body.set_facecolor(color)
        body.set_edgecolor(color)
        body.set_alpha(0.70)
    if "cmeans" in violin:
        violin["cmeans"].set_color("#333333")
        violin["cmeans"].set_linewidth(1.2)
    if "cmedians" in violin:
        violin["cmedians"].set_color("#111111")
        violin["cmedians"].set_linewidth(1.3)

    ax.set_title(title, fontsize=11, pad=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_xticks([])
    ax.grid(axis="y", color="#d9d9d9", linewidth=0.8)
    ax.set_facecolor("white")
    for spine in ["top", "right", "bottom"]:
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("#888888")
    fig.patch.set_facecolor("white")
    fig.tight_layout()

    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    encoded = base64.b64encode(buffer.getvalue()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def new_report_markup(sections):
    return f"""
<style>
  .report-section {{ margin-bottom: 1.5rem; }}
  .section-heading {{
    display: flex;
    align-items: center;
    flex-wrap: wrap;
    gap: 0.75rem;
    margin-bottom: 1.25rem;
  }}
  .section-heading h2 {{
    margin: 0;
    font-size: 1.35rem;
    font-weight: 650;
    color: var(--text-primary);
  }}
  .help-button {{
    display: inline-flex;
    align-items: center;
    justify-content: center;
    width: 28px;
    height: 28px;
    padding: 0;
    border: 0;
    border-radius: 50%;
    color: #fff;
    background: var(--gradient-bg);
    box-shadow: 0 2px 8px rgba(32, 85, 138, 0.25);
    font-weight: 700;
  }}
  .help-panel {{
    display: none;
    flex-basis: 100%;
    padding: 1rem 1.1rem;
    border: 1px solid var(--border-color);
    border-radius: 10px;
    background: var(--bg-surface-hover);
    color: var(--text-secondary);
  }}
  .help-list {{ margin: 0; }}
  .help-list dt {{
    margin-top: 0.7rem;
    color: var(--text-primary);
    font-weight: 650;
  }}
  .help-list dt:first-child {{ margin-top: 0; }}
  .help-list dd {{ margin: 0.15rem 0 0; }}
  .summary-cards {{ margin-bottom: 1.5rem; }}
  .report-grid {{
    display: grid;
    gap: 1rem;
    align-items: stretch;
  }}
  .report-grid.two {{ grid-template-columns: repeat(2, minmax(0, 1fr)); }}
  .report-grid.three {{ grid-template-columns: repeat(3, minmax(0, 1fr)); }}
  .summary-detail-grid {{
    display: grid;
    grid-template-columns: minmax(260px, 0.8fr) minmax(0, 1.7fr);
    gap: 1rem;
    align-items: stretch;
  }}
  .summary-violin-grid {{
    display: grid;
    grid-template-columns: repeat(3, minmax(0, 1fr));
    gap: 0.75rem;
  }}
  .plot-panel,
  .table-panel {{
    box-sizing: border-box;
    min-width: 0;
    max-width: 100%;
    padding: 1rem;
    border: 1px solid var(--border-color);
    border-radius: 10px;
    background: var(--bg-surface);
    box-shadow: var(--shadow-sm);
  }}
  .plot-panel {{ overflow: hidden; }}
  .plot-title {{
    margin: 0 0 0.75rem;
    padding: 0;
    border: 0;
    color: var(--text-primary);
    font-family: Inter, Arial, sans-serif;
    font-size: 1rem;
    font-weight: 600;
    text-align: center;
  }}
  .stats-table-container {{ width: 100%; overflow-x: auto; }}
  .stats-table {{
    width: 100%;
    border-collapse: collapse;
    table-layout: fixed;
  }}
  .stats-table th,
  .stats-table td {{
    padding: 0.75rem 0.65rem;
    text-align: left;
    vertical-align: middle;
  }}
  .stats-table thead tr,
  .stats-table tbody tr {{ border-bottom: 1px solid var(--border-color); }}
  .stats-table th {{
    color: var(--text-secondary);
    font-size: 0.75rem;
    text-transform: uppercase;
    letter-spacing: 0.04em;
  }}
  .stats-table tbody tr:last-child {{ border-bottom: 0; }}
  .stats-table .metric-value {{
    text-align: right;
    color: var(--primary-dark);
    font-weight: 650;
    white-space: nowrap;
  }}
  .stats-table td.metric-value {{ display: table-cell; }}
  .sequencing-table .metric-column {{ width: 50%; }}
  .sequencing-table .read-count-column {{ width: 30%; }}
  .sequencing-table .percent-column {{ width: 20%; }}
  .sequencing-table th:nth-child(2),
  .sequencing-table th:nth-child(3) {{ text-align: right; white-space: nowrap; }}
  .sequencing-table .read-count-value,
  .sequencing-table .percent-value {{ white-space: nowrap; }}
  .dynamic-note {{
    margin-top: 0.5rem;
    color: var(--text-muted);
    font-size: 0.8rem;
  }}
  .barnyard-summary-table .stats-table td,
  .barnyard-summary-table .stats-table th {{
    padding-top: 0.65rem;
    padding-bottom: 0.65rem;
  }}
  .barnyard-plot-panel {{ min-width: 0; }}
  .barnyard-empty-state {{
    min-height: 320px;
    display: flex;
    align-items: center;
    justify-content: center;
    padding: 1.5rem;
    color: var(--text-muted);
    text-align: center;
  }}
  .dynamic-plot,
  .dynamic-plot-wide,
  .dynamic-plot-lg {{
    box-sizing: border-box;
    width: 100%;
    max-width: 100%;
    min-width: 0;
    overflow: hidden;
  }}
  .dynamic-plot .plot-container,
  .dynamic-plot-wide .plot-container,
  .dynamic-plot-lg .plot-container,
  .dynamic-plot .svg-container,
  .dynamic-plot-wide .svg-container,
  .dynamic-plot-lg .svg-container {{ max-width: 100%; }}
  .dynamic-plot {{ height: 360px; }}
  .dynamic-plot-wide,
  .dynamic-plot-lg {{ height: 400px; }}
  .violin-wrap {{
    display: flex;
    align-items: center;
    justify-content: center;
    min-height: 390px;
    max-width: 100%;
    overflow: hidden;
  }}
  .violin-image {{
    display: block;
    width: 100%;
    max-width: 100%;
    height: 390px;
    object-fit: contain;
  }}
  .summary-violin-plot {{
    width: 100%;
    height: 390px;
    min-width: 0;
  }}
  .empty-plot {{ color: var(--text-muted); }}
  @media (max-width: 900px) {{
    .report-grid.two,
    .report-grid.three,
    .summary-detail-grid,
    .summary-violin-grid {{ grid-template-columns: 1fr; }}
  }}
</style>
<div class="tab-content">
  <div class="tab-pane active" id="summary-tab">
    <section class="summary-section report-section" data-library-content="gene-expression">
      {sections["sample_title"]}
      <div class="summary-cards">{sections["summary_cards"]}</div>
    </section>

    <section class="section-card report-section" data-library-content="gene-expression">
      {sections["read_qc_bar"]}
      <div class="report-grid three">
        <div class="plot-panel"><div id="read-qc-quality" class="dynamic-plot"></div></div>
        <div class="plot-panel"><div id="read-qc-length" class="dynamic-plot"></div></div>
        <div class="plot-panel"><div id="read-qc-yield" class="dynamic-plot"></div></div>
      </div>
    </section>
  </div>

  <div class="tab-pane" id="cells-tab">
    <section class="section-card report-section" data-library-content="gene-expression">
      {sections["summary_bar"]}
      <div class="summary-detail-grid">
        <div class="table-panel">{sections["summary_rows"]}</div>
        <div class="summary-violin-grid">{sections["summary_violin_plots"]}</div>
      </div>
    </section>

    <section class="section-card report-section" data-library-content="gene-expression">
      {sections["beads_bar"]}
      <div class="report-grid two">
        <div class="plot-panel">
          <h4 class="plot-title">RNA Barcode Rank Plot</h4>
          <div id="barcode-rank" class="dynamic-plot"></div>
          <div class="dynamic-note" id="barcode-rank-note"></div>
        </div>
        <div class="plot-panel">
          <h4 class="plot-title">Bead Count Distribution</h4>
          <div id="beads-per-droplet" class="dynamic-plot"></div>
          <div class="dynamic-note" id="beads-note"></div>
        </div>
      </div>
    </section>

    {sections.get("barnyard_section", "")}

    <section class="section-card report-section" data-library-content="gene-expression">
      {sections["rna_cluster_bar"]}
      <div class="report-grid two">
        <div class="plot-panel">
          <h4 class="plot-title">RNA Cluster Assignment</h4>
          <div id="rna-cluster-assignment" class="dynamic-plot-wide"></div>
        </div>
        <div class="plot-panel">
          <h4 class="plot-title">UMI Count</h4>
          <div id="rna-umi-counts" class="dynamic-plot-wide"></div>
        </div>
      </div>
    </section>
  </div>

  <div class="tab-pane" id="library-tab">
    <section class="section-card report-section" data-library-content="gene-expression">
      {sections["sequencing_bar"]}
      <div class="report-grid two">
        <div class="table-panel" id="ReadSummaryTable">{sections["read_summary_rows"]}</div>
        <div class="plot-panel"><div id="read-assignment-plot" class="dynamic-plot-wide"></div></div>
      </div>
    </section>

    <section class="section-card report-section" data-library-content="gene-expression">
      {sections["mapping_bar"]}
      <div class="table-panel" id="MappingTable">{sections["mapping_rows"]}</div>
    </section>

    <section class="section-card report-section" data-library-content="gene-expression">
      {sections["saturation_bar"]}
      <div class="report-grid three">
        <div class="plot-panel"><div id="saturation-genes" class="dynamic-plot-lg"></div></div>
        <div class="plot-panel"><div id="saturation-umis" class="dynamic-plot-lg"></div></div>
        <div class="plot-panel"><div id="saturation-rate" class="dynamic-plot-lg"></div></div>
      </div>
    </section>
  </div>
</div>
"""


def barnyard_report_script(payload):
    barnyard = payload.get("barnyard")
    if not barnyard or not barnyard.get("displayed", 0):
        return ""
    return """
  const barnyard = payload.barnyard;
  const barnyardTraceType = barnyard.displayed > 2000 ? "scattergl" : "scatter";
  const barnyardAssignments = ["human_singlet", "mouse_singlet", "mixed", "unclassified"];
  const barnyardClassLabels = {
    human_singlet: "Human singlet",
    mouse_singlet: "Mouse singlet",
    mixed: "Mixed",
    unclassified: "Unclassified"
  };
  const barnyardClassColors = {
    human_singlet: "#2E7D32",
    mouse_singlet: "#1565C0",
    mixed: "#C62828",
    unclassified: "#757575"
  };
  const barnyardTraces = barnyardAssignments.map(assignment => {
    const trace = barnyard.traces[assignment] || {cell_id: [], umi_human: [], umi_mouse: [], umi_total: []};
    const fullCount = barnyard.classCounts[assignment] || 0;
    return {
      x: trace.umi_human,
      y: trace.umi_mouse,
      customdata: trace.cell_id.map((cellId, index) => [cellId, barnyardClassLabels[assignment], trace.umi_total[index]]),
      type: barnyardTraceType,
      mode: "markers",
      name: `${barnyardClassLabels[assignment]} (${fullCount.toLocaleString()})`,
      marker: {color: barnyardClassColors[assignment], size: 5, opacity: 0.68},
      hovertemplate: "Cell ID: %{customdata[0]}<br>Assignment: %{customdata[1]}<br>Human UMI: %{x:,.0f}<br>Mouse UMI: %{y:,.0f}<br>Total UMI: %{customdata[2]:,.0f}<extra></extra>"
    };
  });
  plotIf("barnyard-umi", barnyardTraces, {
    height: 420,
    margin: {l: 65, r: 25, t: 20, b: 60},
    hovermode: "closest",
    dragmode: "zoom",
    xaxis: {
      title: {text: "Human UMI", font: {size: 12}},
      rangemode: "nonnegative",
      gridcolor: "lightgray",
      zeroline: true,
      automargin: true
    },
    yaxis: {
      title: {text: "Mouse UMI", font: {size: 12}},
      rangemode: "nonnegative",
      gridcolor: "lightgray",
      zeroline: true,
      automargin: true
    },
    legend: {
      title: {text: "Assignment"},
      font: {size: 10},
      itemsizing: "constant",
      orientation: "v",
      x: 0.98,
      y: 0.98,
      xanchor: "right",
      yanchor: "top",
      bgcolor: "rgba(255,255,255,0.88)",
      bordercolor: "rgba(44,62,80,0.25)",
      borderwidth: 1
    }
  }, {
    responsive: true,
    displaylogo: false,
    displayModeBar: true,
    scrollZoom: true,
    modeBarButtonsToRemove: ["toImage", "sendDataToCloud"]
  });
  const barnyardSamplingNote = document.getElementById("barnyard-sampling-note");
  if (barnyardSamplingNote && barnyard.sampled) {
    barnyardSamplingNote.textContent = `Showing ${barnyard.displayed.toLocaleString()} of ${barnyard.totalValid.toLocaleString()} valid unique cells using deterministic per-class display sampling.`;
  }
"""


def build_html(args, payload, sections):
    barnyard_script = barnyard_report_script(payload)
    report_body = new_report_markup(sections) + f"""
<script>
  const payload = {script_safe_json(payload)};
  const plotConfig = {{
    responsive: true,
    displaylogo: false,
    displayModeBar: false
  }};
  const baseTemplate = {{
    layout: {{
      colorway: ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"],
      font: {{family: "Arial, sans-serif", color: "#2c3e50"}},
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(0,0,0,0)",
      hoverlabel: {{align: "left"}},
      title: {{x: 0.05}}
    }}
  }};
  const baseLayout = {{
    autosize: true,
    margin: {{t: 60, l: 60, r: 20, b: 55}},
    template: baseTemplate,
    plot_bgcolor: "rgba(0,0,0,0)",
    paper_bgcolor: "rgba(0,0,0,0)",
    font: {{family: "Arial, sans-serif", color: "#2c3e50"}},
    xaxis: {{
      automargin: true,
      gridcolor: "lightgray",
      zeroline: true,
      zerolinecolor: "gray",
      zerolinewidth: 1,
      title: {{font: {{size: 12}}, standoff: 15}}
    }},
    yaxis: {{
      automargin: true,
      gridcolor: "lightgray",
      zeroline: true,
      zerolinecolor: "gray",
      zerolinewidth: 1,
      title: {{font: {{size: 12}}, standoff: 15}}
    }},
    legend: {{tracegroupgap: 0}}
  }};
  function show(id) {{
    const el = document.getElementById(id);
    if (!el) return;
    el.style.display = el.style.display === "none" || !el.style.display ? "block" : "none";
  }}
  document.querySelectorAll('.nav-tabs a').forEach(link => {{
    link.addEventListener('click', ev => {{
      ev.preventDefault();
      document.querySelectorAll('.nav-tabs li').forEach(li => li.classList.remove('active'));
      link.parentElement.classList.add('active');
      document.querySelectorAll('#myTabContent .tab-pane').forEach(p => p.classList.remove('active', 'in'));
      const target = document.querySelector(link.getAttribute('href'));
      if (target) target.classList.add('active', 'in');
      setTimeout(() => window.dispatchEvent(new Event('resize')), 50);
    }});
  }});
  function plotIf(id, data, layout, config) {{
    const el = document.getElementById(id);
    if (!el) return;
    const resolvedConfig = Object.assign({{}}, plotConfig, config || {{}});
    Plotly.newPlot(id, data, Object.assign({{}}, baseLayout, layout || {{}}), resolvedConfig);
    requestAnimationFrame(() => {{
      if (el.offsetParent !== null) Plotly.Plots.resize(el);
    }});
  }}
  function resizeReportPlots(container) {{
    (container || document).querySelectorAll(".js-plotly-plot").forEach(plot => {{
      if (plot.offsetParent !== null) Plotly.Plots.resize(plot);
    }});
  }}
  window.addEventListener("resize", () => resizeReportPlots(document));
  if (window.ResizeObserver) {{
    const reportPlotObserver = new ResizeObserver(entries => {{
      entries.forEach(entry => resizeReportPlots(entry.target));
    }});
    document.querySelectorAll(".plot-panel").forEach(panel => reportPlotObserver.observe(panel));
  }}
  function pbmcLineLayout(title, xTitle, yTitle, extra) {{
    return Object.assign({{
      title: {{text: title, font: {{size: 14}}}},
      xaxis: {{title: {{text: xTitle, font: {{size: 12}}}}}},
      yaxis: {{title: {{text: yTitle}}}}
    }}, extra || {{}});
  }}
  plotIf("read-assignment-plot", [{{
    x: payload.readSummary.labels,
    y: payload.readSummary.counts,
    type: "bar",
    marker: {{color: "rgba(15, 76, 129, 0.8)", line: {{color: "#0F4C81", width: 1}}}},
    hovertemplate: "%{{x}}<br>%{{y:,.0f}} reads<extra></extra>"
  }}], {{
    height: 420,
    title: {{text: "Read assignment summary", font: {{size: 14}}}},
    xaxis: {{title: {{text: ""}}}},
    yaxis: {{title: {{text: "Reads"}}}}
  }});

  const rq = payload.readQc || {{}};
  if (rq.quality) {{
    plotIf("read-qc-quality", [{{x: rq.quality.bins || [], y: rq.quality.counts || [], type: "bar", marker: {{color: "rgba(15, 76, 129, 0.8)", line: {{color: "#0F4C81", width: 1}}}}}}], pbmcLineLayout("Read quality", "Mean Q score", "Reads", {{
      height: 360
    }}));
  }}
  if (rq.length) {{
    plotIf("read-qc-length", [{{x: rq.length.bins_kb || [], y: rq.length.counts || [], type: "bar", marker: {{color: "rgba(46, 124, 188, 0.8)", line: {{color: "#2E7CBC", width: 1}}}}}}], pbmcLineLayout("Read length", "Read length (kb)", "Reads", {{
      height: 360
    }}));
  }}
  if (rq.yield_above_length) {{
    plotIf("read-qc-yield", [{{x: rq.yield_above_length.x_kb || [], y: rq.yield_above_length.y_gb || [], type: "scatter", mode: "lines", line: {{color: "#1358A2", width: 3}}}}], pbmcLineLayout("Base yield above read length", "Read length cutoff (kb)", "Yield above cutoff (Gb)", {{
      height: 360
    }}));
  }}

  const perCell = payload.perCell || {{}};
  function plotTemplateViolin(id, values, label) {{
    const cleanValues = (values || []).filter(value => Number.isFinite(Number(value))).map(Number);
    plotIf(id, [{{
      type: "violin",
      x: cleanValues.map(() => label),
      y: cleanValues,
      box: {{
        visible: true,
        fillcolor: "rgba(255,255,255,0.8)",
        line: {{color: "#2980b9", width: 1}}
      }},
      meanline: {{visible: true, color: "#2980b9", width: 2}},
      points: false,
      fillcolor: "rgba(41, 128, 185, 0.8)",
      line: {{color: "#2980b9"}},
      hoveron: "violins",
      name: "",
      hovertemplate: "%{{y:,.0f}}<extra></extra>"
    }}], {{
      height: 390,
      margin: {{l: 25, r: 0, t: 40, b: 25}},
      plot_bgcolor: "white",
      paper_bgcolor: "white",
      hovermode: "closest",
      showlegend: false,
      font: {{family: "Arial, sans-serif", size: 9, color: "#2c3e50"}},
      xaxis: {{
        title: {{text: ""}},
        showgrid: false,
        zeroline: false,
        showticklabels: true,
        tickfont: {{size: 12, color: "#2c3e50"}}
      }},
      yaxis: {{
        title: {{text: ""}},
        showgrid: true,
        gridcolor: "rgba(0,0,0,0.1)",
        zeroline: false
      }},
      annotations: [{{
        text: "Count",
        x: 0,
        y: 1,
        yshift: 30,
        xref: "paper",
        yref: "paper",
        showarrow: false,
        xanchor: "left",
        yanchor: "top",
        font: {{family: "Inter", size: 12, color: "#2c3e50"}}
      }}]
    }});
  }}
  plotTemplateViolin("violin-reads", perCell.reads, "Reads");
  plotTemplateViolin("violin-umis", perCell.umis, "UMIs");
  plotTemplateViolin("violin-genes", perCell.genes, "Genes");

  const rankData = [];
  if (payload.barcodeRank5p) {{
    const rankMask = payload.barcodeRank5p.is_true || [];
    rankData.push({{
      x: payload.barcodeRank5p.rank,
      y: payload.barcodeRank5p.count.map((value, index) => rankMask[index] ? value : null),
      type: "scattergl",
      mode: "lines",
      name: "TRUE",
      connectgaps: false,
      line: {{color: "#1358A2", width: 3, simplify: false}},
      hovertemplate: "TRUE<br>Rank: %{{x:,.0f}}<br>Read counts: %{{y:,.0f}}<extra></extra>"
    }});
    rankData.push({{
      x: payload.barcodeRank5p.rank,
      y: payload.barcodeRank5p.count.map((value, index) => rankMask[index] ? null : value),
      type: "scattergl",
      mode: "lines",
      name: "NOISE",
      connectgaps: false,
      line: {{color: "#DDDDDD", width: 3, simplify: false}},
      hovertemplate: "NOISE<br>Rank: %{{x:,.0f}}<br>Read counts: %{{y:,.0f}}<extra></extra>"
    }});
  }}
  plotIf("barcode-rank", rankData, {{
    height: 360,
    margin: {{l: 30, r: 20, t: 36, b: 40}},
    plot_bgcolor: "rgba(0,0,0,0)",
    paper_bgcolor: "white",
    hovermode: "closest",
    xaxis: {{
      type: "log",
      title: {{text: "Barcode in Rank-descending Order", standoff: 12, font: {{size: 12, color: "#2c3e50"}}}},
      tickfont: {{size: 10, color: "#2c3e50"}},
      color: "black",
      showline: true,
      ticks: "outside",
      ticklen: 4,
      tickwidth: 1,
      tickcolor: "black",
      showgrid: true,
      gridcolor: "lightgrey",
      linewidth: 1,
      fixedrange: true,
      automargin: false,
      linecolor: "black"
    }},
    yaxis: {{
      type: "log",
      title: {{text: ""}},
      tickfont: {{size: 10, color: "#2c3e50"}},
      color: "black",
      showline: true,
      ticks: "outside",
      ticklen: 4,
      tickwidth: 1,
      showgrid: true,
      gridcolor: "lightgrey",
      linewidth: 1,
      fixedrange: true,
      automargin: false,
      linecolor: "black"
    }},
    annotations: [{{
      text: "Read Counts",
      x: 0,
      y: 1,
      yshift: 30,
      xref: "paper",
      yref: "paper",
      showarrow: false,
      xanchor: "left",
      yanchor: "top",
      font: {{family: "Inter", size: 12, color: "#2c3e50"}}
    }}],
    showlegend: true,
    legend: {{
      orientation: "h",
      x: 0.9,
      y: 1,
      xanchor: "center",
      yanchor: "bottom",
      itemsizing: "constant",
      font: {{family: "Arial", size: 10, color: "black"}},
      bgcolor: "rgba(255,255,255,0)",
      borderwidth: 0
    }}
  }});
  const rankNote = document.getElementById("barcode-rank-note");
  if (rankNote) {{
    const truePoints = payload.barcodeRank5p ? payload.barcodeRank5p.true_points : 0;
    const noisePoints = payload.barcodeRank5p ? payload.barcodeRank5p.noise_points : 0;
    rankNote.textContent = `5' barcode ranks: TRUE ${{truePoints.toLocaleString()}}, NOISE ${{noisePoints.toLocaleString()}}.`;
  }}

  const beadPalette = [
    "rgba(15, 76, 129, 0.8)",
    "rgba(13, 96, 158, 0.8)",
    "rgba(22, 110, 172, 0.8)",
    "rgba(46, 124, 188, 0.8)",
    "rgba(69, 138, 202, 0.8)",
    "rgba(93, 152, 216, 0.8)",
    "rgba(117, 166, 230, 0.8)",
    "rgba(140, 180, 243, 0.8)",
    "rgba(163, 194, 255, 0.8)"
  ];
  const beadData = (payload.beads.x || []).map((x, i) => ({{
    x: [x],
    y: [payload.beads.y[i] || 0],
    type: "bar",
    width: 0.8,
    name: `${{x}} ${{payload.beads.y[i] || 0}}`,
    marker: {{
      color: beadPalette[i % beadPalette.length],
      line: {{color: beadPalette[i % beadPalette.length].replace("0.8", "1"), width: 1}}
    }},
    hovertemplate: `<b>${{x}} beads per droplet</b><br>Count: %{{y:,.0f}}<br><extra></extra>`
  }}));
  const beadMax = Math.max(1, ...(payload.beads.x || []).map(Number));
  plotIf("beads-per-droplet", beadData, {{
    height: 360,
    margin: {{l: 30, r: 20, t: 36, b: 40}},
    xaxis: {{
      title: {{text: "Number of beads per droplet", font: {{family: "Inter", size: 12, color: "#2c3e50"}}}},
      showgrid: false,
      zeroline: false,
      tickfont: {{size: 10, color: "#2c3e50"}},
      range: [0.3, beadMax + 0.7],
      dtick: 1,
      tick0: 1
    }},
    yaxis: {{
      title: {{text: ""}},
      showgrid: true,
      gridcolor: "rgba(0,0,0,0.1)",
      zeroline: false,
      tickfont: {{size: 10, color: "#2c3e50"}}
    }},
    plot_bgcolor: "rgba(0,0,0,0)",
    paper_bgcolor: "rgba(0,0,0,0)",
    annotations: [{{
      text: "Count",
      xref: "paper",
      yref: "paper",
      x: 0,
      y: 1,
      showarrow: false,
      font: {{family: "Inter", size: 14, color: "#2c3e50"}},
      xanchor: "left",
      yanchor: "top",
      yshift: 30
    }}],
    showlegend: true,
    legend: {{
      x: 0.98,
      y: 0.98,
      xanchor: "right",
      yanchor: "top",
      bgcolor: "rgba(255,255,255,0.9)",
      bordercolor: "rgba(0,0,0,0.2)",
      borderwidth: 1,
      font: {{size: 10}}
    }},
    font: {{family: "Arial, sans-serif"}},
    bargap: 0.2,
    bargroupgap: 0.1
  }});
  const beadsNote = document.getElementById("beads-note");
  if (beadsNote) beadsNote.textContent = `${{(payload.beads.n_cells || 0).toLocaleString()}} cells summarized from read_assigned_cell.csv.`;

{barnyard_script}

  const rnaCells = payload.rnaCluster || [];
  const rnaTraceType = rnaCells.length > 2000 ? "scattergl" : "scatter";
  const clusterLabels = [...new Set(rnaCells.map(row => row.leiden))].sort((a, b) => {{
    if (a === "unassigned") return 1;
    if (b === "unassigned") return -1;
    return String(a).localeCompare(String(b), undefined, {{numeric: true}});
  }});
  const clusterPalette = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#13c2c2", "#bdb76b", "#6a5acd", "#20b2aa", "#ff6347"];
  const clusterTraces = clusterLabels.map((label, index) => {{
    const rows = rnaCells.filter(row => row.leiden === label);
    return {{
      x: rows.map(row => row.UMAP_1),
      y: rows.map(row => row.UMAP_2),
      text: rows.map(row => row.cell),
      customdata: rows.map(row => row.total_counts),
      type: rnaTraceType,
      mode: "markers",
      name: label === "unassigned" ? "Unassigned" : `Cluster ${{label}}`,
      marker: {{
        color: label === "unassigned" ? "#B7BDC5" : clusterPalette[index % clusterPalette.length],
        size: 3,
        opacity: 0.7
      }},
      hovertemplate: "Cell: %{{text}}<br>Cluster: " + label + "<br>UMIs: %{{customdata:,.0f}}<extra></extra>"
    }};
  }});
  const umapLayout = {{
    height: 400,
    margin: {{l: 40, r: 20, t: 20, b: 60}},
    plot_bgcolor: "rgba(0,0,0,0)",
    paper_bgcolor: "rgba(0,0,0,0)",
    font: {{family: "Arial, sans-serif"}},
    xaxis: {{
      title: {{text: "UMAP_1", font: {{size: 12, color: "#2c3e50"}}}},
      showgrid: true,
      gridcolor: "lightgray",
      zeroline: false,
      tickfont: {{size: 10, color: "#2c3e50"}}
    }},
    yaxis: {{
      title: {{text: "UMAP_2", font: {{size: 12, color: "#2c3e50"}}}},
      showgrid: true,
      gridcolor: "lightgray",
      zeroline: false,
      tickfont: {{size: 10, color: "#2c3e50"}}
    }}
  }};
  plotIf("rna-cluster-assignment", clusterTraces, Object.assign({{}}, umapLayout, {{
    showlegend: true,
    legend: {{
      title: {{text: "Cluster", font: {{size: 12}}}},
      font: {{size: 10, family: "Arial"}},
      itemsizing: "constant",
      x: 1.02,
      y: 1,
      xanchor: "left",
      yanchor: "top"
    }}
  }}));
  plotIf("rna-umi-counts", [{{
    x: rnaCells.map(row => row.UMAP_1),
    y: rnaCells.map(row => row.UMAP_2),
    text: rnaCells.map(row => row.cell),
    customdata: rnaCells.map(row => row.leiden),
    type: rnaTraceType,
    mode: "markers",
    showlegend: false,
    marker: {{
      color: rnaCells.map(row => row.total_counts),
      colorscale: "Viridis",
      size: 3,
      opacity: 0.8,
      colorbar: {{
        title: {{text: "nUMI", font: {{size: 12}}}},
        thickness: 15,
        len: 0.8,
        x: 1.02,
        tickfont: {{size: 10}}
      }}
    }},
    hovertemplate: "Cell: %{{text}}<br>Cluster: %{{customdata}}<br>UMIs: %{{marker.color:,.0f}}<extra></extra>"
  }}], Object.assign({{}}, umapLayout, {{
    margin: {{l: 40, r: 30, t: 20, b: 60}},
    showlegend: false
  }}));

  const sat = payload.saturation || {{}};
  plotIf("saturation-genes", [{{x: sat.reads_per_cell || [], y: sat.genes_per_cell || [], type: "scatter", mode: "lines", line: {{color: "#1358A2", width: 3}}, showlegend: false, hovertemplate: "x=%{{x}}<br>y=%{{y}}<extra></extra>"}}], {{
    height: 400,
    title: {{text: "Median Genes per Cell", font: {{size: 14}}}},
    xaxis: {{title: {{text: "Mean Reads per Cell", font: {{size: 12}}}}, tickformat: "~s", hoverformat: ",.0f"}},
    yaxis: {{title: {{text: "Median Genes per Cell", standoff: 0}}}}
  }});
  plotIf("saturation-umis", [{{x: sat.reads_per_cell || [], y: sat.umis_per_cell || [], type: "scatter", mode: "lines", line: {{color: "#1358A2", width: 3}}, showlegend: false, hovertemplate: "x=%{{x}}<br>y=%{{y}}<extra></extra>"}}], {{
    height: 400,
    title: {{text: "Median UMI counts per cell", font: {{size: 14}}}},
    xaxis: {{title: {{text: "Mean Reads per Cell", font: {{size: 12}}}}, tickformat: "~s", hoverformat: ",.0f"}},
    yaxis: {{title: {{text: "Median UMI Counts per Cell"}}}}
  }});
  plotIf("saturation-rate", [{{x: sat.reads_per_cell || [], y: (sat.saturation || []).map(v => v * 100), type: "scatter", mode: "lines", line: {{color: "#1358A2", width: 3}}, showlegend: false, hovertemplate: "x=%{{x}}<br>y=%{{y}}<extra></extra>"}}], {{
    height: 400,
    title: {{text: "Sequencing saturation", font: {{size: 14}}}},
    xaxis: {{title: {{text: "Mean Reads per Cell", font: {{size: 12}}}}, tickformat: "~s", hoverformat: ",.0f"}},
    yaxis: {{title: {{text: "Sequencing Saturation"}}, range: [0, 100]}}
  }});
</script>
"""
    template = load_report_template()
    return (
        template
        .replace("__SAMPLE__", html.escape(args.sample_id))
        .replace("__PLOTLY_LOADER__", embedded_plotly_loader())
        .replace("__REPORT_BODY__", report_body)
    )


def main():
    args = parse_args()
    report_df = read_tsv(args.report_metrics_tsv)
    qc_df = read_tsv(args.rna_qc_metrics_tsv)
    saturation_df = read_tsv(args.saturation_tsv)
    params = read_parameters(args.parameters_tsv)
    read_qc = safe_read_json(args.read_qc_json)
    rna_clusters = rna_cluster_payload(args.rna_cluster_tsv)
    barnyard = barnyard_payload(
        args.barnyard_summary_tsv,
        args.barnyard_per_cell_tsv,
    )

    read_summary_rows, total_reads = build_read_summary(report_df, qc_df, args.skip_glycine, args.glycine_stats)
    per_cell = per_cell_payload(args.per_cell_qc_tsv)
    summary_violin_rows = summary_violin_cards([
        {
            "title": "Reads per cell",
            "plot_id": "violin-reads",
        },
        {
            "title": "UMIs per cell",
            "plot_id": "violin-umis",
        },
        {
            "title": "Genes per cell",
            "plot_id": "violin-genes",
        },
    ])

    summary = summary_rows(report_df, qc_df)
    summary_card_rows = []
    for source_label, display_label in [
        ("Estimated cells", "Estimated number of cells"),
        ("UMIs per cell (median)", "Median UMI counts per cell"),
        ("Genes per cell (median)", "Median genes per cell"),
        ("Reads per cell (mean)", "Mean reads per cell"),
    ]:
        value = next((row["Value"] for row in summary if row["Metric"] == source_label), "NA")
        summary_card_rows.append({"label": display_label, "value": value})

    payload = {
        "readSummary": {
            "labels": [row["Metric"] for row in read_summary_rows if row["Metric"] not in {"Raw reads", "Clean reads"}],
            "counts": [value_to_float(row["Read count"]) or 0 for row in read_summary_rows if row["Metric"] not in {"Raw reads", "Clean reads"}],
        },
        "readQc": read_qc,
        "perCell": per_cell,
        "saturation": dataframe_payload(saturation_df) if not saturation_df.empty else {},
        "barcodeRank5p": barcode_rank_payload(args.barcode_counts_5p_tsv, args.whitelist_5p),
        "beads": beads_per_droplet_payload(args.read_assigned_cell),
        "rnaCluster": rna_clusters,
    }
    if barnyard is not None:
        payload["barnyard"] = barnyard

    sections = {
        "sample_title": pbmc_title_block(
            "Summary",
            "sample-information",
            help_dl([
                ("Estimated number of cells", "The number of barcodes identified as real cells in the sequencing data after barcode merging and cell calling."),
                ("Median UMI counts per cell", "The median number of unique molecular identifiers detected per cell among all identified cells."),
                ("Median genes per cell", "The median number of unique genes with detectable expression in each cell."),
                ("Mean reads per cell", "The average number of sequencing reads per cell, calculated by dividing the total reads assigned to cells by the number of identified cells."),
            ]),
        ),
        "summary_cards": pbmc_summary_cards(summary_card_rows),
        "summary_bar": pbmc_section_bar(
            "Summary",
            "sumary-detail",
            help_dl([
                ("Left", "Core run-level summary values calculated from the final cell set."),
                ("Right", "Per-cell read, UMI, and gene distributions drawn from per-cell QC outputs."),
            ]),
            width_px=190,
        ),
        "summary_rows": pbmc_metric_rows(summary),
        "summary_violin_plots": summary_violin_rows,
        "read_qc_bar": pbmc_section_bar(
            "Read QC",
            "read-qc-detail",
            help_dl([
                ("Read quality", "Distribution of per-read mean base quality scores."),
                ("Read length", "Distribution of full-length read lengths."),
                ("Base yield above read length", "Total retained bases after applying each minimum read-length cutoff."),
            ]),
            width_px=180,
        ),
        "beads_bar": pbmc_section_bar(
            "Beads to cells",
            "bead-detail",
            help_dl([
                ("Left", "Barcode rank plot showing abundance-ranked corrected 5' barcodes."),
                ("Right", "Distribution of the number of unique corrected barcodes associated with each final cell."),
            ]),
            width_px=230,
        ),
        "rna_cluster_bar": pbmc_section_bar(
            "RNA Cluster Analysis",
            "rna-cluster-detail",
            help_dl([
                ("RNA clusters", "Leiden communities calculated from normalized gene-expression profiles and displayed on a Scanpy UMAP embedding."),
                ("UMI count", "The same UMAP embedding colored by each cell's raw total UMI count."),
                ("Unassigned", "Final cells with zero detected gene-expression UMIs are retained in the report but excluded from Scanpy model fitting."),
            ]),
            width_px=300,
        ),
        "sequencing_bar": pbmc_section_bar(
            "Sequencing",
            "sequence-detail",
            help_dl([
                ("Raw reads", "Total reads in --fastq when Glycine runs, or in --full-length-fastq when Glycine is skipped. This is the denominator for every percentage in the table."),
                ("Clean reads", "Input reads retained after glycine filtering when enabled; otherwise equal to full-length reads."),
                ("Full length", "Reads retained as full-length cDNA reads before downstream cell assignment."),
                ("Barcode-valid / Cell-assigned / Gene assigned / Transcript assigned", "Read counts carried through barcode correction, cell assignment, gene assignment, and transcript assignment."),
            ]),
            width_px=180,
        ),
        "read_summary_rows": pbmc_metric_rows_3col(read_summary_rows),
        "mapping_bar": pbmc_section_bar(
            "Mapping & Annotation",
            "mapping-detail",
            help_dl([
                ("Aligned BAM reads / total reads", "Fraction of total reads that are present in the aligned BAM output."),
                ("Unmapped / total reads", "Fraction of total reads that remain unmapped."),
                ("Unique genes / Unique isoforms", "Numbers of unique genes and isoforms detected in the final outputs."),
            ]),
            width_px=320,
        ),
        "mapping_rows": pbmc_metric_rows(mapping_rows(report_df, qc_df, total_reads)),
        "saturation_bar": pbmc_section_bar(
            "Saturation",
            "saturation-detail",
            help_dl([
                ("Median genes per cell", "Downsampling preview of median genes per cell using the same known-gene definition as the summary metric."),
                ("Median UMI counts per cell", "Downsampling preview of median UMI counts per cell."),
                ("Sequencing saturation", "Estimated saturation as sequencing depth increases."),
            ]),
            width_px=190,
        ),
    }
    if barnyard is not None:
        sections["barnyard_section"] = barnyard_report_section(barnyard)
    html_text = build_html(args, payload, sections)
    Path(args.output_html).write_text(html_text, encoding="utf-8")
    print(f"Wrote report: {args.output_html}")


if __name__ == "__main__":
    main()
