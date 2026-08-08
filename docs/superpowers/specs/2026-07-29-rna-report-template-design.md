# RNA Report Template Migration Design

## Goal

Use the visual shell of `report_new_2/P2026042903_multi_report.html` for the Rust pipeline report while retaining every metric and plot currently produced by `scripts/build_report.py`.

## Boundaries

- The generated report contains RNA content only.
- VDJ-T, VDJ-B, ATAC, marker, cluster, and cell-annotation sample sections are excluded.
- No sample values or chart arrays from `P2026042903_multi_report.html` enter generated reports.
- All displayed values and plots come from the current pipeline inputs consumed by `scripts/build_report.py`.
- The `run_all.sh` report invocation and command-line interface remain unchanged.

## Structure

`scripts/report_template.html` becomes a compact, offline-capable shell based on the new report's typography, colors, sidebar, sample header, cards, responsive layout, and footer. It exposes only `__SAMPLE__` and `__REPORT_BODY__` placeholders.

`scripts/build_report.py` retains its existing data readers, metric definitions, plot payload construction, and display order. Its generated body is adapted to the new shell's Summary, Cells, and Library navigation:

1. Summary: Sample information, Summary, Read QC.
2. Cells: Beads to cells.
3. Library: Sequencing, Mapping & Annotation, Saturation.

## Verification

Automated tests check that the template has the required placeholders and RNA navigation, contains no example sample identifier or VDJ/ATAC content, and that generated HTML replaces placeholders while retaining every current report section and plot host.
