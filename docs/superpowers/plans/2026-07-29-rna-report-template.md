# RNA Report Template Migration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Migrate the existing Rust pipeline RNA report content into the visual shell of `P2026042903_multi_report.html` without carrying over example or non-RNA data.

**Architecture:** Keep data extraction and chart payload generation in `scripts/build_report.py`. Replace `scripts/report_template.html` with a compact RNA-only shell derived from the new template and adapt the generated report body to its three-page navigation.

**Tech Stack:** Python 3, HTML/CSS, vanilla JavaScript, embedded Plotly.js, unittest.

---

### Task 1: Add report contract tests

**Files:**
- Create: `tests/test_build_report_template.py`
- Test: `scripts/build_report.py`
- Test: `scripts/report_template.html`

- [ ] Write tests for placeholders, RNA-only navigation, absence of example data, section preservation, and complete placeholder replacement.
- [ ] Run the tests and verify they fail against the existing PBMC template.

### Task 2: Create the RNA-only new template shell

**Files:**
- Modify: `scripts/report_template.html`

- [ ] Preserve the new template's favicon, embedded Plotly runtime, colors, typography, sidebar, sample header, cards, responsive behavior, and footer.
- [ ] Remove all example content, VDJ-T, VDJ-B, ATAC, cluster, marker, and annotation code.
- [ ] Add `__SAMPLE__` and `__REPORT_BODY__` placeholders.
- [ ] Run the template contract tests.

### Task 3: Adapt the report body to the new navigation

**Files:**
- Modify: `scripts/build_report.py`

- [ ] Keep all existing report input readers, metric calculations, section text, and plot payloads.
- [ ] Render existing sections under Summary, Cells, and Library panes in their current order.
- [ ] Replace Bootstrap-specific tab handling with the new template's button navigation.
- [ ] Run the report contract tests and Python syntax compilation.

### Task 4: Generate and inspect a fixture report

**Files:**
- Create during test only: temporary TSV/JSON/PNG/HTML fixtures.

- [ ] Generate a complete report using synthetic pipeline outputs.
- [ ] Verify all expected section titles and plot hosts are present.
- [ ] Verify example data, VDJ/ATAC terms, and unresolved placeholders are absent.
- [ ] Validate HTML structure and report generation exit status.
