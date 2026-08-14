# Flora Peak Memory Optimization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reduce Flora barcode/cell-assignment peak memory while preserving all analysis results and output row ordering from `main@0857cfa`.

**Architecture:** Keep the current analysis model, but remove full-read temporary generations. Write correction maps incrementally, consume records when converting between stages, retain summary scalars before releasing source vectors, and perform cell assignment as a count pass followed by a streaming output pass. Disk-backed chunking remains a second-stage fallback only if large-data profiling still exceeds the memory target.

**Tech Stack:** Rust 2021, Rayon, `csv`, `rustc-hash`, existing workflow fixtures.

---

### Task 1: Stream correction-map output

**Files:**
- Modify: `src/pipeline.rs`

- [x] Add a test proving correction-map rows and ordering are unchanged.
- [x] Verify the test fails against a writer API that does not yet exist.
- [x] Serialize each row directly to `csv::Writer` without collecting `Vec<CorrectionMapRow>`.
- [x] Run focused tests and commit.

### Task 2: Consume intermediate read generations

**Files:**
- Modify: `src/pipeline.rs`

- [x] Add tests for order-preserving corrected-to-clean conversion and in-place filtering.
- [x] Convert `corrected` into `CleanRead` by ownership movement rather than cloning.
- [x] Save summary counts before releasing whitelists, caches, raw counts, and `putative`.
- [x] Make pair filtering consume its input and avoid `cleaned`/`paired_final` duplicate generations.
- [x] Run focused tests and commit.

### Task 3: Stream final cell assignment

**Files:**
- Modify: `src/pipeline.rs`

- [x] Add a test that compares streamed assignment output and statistics with legacy ordering semantics.
- [x] Extract a per-read cell-resolution helper.
- [x] First pass: count assigned reads per cell and resolve the minimum-read cell set.
- [x] Second pass: serialize retained reads directly and optionally collect read IDs only when FASTQ outputs are requested.
- [x] Remove the full `Vec<AssignedRead>` and release `clean` before optional FASTQ scans.
- [x] Run focused tests and commit.

### Task 4: Regression and memory-shape verification

**Files:**
- Modify: `tests/` only if additional regression coverage is needed.

- [x] Run the complete Rust test suite with a local test-only embedded-asset compiler; retain the formal Python 3.11 requirement for release packaging.
- [x] Compare fixture upstream CSV/TSV outputs with the `main` oracle, including correction maps and read assignments.
- [x] Run `cargo fmt --check` and `git diff --check`.
- [x] Review remaining read-count-proportional containers and document the optional-output exception.
- [x] Commit verification adjustments.

### Task 5: Conditional disk-backed phase

- [ ] Re-profile the 149G dataset with one-minute PSS/RSS and stage markers.
- [ ] Implement disk-backed chunking only if Steps 1-4 do not meet the operational memory target.
