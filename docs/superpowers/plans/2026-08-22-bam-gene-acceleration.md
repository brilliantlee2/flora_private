# BAM and Gene-stage Acceleration Plan

## 1. Detailed timing

- Add tests for stable timing labels where practical.
- Instrument CB/UR tagging, gene assignment, GN tagging, and UMI clustering.
- Report load/read/cluster/write/index phases to stderr.
- Run focused tests and commit.

## 2. Controlled BAM threads

- Add failing CLI/config tests for bounded thread normalization.
- Apply reader/writer/index threads to BAM tools.
- Preserve output bytes at the record/tag level across thread counts.
- Run focused and full tests, then commit.

## 3. Deterministic parallel UMI clustering

- Add a regression fixture with multiple cell-gene groups and ties.
- Compare one-thread and multi-thread correction maps exactly.
- Use a local Rayon pool controlled by `--threads`; keep chromosome order and
  BAM writing serial.
- Run focused and full tests, then commit.

## 4. Direct BAM gene tagging

- Add CIGAR-to-unsplit-BED regression tests.
- Implement direct BAM record conversion and fused CB/UR/GN tagging.
- Add `legacy`, `rust`, and `validate` workflow modes with a rollback path.
- Compare assignment TSV rows and BAM tags/order on synthetic fixtures.
- Run shell syntax, focused, and full tests, then commit.

## 5. Barcode correction follow-up

- Profile correction subphases before changing algorithms.
- Add exact correction-table fixtures and thread-count equivalence tests.
- Optimize only allocations/caches that leave selected barcodes unchanged.
- Run full tests and commit separately.

## 6. Release verification

- Run `cargo fmt --check`, `cargo test`, and `bash -n` for both workflows.
- Review the diff for generated files or unrelated changes.
- Record commits with date, reason, validation commands, and residual real-data
  validation requirements.
