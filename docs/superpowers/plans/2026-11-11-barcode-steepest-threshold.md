# Barcode Whitelist Threshold Implementation Plan

Historical implementation record for strict commit 6782643. The user-approved
relaxation in `2026-11-11-whitelist-main-window-relaxation.md` supersedes the
strict selection guards described below; current behavior is documented in
`docs/barcode-steepest-threshold.md`.

**Goal:** Add an opt-in steepest-descent threshold for pre-merge 3p/5p observed whitelists, preserving default main behavior and all downstream correction/pair/graph rules.

**Baseline:** main 291de65; branch codex/1111-barcode-steepest-threshold. Rust baseline: 88 library tests passed. Python baseline uses the existing Python 3.11 test environment and a pysam-backed samtools dispatcher, not the system Python without dependencies.

**Design:** One Rust implementation serves both pipeline and `flora barcode-thresholds` offline preview. Default `--whitelist-method legacy` retains Q95(top exp_cells)/20 and strict `>` comparison. Opt-in `steepest` uses inclusive `>=` counts at the selected rank, log-log centered secants, half-windows 50/100/200/500 by default, both endpoints >100, cross-window rank spread <=0.10, peak/background >=3, count drop >=1.5 and broad-tie rejection. All boundaries are numerical heuristics, not cell probabilities or FDR. Both sides must pass; diagnostics are written before an explicit failure. Never silently fall back.

**Scope:** No changes to barcode extraction, quality filtering, theoretical whitelist/orientation, correction, pair support, UMI pruning, graph membership, final cell calling, alignment or assembly. No additional environment dependencies. Not based on the alignment/assembly performance branch.

- [x] Add failing tests for threshold edge cases, legacy equivalence, NGS-reference parity, preview validation and actual plot thresholds.
- [x] Implement shared Rust threshold module, pipeline diagnostics and an offline preview CLI with safe output handling.
- [x] Wire both workflow runners, embedded CLI, parameter records and plotting to authoritative diagnostics. Preserve legacy Python fallback behavior or fail explicitly for unsupported new-mode execution.
- [x] Document opt-in full-run and offline preview commands in English/Chinese; keep old mode as default.
- [x] Run focused/full Python and Rust tests, shell checks, synthetic pipeline parity and independent code review; fix findings and record actual limitations.

Validation will explicitly cover flat/short/zero curves, power laws, low-count tails, multi-transition instability, ties at the cutoff, invalid options, side asymmetry and empty whitelist intersections. Offline preview and production must produce identical thresholds and observed whitelist files for identical counts. Real-data biological validation remains pending on the Linux server.

Design review clarified that rejection guards apply only to steepest, the theoretical-whitelist intersection and orientation must be identical in preview/production, and failed reruns must not mix with old results. Steepest therefore refuses existing full-workflow output directories/nonempty standalone outputs. Legacy preserves empty-side handling and exp_cells=0 normalization.

Local differential evidence: 10 independently evaluated curves match the NGS reference's thresholds, main ranks and accept/reject decisions. Default legacy mode matches original-main executable outputs byte-for-byte for whitelists, read/cell assignment, retained pairs, clean read table, cell read counts and assignment statistics. A 3,968-read synthetic FASTQ produces byte-identical preview/production threshold reports and whitelists. No real server FASTQ was run locally.

Independent review found three boundary defects, reproduced by regression tests and fixed: concurrent empty-output reuse is prevented by atomic directory/marker creation; previous Rust reports are invalidated before a Python fallback; plotting whitelist suffix/case/unknown-base normalization matches Rust. Both authoritative-report plots were actually rendered and inspected locally.

Final local verification: 113 Rust tests passed across all targets; 172 Python tests run with 3 pre-existing explicit skips (external assembly toolchain, opt-in workflow parity and reference executable setup). The 19 new threshold tests ran without skips, including original-main executable output parity and a concurrent FIFO/input reservation check. Both shell runners and the release-build script passed syntax checks. Independent re-review found no remaining actionable issues in the reviewed changes. Linux real-data validation is still pending.
