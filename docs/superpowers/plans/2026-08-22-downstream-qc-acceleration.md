# Downstream and QC Acceleration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reduce Flora Steps 3, 4, and 6 wall time without changing analytical results, while fixing Glycine rescued-read slice panics.

**Architecture:** Add measured stage boundaries first, then make four independently testable changes: safe Glycine rescue ranges, reusable lossless Read-QC summaries, compact Rust saturation, and a combined cell-table/gene-matrix BAM scan. Keep existing Python and standalone Rust commands as compatibility fallbacks.

**Tech Stack:** Rust, rust-htslib, csv/serde, Bash workflow drivers, Python/pandas reference scripts, pytest, Cargo tests.

---

Before running Cargo commands, select a Python 3.11 interpreter that provides
headers for `pyo3`:

```bash
export PYTHON_BIN="${PYTHON_BIN:-$(command -v python3.11 || command -v python3)}"
"${PYTHON_BIN}" --version
```

On the maintainer Mac worktree, the already validated compatibility wrapper is
`/tmp/flora-python311-build-wrapper`; set `PYTHON_BIN` to that path when a native
Python 3.11 is unavailable.

### Task 1: Commit the approved design and establish the baseline

**Files:**
- Add: `docs/superpowers/specs/2026-08-22-downstream-qc-acceleration-design.md`
- Add: `docs/superpowers/plans/2026-08-22-downstream-qc-acceleration.md`

- [ ] **Step 1: Run the unchanged baseline tests**

Run: `cargo test --locked`

Expected: all existing Rust tests pass.

- [ ] **Step 2: Run workflow contract tests**

Run: `python3 -m unittest tests.test_workflow_contract tests.test_workflow_parity tests.test_release_layout`

Expected: all available tests pass; unavailable local Python packages are recorded rather than hidden.

- [ ] **Step 3: Commit the reviewed design and plan**

```bash
git add docs/superpowers/specs/2026-08-22-downstream-qc-acceleration-design.md \
        docs/superpowers/plans/2026-08-22-downstream-qc-acceleration.md
git commit -m "docs: plan downstream and QC acceleration"
```

### Task 2: Add stage timing without changing execution

**Files:**
- Modify: `run_all.sh`
- Modify: `run_all_mixed_species.sh`
- Modify: `tests/test_workflow_contract.py`

- [ ] **Step 1: Add failing contract assertions**

Assert both workflows log `Stage <name> (rust|python) elapsed: <n>s` from `run_stage`, preserve the invoked command's exit status through `tee`, and emit named timings for Glycine, the complete `minimap2 | samtools view | samtools sort` pipeline, `bedtools bamtobed`, `samtools index`, and RNA clustering.

- [ ] **Step 2: Run the focused test and confirm failure**

Run: `python3 -m unittest tests.test_workflow_contract`

Expected: FAIL because implementation/timing markers do not yet exist.

- [ ] **Step 3: Implement the timing wrapper**

Capture start/end epoch seconds, select `stage_impl=rust|python`, execute the same argument vector, log timing even on failure, and return the original status. Add named timing around direct external commands without altering their pipelines.

- [ ] **Step 4: Run focused and full contract tests**

Run: `python3 -m unittest tests.test_workflow_contract tests.test_workflow_parity`

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add run_all.sh run_all_mixed_species.sh tests/test_workflow_contract.py
git commit -m "perf: add downstream stage timing"
```

### Task 3: Reject invalid Glycine rescue slices

**Files:**
- Modify: `src/glycine/identifier.rs`

- [ ] **Step 1: Add failing range tests without changing production code**

Add tests through the existing identifier entry point for valid plus/minus rescue candidates and inverted, equal, underflowing, and out-of-read coordinates corresponding to the observed panics.

- [ ] **Step 2: Confirm the tests fail before implementation**

Run: `cargo test glycine::identifier::tests --locked`

Expected: FAIL because safe range construction is absent.

- [ ] **Step 3: Implement checked rescue ranges**

Extract a helper returning `Option<Range<usize>>`; use `checked_add`/`checked_sub`, require `start < end && end <= read_len`, and only construct rescued sequence/quality records for valid windows. Invalid candidates retain the existing non-rescued classification path and never panic.

- [ ] **Step 4: Verify Glycine and full Rust tests**

Run: `cargo test glycine --locked`

Run: `cargo test --locked`

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/glycine/identifier.rs
git commit -m "fix: guard Glycine rescued-read ranges"
```

### Task 4: Reuse lossless Read-QC summaries

**Files:**
- Add: `src/read_qc.rs`
- Modify: `src/lib.rs`
- Modify: `src/bin/read_qc_summary.rs`
- Modify: `src/pipeline.rs`
- Modify: `src/glycine/mod.rs`
- Modify: `src/glycine/batch.rs`
- Modify: `run_all.sh`
- Modify: `run_all_mixed_species.sh`
- Modify: `tests/test_workflow_contract.py`
- Modify: `tests/test_workflow_parity.py`

- [ ] **Step 1: Add one failing single-versus-merged accumulator test**

Test that observing a fixed record set at once and merging two chunk accumulators produce identical payloads, including empty reads, trimmed quality means, medians, 99.5% length quantile, histogram edges/counts, N50, and yield curves. Assert the accumulator stores a length frequency map and one `f32` quality per non-empty read, and that merge consumes rather than clones the source.

Freeze the current `read_qc_summary` JSON for the same FASTQ before moving the
implementation, then require exact integer/array parity and `1e-6` floating
parity. Add a 256-entry thread-safe static Phred error-probability lookup test
against the existing per-base `powf` formula for every byte value.

- [ ] **Step 2: Run the focused accumulator test red, implement, then pass**

Run: `cargo test read_qc::tests::merged_payload_matches_single_pass --locked`

Expected before implementation: FAIL. Implement only `ReadQcAccumulator` and its merge; rerun and expect PASS.

- [ ] **Step 3: Add failing atomic publication and cache-validation tests**

Use two-plus FASTQ chunks with deliberately non-lexical names (`chunk_2`, `chunk_10`). For explicit `--fastq` and `--fastq-dir`, require exact merged FASTQ record order, summed `identifying_statistic.txt` counts, one final JSON/count pair, per-worker temporary paths, and atomic rename. Add independent cache rejection cases for schema version, sample ID, source size/mtime fingerprint, read count, malformed JSON, and missing file. Require failed workers to leave no reusable final summary.

- [ ] **Step 4: Run each publication/cache case red, implement minimally, then pass**

Run: `cargo test glycine::batch::tests::merged_read_qc_preserves_natural_chunk_order --locked`

Run: `cargo test read_qc::tests::cache_validation --locked`

Expected before each corresponding implementation: FAIL; expected afterward: PASS.

- [ ] **Step 5: Accumulate in existing FASTQ passes**

For `--skip-glycine`, observe records in the Flora upstream FASTQ pass. For Glycine, summarize emitted full-length-plus-rescued records per chunk and merge in natural chunk order. Do not add another FASTQ scan.

- [ ] **Step 6: Add and pass workflow reuse/fallback cases one at a time**

Add separate contract cases for single/mixed species, light/full output, `--skip-glycine`, valid reuse, and every invalid-cache fallback. After each new assertion fails for the expected reason, implement only that path and rerun the focused test.

- [ ] **Step 7: Use validated summaries in Step 6**

Copy/reuse valid precomputed artifacts. Log the reason and call `read_qc_summary` only when absent or invalid. Preserve the standalone command and add `--output-fastq-count` compatibility.

- [ ] **Step 8: Run focused and full tests**

Run: `cargo test read_qc --locked`

Run: `python3 -m unittest tests.test_workflow_contract tests.test_workflow_parity`

Run: `cargo test --locked`

Expected: PASS.

- [ ] **Step 9: Commit**

```bash
git add src/read_qc.rs src/lib.rs src/bin/read_qc_summary.rs src/pipeline.rs \
        src/glycine/mod.rs src/glycine/batch.rs run_all.sh \
        run_all_mixed_species.sh tests/test_workflow_contract.py \
        tests/test_workflow_parity.py
git commit -m "perf: reuse upstream read QC summaries"
```

### Task 5: Align and reduce RNA-QC scans

**Files:**
- Modify: `src/bin/rna_qc_metrics.rs`
- Modify: `run_all.sh`
- Modify: `run_all_mixed_species.sh`

- [ ] **Step 1: Add failing CLI/count-source tests**

Test `--glycine-stats`, `--fastq-count-file`, and `--threads`; test raw/full-length same-file detection and Glycine `Total`/`Full-length` parsing.

- [ ] **Step 2: Confirm failure**

Run: `cargo test stage_rna_qc_metrics --locked`

Expected: FAIL on unsupported arguments/count reuse.

- [ ] **Step 3: Implement exact count reuse**

Prefer validated count files, use Glycine statistics for raw counts, avoid counting the same canonical FASTQ twice, and apply bounded htslib reader threads. Keep compatibility scans when artifacts are absent.

- [ ] **Step 4: Add a failing interned-summary parity test**

On a fixture containing repeated read IDs, cells, genes, UMIs, mitochondrial
genes, placeholders, and duplicate alignments, compare every metric from the
current string-set implementation with an interned-ID implementation. Assert
that per-cell structures store integer IDs and that mapped-read membership uses
one owned read-ID set rather than duplicate global string sets.

- [ ] **Step 5: Implement bounded interned summaries and pass the focused test**

Intern cell/gene/UMI text once while streaming `cell_umi_gene.tsv`; derive global
counts from per-cell/interned sets when possible. Store final-cell read IDs once
and update BAM counters without creating duplicate owned qname strings. Preserve
all existing output values and ordering.

- [ ] **Step 6: Verify**

Run: `cargo test stage_rna_qc_metrics --locked`

Run: `cargo test --locked`

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add src/bin/rna_qc_metrics.rs run_all.sh run_all_mixed_species.sh
git commit -m "perf: avoid duplicate RNA QC input scans"
```

### Task 6: Replace pandas saturation with compact Rust computation

**Files:**
- Add: `src/bin/saturation.rs`
- Modify: `src/lib.rs`
- Modify: `src/main.rs`
- Modify: `src/workflow_runtime.rs`
- Modify: `scripts/Saturation.py`
- Modify: `run_all.sh`
- Modify: `run_all_mixed_species.sh`
- Modify: `tests/test_workflow_parity.py`
- Add: `tests/fixtures/saturation/reference_input.tsv`
- Add: `tests/fixtures/saturation/reference_output.tsv`

- [ ] **Step 1: Add and pass the MT19937 permutation test**

Add the fixed NumPy seed-42 permutation and selected-prefix expectations. Confirm the missing stage fails, implement only MT19937/permutation helpers, and rerun to PASS.

- [ ] **Step 2: Add and pass half-even size and fixed-width row tests**

Add rounding-boundary cases and a compile-time `size_of::<Row>()` bound proving compact rows contain only fixed-width IDs/flags, not `String`. Run focused red/green cycles.

- [ ] **Step 3: Add and pass metric-equivalence cases**

Add placeholder-gene, duplicate-row, underscore-collision, and zero-known-gene fixtures. Generate `reference_output.tsv` once with the unchanged pandas script and compare every Rust field with exact integers and `1e-12` float tolerance.

- [ ] **Step 4: Implement compact IDs and cumulative fractions**

Intern repeated strings to `u32` IDs, preserve the existing concatenated molecule collision semantics, reproduce MT19937 permutation order, and compute all increasing fractions in one shuffled pass.

- [ ] **Step 5: Retain Python plotting only**

Add `--plot-existing-tsv` to `scripts/Saturation.py`. The workflow runs the Rust metric stage, then invokes Python plotting only if the Rust stage did not create a PNG.

- [ ] **Step 6: Update and verify stage-selection contracts**

Update `tests/test_workflow_parity.py` so optimized workflows require Rust `saturation` while preserving Python plotting. Confirm the test fails before workflow/runtime registration and passes afterward.

- [ ] **Step 7: Verify Rust and workflow contracts**

Run: `cargo test stage_saturation --locked`

Run: `python3 -m unittest tests.test_workflow_contract tests.test_workflow_parity`

Expected: PASS.

- [ ] **Step 8: Commit**

```bash
git add src/bin/saturation.rs src/lib.rs src/main.rs src/workflow_runtime.rs \
        scripts/Saturation.py run_all.sh run_all_mixed_species.sh \
        tests/test_workflow_parity.py tests/fixtures/saturation/reference_input.tsv \
        tests/fixtures/saturation/reference_output.tsv
git commit -m "perf: compute saturation in compact Rust stage"
```

### Task 7: Generate cell table and gene matrix in one BAM scan

**Files:**
- Modify: `src/bin/cell_umi_gene_table.rs`
- Modify: `run_all.sh`
- Modify: `run_all_mixed_species.sh`
- Modify: `tests/test_gene_expression_cell_parity.py`
- Modify: `tests/test_workflow_parity.py`

- [ ] **Step 1: Add failing combined-output parity test**

Run the standalone table and matrix stages and the combined stage on the same tagged BAM fixture. Compare table bytes and matrix content, including duplicate UMIs and cells containing only genomic placeholder genes.

- [ ] **Step 2: Confirm failure**

Run: `python3 -m unittest tests.test_gene_expression_cell_parity`

Expected: FAIL because `--gene-expression-output` is unsupported.

- [ ] **Step 3: Implement optional matrix output**

During the existing BAM pass, write each table row and accumulate the same unique `(gene, cell, UMI)` sets used by `gene_expression`. Preserve sorted axes and placeholder-only cell columns.

- [ ] **Step 4: Update stage selection, then remove only the redundant invocation**

First update `tests/test_workflow_parity.py` to expect `cell_umi_gene_table` with combined output and no workflow-level `gene_expression`; verify red. Then pass `--gene-expression-output` from both workflows. Retain standalone `gene_expression` for users and compatibility tests.

- [ ] **Step 5: Verify parity and all tests**

Run: `python3 -m unittest tests.test_gene_expression_cell_parity tests.test_workflow_contract tests.test_workflow_parity`

Run: `cargo test --locked`

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/bin/cell_umi_gene_table.rs run_all.sh run_all_mixed_species.sh \
        tests/test_gene_expression_cell_parity.py tests/test_workflow_parity.py
git commit -m "perf: build gene matrix during cell table scan"
```

### Task 8: Final parity and benchmark gate

**Files:**
- Modify: `docs/superpowers/plans/2026-08-22-downstream-qc-acceleration.md` (record results)
- Add: `tests/benchmark_downstream_acceleration.py`

- [ ] **Step 1: Run formatting and all local tests**

Run: `cargo fmt --check`

Run: `cargo test --locked`

Run: `python3 -m unittest discover -s tests -p 'test_*.py'`

Expected: PASS, with any environment-only skips documented.

- [ ] **Step 2: Build the release artifact**

Run: `cargo build --release --locked`

Expected: release build succeeds.

- [ ] **Step 3: Execute all fixture workflow modes**

Run:

```bash
python3 -m unittest \
  tests.test_workflow_parity.WorkflowOracleFixtureTests \
  tests.test_workflow_contract
```

The fixture harness must execute single and mixed workflows for multi-file
`--fastq`, `--fastq-dir`, light/full output, `--skip-glycine`, malformed input,
forced worker failure, and stale-cache scenarios. Expected: PASS.

- [ ] **Step 4: Compare baseline and optimized real-sample outputs**

On the Linux benchmark host, build one binary from `main` and one from the
optimization branch, then run identical commands under the existing memory
monitor:

```bash
git worktree add /tmp/flora-main-benchmark main
(cd /tmp/flora-main-benchmark && cargo build --release --locked)
cargo build --release --locked

MEMTIME=/home/liyy/2.project/C4_V3/Strint2.6/TroubleShootingScripts/MemMonitor/memtime.sh
FASTQ_DIR=/home/liyy/1.data/C4_Cyclone/260F401343011
BC=/home/liyy/2.project/C4_V3/Strint2.6/BC_1536.txt
REF=/home/liyy/1.data/REF/GRCH38

bash "$MEMTIME" /data/flora-benchmark/baseline/memory 1 \
  /tmp/flora-main-benchmark/target/release/flora \
  --fastq-dir "$FASTQ_DIR" --barcode-list-10bp "$BC" --ref-dir "$REF" \
  --out-dir /data/flora-benchmark/baseline/sample_output --sample-id sample \
  > /data/flora-benchmark/baseline/workflow.log 2>&1

bash "$MEMTIME" /data/flora-benchmark/optimized/memory 1 \
  ./target/release/flora \
  --fastq-dir "$FASTQ_DIR" --barcode-list-10bp "$BC" --ref-dir "$REF" \
  --out-dir /data/flora-benchmark/optimized/sample_output --sample-id sample \
  > /data/flora-benchmark/optimized/workflow.log 2>&1

! grep -E "panicked at|Glycine worker panicked" \
  /data/flora-benchmark/optimized/workflow.log
```

Repeat the optimized fixture-sized command once with explicit
`--fastq chunk_2.fastq.gz chunk_10.fastq.gz` to verify that input form as well.
Then run the comparison tool:

```bash
python3 tests/benchmark_downstream_acceleration.py compare \
  --baseline-output /data/flora-benchmark/baseline/sample_output \
  --optimized-output /data/flora-benchmark/optimized/sample_output \
  --sample-id sample \
  --report /data/flora-benchmark/output_parity.json
```

The tool compares sorted assignment tables, matrices, QC JSON/TSV metrics,
saturation TSV, report inputs, and Glycine statistics. Expected: exit 0 and
`"parity": true`.

- [ ] **Step 5: Apply the performance gate**

Run:

```bash
python3 tests/benchmark_downstream_acceleration.py performance \
  --baseline-log /data/flora-benchmark/baseline/workflow.log \
  --optimized-log /data/flora-benchmark/optimized/workflow.log \
  --baseline-mem /data/flora-benchmark/baseline/mem_rss.csv \
  --optimized-mem /data/flora-benchmark/optimized/mem_rss.csv \
  --report /data/flora-benchmark/performance_gate.json
```

Expected: exit 0 when each replaced stage improves at least 20% or combined
Step 4/6 improves at least 15%, and optimized peak PSS is no more than 105% of
baseline.

- [ ] **Step 6: Commit verification tooling and records**

```bash
git add tests/benchmark_downstream_acceleration.py \
        docs/superpowers/plans/2026-08-22-downstream-qc-acceleration.md
git commit -m "test: record downstream acceleration verification"
```
