# Downstream and QC Acceleration Design

Date: 2026-08-22
Branch: `codex/0822-downstream-qc-acceleration`
Baseline: `main` at `bcf686c`

## Motivation

A 179,779,793-read run completed in 84,001 seconds (23.33 hours). The measured stage times were:

| Stage | Seconds | Hours | Share |
| --- | ---: | ---: | ---: |
| Glycine | 27,337 | 7.59 | 32.5% |
| Alignment and sorting | 9,520 | 2.64 | 11.3% |
| Barcode/cell assignment | 6,175 | 1.72 | 7.4% |
| BAM barcode/UMI tagging | 12,336 | 3.43 | 14.7% |
| Gene tagging, UMI clustering, and matrices | 14,344 | 3.98 | 17.1% |
| Isoform assignment | 2,448 | 0.68 | 2.9% |
| RNA QC, saturation, and report | 11,806 | 3.28 | 14.1% |

The next optimization target is the repeated full-file work in Steps 3, 4, and 6. The barcode correction substage remains a separate later target because it already has detailed internal timings and a different data model.

The same run also reported repeated slice-range panics in Glycine. Those panics can drop worker results, so they must be fixed before performance comparisons are considered valid.

## Behavioral Invariants

The optimized workflow must preserve:

- FASTQ classification and `identifying_statistic.txt` counts.
- Read-to-cell, barcode, UMI, gene, and transcript assignments.
- Existing TSV/CSV/JSON schemas, row semantics, and report metric definitions.
- Gene-expression and isoform-expression matrix values and axes.
- Saturation fractions and deterministic sampling seed/selection semantics.
- Single-species and mixed-species workflow behavior.
- Light-output and full-output modes.

Any deliberate byte-order difference must be documented and validated as semantically equivalent. Exact byte equality is preferred where the existing order is part of the workflow contract.

## Phase 0: Instrumentation and Correctness

### Per-stage timing

Extend `run_stage` to log the implementation and elapsed wall time:

```text
[run_all] Stage <name> (rust|python) elapsed: <seconds>s
```

Add equivalent timing around direct external commands that dominate a stage, including `bedtools bamtobed`, `samtools index`, and RNA clustering. This allows future optimization decisions to use evidence rather than aggregate Step timings.

### Glycine slice safety

Reproduce the invalid ranges around the rescued-read slices in `identifier.rs`. Add unit tests for overlapping/inverted primer and poly(A/T) coordinates. Invalid candidate windows must be rejected without panicking, matching the prior intended `irrescuable`/non-rescued behavior rather than silently losing a read.

## Phase 1: Reuse Read QC From Existing FASTQ Scans

Introduce a mergeable `ReadQcAccumulator` that records:

- read count and total bases;
- read-length distribution required by the report;
- per-read mean quality using a static Phred error-probability lookup;
- data needed for median, histogram, N50, and yield-above-length curves.

The lossless accumulator representation is a read-length frequency map plus one
`f32` mean-quality value per emitted read. This matches the current Rust
`read_qc_summary` numeric path exactly while reducing the length side from one
integer per read to one entry per observed length. At 180 million emitted reads,
the quality vector is approximately 687 MiB; chunk workers publish summaries one
at a time so their vectors are moved into the merger rather than cloned.

Chunk summaries are written to temporary paths and atomically renamed only after
the worker succeeds. The reusable final JSON is accepted only when its schema
version, source FASTQ size/mtime fingerprint, read count, and sample identifier
match the current run; otherwise Step 6 logs the mismatch and performs the
compatibility scan. Merge tests require exact integer fields and JSON array
contents, with `1e-6` tolerance only for floating-point serialization.

When Glycine runs, each FASTQ worker accumulates QC for emitted full-length-plus-rescued records and writes a chunk summary. Batch merging combines those summaries into the final read-QC JSON and count file.

When `--skip-glycine` is used, the existing Flora upstream FASTQ scan accumulates the same summary. Step 6 reuses the precomputed files. It scans the FASTQ only as a compatibility fallback when precomputed data are unavailable.

`rna_qc_metrics` must accept precomputed read counts and Glycine statistics so it does not decompress the same full-length FASTQ twice.

## Phase 2: Compact Rust Saturation

Replace the pandas implementation that materializes the full `cell_umi_gene.tsv` and performs thirteen repeated samples/group-bys.

The Rust implementation will:

- intern repeated read, gene, cell, UMI, and molecule strings into compact integer IDs;
- reproduce NumPy `RandomState(42)` MT19937 permutation semantics used by
  `DataFrame.sample(..., random_state=42, replace=False)`;
- use pandas/NumPy half-even sample-size rounding for every configured fraction;
- shuffle the compact rows once and update cumulative statistics at the increasing
  fraction boundaries in one pass;
- compute reads, unique UMI observations, known genes per cell, and saturation with the current definitions;
- write the same saturation TSV schema;
- leave plotting to the lightweight Python plotting path when needed.

The implementation must be tested against the current Python reference on fixed
fixtures, including exact selected row indices for multiple input sizes and
rounding boundaries, zero-known-gene cells, and duplicate
read/gene/barcode/UMI rows. It must also assert that the compact representation
uses fixed-width IDs rather than retaining duplicate input strings.

## Phase 3: One BAM Scan for Cell Table and Gene Matrix

Extend the Rust `cell_umi_gene_table` stage with an optional gene-expression output. During one pass over the final tagged BAM it will:

- stream `read_id`, `GN`, `CB`, and `UB` to `cell_umi_gene.tsv`;
- accumulate unique UMIs per known gene and cell;
- retain cells that only have genomic placeholder assignments;
- write the existing dense gene-expression TSV after the scan.

The standalone `gene_expression` command remains available for compatibility, but the complete workflows no longer invoke it separately.

## Phase 4: RNA QC Memory and Scan Reduction

Align the Rust `rna_qc_metrics` CLI with the workflow arguments, including Glycine statistics, precomputed FASTQ counts, and bounded worker settings.

Reduce repeated string allocation and duplicate global sets where counts can be derived during the existing table/BAM passes. The first implementation remains exact and conservative; more invasive BAM/gene-tagging fusion is deferred until instrumentation identifies the dominant remaining substage.

## Alternatives Considered

### Instrument only

Lowest risk, but it leaves known repeated scans and the existing three-hour QC stage untouched. Instrumentation is included, but not used as the stopping point.

### Fuse Steps 3 and 4 immediately

Combining barcode tagging, gene assignment, GN tagging, UMI clustering, and matrix generation could remove several BAM/TSV intermediates. It has the largest theoretical gain but also the largest result-parity risk. It is deferred until the incremental changes provide stage-level evidence.

### Selected approach

Implement correctness and instrumentation first, then independently verifiable scan elimination and Rust streaming stages. This preserves rollback points and makes each measured gain attributable to one change.

## Verification

- Unit tests must fail before each production change and pass afterward.
- Run the full Rust test suite and relevant Python contract tests.
- Compare optimized and baseline outputs on the same fixture and a real sample.
- Exercise both explicit multi-file `--fastq` input and `--fastq-dir`, preserving
  natural input ordering and merged FASTQ/statistics/QC/count artifacts.
- Exercise single- and mixed-species workflows, light/full output,
  `--skip-glycine`, worker failure cleanup, and missing, malformed, stale, or
  mismatched precomputed Read-QC fallback.
- For real data, compare key counts, sorted table content where row order is not contractual, matrix checksums/content, QC metrics, saturation TSV, and HTML report inputs.
- Record per-stage wall times and peak RSS/PSS with identical inputs and parameters.
- On the fixed benchmark, each replaced standalone stage must be at least 20%
  faster, or the combined Step 4/6 wall time must improve by at least 15%.
  Peak PSS must not increase by more than 5%; any exception requires measured
  justification before merge.
- Do not merge until Glycine runs without worker panics and all result-parity checks pass.
