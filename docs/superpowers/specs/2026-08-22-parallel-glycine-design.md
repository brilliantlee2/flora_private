# Parallel Glycine Design

## Goal

Allow the full Flora workflows to accept either multiple FASTQ paths after one
`--fastq` option or a non-recursive `--fastq-dir`, run Glycine on the inputs
concurrently, and expose the same single sample-level Glycine outputs consumed
by the rest of Flora.

## Interface

- `--fastq sample_1.fastq.gz sample_2.fastq.gz ...`
- `--fastq-dir /path/to/chunks`
- `--glycine-jobs 10` controls the maximum number of concurrent Glycine jobs.
- `--glycine-threads 64` is the total Glycine thread budget.
- `--fastq` and `--fastq-dir` may not be combined.
- Directory discovery is non-recursive and accepts `.fastq.gz` and `.fq.gz`.
- Inputs are naturally sorted and duplicate canonical paths are rejected.

## Execution

Each input is processed in an isolated chunk directory. At most
`min(glycine_jobs, input_count)` child Glycine processes run concurrently. The
total thread budget is divided across active jobs, with at least one thread per
job. A failed child terminates the batch and prevents publication of partial
sample-level outputs.

Single-input workflows retain the same output contract and analysis semantics.
Multi-input workflows concatenate the independently compressed full-length
FASTQ members in deterministic input order without decompression or
recompression.

## Statistics

Every chunk must first produce its own `identifying_statistic.txt`. The merger
parses all three sections structurally, sums every count and base total, and
recomputes every percentage from the merged denominators. Percentages are never
averaged. The final file keeps the current Glycine text format and filename.

The per-chunk `full_length_fastq_count.txt` values are summed. Per-chunk
`read_qc_summary.json` data are merged by additive histograms/counters where the
schema permits exact reconstruction; otherwise the final merged FASTQ is
streamed once to regenerate the sample-level QC artifact rather than averaging
derived values.

## Cleanup And Diagnostics

Chunk outputs live below the Glycine output directory during execution. They
are removed after a successful merge unless intermediate retention is
requested. Logs identify each input, assigned thread count, elapsed time, and
failure status.

## Verification

Tests cover argument parsing, non-recursive natural ordering, bounded
concurrency/thread allocation, concatenated FASTQ readability, exact statistics
aggregation, failure handling, single-input compatibility, and both single- and
mixed-species runner contracts.
