# BAM and Gene-stage Acceleration Design

## Goal

Reduce the wall time of Flora Steps 3 and 4 while preserving the existing
analysis results, output ordering, BAM tags, assignment categories, matrices,
and default output files.

## Compatibility contract

The optimized implementation must preserve:

- CB/CR/UR/C5/C3/U5/U3 tag values and the set/order of retained alignments.
- The current `bedtools bamtobed` BED6 semantics: one unsplit interval per
  mapped BAM record, zero-based start, reference-consuming end (including CIGAR
  `D` and `N`), MAPQ score, read name, and alignment strand.
- Gene assignment status, score, and gene label for every retained alignment.
- GN and UB tags, including `NA` regional fallback naming.
- `cell_gene_max_reads` truncation behavior and encounter order.
- Directional UMI clustering results independent of thread count.
- Final BAM record order, indexes, cell-UMI-gene table, and gene matrix.

## Rollout order

1. Add internal timings to identify BAM read, clustering, write, and index time.
2. Add bounded HTSlib reader/writer/index threads without changing records.
3. Parallelize independent cell-gene UMI groups with deterministic collection.
4. Add a fused Rust CB/UR/GN tagging path that directly converts BAM records to
   the same unsplit BED6 model, while retaining the legacy path and a validation
   mode for exact comparison on real data.
5. Optimize barcode correction only after downstream equivalence is proven.

## Safety and validation

- Every behavior change starts with a failing regression test.
- Thread-count equivalence tests compare complete correction maps and BAM tags.
- CIGAR conversion tests cover match, soft clipping, insertion, deletion,
  skipped reference, reverse strand, and unmapped records.
- The legacy path remains selectable for rollback.
- Validation mode generates legacy and fused results and fails on the first
  assignment or BAM-tag mismatch.
- A real-data release is not considered validated until final matrices and BAM
  tag streams are compared, not merely summary counts.

## Memory bounds

Parallelism is limited to UMI groups inside one chromosome. Chromosomes are not
processed concurrently, and group input maps are shared by reference. This
avoids multiplying the chromosome-scale read-state memory by the thread count.

