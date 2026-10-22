# ONT feature alignment (2026-10-22)

Base: `codex/1013-ont-prealignment-trimming` (`b931d16`).
Branch: `codex/1022-ont-feature-alignment`.
Reference: latest stable wf-single-cell v3.3.4, checked using GitHub Releases
on 2026-10-22. Container tag in that release:
`ontresearch/wf-single-cell:sha69b55755e85ed2abf209bf19d9f72d7b5dd5f807`.

## Design and scope

Flora retains Glycine, dual-end barcode extraction, UMI extraction, cell
assignment and cDNA streaming. The default feature backend becomes `ont`;
`--feature-backend legacy` reproduces the pre-existing annotation backend.

Genome mapping uses splice/forward/MD, junction hints, 10M batches, a 100m
allocation cap, and secondary-record filtering, as in ONT preprocess.nf.
Use a single genome index (-I 16G), as in ONT common.nf. Supplementary
alignments remain in the genomic BAM.

The ONT backend assembles and assigns independently per annotated chromosome:

1. StringTie `-L -c 2`, guided by the reference GTF, using all mapped cDNA
   alignments rather than only reads assigned to final cells.
2. gffread extracts the sample transcriptome; gffcompare maps assembled
   transcripts to reference transcripts and gene names.
3. Re-extract primary mapped cDNA reads from genomic BAM and remap using
   minimap2 `map-ont --end-bonus 10 -p 0.9 -N 3 --cap-kalloc 100m
   --cap-sw-mem 50m`; keep secondary hits and discard unmapped/supplementary.
4. Joint feature assignment follows ONT v3.3.4 score/coverage rules. Genomic
   MAPQ >=30 is required for gene assignment; transcript assignment has its
   own transcriptome ambiguity rules and is not gated by genomic MAPQ.
   Classes i/y/p/s are excluded from transcript counting. Unmatched novel
   transcripts are exported in the assembly, not promoted into count matrices.
   Join reference annotations before determining unique/multiple candidates;
   a unique candidate needs transcript MAPQ >0 and has no coverage cutoff.
   Multiple candidates rank by AS, query coverage and transcript coverage.
   Distinct AS/query coverage requires both coverages >=0.4; otherwise use
   distinct transcript coverage >0.8. Excluded classes still participate in
   ranking, and keep their gene label.
5. Feed backend gene labels into Flora's existing directional UMI correction.
   Reuse corrected UMIs for transcript counts. Keep existing TSV/QC/report
   interfaces, and additionally export raw and processed MEX matrices.

Every chromosome uses a unique StringTie label. Feature lookup uses
`(chromosome, read_id)`; the filtered counting BAM contains primary mapped
records only, while the original genome BAM retains supplementary records.
Duplicate primary read names in the tagged population are rejected. Translate
ONT's missing-feature sentinel to Flora's `NA` for positional UMI grouping;
exclude generated region labels from gene matrices. Transcript-only reads are
retained independently of gene counting. This intentionally avoids the upstream
gene-matrix builder's in-place removal of unknown genes before transcript counts.

Additional differences: processed MEX recognizes both MT- and mt- identifiers;
the report and existing TSVs retain Flora's cell population, not the extra MEX
filters. ONT UMI mode has no per-cell/gene read cap by default (explicit caps
remain available), but retains Flora's extracted UMI representation and its
existing directional edit-distance correction engine.

Implementation is independent Flora code; ONT source is used as a behavioral
reference, not bundled. Exact numerical equivalence with an ONT run is not
claimed: library preparation, cell assignment, UMI representation, upstream
read population and executable versions can differ. External versions are
recorded per run rather than inferred from the ONT container tag.

## Resource and failure handling

Process chromosomes sequentially to bound assembly/remapping memory. Use a
temporary SQLite index for read assignments and genome MAPQ rather than a
whole-sample in-memory dictionary. Remove per-chromosome remapping BAM/FASTQ
on success unless --save-intermediate or --full-output is supplied. Preserve
assembled transcriptome FASTA/GTF and assignment diagnostics. Fail explicitly
on missing tools or external process failures; never silently fall back.
Saved intermediates have a unique per-run directory under `intermediates/`.
Generated genome indexes live in the per-run temporary directory, so reruns
cannot accidentally use an earlier reference. Sparse matrix thresholds use
simultaneous cell/feature masks, then mitochondrial filtering, normalization
to 10,000 and log10(1+x). Empty processed results are written as empty matrices.

## Implementation checklist

- [x] Add ONT candidate selection tests, including MAPQ independence,
  coverage boundaries, ambiguity, forbidden classes and novel transcripts.
- [x] Implement per-chromosome assembly/remapping and BAM tagging.
- [x] Integrate both single/mixed runners and embedded Python assets.
- [x] Export sparse raw/processed matrices with explicit preprocessing.
- [x] Update dependencies, help, release documentation and comparison commands.
- [x] Run tests, external-tool smoke test when tools are available, and review.

## Verification record (2026-10-22)

All Rust targets passed. 68 focused Python/contract/release/oracle checks
passed. A real-tool synthetic test passed assembly, annotation, transcript
remapping, gene/transcript assignment, Rust UMI correction and TSV matrix
generation, including saved-intermediate reruns. Local smoke tool versions:
StringTie 2.2.2, gffread 0.12.7, gffcompare 0.12.6, minimap2 2.31-r1302;
samtools 1.24 was invoked through pysam's command dispatcher on macOS.
Release compilation, embedded single/mixed help, shell syntax and formatting
checks passed. Code review fixes cover per-run reference indexes, unique
intermediate directories, subprocess cleanup and prefixed chromosome region
labels. Large Linux datasets, full HTML/report generation and a direct ONT
Nextflow numerical comparison have not been tested in this change.

Changes are left uncommitted for the repository owner to review and publish.

Sources: [release](https://github.com/epi2me-labs/wf-single-cell/releases/tag/v3.3.4),
[assembly/mapping](https://github.com/epi2me-labs/wf-single-cell/blob/v3.3.4/subworkflows/assign_features.nf),
[assignment](https://github.com/epi2me-labs/wf-single-cell/blob/v3.3.4/bin/workflow_glue/assign_features.py).
