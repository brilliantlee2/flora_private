# ONT alignment toolchain comparison

This branch provides a controlled dependency profile for comparing Flora with
the alignment toolchain used by the `wf-single-cell` 3.3 generation. It changes
software versions only. Flora's alignment command, filtering, barcode/UMI
logic, gene assignment, and transcript assignment are unchanged.

## Dependency profile

| Component | Flora main before this branch | This comparison branch | Role |
| --- | --- | --- | --- |
| minimap2 | `>=2.31,<2.32` | `>=2.28,<2.29` | Produces alignments and MAPQ |
| samtools | `>=1.24,<1.25` | `>=1.21,<1.22` | Converts, filters, sorts, and indexes BAM |
| bedtools | `>=2.31,<2.32` | `>=2.31,<2.32` | BED projection used by transcript assignment |
| StringTie | Not used | Not added | ONT sample-transcriptome assembly |
| gffcompare | Not used | Not added | ONT transcript annotation comparison |

`wf-single-cell` pins a complete container image rather than declaring each
binary version in the workflow repository. Its 3.3.1 changelog explicitly
records StringTie 2.2.3. The minimap2 2.28 and samtools 1.21 pins in this branch
form a reproducible ONT-generation comparison profile; they must not be
described as a byte-for-byte reconstruction of the ONT container.

## Why minimap2 is the important variable

Minimap2 2.29 changed MAPQ calculation for spliced alignments. Flora maps RNA
reads with the `splice` preset, so a 2.28 versus 2.31 comparison can legitimately
move reads across the MAPQ threshold even when every Flora argument remains the
same.

Samtools does not calculate minimap2 MAPQ. In Flora it filters, sorts, and
indexes alignment records. Bedtools also does not calculate MAPQ. Their versions
are pinned here to keep the comparison environment reproducible, not because
they are expected to rescue low-MAPQ reads.

StringTie and gffcompare are intentionally absent. Adding the programs without
rewiring Flora's gene/transcript algorithms would not change any result and
would make the experiment harder to interpret.

## Build and verify

Create a fresh environment. Do not reuse an environment that already contains
the main-branch toolchain:

```bash
conda env remove -n flora-ont-toolchain
conda env create -n flora-ont-toolchain -f environment.yml
conda activate flora-ont-toolchain

minimap2 --version
samtools --version | head -n 2
bedtools --version
```

Expected version families are:

```text
minimap2 2.28
samtools 1.21
bedtools 2.31.x
```

Build the release from this same activated environment:

```bash
PYTHON_BIN="$CONDA_PREFIX/bin/python" \
  bash packaging/build_binary_release.sh
```

## Comparison requirements

For an interpretable result, run main and this branch with the same FASTQ,
reference FASTA, GTF, barcode list, Flora arguments, thread counts, and Conda
channels. Use different output directories. Compare at least:

- the exact MAPQ distribution;
- counts of `Assigned`, `Unassigned_mapq`, `Unassigned_no_features`, and
  `Unassigned_ambiguous`;
- read IDs moving from MAPQ below 30 to MAPQ 30 or higher, and vice versa;
- gene and transcript expression matrices.

A version-only comparison does not reproduce ONT preprocessing or ONT's
StringTie/gffcompare assignment model. It answers the narrower question of
whether the older minimap2 generation changes Flora's alignment and MAPQ loss.
