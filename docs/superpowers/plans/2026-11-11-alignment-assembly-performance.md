# Alignment and Assembly Performance Plan

**Goal:** Address the measured small-batch minimap2 pipeline and suspected StringTie same-start input lookup cost without filtering reads or changing assembly thresholds.

**Baseline:** Local main 291de65. Isolated branch codex/1111-alignment-assembly-performance. Python baseline: 153 tests, 3 skipped with pysam's test-only samtools dispatcher.

**Design:** Keep splice preset, k15/w5, junction BED, MAPQ and secondary/supplementary rules. Expose batch, I/O and gzip supplier controls with reversible defaults. Add explicit StringTie executable selection without changing the active Conda installation. An optional source-pinned StringTie 2.2.2 lookup accelerator preserves comparisons/model thresholds and remains opt-in until real-data parity is established. Do not split loci or raise -j/-c by default.

- [x] Inspect main, prior sample benchmark, official issue/source evidence; establish clean baseline.
- [x] Add failing executable pipeline tests for plain/gzip input, failure propagation, controls and CLI contracts.
- [x] Implement in both single/mixed workflows and normalize the explicit StringTie path for embedded runtime.
- [x] Build/review optional upstream lookup accelerator and compare stock/patched synthetic assemblies.
- [x] Record official sources, changes that affect sensitivity, limitations, and validation commands for 260F800579011/sample_dir.
- [x] Run focused/full tests, Rust tests and shell checks; leave main unchanged and report remaining server validation.

Validation: run serial controls on identical cDNA/chr2/chr22 inputs, compare decoded alignment records and assembly/assignment/matrix output, measure tool wall/CPU/RSS. Batch sizes and profiler/build variants must not be credited with speedups from concurrent load differences. No server FASTQ/BAM is accessible on this local Mac; full-run claims require the user's Linux validation.

Local verification completed: 170 Python tests run (3 explicitly skipped for unavailable whole-toolchain/workflow parity setups), 109 Rust tests passed across all targets, both shell runners passed independent syntax checks, native minimap2 2.31 pipeline parity passed. Both source-pinned StringTie release builds and actual processRead/assembly/abundance/merge differential tests ran on macOS. Production path review found no remaining introduced regressions; the selected executable preflight and standalone ./stringtie resolution issues were fixed with regression tests. Real Linux/Conda build and 260F800579011 validation remain pending, not marked complete.
