# Operation record: isolated ONT hotspot fixes

Time: 2026-11-06T20:46:31+08:00 (Asia/Shanghai).
Reason: investigate human/mouse runtime differences without changing the tested
`codex/1028-ont-v330-performance` branch or biological selection rules.

## Branch and changes

- Created `codex/1106-ont-hotspot-fix` from
  `9025ff66c219b19696b3e64e3759bb881c26eb68` in a linked worktree below the
  project's `Flora_private/.worktrees/` directory. The GitHub main worktree stays
  on `codex/1028-ont-v330-performance` at that same commit, with no local changes.
- `src/umi_cluster.rs`: exact active-trie traversal for large groups instead of
  materialized candidate pairs/adjacency; retain small exhaustive clustering.
- `src/bin/cluster_umis_allbam.rs`: bounded top-group and large-group diagnostics,
  borrowing existing cell labels rather than duplicating them for every group.
- `scripts/ont_features.py`: bounded completion-driven worker refill and
  StringTie chromosome/options/work-directory start messages.
- Added exhaustive, scalability, scheduler and decoded BAM/matrix parity tests;
  linked focused upstream-evidence documentation from both READMEs.

No dependency/version, default threshold, read cap, filtering, gene/transcript
assignment, output schema or report metric is changed. No automatic StringTie
`-j 3`, upgrade, retry with altered parameters or immunoglobulin filtering.

## Verification

- TDD: missing trie/top-group functionality and blocked worker replenishment
  were observed as failing tests before the production implementation.
- `cargo test --locked --offline --all-targets`: 115 passed, 0 failed;
  1 synthetic timing test ignored in the normal suite and run separately.
- `python -m unittest discover -s tests -v`: 161 tests, 158 passed, 3 skipped.
  Skips: external ONT assembly toolchain, original ONT source oracle, and Scanpy.
- The full Python suite initially failed because system samtools was absent.
  Reran successfully using an existing pysam-embedded samtools dispatcher on
  PATH for frozen BAM fixture comparisons; this is not Linux samtools 1.24.
- New BAM/matrix integration passed against exhaustive Python UMI corrections
  with debug and release binaries, including one/four-thread record parity.
- `cargo build --release --locked --offline --bin flora --bin cluster_umis_allbam
  --bin cell_umi_gene_table` passed; executable reports `flora 0.1.2`.
- Shell syntax, changed-Rust formatting and `git diff --check` passed.
- Independent specification and code-quality reviews found no blocking issue.
- Optimized shared-prefix benchmark matched four original-indexed UMI maps;
  timings and limitations are in `docs/ont-hotspot-fix-20261106.md`.
- Builds used Python 3.11 via `PYTHON_BIN`; the host default is Python 3.9.
  Existing dependency warning: `buf_redux 0.8.4` future Rust incompatibility.

Real server FASTQ/BAM end-to-end performance, peak RSS and actual StringTie
chr2/chr22 bundle causes were not measured locally. The original production
branch remains unchanged. Record these changes only on the new branch in Git;
keep the branch/worktree for user review. Do not push, merge or create a release
automatically.
