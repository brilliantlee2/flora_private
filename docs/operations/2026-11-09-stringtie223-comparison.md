# Operation record: isolated StringTie 2.2.3 comparison

Date: 2026-11-09, Asia/Shanghai.
Repository: `/Users/liyy/Documents/GitHub/flora_private`.
Branch: `codex/1109-stringtie-2.2.3`.
Base: local `main`, `291de65fa6b38a3abadd012095076124bc8d9e25`.

## Boundaries

The user requested a new branch to compare StringTie 2.2.3 with the current
main workflow. An external linked worktree isolates these changes; the primary
checkout remains clean on main. No merge, remote push or existing Conda
environment update was performed.

Only the StringTie environment pin, both workflow version guards and the Python
backend guard are changed. Parameters TSV and tool-version JSON identify the
required assembly version while retaining honest historical feature-selection
provenance. Bilingual READMEs and the comparison guide explain separate Conda
environments, clones, executables, output directories and timing controls.

Default assembly options, mapper invocations, feature ranking, UMI correction,
matrix/QC/report algorithms and Cargo version/lockfile are unchanged. The
separate UMI hotspot optimization branch is not merged. StringTie itself may
change assemblies and consequently downstream counts; parity of biological
outputs or a real-data speedup is not claimed.

## Verification

- Version regression tests were first run against the original 2.2.2 checks
  and failed as expected, then passed after the implementation update.
- Rust baseline and final: 108 passed across 16 test targets, locked/offline.
- Fresh CPython 3.11 main baseline: 153 total, 150 passed, 3 skipped.
- Branch Python suite: 156 total, 153 passed, 3 skipped.
- Focused ONT performance/BAM tests: 14 passed, including rejection of an old
  assembler before assembly and actual/required-version provenance.
- Skips cover the external assembly toolchain, original ONT source oracle and
  Scanpy-dependent integration. Samtools fixture inspection uses pysam's native
  command dispatcher in a private test-only wrapper, not an independently
  installed samtools executable.
- Release build succeeded on macOS with CPython 3.11 embedded assets. Both the
  single- and mixed-species executable help show the 2.2.3 requirement.
- Both workflow shells and both packaging shells pass separate `bash -n` checks;
  `git diff --check` passed.
- AST comparison verifies all 23 feature functions/classes outside `run()` are
  unchanged. Each environment file differs only in its StringTie pin. Rust,
  Cargo files, runtime manifest and historical source-hash manifest are unchanged.

The old temporary Python cache lacked source modules; a new private environment
under `/private/tmp` was used instead. No production dependency was changed to
work around that test-host issue. Test logs and build artifacts remain in
`/private/tmp/flora-stringtie223-*` and `/private/tmp/flora-main-audit-target`.

Not performed: Linux release packaging/ABI checks, SIF build, real StringTie
2.2.3 assembly smoke test, or the user's GRCh38 FASTQ/BAM performance comparison.
Those require the Linux toolchain and sample inputs on the server.
