# Flora Single-Binary Workflow Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `flora` run the complete single-species workflow, make `flora mixed` run the complete mixed-species workflow, retain `flora glycine` and `flora analyze`, and publish all runtime code in one stripped `flora` executable without changing scientific results.

**Architecture:** Replace both Bash orchestrators with typed Rust workflow modules and a shared stage engine. Invoke Glycine and existing Rust stages through library APIs, supervise external tools through an argument-safe process runner, and embed deterministic CPython 3.11 bytecode plus report assets into the executable for verified temporary extraction at runtime.

**Tech Stack:** Rust 2021, Clap 4, `anyhow`, `serde`, `sha2`, `zip`, Unix process groups/signals, CPython 3.11, minimap2, samtools, bedtools, Bash release tooling, Rust tests, and Python `unittest`.

**Design spec:** `docs/superpowers/specs/2026-08-07-flora-single-binary-workflow-design.md`

---

## File Structure

```text
src/cli.rs                    four public command contracts and dispatch
src/workflow/config.rs        normalized single/mixed configuration
src/workflow/paths.rs         references, outputs, and stage paths
src/workflow/command.rs       supervised commands and pipelines
src/workflow/logging.rs       stage logs and progress records
src/workflow/runtime.rs       embedded Python asset extraction/validation
src/workflow/stages.rs        shared ordered workflow stages
src/workflow/species.rs       single/mixed behavior hooks
src/workflow/mod.rs           full-workflow entrypoints
src/stages/*.rs               reusable Rust implementations from src/bin
build.rs                      deterministic embedded-runtime bundle build
packaging/compile_embedded_runtime.py
tests/fixtures/cli/*.json
tests/fixtures/workflows/*
tests/test_workflow_contract.py
tests/test_stage_parity.py
tests/test_single_binary_release.py
```

### Task 1: Freeze All Four CLI Contracts

**Files:**
- Create: `tests/fixtures/cli/single.json`
- Create: `tests/fixtures/cli/mixed.json`
- Create: `tests/fixtures/cli/glycine.json`
- Create: `tests/fixtures/cli/analyze.json`
- Create: `tests/test_workflow_contract.py`
- Reference: `run_all.sh`
- Reference: `run_all_mixed_species.sh`
- Reference: `src/main.rs`
- Reference: `src/glycine/args.rs`

- [ ] Record every option, alias, type, default, conflict, conditional requirement, and accepted/rejected example for all four commands.
- [ ] Record `--glycine-bin-dir` as rejected by `flora` and `flora mixed`; do not silently ignore it.
- [ ] Record `--python` plus `FLORA_PYTHON` precedence and `--singlet-threshold` mixed-species default, all float spellings accepted by legacy Python parsing, rejected non-float examples, and barnyard-only propagation.
- [ ] Add tests asserting `flora glycine` has no reference/alignment/Python requirements and `flora analyze` has no Glycine/alignment/Python requirements.
- [ ] Run `python3 -m unittest tests.test_workflow_contract -v`; expect failures because the new CLI does not exist yet.
- [ ] Commit with `git commit -m "test: freeze Flora workflow CLI contracts"`.

### Task 2: Capture Legacy Workflow Oracles

**Files:**
- Create: `tests/fixtures/workflows/single/`
- Create: `tests/fixtures/workflows/mixed/`
- Create: `tests/fixtures/workflows/tool_versions.tsv`
- Create: `tests/artifact_manifest_single.json`
- Create: `tests/artifact_manifest_mixed.json`
- Create: `tests/compare_workflow_outputs.py`
- Create: `tests/test_workflow_parity.py`
- Create: `tests/fixtures/workflows/generate_legacy_oracles.sh`

- [ ] Build deterministic small FASTQ, 10 bp whitelist, genome FASTA, GTF, junction BED, and chromosome-size fixtures for both modes.
- [ ] Record exact Python, package, minimap2, samtools, bedtools, libc, and OS versions and refuse oracle regeneration under a different environment.
- [ ] Make `generate_legacy_oracles.sh` run `cargo build --release --locked --bins`, fail if any expected Rust stage binary is absent, run both legacy workflows, and assert their logs selected every available release Rust stage rather than a Python fallback.
- [ ] Run unchanged `run_all.sh` and `run_all_mixed_species.sh` and classify every output as byte-exact, parsed-exact, numeric-tolerance, canonicalized HTML/log, or intentionally absent for base, mixed, skip-Glycine, skip-isoform, light/full, upstream-only, stale/malformed, and forced-failure cases.
- [ ] Include barcode statistics, cell assignments, BAM tags/index validity, gene/isoform matrices, Scanpy tables, QC, saturation, barnyard output, reports, logs, and symlinks.
- [ ] Implement manifest tests that reject every unclassified missing or extra path and verify all declared comparison rules can parse their oracle artifacts.
- [ ] Freeze oracle artifacts and manifests after generation. Changes require an explicit regeneration command under the exact pinned environment and a dedicated review; new-workflow failures must never update the oracle.
- [ ] Run `python3 -m unittest tests.test_workflow_contract tests.test_workflow_parity -v`; expect PASS for the frozen legacy oracle.
- [ ] Commit with `git commit -m "test: capture legacy Flora workflow oracles"`.

### Task 3: Introduce the New Top-Level CLI

**Files:**
- Create: `src/cli.rs`
- Modify: `src/main.rs`
- Modify: `src/lib.rs`
- Test: `tests/test_workflow_contract.py`

- [ ] Add failing tests for `flora --help`, `flora mixed --help`, `flora glycine --help`, and `flora analyze --help`.
- [ ] Define `SingleArgs`, `MixedArgs`, and explicit dispatch for `mixed`, `glycine`, and `analyze`; all other option-led invocations parse as the single workflow.
- [ ] Parse `--python` for both full workflows and `--singlet-threshold` only for mixed mode, including default/custom/invalid tests from the fixtures.
- [ ] Remove the current fallback that lets the word `run` become an analyze FASTQ positional argument.
- [ ] Make `flora --help` show full single-workflow options and the three explicit subcommands.
- [ ] Keep `flora glycine` and `flora analyze` argument behavior unchanged and ensure `flora run` fails rather than aliasing another command.
- [ ] Run `cargo test` and `python3 -m unittest tests.test_workflow_contract -v`; expect PASS for CLI tests.
- [ ] Commit with `git commit -m "feat: define Flora workflow CLI"`.

### Task 4: Normalize Workflow Configuration and Paths

**Files:**
- Create: `src/workflow/mod.rs`
- Create: `src/workflow/config.rs`
- Create: `src/workflow/paths.rs`
- Create: `src/workflow/species.rs`
- Modify: `src/lib.rs`

- [ ] Write failing Rust tests for defaults (`threads`, `cluster-threads`, `max-ed=2`, `top1-alpha=0.1`), aliases, output modes, and skip conditions.
- [ ] Preserve legacy `--singlet-threshold` float parsing exactly: accept finite, nonfinite, and out-of-range floats accepted by Python `float`, and reject only strings legacy parsing rejects.
- [ ] Define `SpeciesMode`, `WorkflowArgs`, `WorkflowConfig`, `ReferencePaths`, `OutputPaths`, and `StagePlan` with no process-global state.
- [ ] Implement `--ref-dir` inference for `genome.fa`, `genes.gtf`, `genes.bed`, and `chrom_sizes.tsv`, with explicit arguments taking precedence.
- [ ] Implement normalization and pure input/path validation only; defer stage-derived executable/package preflight until the complete stage plan and runtime exist.
- [ ] Normalize interpreter selection as `--python` then `FLORA_PYTHON` then `python3`, without probing the interpreter in `flora glycine` or `flora analyze`.
- [ ] Validate `--skip-glycine` requires `--full-length-fastq`; otherwise require raw `--fastq`, TSO, and RTP inputs exactly as the legacy contract specifies.
- [ ] Run `cargo test workflow::config` and `cargo test workflow::paths` separately; expect PASS.
- [ ] Commit with `git commit -m "feat: add typed workflow configuration"`.

### Task 5: Build Stage Logging and External Process Supervision

**Files:**
- Create: `src/workflow/logging.rs`
- Create: `src/workflow/command.rs`
- Modify: `Cargo.toml`
- Modify: `Cargo.lock`

- [ ] Add `signal-hook` and `nix` (process/signal features) dependencies and failing tests for paths with spaces, argument preservation, stdout/stderr teeing, append/truncate modes, and large stderr.
- [ ] Add explicit failures for each member of `minimap2 -> samtools view -> samtools sort`, early downstream exit with upstream SIGPIPE, deterministic primary-error selection, hung children, cancellation, and no surviving descendants.
- [ ] Add Unix process-group tests for SIGINT/SIGTERM forwarding, TERM-to-KILL timeout, child reaping, and no surviving descendants.
- [ ] Implement `CommandSpec`, `StageLog`, and `ProcessSupervisor` using argument vectors rather than shell strings.
- [ ] Implement `minimap2 -> samtools view -> samtools sort` with all three exit statuses checked and deterministic primary-error reporting.
- [ ] Ensure unused pipe handles close immediately and stderr is drained concurrently with bounded memory.
- [ ] Run `cargo test workflow::command` and `cargo test workflow::logging` separately; expect PASS.
- [ ] Commit with `git commit -m "feat: supervise Flora workflow processes"`.

### Task 6: Extract Whitelist, Tagging, and Assignment Stage APIs

**Files:**
- Create: `src/stages/mod.rs`
- Create: `src/stages/generate_whitelists.rs`
- Create: `src/stages/prepare_read_tags.rs`
- Create: `src/stages/add_cb_ur_tags.rs`
- Create: `src/stages/assign_genes.rs`
- Create: `src/stages/add_gene_tags.rs`
- Modify: `src/bin/generate_26bp_whitelists.rs`
- Modify: `src/bin/prepare_read_tags.rs`
- Modify: `src/bin/add_cb_ur_tags.rs`
- Modify: `src/bin/assign_genes.rs`
- Modify: `src/bin/add_gene_tags.rs`
- Modify: `src/lib.rs`
- Create: `tests/test_stage_parity.py`

- [ ] Add wrapper-versus-library fixture tests for each listed stage before moving code.
- [ ] Preserve the exact binary-to-module mapping: `generate_26bp_whitelists` to `generate_whitelists`, `prepare_read_tags` to `prepare_read_tags`, `add_cb_ur_tags` to `add_cb_ur_tags`, `assign_genes` to `assign_genes`, and `add_gene_tags` to `add_gene_tags`.
- [ ] Move computation into typed `run(config, writers) -> Result<()>` functions; leave each private binary as a thin Clap wrapper.
- [ ] Prohibit `process::exit`, working-directory/environment mutation, and unmanaged global stdout/stderr inside library stages.
- [ ] Compare files and normalized diagnostics against the legacy oracle after each extraction.
- [ ] Run `cargo test` and `python3 -m unittest tests.test_stage_parity -v`; expect PASS.
- [ ] Commit with `git commit -m "refactor: expose tagging and assignment stage APIs"`.

### Task 7: Extract Expression, UMI, and QC Stage APIs

**Files:**
- Create: `src/stages/cluster_umis.rs`
- Create: `src/stages/cell_umi_gene.rs`
- Create: `src/stages/gene_expression.rs`
- Create: `src/stages/assign_transcripts.rs`
- Create: `src/stages/isoform_expression.rs`
- Create: `src/stages/rna_qc.rs`
- Create: `src/stages/read_qc.rs`
- Modify: `src/bin/cluster_umis_allbam.rs`
- Modify: `src/bin/cell_umi_gene_table.rs`
- Modify: `src/bin/gene_expression.rs`
- Modify: `src/bin/assign_transcripts.rs`
- Modify: `src/bin/isoform_expression.rs`
- Modify: `src/bin/rna_qc_metrics.rs`
- Modify: `src/bin/read_qc_summary.rs`
- Modify: `src/lib.rs`
- Modify: `Cargo.toml`
- Modify: `tests/test_stage_parity.py`

- [ ] Add failing wrapper-versus-library parity tests for each remaining Rust stage.
- [ ] Preserve the exact mappings for `cluster_umis_allbam`, `cell_umi_gene_table`, `gene_expression`, `assign_transcripts`, `isoform_expression`, `rna_qc_metrics`, and `read_qc_summary` to their named `src/stages` modules.
- [ ] Move implementations into typed library APIs while preserving streaming, row order, thresholds, tag semantics, and memory behavior.
- [ ] Use stage-local thread pools/resources and require release `panic = "unwind"` for the workflow panic boundary.
- [ ] Add a Cargo/profile contract test that fails if `[profile.release] panic = "unwind"` is removed or changed to abort.
- [ ] Compare all matrices, assignments, QC tables, BAM tags, and diagnostics against the oracle.
- [ ] Run `cargo test` and `python3 -m unittest tests.test_stage_parity -v`; expect PASS.
- [ ] Commit with `git commit -m "refactor: expose expression and QC stage APIs"`.

### Task 8: Build the Embedded Python Runtime

**Files:**
- Create: `build.rs`
- Create: `packaging/compile_embedded_runtime.py`
- Create: `src/workflow/runtime.rs`
- Modify: `Cargo.toml`
- Modify: `Cargo.lock`
- Modify: `src/workflow/mod.rs`
- Create: `tests/test_embedded_runtime.py`

- [ ] Add failing tests for deterministic CPython 3.11 bytecode names, manifest schema, hashes, sizes, compression, ABI, dependency constraints, and required entrypoints.
- [ ] Generate one deterministic bundle containing `generate_knee_plots.py`, `rna_cluster_analysis.py`, `rna_violin_plot.py`, `Saturation.py`, `barnyard_qc.py`, and `build_report.py` bytecode plus `report_template.html` and `plotly-2.26.0.min.js`; embed it with `include_bytes!`. Do not embed the unused `rna_qc_metrics_mixed.py` fallback.
- [ ] Define manifest resource records with logical path, type, entrypoint, SHA-256, compression, compressed size, and uncompressed size.
- [ ] Atomically create a random mode-`0700` runtime directory; reject absolute, traversal, duplicate, symlink, pre-existing, missing-entrypoint, unsupported-compression, oversized, wrong-size, and wrong-hash resources.
- [ ] Execute Python with `-I`, clear `PYTHONHOME/PYTHONPATH`, disable user site, and expose only the extracted entrypoint while preserving declared Conda site-packages.
- [ ] Remove only the current process-owned runtime directory on normal/error exit; do not perform unsafe stale-directory cleanup.
- [ ] Test hostile shadow modules, concurrent runs, interrupted extraction, cleanup failure, corrupt bundles, ABI/package mismatch, package-relative imports, and `build_report.py` resolving its template/Plotly through `__file__`.
- [ ] Run `cargo test workflow::runtime` and `python3 -m unittest tests.test_embedded_runtime -v`; expect PASS.
- [ ] Commit with `git commit -m "feat: embed Flora Python runtime assets"`.

### Task 9: Port Glycine, Whitelist, Alignment, and Analyze Stages

**Files:**
- Create: `src/workflow/stages.rs`
- Create: `src/workflow/provenance.rs`
- Modify: `src/workflow/mod.rs`
- Create: `tests/test_workflow_stage_plan.py`

- [ ] Add failing stage-plan snapshots for optional embedded Glycine, whitelist generation, FASTQ selection, alignment pipeline, BAM index, analyze, knee plotting, logs, provenance, interpreter selection, and mixed singlet threshold propagation.
- [ ] Resolve the complete stage plan and derive executable, Python ABI/package, input, and reference preflight from that plan; test `--python > FLORA_PYTHON > python3`, missing interpreter, ABI mismatch, and no Python lookup for `flora glycine`/`flora analyze`.
- [ ] Invoke Glycine directly through `flora::glycine`; never resolve a binary path or spawn `flora glycine` from the full workflow.
- [ ] Invoke whitelist and analyze logic through Rust library APIs, alignment through `ProcessSupervisor`, and knee plotting through the embedded runtime.
- [ ] Preserve working directories, filenames, parameter TSV fields, light/full/debug/intermediate behavior, and `--upstream-only` outputs.
- [ ] Compare Glycine statistics, barcode validity, pair threshold, core cells/barcodes, assignments, knee outputs, and BAM validity to the oracle.
- [ ] Run `cargo test workflow` and `python3 -m unittest tests.test_workflow_stage_plan -v`; expect PASS.
- [ ] Commit with `git commit -m "feat: port Flora upstream workflow stages"`.

### Task 10: Port BAM, Gene, Isoform, and Matrix Stages

**Files:**
- Modify: `src/workflow/stages.rs`
- Modify: `src/workflow/mod.rs`
- Test: `tests/test_workflow_stage_plan.py`

- [ ] Add failing snapshots for read-tag preparation, CB/UR and gene tags, BED conversion, UMI clustering, cell table, gene matrix, transcript assignment, and isoform matrix.
- [ ] Call Rust stages through their library APIs; call samtools/bedtools through `ProcessSupervisor`.
- [ ] Preserve exact BAM inputs, MAPQ thresholds, file moves/symlinks, index points, optional isoform behavior, log modes, and output names.
- [ ] Compare BAM records/tags, barcode-to-cell CSV, cell/UMI/gene table, and gene/isoform matrices to the oracle.
- [ ] Run `cargo test workflow` and `python3 -m unittest tests.test_workflow_stage_plan -v`; expect PASS.
- [ ] Commit with `git commit -m "feat: port Flora expression workflow stages"`.

### Task 11: Port Python QC, Clustering, Saturation, and Reporting

**Files:**
- Modify: `src/workflow/stages.rs`
- Modify: `src/workflow/mod.rs`
- Modify: `tests/test_workflow_stage_plan.py`
- Test: `tests/test_embedded_runtime.py`

- [ ] Add failing snapshots for knee plots, RNA clustering, violin fallback, saturation, mixed RNA QC, barnyard QC, read QC, and report arguments.
- [ ] Invoke every Python-only stage through the embedded runtime, never an installed script path.
- [ ] Pass the legacy-parsed `--singlet-threshold` exactly to `barnyard_qc.py` only and record it in the parameter TSV; do not pass it to unrelated mixed-QC stages.
- [ ] Run mixed RNA QC through the extracted Rust `rna_qc_metrics` API with mixed mode enabled and compare it to an oracle whose logs prove the Rust stage binary was selected.
- [ ] Preserve Glycine/raw/full-length read accounting, `--skip-isoform`, optional fallback images, report payload, and final report path.
- [ ] Verify report generation succeeds when the source tree and external template/Plotly files are absent.
- [ ] Compare tables, JSON, PNG dimensions/pixel hashes, saturation payload, barnyard results, and canonicalized HTML to the oracle.
- [ ] Run all workflow/runtime tests; expect PASS.
- [ ] Commit with `git commit -m "feat: port Flora Python workflow stages"`.

### Task 12: Complete Shared Single/Mixed Workflow Execution

**Files:**
- Modify: `src/workflow/species.rs`
- Modify: `src/workflow/stages.rs`
- Modify: `src/workflow/mod.rs`
- Modify: `src/main.rs`
- Test: `tests/test_workflow_stage_plan.py`

- [ ] Add failing end-to-end stage-plan tests for single, mixed, skip-Glycine, skip-isoform, upstream-only, light-output, full-output, and debug/intermediate combinations.
- [ ] Implement one shared stage sequence with explicit single/mixed hooks; do not duplicate the two legacy scripts in Rust.
- [ ] Add stage panic boundaries, stage-labelled errors, log paths, child cleanup, and nonzero exits.
- [ ] Add a release-only injected-panic fixture that verifies stage-labelled stderr, stage log creation, nonzero exit, child reaping, and embedded-runtime cleanup.
- [ ] Ensure `flora glycine` and `flora analyze` produce no downstream artifacts and do not preflight downstream dependencies.
- [ ] Run the new development binary against both fixed fixtures and compare every artifact to the legacy oracle.
- [ ] Run `cargo test` and all Python tests; expect PASS.
- [ ] Commit with `git commit -m "feat: complete Flora single and mixed workflows"`.

### Task 13: Replace the Linux Release Packager

**Files:**
- Modify: `packaging/build_binary_release.sh`
- Modify: `packaging/refresh_binary_release_metadata.sh`
- Modify: `tests/test_release_layout.py`
- Create: `tests/test_single_binary_release.py`
- Create: `packaging/containers/glibc217/Dockerfile`
- Create: `packaging/containers/glibc217/run_release_fixtures.sh`
- Create: `docs/release-size-approvals.md`
- Modify: `.gitignore`

- [ ] Add failing tests requiring the exact archive allowlist and `flora` as the only regular file with an executable bit.
- [ ] Build only `cargo build --release --locked --target x86_64-unknown-linux-gnu --bin flora` with baseline x86-64 and path remapping.
- [ ] Stage only `flora`, environment/requirements files, public READMEs, notices, and enumerated licenses; do not stage scripts, bytecode, templates, JavaScript, Cargo files, or helper binaries.
- [ ] Change the release gate to require `flora --help`, `flora mixed --help`, `flora glycine --help`, and `flora analyze --help`; remove obsolete `run/run-mixed` checks.
- [ ] Reject symlinks, traversal, unlisted entries, extra executable regular files, private paths, debug symbols, xattrs, and GLIBC symbols newer than 2.17.
- [ ] Scan archive names/modes, ELF sections/strings, decoded embedded manifest, Python code-object filenames, extracted HTML/JavaScript, and packaged documentation for private absolute paths and legacy source-tree names.
- [ ] Add a digest-pinned `quay.io/pypa/manylinux2014_x86_64`-based test image, record the digest in the Dockerfile, inspect ISA requirements, and run `docker build -t flora-glibc217-test -f packaging/containers/glibc217/Dockerfile .`.
- [ ] Generate and verify SHA256, then run `docker run --rm -v "$PWD/dist:/dist:ro" flora-glibc217-test /opt/flora/run_release_fixtures.sh /dist/Flora-<version>-linux-x86_64.tar.gz` to execute functional fixture runs of all four commands, including complete single/mixed reports, without source, legacy scripts, helper binaries, external Glycine, or external report assets.
- [ ] Run `python3 -m unittest tests.test_release_layout tests.test_single_binary_release -v`; expect PASS.
- [ ] Commit with `git commit -m "build: package Flora as one runtime executable"`.

### Task 14: Update Public and Private Documentation

**Files:**
- Modify: `README.md`
- Modify: `README_zh-CN.md`
- Modify: `docs/repository-templates/public/README.md`
- Modify: `docs/repository-templates/public/README_zh-CN.md`
- Modify: `docs/repository-templates/private/README.md`
- Modify: `docs/repository-templates/private/README_zh-CN.md`

- [ ] Replace every `run_all*.sh`, `flora run`, and `flora run-mixed` user example with `flora` or `flora mixed` as appropriate.
- [ ] Document that Glycine is internal and `--glycine-bin-dir` no longer exists; retain `flora glycine` and `flora analyze` advanced examples.
- [ ] Document external tools, CPython 3.11 packages, reference files, HPC examples, checksum verification, supported glibc/CPU baseline, and output/resource recommendations.
- [ ] Add documentation tests that reject stale command names and source-tree installation assumptions.
- [ ] Commit with `git commit -m "docs: document the single-binary Flora CLI"`.

### Task 15: Final Scientific and Release Verification

**Files:**
- Test (frozen): `tests/artifact_manifest_single.json`
- Test (frozen): `tests/artifact_manifest_mixed.json`
- Modify only for comparator bugs, never new-output mismatches: `tests/compare_workflow_outputs.py`
- Modify only after explicit approval: `docs/release-size-approvals.md`

- [ ] Run clean legacy/new single workflows and compare every classified artifact.
- [ ] Run clean legacy/new mixed workflows and compare every classified artifact.
- [ ] Repeat skip/light/full/upstream combinations and interrupt/restart both workflows at every stage with stale zero-byte, truncated, malformed, and mismatched-parameter files; verify none is treated as a completion marker.
- [ ] Run `cargo test --all-targets` and `python3 -m unittest discover -s tests -v`; expect all PASS.
- [ ] Build on Linux, verify the exact allowlist, sole executable, no private/source paths, baseline GLIBC/ISA, no xattrs, valid SHA256, and successful extracted-package workflows.
- [ ] Record final archive size; fail above 28 MB unless `docs/release-size-approvals.md` contains a reviewed entry identifying the measured dependency/resource change and approved new ceiling.
- [ ] Commit with `git commit -m "test: verify Flora single-binary release parity"`.
