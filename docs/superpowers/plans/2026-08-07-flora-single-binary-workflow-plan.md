# Flora Single-Binary Workflow Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace public shell workflow entrypoints and standalone stage executables with `flora run` and `flora run-mixed` in one stripped release executable while preserving scientific outputs.

**Architecture:** Port orchestration into typed Rust workflow modules. Extract each existing Rust stage into a library `run(config, writers) -> Result<()>` API called directly by the workflow; private standalone binaries remain thin wrappers for regression testing. Keep Python-only stages as deterministic CPython 3.11 bytecode and package from an explicit allowlist.

**Tech Stack:** Rust 2021, Clap 4, `std::process::Command`, Unix process groups/signals, CPython 3.11 bytecode, Bash packaging, Rust tests and Python `unittest`.

**Design spec:** `docs/superpowers/specs/2026-08-07-flora-single-binary-workflow-design.md`

---

### Task 1: Capture the Legacy Oracle Before Refactoring

**Files:**
- Create: `tests/fixtures/run_cli_contract.json`
- Create: `tests/fixtures/run_mixed_cli_contract.json`
- Create: `tests/fixtures/workflows/single/input.fastq.gz`
- Create: `tests/fixtures/workflows/single/BC_test.txt`
- Create: `tests/fixtures/workflows/single/ref/{genome.fa,genes.bed,chrom_sizes.tsv,genes.gtf}`
- Create: `tests/fixtures/workflows/mixed/input.fastq.gz`
- Create: `tests/fixtures/workflows/mixed/BC_test.txt`
- Create: `tests/fixtures/workflows/mixed/ref/{genome.fa,genes.bed,chrom_sizes.tsv,genes.gtf}`
- Create: `tests/fixtures/workflows/environment.lock.yml`
- Create: `tests/fixtures/workflows/tool_versions.tsv`
- Create: `tests/artifact_manifest_single.json`
- Create: `tests/artifact_manifest_mixed.json`
- Create: `tests/compare_workflow_outputs.py`
- Create: `tests/test_workflow_contract.py`

- [ ] Record every option, alias, default, conditional requirement, accepted invocation, rejected invocation, stdout/stderr destination, and exit status from both legacy scripts.
- [ ] Create deterministic synthetic single/mixed FASTQ, whitelist, genome,
  junction BED, chromosome sizes, and GTF files at the exact paths listed above.
- [ ] Pin CPython and every Python package in `environment.lock.yml`; record exact
  Flora build, minimap2, samtools, bedtools, libc, and OS versions in
  `tool_versions.tsv`. Refuse oracle regeneration when versions differ.
- [ ] Run the unchanged scripts on those fixed fixtures in the locked environment and check in manifests classifying every output path by the spec's byte/structured/tolerance rule.
- [ ] Capture per-stage command vectors, working directories, log modes, symlink targets, and representative stage diagnostics before changing code.
- [ ] Implement a comparator that rejects unclassified missing or extra paths.
- [ ] Run `python -m unittest tests.test_workflow_contract -v`; expect PASS against the legacy fixtures.

### Task 2: Lock the New CLI and Release Contracts

**Files:**
- Create: `tests/test_single_binary_release.py`
- Modify: `tests/test_release_layout.py`

- [ ] Add failing tests requiring source definitions for `flora run --help` and `flora run-mixed --help` matching the normalized legacy contracts.
- [ ] Add failing tests requiring one root-level `flora`, no `run_all*.sh`, no additional ELF, no source extensions, and all approved assets including optional `BC_1536.txt` when present.
- [ ] Run `python -m unittest tests.test_release_layout tests.test_single_binary_release -v`; expect the new layout assertions to FAIL.

### Task 3: Add Typed Workflow Configuration and Paths

**Files:**
- Create: `src/workflow/mod.rs`
- Create: `src/workflow/config.rs`
- Create: `src/workflow/paths.rs`
- Modify: `src/lib.rs`
- Modify: `src/main.rs`

- [ ] Write failing Rust tests for all defaults, aliases, flags, mode-specific options, and conditional requirements from the contract fixtures.
- [ ] Define `RunArgs`, `SpeciesMode`, `WorkflowConfig`, and normalized output flags with Clap.
- [ ] Implement reference inference/override precedence, output path normalization, bundled `BC_1536.txt`, installation root, and interpreter selection.
- [ ] Resolve the complete stage plan first, then preflight only the external
  tools, input files, references, and assets required by that plan; add tests
  proving skipped stages do not require their tools and required stages fail
  before expensive work.
- [ ] Dispatch `run` and `run-mixed` without changing `glycine`, `analyze`, version, or existing top-level behavior.
- [ ] Run `cargo test workflow::config` and `cargo test workflow::paths`; expect PASS.

### Task 4: Build the External Process and Logging Runtime

**Files:**
- Create: `src/workflow/command.rs`
- Create: `src/workflow/logging.rs`
- Modify: `Cargo.toml`
- Modify: `Cargo.lock`

- [ ] Write failing tests for argument-safe execution, tee logging, append/truncate modes, large stderr, paths with spaces, and nonzero exit propagation.
- [ ] Implement bounded concurrent stream draining and exact stage log destinations.
- [ ] Add Unix process-group creation, signal forwarding, TERM/KILL timeout, and complete reaping.
- [ ] Write failing tests for each position of `minimap2 -> samtools view -> samtools sort`, downstream SIGPIPE, hung children, cancellation, and no survivors.
- [ ] Implement the streaming pipeline and deterministic primary-error selection.
- [ ] Run `cargo test workflow::command`; expect PASS, then run `cargo test`.

### Task 5: Extract Rust Stage Library APIs

**Files:**
- Create: `src/stages/mod.rs`
- Create one focused library module under `src/stages/` for every current `src/bin/*.rs`
- Modify every `src/bin/*.rs` into a thin parse-and-call wrapper
- Modify: `src/lib.rs`
- Create: `tests/test_stage_parity.py`

- [ ] For one stage at a time, add a failing wrapper-versus-library fixture comparison before moving its logic.
- [ ] Move computation into typed `run(config, stdout, stderr) -> Result<()>` APIs without `process::exit`, global CWD/environment mutation, global thread pools, or unmanaged output.
- [ ] Keep private standalone binaries as thin wrappers using the same API.
- [ ] Add workflow panic boundaries and ensure stage-local thread pools/resources drop before the next stage.
- [ ] Repeat red/green parity checks for all twelve Rust stages; compare files and diagnostics to the Task 1 oracle.
- [ ] Run `cargo test` and `python -m unittest tests.test_stage_parity -v`; expect PASS.

### Task 6: Build Deterministic Python Runtime Assets

**Files:**
- Create: `runtime_manifest.json`
- Create: `packaging/compile_python_assets.py`
- Create: `src/workflow/python.rs`
- Modify: `tests/test_single_binary_release.py`

- [ ] Add failing tests for CPython 3.11, manifest schema, import/version constraints, hash-based bytecode, deterministic public `co_filename`, and report asset lookup.
- [ ] Compile each Python-only stage with a public `flora/scripts/<name>.py` code filename.
- [ ] Implement manifest parsing plus `--python`, `FLORA_PYTHON`, and `python3` precedence.
- [ ] Implement stage-plan-dependent interpreter/package preflight in this module only.
- [ ] Extract a test runtime tree and execute every packaged Python `--help`/smoke path.

### Task 7: Port Whitelist, Glycine, and Alignment Stages

**Files:**
- Create: `src/workflow/stages.rs`
- Create: `src/workflow/provenance.rs`
- Modify: `src/workflow/mod.rs`
- Create: `tests/test_workflow_stage_plan.py`

- [ ] Add failing snapshots for whitelist generation, optional Glycine, FASTQ selection, minimap2/view/sort, BAM index, parameters, provenance, working directories, and log modes.
- [ ] Implement these stages with direct Rust library calls for internal stages and the supervised runner for external tools.
- [ ] Compare stage outputs and normalized logs to the legacy oracle before proceeding.
- [ ] Test `--skip-glycine`, explicit full-length FASTQ, inferred references, explicit overrides, and paths containing spaces.

### Task 8: Port Barcode Assignment and Knee Plot Stages

**Files:**
- Modify: `src/workflow/stages.rs`
- Modify: `src/workflow/mod.rs`
- Test: `tests/test_workflow_stage_plan.py`

- [ ] Add failing snapshots for every `analyze` option, optional pair minimum, whitelist orientation, debug/intermediate/light outputs, and knee plot arguments.
- [ ] Invoke existing barcode analysis directly through its typed library API.
- [ ] Invoke Python knee plotting through the Python runtime module.
- [ ] Compare barcode validity, merged/filtered/clean/assigned rows, pair threshold, core cells/barcodes, assignments, knee data, and plots to the oracle.
- [ ] Verify `--upstream-only` exits with the legacy artifact set and status.

### Task 9: Port BAM Tagging, Gene, Matrix, and RNA Cluster Stages

**Files:**
- Modify: `src/workflow/stages.rs`
- Test: `tests/test_workflow_stage_plan.py`

- [ ] Add failing snapshots for prepare tags, CB/UR tags, BED conversion, gene assignment/tagging, BAM index, UMI clustering, cell table, gene matrix, and RNA clustering.
- [ ] Implement Rust stages through typed library calls and bedtools/samtools through the supervised runner.
- [ ] Preserve exact file moves, symlinks, arguments, working directories, and log append behavior.
- [ ] Compare BAM tags/index validity, assignments, barcode-to-cell mappings, matrices, and RNA cluster tables to the oracle.

### Task 10: Port Isoform, QC, Saturation, and Report Stages

**Files:**
- Modify: `src/workflow/stages.rs`
- Test: `tests/test_workflow_stage_plan.py`

- [ ] Add failing snapshots for transcript assignment, isoform matrix, RNA QC, violin fallback, file renames, saturation, read QC, and report arguments.
- [ ] Implement Rust and Python stages through their established typed/runtime interfaces.
- [ ] Preserve `--skip-isoform`, Glycine statistics handling, raw/full-length read accounting, optional plot fallback, and final report location.
- [ ] Compare matrices, QC TSV/JSON/PNG, saturation, report payloads, logs, and empty/absent optional artifacts to the oracle.

### Task 11: Add Mixed-Species Hooks

**Files:**
- Create: `src/workflow/single_species.rs`
- Create: `src/workflow/mixed_species.rs`
- Modify: `src/workflow/stages.rs`
- Test: `tests/test_workflow_stage_plan.py`

- [ ] Add failing tests for `--singlet-threshold`, mixed RNA QC, barnyard QC, and mixed report inputs.
- [ ] Implement mode hooks without duplicating the shared stage sequence.
- [ ] Compare the complete mixed artifact manifest to the legacy oracle.

### Task 12: Replace and Harden Linux Packaging

**Files:**
- Modify: `packaging/build_binary_release.sh`
- Modify: `packaging/refresh_binary_release_metadata.sh`
- Modify: `.gitignore`
- Modify: `tests/test_release_layout.py`
- Modify: `tests/test_single_binary_release.py`

- [ ] Require a Linux x86-64 host/target, `x86_64-unknown-linux-gnu`, baseline CPU flags, Rust path remapping, and a glibc symbol ceiling no newer than 2.17.
- [ ] Build only `cargo build --release --locked --target x86_64-unknown-linux-gnu --bin flora` and stage `target/x86_64-unknown-linux-gnu/release/flora` as `<release-root>/flora`.
- [ ] Copy through an explicit allowlist: deterministic Python assets, report template, Plotly, runtime manifests, runtime environment, requirements, public READMEs, notices/licenses, and `BC_1536.txt` when present.
- [ ] Do not copy shell workflows, standalone bins, `main.pyc`, `args_parser.pyc`, or `utils.pyc` unless an import test proves a runtime dependency and the asset is explicitly approved.
- [ ] Reject source files, extra ELF, symlinks, traversal paths, debug/unstripped binaries, private paths, legacy source names, disallowed dynamic dependencies, and extended attributes.
- [ ] Generate archive and SHA256, then extract and validate the completed archive.
- [ ] In metadata refresh, validate names/types before extraction, build into a new path, fully validate the replacement, and only then atomically replace the original.
- [ ] Run `python -m unittest tests.test_release_layout tests.test_single_binary_release -v`; expect PASS.

### Task 13: Update Documentation and Run Final Parity

**Files:**
- Modify: `README.md`
- Modify: `README_zh-CN.md`
- Modify public/private README templates
- Modify: `tests/compare_workflow_outputs.py`

- [ ] Replace public Bash examples with `./flora run` and `./flora run-mixed`.
- [ ] Document one-command Linux packaging, CPython 3.11 runtime, external dependencies, checksum and compatibility verification, and private legacy regression use.
- [ ] Run clean single/mixed workflows for all skip/light/full/upstream combinations.
- [ ] Interrupt once during every stage, then restart in the same output
  directory and compare with a clean run. Repeat with stale zero-byte, truncated,
  malformed, and mismatched-parameter outputs at each stage boundary; verify
  none is treated as an implicit completion marker.
- [ ] Compare every artifact to the pre-refactor oracle according to the manifest.
- [ ] Run `cargo test` and all Python tests.
- [ ] Build on Linux and verify exactly one stripped ELF, no source/private paths/xattrs, glibc/CPU baseline, valid SHA256, extracted command smoke tests, and compressed size no larger than 28 MB.
