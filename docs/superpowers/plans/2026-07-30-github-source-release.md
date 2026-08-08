# GitHub Source Release Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Prepare a GitHub source repository that can be cloned or extracted, configured with documented dependencies, and built with `cargo build --release`.

**Architecture:** Publish the complete Rust/Python workflow rather than the `scripts` directory alone. Keep the patched `rust-htslib` source at the path referenced by `Cargo.toml`, provide Conda and pip dependency manifests, and use two linked README files for English and Chinese documentation.

**Tech Stack:** Rust 2021/Cargo, Python 3.11, Conda, minimap2, samtools, bedtools, Bash.

---

### Task 1: Repository packaging rules

**Files:**
- Create: `.gitignore`
- Create: `tests/test_release_layout.py`

- [ ] Add a test that requires the release manifests, bilingual README links, and vendored `rust-htslib/Cargo.toml`.
- [ ] Run the release-layout test and verify it fails before the files are added.
- [ ] Add ignore rules for build products, analysis outputs, large sequencing/reference files, macOS metadata, and template-development directories.

### Task 2: Reproducible dependency manifests

**Files:**
- Create: `environment.yml`
- Create: `requirements.txt`

- [ ] Pin Python to 3.11 and install Rust through the `rust` Conda package, which includes Cargo.
- [ ] Include samtools, minimap2, bedtools, build tools, and Python libraries used by both the Rust workflow and Python fallbacks.
- [ ] Keep pip-only packages in `requirements.txt` and install them from `environment.yml`.

### Task 3: Buildable vendored dependency

**Files:**
- Create from archive: `vendor/rust-htslib/**`
- Retain: `Cargo.toml`
- Retain: `Cargo.lock`

- [ ] Extract `vendor.zip` so `[patch.crates-io] rust-htslib = { path = "vendor/rust-htslib" }` resolves after cloning.
- [ ] Ensure archive metadata is excluded by `.gitignore`.
- [ ] Verify Cargo metadata and release compilation from the repository root.

### Task 4: Bilingual installation and usage documentation

**Files:**
- Create: `README.md`
- Create: `README_zh-CN.md`

- [ ] Add English/Chinese language links at the top of both files.
- [ ] Document repository contents, hardware expectations, Conda setup, Cargo compilation, external Glycine setup, and reference preparation.
- [ ] Add complete `run_all.sh`, `--skip-glycine`, mixed-species, and SGE/qsub examples.
- [ ] List required and excluded GitHub files and explain that FASTQ/reference/output data must not be committed.

### Task 5: Verification

**Files:**
- Test: `tests/test_release_layout.py`
- Test: `tests/test_build_report_template.py`

- [ ] Run Python syntax and unit tests.
- [ ] Run `bash -n run_all.sh run_all_mixed_species.sh`.
- [ ] Run `cargo build --release` in a clean copied source tree.
- [ ] Confirm the expected Rust binaries are present under `target/release/`.
