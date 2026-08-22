# Parallel Glycine Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add bounded concurrent Glycine processing for multiple FASTQs while preserving one sample-level output contract.

**Architecture:** The embedded workflow runners collect and validate input paths, then call a Rust batch orchestration module. The module launches existing single-file Glycine workers, merges gzip members deterministically, and structurally aggregates statistics before downstream Flora stages continue.

**Tech Stack:** Rust, Clap, standard process/thread APIs, Bash workflow wrappers, Python contract tests.

---

### Task 1: Freeze the command contract

**Files:**
- Modify: `tests/test_release_layout.py`
- Modify: `tests/fixtures/cli/single.json`
- Modify: `tests/fixtures/cli/mixed.json`
- Modify: `tests/test_workflow_contract.py`

- [ ] Add failing tests for multi-value `--fastq`, `--fastq-dir`, and default job/thread controls.
- [ ] Run targeted tests and confirm failure reflects missing behavior.
- [ ] Update both embedded runners consistently.
- [ ] Run targeted tests and commit.

### Task 2: Implement exact output aggregation

**Files:**
- Create: `src/glycine/batch.rs`
- Modify: `src/glycine/mod.rs`

- [ ] Add failing Rust tests for natural sorting, thread allocation, gzip concatenation, and statistics aggregation.
- [ ] Run tests and confirm the expected failures.
- [ ] Implement input discovery and structured statistics aggregation.
- [ ] Run tests and commit.

### Task 3: Implement bounded concurrent execution

**Files:**
- Modify: `src/glycine/batch.rs`
- Modify: `src/main.rs`
- Modify: `run_all.sh`
- Modify: `run_all_mixed_species.sh`

- [ ] Add a failing integration test using small FASTQ fixtures.
- [ ] Implement worker scheduling, chunk directories, deterministic publication, and cleanup.
- [ ] Connect both full workflows to the batch subcommand.
- [ ] Run integration and contract tests and commit.

### Task 4: Document and verify release behavior

**Files:**
- Modify: `README.md`
- Modify: `README_zh-CN.md`
- Modify: `docs/repository-templates/private/README.md`
- Modify: `docs/repository-templates/private/README_zh-CN.md`

- [ ] Document multi-file and directory examples plus resource controls.
- [ ] Run formatting, Python tests, Rust tests, and release build.
- [ ] Inspect branch diff for unrelated changes and commit.
