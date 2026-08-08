# Scanpy RNA Cluster Analysis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add sparse Scanpy RNA clustering to both workflows and render Leiden and raw-UMI UMAPs after Beads to cells in the self-contained HTML report.

**Architecture:** A focused Python command converts the wide gene-by-cell TSV into a sparse AnnData object, runs a pinned deterministic Scanpy workflow, and writes a compact per-cell TSV. Both shell workflows invoke that command and pass its absolute output path to the existing report builder, which validates the TSV and renders responsive Plotly traces.

**Tech Stack:** Python 3.14, Scanpy 1.12, NumPy, pandas, SciPy CSR, igraph/Leiden, Plotly embedded in HTML, Bash, Rust/rust-htslib, unittest, Cargo tests.

---

### Task 1: Lock The Python 3.14 Analysis Environment

**Files:**
- Modify: `tests/test_release_layout.py`
- Modify: `environment.yml`
- Modify: `requirements.txt`
- Modify: `README.md`
- Modify: `README_zh-CN.md`

- [ ] **Step 1: Write the failing release-layout assertions**

Update the environment contract test to require:

```python
self.assertIn("python=3.14", environment)
self.assertIn("scanpy>=1.12,<1.13", environment)
self.assertIn("numba>=0.65", environment)
self.assertIn("python-igraph", environment)
self.assertIn("leidenalg", environment)
```

Assert the pip requirements contain `scanpy>=1.12,<1.13`, `numba>=0.65`,
`igraph`, and `leidenalg`, and both README files mention standard CPython 3.14
and the non-default status of free-threading 3.14t.

- [ ] **Step 2: Run the release-layout test and verify RED**

Run:

```bash
python -m unittest tests.test_release_layout -v
```

Expected: FAIL because the current environment still uses Python 3.11 and has
no Scanpy clustering dependencies.

- [ ] **Step 3: Update environment and documentation**

Use standard `python=3.14`, `scanpy>=1.12,<1.13`, and `numba>=0.65`.
Add conda `python-igraph` and `leidenalg`; add pip equivalents to
`requirements.txt`. Update both README languages consistently.

- [ ] **Step 4: Verify GREEN and create an isolated test environment**

Run the release-layout test, then create a standard (non-free-threading)
CPython 3.14 environment outside the source tree:

```bash
conda env create --dry-run -f environment.yml
conda env create -p /private/tmp/strint-fl-py314 -f environment.yml
conda run -p /private/tmp/strint-fl-py314 python -c \
  'import sys, scanpy, numba, igraph, leidenalg; assert sys.version_info[:2] == (3, 14)'
```

Expected: tests PASS, Conda resolves the environment, and all clustering imports
succeed under standard CPython 3.14. Use
`conda run -p /private/tmp/strint-fl-py314` for every Python test in Tasks 3-7.

### Task 2: Make Rust And Python Gene Matrices Preserve The Same Cells

**Files:**
- Modify: `src/bin/gene_expression.rs`
- Create: `tests/test_gene_expression_cell_parity.py`
- Test: unit tests inside `src/bin/gene_expression.rs`

- [ ] **Step 1: Add a passing Python characterization test**

For Python, feed `read_bam_entries` fake alignments containing one normal cell
and one placeholder-only cell, then call `write_matrix` and assert both cell
columns exist.

- [ ] **Step 2: Run the Python characterization test**

```bash
conda run -p /private/tmp/strint-fl-py314 \
  python -m unittest tests.test_gene_expression_cell_parity -v
```

Expected: PASS, documenting the existing Python behavior.

- [ ] **Step 3: Add a failing Rust test**

Define the wished-for record collector API and exercise collection before
placeholder filtering:

```rust
collect_gene_observation(
    "chr1_1000_2000", "CELL_B", "UMI_B", &mut matrix, &mut cells
);
```

Then write a matrix with expression entries only for `CELL_A`. Assert the header
includes both cells and `CELL_B` receives zero for every retained gene.

- [ ] **Step 4: Run the focused Rust test and verify RED**

```bash
cargo test --locked --bin gene_expression
```

Expected: FAIL because columns currently come only from matrix keys.

- [ ] **Step 5: Track all CB tags independently**

Implement `collect_gene_observation` so every valid CB is inserted before
genomic-placeholder GN filtering. Pass the complete sorted cell set into the
writer rather than deriving columns from retained expression entries. Keep the
existing Python behavior unchanged.

- [ ] **Step 6: Verify GREEN**

```bash
conda run -p /private/tmp/strint-fl-py314 \
  python -m unittest tests.test_gene_expression_cell_parity -v
cargo test --locked --bin gene_expression
cargo test --locked
```

### Task 3: Build And Validate The Sparse TSV Loader

**Files:**
- Create: `scripts/rna_cluster_analysis.py`
- Create: `tests/test_rna_cluster_analysis.py`

- [ ] **Step 1: Write failing loader tests**

Test a small gene-by-cell TSV and assert:

```python
loaded.cells == ["CELL_B", "CELL_A"]
loaded.genes == ["G1", "G2"]
loaded.matrix.shape == (2, 2)  # cells x genes
scipy.sparse.isspmatrix_csr(loaded.matrix)
loaded.total_counts.tolist() == [3, 5]
```

Add separate failures for a first column not exactly named `gene`, no cell
columns, blank/duplicate cell headers, blank genes,
nonnumeric/fractional/negative/NaN/Inf/greater-than-uint32 values, plus a valid
fixture whose per-cell sum exceeds `2**32 - 1`. Assert all-zero genes are
removed. Unit-test a checked uint64 accumulation helper and its overflow error.

- [ ] **Step 2: Run the focused tests and verify RED**

```bash
conda run -p /private/tmp/strint-fl-py314 \
  python -m unittest tests.test_rna_cluster_analysis.SparseLoaderTests -v
```

Expected: import/file-not-found failure because the command does not exist.

- [ ] **Step 3: Implement header validation and chunked CSR loading**

Implement:

```python
@dataclass
class LoadedExpression:
    cells: list[str]
    genes: list[str]
    matrix: scipy.sparse.csr_matrix
    total_counts: numpy.ndarray
```

Validate the raw tab-delimited header using `csv.reader`. Read 256 gene rows per
pandas chunk as strings, validate values before conversion, convert numeric
blocks to `uint32` CSR, stack, transpose, convert to CSR, and compute raw totals
using checked uint64 addition. Drop all-zero genes before block storage.

- [ ] **Step 4: Verify GREEN**

```bash
conda run -p /private/tmp/strint-fl-py314 \
  python -m unittest tests.test_rna_cluster_analysis.SparseLoaderTests -v
```

### Task 4: Implement Deterministic Scanpy And Fallback Analysis

**Files:**
- Modify: `scripts/rna_cluster_analysis.py`
- Modify: `tests/test_rna_cluster_analysis.py`

- [ ] **Step 1: Write failing fallback tests**

Cover one, two, and three informative cells, one gene, identical profiles,
all-zero cells, and a mixture of informative and zero-count cells. Assert:

- Input cell order is preserved.
- Every cell appears exactly once.
- Informative fallback rows use cluster `0` and status `fallback`.
- Zero-count rows use cluster/status `unassigned`.
- Analysis records and coordinates are deterministic.

- [ ] **Step 2: Verify fallback tests fail**

```bash
conda run -p /private/tmp/strint-fl-py314 \
  python -m unittest tests.test_rna_cluster_analysis.FallbackAnalysisTests -v
```

Confirm failures are due to missing analysis functions.

- [ ] **Step 3: Implement explicit fallback behavior**

Use deterministic horizontal coordinates. For mixed informative/unassigned
sets, place unassigned points below/outside the informative coordinate range.
Do not catch dependency or arbitrary Scanpy exceptions.

- [ ] **Step 4: Write a failing normal-path integration test**

Create a ranked synthetic matrix with at least six informative cells, two
zero-count cells, and enough varying genes. Execute the actual Scanpy pipeline
twice and assert:

- informative statuses are `scanpy` and zero-count statuses are `unassigned`;
- all coordinates are finite;
- cell order and uint64 totals match input;
- repeated analysis records and coordinates are identical.

- [ ] **Step 5: Verify normal-path test fails**

```bash
conda run -p /private/tmp/strint-fl-py314 \
  python -m unittest tests.test_rna_cluster_analysis.ScanpyAnalysisTests -v
```

Expected: FAIL because the Scanpy path is not implemented.

- [ ] **Step 6: Implement the normal Scanpy path**

Use:

```python
sc.pp.filter_genes(adata, min_cells=1)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(
    adata, flavor="seurat", n_top_genes=min(2000, adata.n_vars)
)
adata = adata[:, adata.var["highly_variable"]].copy()
n_hvg = adata.n_vars
n_comps = min(50, adata.n_obs - 1, n_hvg - 1)
sc.pp.pca(
    adata,
    n_comps=n_comps,
    zero_center=True,
    svd_solver="arpack",
    random_state=0,
)
sc.pp.neighbors(
    adata,
    n_neighbors=min(15, adata.n_obs - 1),
    use_rep="X_pca",
    random_state=0,
)
sc.tl.leiden(
    adata,
    resolution=1,
    random_state=0,
    flavor="igraph",
    n_iterations=2,
    directed=False,
)
sc.tl.umap(adata, random_state=0)
```

Construct AnnData from only `total_counts > 0` cells and merge unassigned rows
back in input order afterward. Assert SciPy sparse type after normalization,
log1p, and HVG subsetting. Check explicit fallback preconditions before Scanpy
calls. After neighbors, fallback only when connectivities have zero edges.
After UMAP, fallback only for wrong-shaped or non-finite coordinates. Let
unexpected exceptions propagate.

- [ ] **Step 7: Verify all cluster-analysis tests pass**

Before the final verification, add injected tests for a zero-edge connectivity
matrix and malformed/non-finite UMAP output. Verify these tests fail against the
new normal path, implement only the two explicit post-Scanpy fallback checks,
then rerun them to GREEN.

```bash
conda run -p /private/tmp/strint-fl-py314 \
  python -m unittest tests.test_rna_cluster_analysis -v
conda run -p /private/tmp/strint-fl-py314 \
  python -m unittest tests.test_rna_cluster_analysis -v
```

Compare deterministic analysis records in the tests.

### Task 5: Add The Command-Line And TSV Serialization Contract

**Files:**
- Modify: `scripts/rna_cluster_analysis.py`
- Modify: `tests/test_rna_cluster_analysis.py`

- [ ] **Step 1: Write failing subprocess tests**

Invoke the script with `--input` and `--output`. Assert successful output has
exactly:

```text
cell	UMAP_1	UMAP_2	leiden	total_counts	status
```

Assert row order matches input cell columns. Run an invalid matrix and assert a
nonzero exit code and actionable stderr. Run the valid command twice and assert
the complete TSV output bytes are identical.

- [ ] **Step 2: Run and verify RED**

```bash
conda run -p /private/tmp/strint-fl-py314 \
  python -m unittest tests.test_rna_cluster_analysis.CliTests -v
```

- [ ] **Step 3: Implement `parse_args`, `main`, and atomic TSV writing**

Require `--input` and `--output`. Write to a sibling temporary file and replace
the destination only after all rows serialize successfully, preventing a
partial TSV after failure.

- [ ] **Step 4: Verify GREEN**

Run all cluster-analysis tests under the Python 3.14 environment.

### Task 6: Integrate The Stage Into Both Workflows

**Files:**
- Create: `tests/test_rna_cluster_workflow.py`
- Modify: `run_all.sh`
- Modify: `run_all_mixed_species.sh`

- [ ] **Step 1: Write failing shell-contract tests**

Assert both scripts:

- require `rna_cluster_analysis.py`;
- invoke it after gene expression;
- write `${MATRIX_DIR}/${SAMPLE_ID}.rna_cluster.tsv`;
- pass absolute `--rna-cluster-tsv` to `build_report.py`;
- list the TSV in key outputs.

- [ ] **Step 2: Run and verify RED**

```bash
conda run -p /private/tmp/strint-fl-py314 \
  python -m unittest tests.test_rna_cluster_workflow -v
```

- [ ] **Step 3: Modify both workflows**

Run the Python command inside `MATRIX_DIR` immediately after gene-expression
generation. Add the absolute matrix output to `BUILD_REPORT_ARGS`, required
script checks, parameter/output summaries, and existing QC log routing.

- [ ] **Step 4: Verify GREEN and Bash syntax**

```bash
conda run -p /private/tmp/strint-fl-py314 \
  python -m unittest tests.test_rna_cluster_workflow -v
bash -n run_all.sh run_all_mixed_species.sh
```

### Task 7: Add RNA Cluster Analysis To The Cells Report

**Files:**
- Modify: `tests/test_build_report_template.py`
- Modify: `scripts/build_report.py`

- [ ] **Step 1: Write failing report loader/payload tests**

Extend report fixtures with a valid cluster TSV. Assert:

- `--rna-cluster-tsv` is required;
- exact columns and values are validated;
- `rnaCluster` payload preserves cell order and totals;
- `clusters` is sorted and unique;
- safe JSON escapes `<`, `>`, `&`, U+2028, and U+2029;

- [ ] **Step 2: Run loader/payload tests and verify RED**

```bash
conda run -p /private/tmp/strint-fl-py314 \
  python -m unittest \
  tests.test_build_report_template.ReportTemplateContractTests.test_rna_cluster_payload_validation \
  tests.test_build_report_template.ReportTemplateContractTests.test_inline_json_is_script_safe -v
```

- [ ] **Step 3: Implement the TSV payload loader and verify it**

Add strict column, uniqueness, finite numeric, integer-total, label, and status
validation. Serialize inline JSON through a helper that escapes HTML/script
terminators.

Run the two focused loader/payload tests and confirm GREEN before changing
markup.

- [ ] **Step 4: Add and run failing markup-order tests**

Assert RNA Cluster Analysis occurs after Beads to cells, plot IDs
`rna-cluster-assignment` and `rna-umi-counts` exist, and fallback/unassigned
status is disclosed in visible section text.

```bash
conda run -p /private/tmp/strint-fl-py314 \
  python -m unittest \
  tests.test_build_report_template.ReportTemplateContractTests.test_rna_cluster_section_follows_beads -v
```

- [ ] **Step 5: Add report markup and verify markup GREEN**

Insert a section after Beads to cells with the existing section-bar/help style
and two responsive plot panels.

- [ ] **Step 6: Add and run failing rendering-contract tests**

Assert `scattergl` is used at 2,000 cells and `scatter` below, one trace is
generated per sorted cluster, hover fields include cell/cluster/raw UMI total,
the UMI plot uses Viridis with color bar title `UMI Counts`, and both plots use
equal axes with hidden ticks/grid.

```bash
conda run -p /private/tmp/strint-fl-py314 \
  python -m unittest \
  tests.test_build_report_template.ReportTemplateContractTests.test_rna_cluster_plot_contract -v
```

- [ ] **Step 7: Add Plotly rendering**

Build one trace per sorted Leiden cluster for Cluster Assignment and one
continuous-color trace for UMI Counts. Use equal axis scaling, hidden ticks,
Viridis, responsive dimensions, and informative hover text.

- [ ] **Step 8: Verify GREEN**

```bash
conda run -p /private/tmp/strint-fl-py314 \
  python -m unittest tests.test_build_report_template -v
```

Inspect a generated fixture HTML for section order, payload, and both plot
containers.

### Task 8: Full Verification And Release Readiness

**Files:**
- Modify if needed: `README.md`
- Modify if needed: `README_zh-CN.md`

- [ ] **Step 1: Run Python tests**

```bash
conda run -p /private/tmp/strint-fl-py314 \
  python -m unittest discover -s tests -v
```

- [ ] **Step 2: Run shell and Rust verification**

```bash
bash -n run_all.sh run_all_mixed_species.sh
cargo metadata --locked --no-deps --format-version 1
cargo test --locked
cargo build --release --locked
```

- [ ] **Step 3: Run a clean-source build**

From the repository root:

```bash
release_dir="$(mktemp -d /private/tmp/strint-fl-release.XXXXXX)"
archive="${release_dir}/strint-fl-source.tar"
tar -cf "${archive}" \
  --exclude='*/__pycache__' \
  --exclude='.DS_Store' \
  .gitignore Cargo.toml Cargo.lock README.md README_zh-CN.md \
  environment.yml requirements.txt args_parser.py main.py utils.py \
  run_all.sh run_all_mixed_species.sh src scripts tests docs \
  vendor/rust-htslib
mkdir -p "${release_dir}/source"
tar -xf "${archive}" -C "${release_dir}/source"
cd "${release_dir}/source"
cargo build --release --locked
cargo test --locked
conda run -p /private/tmp/strint-fl-py314 \
  python -m unittest discover -s tests -v
bash -n run_all.sh run_all_mixed_species.sh
for bin in \
  strint-rust add_cb_ur_tags add_gene_tags assign_genes assign_transcripts \
  cell_umi_gene_table cluster_umis_allbam gene_expression \
  generate_26bp_whitelists isoform_expression prepare_read_tags \
  read_qc_summary rna_qc_metrics
do
  test -x "target/release/${bin}"
done
```

The archive command is the release manifest: it includes only named root files,
source/script/test/docs directories, and the patched rust-htslib tree. The
executable loop verifies every Cargo-declared binary.

- [ ] **Step 4: Verify Conda resolution**

```bash
conda env create --dry-run -f environment.yml
```

- [ ] **Step 5: Review final documentation**

Confirm both README languages describe the RNA cluster TSV, report plots,
Python 3.14 environment, Scanpy dependencies, and standard installation/run
commands.
