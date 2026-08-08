# Scanpy RNA Cluster Analysis Design

## Goal

Add an RNA clustering stage to both Strint-FL workflows. The stage consumes the
final gene-by-cell UMI count matrix at
`matrix/<sample>.gene_expression.tsv`, preserves every final called cell, and
adds two UMAP views to the Cells page immediately after Beads to cells:

1. Leiden cluster assignment.
2. Raw UMI count per cell.

The Scanpy stage writes a compact, reusable TSV. The report builder renders the
plots with its bundled Plotly runtime so the final HTML remains self-contained.

## Input And Output Contract

Input:

```text
gene    CELL1    CELL2    ...
GENE1   0        3        ...
GENE2   1        0        ...
```

The first column contains gene identifiers. Remaining columns are final cell
identifiers and values are deduplicated UMI counts.

Output:

`matrix/<sample>.rna_cluster.tsv`

```text
cell    UMAP_1    UMAP_2    leiden    total_counts    status
CELL1   ...       ...       0         1234            scanpy
```

Requirements:

- Exactly one output row per input cell column.
- Input cell order is preserved.
- `total_counts` is calculated from the original unnormalized UMI matrix.
- Leiden labels are strings and use the key `leiden`.
- UMAP coordinates and cluster labels are deterministic for a fixed input.
- No `.h5ad` output is required.

## Sparse Matrix Loading

The expression TSV may be large and must not be loaded as one dense DataFrame.
The loader reads genes in bounded chunks, converts each numeric chunk to a
SciPy CSR matrix, stacks the chunks as gene by cell, and transposes once to a
cell-by-gene CSR matrix for AnnData.

The loader validates:

- A `gene` column and at least one cell column exist.
- Cell identifiers are non-empty and unique.
- Gene identifiers are non-empty.
- Counts are numeric, finite, non-negative integers representable as unsigned
  32-bit values.

Only genes observed in at least one cell are retained. Cells are never filtered.
The raw header is validated with Python's `csv` module before pandas sees it,
because pandas otherwise mangles duplicate or blank column names.

The loader uses 256-gene chunks. Each chunk is parsed as `uint32`, converted
immediately to CSR, and appended to the gene-block list. After `scipy.sparse.vstack`,
the gene-by-cell matrix is transposed once to cell-by-gene CSR and the block list
is released. The Scanpy preprocessing matrix must remain sparse through gene
filtering, normalization, log transformation, and HVG subsetting. Raw
`total_counts` are accumulated as `uint64` with overflow checks rather than in
the matrix's `uint32` storage dtype.

### Cell-set consistency

The Python gene-expression writer records every CB barcode before excluding
genomic-placeholder genes, but the current Rust writer derives cell columns only
from retained `(gene, cell)` entries. Update the Rust writer to track all CB
barcodes independently and write the same sorted cell columns as Python. Add a
Rust regression test for a cell whose reads only have placeholder genes.

Cells with zero retained gene UMIs cannot be clustered from RNA expression.
They remain in the cluster TSV and both report plots with Leiden label
`unassigned`. They are excluded only from PCA, neighbors, Leiden, and UMAP
fitting, then placed at deterministic coordinates outside the informative-cell
embedding. This preserves all final cells without inventing an RNA cluster for
cells having no usable gene signal.

## Scanpy Workflow

The recommended Conda environment uses standard CPython 3.14 and pins
`scanpy>=1.12,<1.13` plus `numba>=0.65`. Explicit arguments below make the result reproducible
across that supported range instead of relying on defaults that may change.
Free-threading CPython 3.14t is not a default supported target because native
extensions such as pysam may re-enable the GIL and have not all declared
free-threading safety.

For normal data with enough informative cells and matrix rank:

1. Store raw total UMI counts in `adata.obs["total_counts"]`.
2. Remove genes with zero counts using `sc.pp.filter_genes(min_cells=1)`.
3. Normalize total counts with `sc.pp.normalize_total`.
4. Apply `sc.pp.log1p`.
5. Select up to 2,000 highly variable genes using
   `sc.pp.highly_variable_genes(flavor="seurat",
   n_top_genes=min(2000, n_vars))`.
6. Subset to the selected features and calculate
   `n_comps=min(50, n_obs - 1, n_hvg - 1)`.
7. Run sparse-compatible PCA with `random_state=0`, `svd_solver="arpack"`,
   `zero_center=True`, and the explicit component count.
8. Build the neighbor graph from `X_pca` with `n_neighbors=min(15, n_obs - 1)`
   `use_rep="X_pca"`, and `random_state=0`.
9. Run Leiden with `resolution=1`, `random_state=0`,
   `flavor="igraph"`, and `n_iterations=2`.
10. Run UMAP with `random_state=0`.

The workflow deliberately omits cell QC filtering, regression, and scaling.
This preserves the final cell set and avoids sparse-to-dense memory expansion.

If there are fewer than four informative cells, fewer than two selected
features, `n_comps < 1`, identical informative profiles, a neighbor graph with
zero edges, or a UMAP result with the wrong shape/non-finite coordinates, use
deterministic fallback coordinates on a horizontal line and assign informative
cells cluster `0`.
This supports rank-deficient matrices and pipeline smoke tests without
pretending that a biological clustering was performed. The output records the
analysis mode in a `status` column (`scanpy`, `fallback`, or `unassigned`) so
the report can disclose fallback results.

Fallback is selected only from these explicit, testable conditions. Dependency
errors, Scanpy API errors, and other unexpected exceptions are not caught and
must stop the pipeline with their original error context.

If no genes remain, the input is invalid and the stage exits with an actionable
error unless all cells are zero-count cells; in that case all rows are emitted
as `unassigned`.

## Workflow Integration

Add `scripts/rna_cluster_analysis.py` and require it in both:

- `run_all.sh`
- `run_all_mixed_species.sh`

Run it once inside `MATRIX_DIR`, after gene expression matrix generation and
before report building:

```bash
python3 scripts/rna_cluster_analysis.py \
  --input <sample>.gene_expression.tsv \
  --output <sample>.rna_cluster.tsv
```

The stage is required. Missing Scanpy dependencies or invalid matrices stop the
pipeline with a clear error. Report building currently runs from `QC_DIR`, so
both workflows pass absolute paths to `build_report.py`:

```text
--rna-cluster-tsv "${MATRIX_DIR}/<sample>.rna_cluster.tsv"
```

The final key-output summary lists the cluster TSV.

## Report Integration

`build_report.py` loads and validates the compact cluster TSV and adds
`payload.rnaCluster` with:

- `cells`
- `umap1`
- `umap2`
- `leiden`
- `totalCounts`
- sorted unique `clusters`
- `status`

The report loader requires the exact output columns, unique non-empty cells,
finite UMAP coordinates, non-empty cluster/status strings, and finite
non-negative integral totals.

The Cells tab gains an RNA Cluster Analysis section after Beads to cells.

Left panel:

- Title: Cluster Assignment
- One scatter trace per Leiden cluster
- Categorical palette matching the report theme
- Cluster legend
- Hover text containing cell ID, cluster, and raw UMI count

Right panel:

- Title: UMI Counts
- One scatter trace
- Marker color is raw `total_counts`
- Continuous Viridis color bar labelled UMI Counts
- Hover text containing cell ID and raw UMI count

Both plots use `scattergl` at 2,000 or more cells and regular `scatter` below
that threshold, responsive Plotly sizing,
equal axis scaling, no grid, hidden tick labels, and the existing report panel
style. The section includes a help button explaining both plots.

Before embedding the payload in an inline script, JSON is serialized with
`ensure_ascii=False` and escapes `<`, `>`, `&`, U+2028, and U+2029. This
prevents cell identifiers from terminating the script while keeping the report
self-contained.

## Dependencies

Add to `environment.yml`:

- `python=3.14`
- `scanpy>=1.12,<1.13`
- `numba>=0.65`
- `python-igraph`
- `leidenalg`

Add to `requirements.txt`:

- `scanpy>=1.12,<1.13`
- `numba>=0.65`
- `igraph`
- `leidenalg`

Scanpy supplies AnnData, SciPy, scikit-learn, UMAP, and related Python
dependencies. Conda remains the recommended installation method.

Update both README languages to recommend standard CPython 3.14 instead of
Python 3.11 and to state that free-threading 3.14t is not the default supported
target. Update `tests/test_release_layout.py` so its environment contract
expects Python 3.14, the Scanpy/Numba ranges, and the graph dependencies.

## Tests

Add tests before implementation for:

1. Sparse loader preserves matrix orientation, cell order, and raw counts.
2. Rust and Python matrix generation both preserve placeholder-only cells.
3. One-, two-, and three-cell fallback preserves all cells.
4. One-gene, identical-profile, disconnected-graph, and zero-count-cell inputs
   produce deterministic valid output without silently dropping cells.
5. Invalid duplicate/blank cells and nonnumeric, fractional, negative, NaN,
   Inf, or overflowing counts fail clearly.
6. Valid `uint32` entries whose per-cell sum exceeds `2**32 - 1` retain the
   correct `uint64` total.
7. A ranked fixture with at least four informative cells executes PCA,
   neighbors, Leiden, and UMAP; returns `status=scanpy`; preserves cell
   order/raw totals; emits finite coordinates; and is byte-identical over two
   runs.
8. Report CLI accepts and validates `--rna-cluster-tsv`, including missing
   columns and malformed values.
9. RNA Cluster Analysis appears after Beads to cells on the Cells tab.
10. Both Plotly containers, payload fields, safe JSON escaping, and the 2,000
   cell `scattergl` threshold are present.
11. Both workflow scripts require and invoke the clustering script using absolute
   report paths.
12. Environment and pip dependency lists include pinned Scanpy clustering
    packages; release-layout tests and both README languages require standard
    Python 3.14 and no longer recommend Python 3.11.

The complete existing Python, report, Bash syntax, Cargo metadata, and Rust test
suites must continue to pass.
