#!/usr/bin/env python3
"""Run a compact Scanpy RNA clustering workflow on a gene-expression TSV."""

from __future__ import annotations

import argparse
import csv
import os
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


UINT32_MAX = int(np.iinfo(np.uint32).max)
UINT64_MAX = np.iinfo(np.uint64).max
COUNT_PATTERN = re.compile(r"^[0-9]+$")
OUTPUT_COLUMNS = [
    "cell",
    "UMAP_1",
    "UMAP_2",
    "leiden",
    "total_counts",
    "status",
]


@dataclass
class LoadedExpression:
    cells: list[str]
    genes: list[str]
    matrix: sparse.csr_matrix
    total_counts: np.ndarray


@dataclass(frozen=True)
class ClusterRecord:
    cell: str
    umap1: float
    umap2: float
    leiden: str
    total_counts: int
    status: str


def checked_add_uint64(current: np.ndarray, increment: np.ndarray) -> np.ndarray:
    """Add uint64 vectors while rejecting overflow."""
    current = np.asarray(current, dtype=np.uint64)
    increment = np.asarray(increment, dtype=np.uint64)
    if current.shape != increment.shape:
        raise ValueError("count vectors have incompatible shapes")
    if np.any(increment > UINT64_MAX - current):
        raise ValueError("total UMI count exceeds uint64 range")
    return current + increment


def _read_header(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        try:
            header = next(csv.reader(handle, delimiter="\t"))
        except StopIteration as exc:
            raise ValueError("gene-expression TSV is empty") from exc

    if len(header) < 2 or header[0] != "gene":
        raise ValueError("gene-expression TSV must start with a 'gene' column")

    cells = header[1:]
    if any(not cell.strip() for cell in cells):
        raise ValueError("gene-expression TSV contains a blank cell name")
    if len(set(cells)) != len(cells):
        raise ValueError("gene-expression TSV contains duplicate cell names")
    return cells


def _parse_count_block(frame: pd.DataFrame, expected_columns: list[str]) -> np.ndarray:
    if frame.columns.tolist() != ["gene", *expected_columns]:
        raise ValueError("gene-expression TSV columns changed while reading")

    values = frame.iloc[:, 1:].to_numpy(dtype=object, copy=False)
    flat = pd.Series(values.reshape(-1), dtype="string", copy=False)
    valid = flat.str.fullmatch(COUNT_PATTERN).fillna(False)
    if not bool(valid.all()):
        token = flat.iloc[int(np.flatnonzero(~valid.to_numpy())[0])]
        raise ValueError(f"invalid UMI count {token!r} in gene-expression TSV")

    numeric = pd.to_numeric(flat, errors="raise")
    if bool((numeric > UINT32_MAX).any()):
        value = int(numeric[numeric > UINT32_MAX].iloc[0])
        raise ValueError(f"UMI count {value} exceeds uint32 range")
    return numeric.to_numpy(dtype=np.uint32).reshape(values.shape)


def load_expression_matrix(
    path: str | Path, chunk_size: int = 256
) -> LoadedExpression:
    """Load a genes-by-cells TSV as a cells-by-genes sparse matrix."""
    input_path = Path(path)
    cells = _read_header(input_path)
    genes: list[str] = []
    blocks: list[sparse.csr_matrix] = []
    totals = np.zeros(len(cells), dtype=np.uint64)

    with pd.read_csv(
        input_path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        chunksize=chunk_size,
    ) as reader:
        for frame in reader:
            gene_names = frame.iloc[:, 0].tolist()
            if any(not str(gene).strip() for gene in gene_names):
                raise ValueError("gene-expression TSV contains a blank gene name")

            counts = _parse_count_block(frame, cells)
            totals = checked_add_uint64(
                totals,
                counts.sum(axis=0, dtype=np.uint64),
            )
            keep = np.any(counts != 0, axis=1)
            if np.any(keep):
                genes.extend(str(gene) for gene in np.asarray(gene_names)[keep])
                blocks.append(sparse.csr_matrix(counts[keep].T, dtype=np.uint32))

    if blocks:
        matrix = sparse.hstack(blocks, format="csr", dtype=np.uint32)
    else:
        matrix = sparse.csr_matrix((len(cells), 0), dtype=np.uint32)
    return LoadedExpression(cells, genes, matrix, totals)


def _centered_positions(count: int, y: float = 0.0) -> np.ndarray:
    if count == 0:
        return np.empty((0, 2), dtype=np.float64)
    x = np.arange(count, dtype=np.float64) - (count - 1) / 2.0
    return np.column_stack((x, np.full(count, y, dtype=np.float64)))


def _assemble_records(
    data: LoadedExpression,
    informative_indices: np.ndarray,
    informative_coordinates: np.ndarray,
    informative_labels: list[str],
    informative_status: str,
) -> list[ClusterRecord]:
    coordinates = np.zeros((len(data.cells), 2), dtype=np.float64)
    labels = np.full(len(data.cells), "unassigned", dtype=object)
    statuses = np.full(len(data.cells), "unassigned", dtype=object)

    if informative_indices.size:
        coordinates[informative_indices] = informative_coordinates
        labels[informative_indices] = informative_labels
        statuses[informative_indices] = informative_status

    unassigned = np.flatnonzero(data.total_counts == 0)
    if unassigned.size:
        if informative_indices.size:
            x_center = float(np.mean(informative_coordinates[:, 0]))
            y_min = float(np.min(informative_coordinates[:, 1]))
            y_span = float(np.ptp(informative_coordinates[:, 1]))
            gap = max(1.0, y_span * 0.15)
            positions = _centered_positions(unassigned.size, y_min - gap)
            positions[:, 0] += x_center
        else:
            positions = _centered_positions(unassigned.size)
        coordinates[unassigned] = positions

    return [
        ClusterRecord(
            cell=cell,
            umap1=float(coordinates[index, 0]),
            umap2=float(coordinates[index, 1]),
            leiden=str(labels[index]),
            total_counts=int(data.total_counts[index]),
            status=str(statuses[index]),
        )
        for index, cell in enumerate(data.cells)
    ]


def _fallback_records(
    data: LoadedExpression,
    informative_indices: np.ndarray,
) -> list[ClusterRecord]:
    coordinates = _centered_positions(informative_indices.size)
    return _assemble_records(
        data,
        informative_indices,
        coordinates,
        ["0"] * informative_indices.size,
        "fallback",
    )


def _profiles_are_identical(matrix: sparse.csr_matrix) -> bool:
    if matrix.shape[0] < 2:
        return True
    first = matrix.getrow(0)
    return all((matrix.getrow(index) != first).nnz == 0 for index in range(1, matrix.shape[0]))


def analyze_expression(data: LoadedExpression) -> list[ClusterRecord]:
    """Cluster informative cells and retain zero-count cells as unassigned."""
    informative = np.flatnonzero(data.total_counts > 0)
    informative_matrix = data.matrix[informative]

    if (
        informative.size < 4
        or informative_matrix.shape[1] < 2
        or _profiles_are_identical(informative_matrix)
    ):
        return _fallback_records(data, informative)

    import anndata as ad
    import scanpy as sc

    adata = ad.AnnData(
        X=informative_matrix.astype(np.float32),
        obs=pd.DataFrame(index=[data.cells[index] for index in informative]),
        var=pd.DataFrame(index=data.genes),
    )
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat",
        n_top_genes=min(2000, adata.n_vars),
    )
    highly_variable = np.asarray(adata.var["highly_variable"], dtype=bool)
    if int(highly_variable.sum()) < 2:
        return _fallback_records(data, informative)

    adata = adata[:, highly_variable].copy()
    n_comps = min(50, adata.n_obs - 1, adata.n_vars - 1)
    if n_comps < 1:
        return _fallback_records(data, informative)

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
    connectivities = adata.obsp.get("connectivities")
    if connectivities is None or connectivities.nnz == 0:
        return _fallback_records(data, informative)

    sc.tl.leiden(
        adata,
        resolution=1.0,
        random_state=0,
        flavor="igraph",
        n_iterations=2,
        directed=False,
        key_added="leiden",
    )
    sc.tl.umap(adata, random_state=0)
    coordinates = np.asarray(adata.obsm["X_umap"], dtype=np.float64)
    if coordinates.shape != (adata.n_obs, 2) or not np.all(np.isfinite(coordinates)):
        return _fallback_records(data, informative)

    labels = adata.obs["leiden"].astype(str).tolist()
    return _assemble_records(
        data,
        informative,
        coordinates,
        labels,
        "scanpy",
    )


def _format_coordinate(value: float) -> str:
    return format(value, ".8g")


def write_records_atomic(records: list[ClusterRecord], output_path: str | Path) -> None:
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_name = handle.name
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(OUTPUT_COLUMNS)
            for record in records:
                writer.writerow(
                    [
                        record.cell,
                        _format_coordinate(record.umap1),
                        _format_coordinate(record.umap2),
                        record.leiden,
                        str(record.total_counts),
                        record.status,
                    ]
                )
        os.replace(temporary_name, destination)
    except Exception:
        if temporary_name is not None:
            try:
                os.unlink(temporary_name)
            except FileNotFoundError:
                pass
        raise


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run Scanpy RNA clustering on a Flora gene-expression TSV."
    )
    parser.add_argument("--input", required=True, help="Input *.gene_expression.tsv")
    parser.add_argument("--output", required=True, help="Output *.rna_cluster.tsv")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data = load_expression_matrix(args.input)
    records = analyze_expression(data)
    write_records_atomic(records, args.output)
    assigned = sum(record.status != "unassigned" for record in records)
    print(
        f"RNA clustering complete: cells={len(records)}, "
        f"analyzed={assigned}, output={args.output}"
    )


if __name__ == "__main__":
    main()
