#!/usr/bin/env python3
"""Strictly compare workflow output trees against frozen legacy oracles."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import subprocess
from pathlib import Path, PurePosixPath


CLASSIFICATIONS = {
    "byte_exact",
    "parsed_exact",
    "numeric_tolerance",
    "canonicalized_html_log",
    "intentionally_absent",
}
REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
ELAPSED_PATTERNS = (
    re.compile(r"(?<=elapsed: )\d+(?:\.\d+)?s"),
    re.compile(r"(?<=Total elapsed: )\d+(?:\.\d+)?s"),
    re.compile(r"(?<=\[Flora\] Total elapsed: )\d+(?:\.\d+)?s"),
)


class ComparisonError(AssertionError):
    pass


def _relative_path(value: str) -> PurePosixPath:
    if not isinstance(value, str) or not value:
        raise ComparisonError("artifact path must be a non-empty string")
    path = PurePosixPath(value)
    if path.is_absolute() or ".." in path.parts or "." in path.parts:
        raise ComparisonError(f"path traversal is forbidden: {value!r}")
    return path


def _inside(root: Path, target: Path) -> bool:
    try:
        target.relative_to(root)
        return True
    except ValueError:
        return False


def _validate_path(root: Path, relative: PurePosixPath) -> Path:
    root = root.resolve()
    path = root.joinpath(*relative.parts)
    if path.is_symlink():
        resolved = path.resolve(strict=False)
        if not _inside(root, resolved):
            raise ComparisonError(f"symlink escapes comparison root: {relative}")
    else:
        parent = path.parent.resolve(strict=False)
        if not _inside(root, parent):
            raise ComparisonError(f"path escapes comparison root: {relative}")
    return path


def sha256_path(path: Path) -> str:
    if path.is_symlink():
        payload = os.readlink(path).encode("utf-8")
    else:
        payload = path.read_bytes()
    return hashlib.sha256(payload).hexdigest()


def _listed_paths(root: Path) -> set[str]:
    if not root.is_dir():
        raise ComparisonError(f"comparison root is not a directory: {root}")
    paths = set()
    for path in root.rglob("*"):
        if path.is_file() or path.is_symlink():
            relative = path.relative_to(root).as_posix()
            _validate_path(root, _relative_path(relative))
            paths.add(relative)
    return paths


def _read_text(path: Path) -> str:
    if path.suffix == ".gz":
        try:
            with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
                return handle.read()
        except (gzip.BadGzipFile, UnicodeDecodeError) as error:
            raise ComparisonError(f"invalid gzip text file: {path.name}") from error
    try:
        return path.read_text(encoding="utf-8")
    except UnicodeDecodeError as error:
        raise ComparisonError(f"invalid UTF-8 text file: {path.name}") from error


def _parse_delimited(path: Path, text: str):
    suffixes = path.suffixes
    effective = suffixes[-2] if suffixes and suffixes[-1] == ".gz" else path.suffix
    delimiter = "\t" if effective in {".tsv", ".bed", ".gtf", ".txt"} else ","
    try:
        return list(csv.reader(text.splitlines(), delimiter=delimiter, strict=True))
    except csv.Error as error:
        raise ComparisonError(f"malformed delimited file: {path.name}") from error


def _parse_exact(path: Path):
    if path.suffix == ".json":
        try:
            return json.loads(_read_text(path), parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value)))
        except (json.JSONDecodeError, ValueError) as error:
            raise ComparisonError(f"malformed JSON: {path.name}") from error
    if path.suffix == ".bam":
        return _bam_records(path, {})
    return _parse_delimited(path, _read_text(path))


def _numeric_token(value: str):
    stripped = value.strip()
    if not stripped:
        return None
    try:
        result = float(stripped)
    except ValueError:
        return None
    if not math.isfinite(result):
        raise ComparisonError(f"non-finite numeric token: {value!r}")
    return result


def _compare_numeric(expected: Path, actual: Path, absolute: float, relative: float):
    expected_rows = _parse_delimited(expected, _read_text(expected))
    actual_rows = _parse_delimited(actual, _read_text(actual))
    if len(expected_rows) != len(actual_rows):
        raise ComparisonError(f"row count differs: {expected.name}")
    for row_index, (expected_row, actual_row) in enumerate(zip(expected_rows, actual_rows)):
        if len(expected_row) != len(actual_row):
            raise ComparisonError(f"column count differs at row {row_index}: {expected.name}")
        for column_index, (left, right) in enumerate(zip(expected_row, actual_row)):
            left_number = _numeric_token(left)
            right_number = _numeric_token(right)
            if left_number is None or right_number is None:
                if left != right:
                    raise ComparisonError(
                        f"text differs at row {row_index}, column {column_index}: {expected.name}"
                    )
            elif not math.isclose(left_number, right_number, rel_tol=relative, abs_tol=absolute):
                raise ComparisonError(
                    f"number differs at row {row_index}, column {column_index}: {expected.name}"
                )


def _canonical_text(path: Path, roots: tuple[Path, Path], approved: list[str]) -> str:
    allowed = {"elapsed_seconds", "workspace_root"}
    unknown = set(approved) - allowed
    if unknown:
        raise ComparisonError(f"unapproved canonicalization rules: {sorted(unknown)}")
    text = _read_text(path).replace("\r\n", "\n")
    if "workspace_root" in approved:
        text = text.replace(str(REPOSITORY_ROOT), "<REPOSITORY_ROOT>")
        for root in roots:
            text = text.replace(str(root.resolve()), "<WORKFLOW_ROOT>")
    if "elapsed_seconds" in approved:
        for pattern in ELAPSED_PATTERNS:
            text = pattern.sub("<ELAPSED>", text)
    return text


def _samtools(metadata: dict) -> str:
    configured = metadata.get("samtools") or os.environ.get("FLORA_SAMTOOLS") or "samtools"
    path = Path(configured)
    if not path.is_absolute() and "/" in configured:
        path = REPOSITORY_ROOT / path
        return str(path)
    return configured


def _bam_records(path: Path, metadata: dict) -> list[str]:
    samtools = _samtools(metadata)
    checks = ([samtools, "quickcheck", "-v", str(path)],)
    if metadata.get("require_index", True):
        checks += ([samtools, "idxstats", str(path)],)
    for command in checks:
        result = subprocess.run(command, text=True, capture_output=True, check=False)
        if result.returncode:
            raise ComparisonError(
                f"BAM validation failed for {path.name}: {result.stderr.strip()}"
            )
    result = subprocess.run(
        [samtools, "view", "-h", "--no-PG", str(path)],
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode:
        raise ComparisonError(f"unable to parse BAM: {path.name}")
    records = []
    required_tags = set(metadata.get("required_tags", []))
    for line in result.stdout.splitlines():
        if line.startswith("@"):
            if not line.startswith("@PG"):
                records.append(line)
            continue
        fields = line.split("\t")
        tags = {field.split(":", 1)[0] for field in fields[11:]}
        missing = required_tags - tags
        if missing:
            raise ComparisonError(f"BAM record is missing tags {sorted(missing)}: {path.name}")
        records.append(line)
    return records


def _compare_artifact(expected: Path, actual: Path, artifact: dict, roots: tuple[Path, Path], tolerance: dict):
    classification = artifact["classification"]
    if expected.is_symlink() or actual.is_symlink():
        if not expected.is_symlink() or not actual.is_symlink():
            raise ComparisonError(f"symlink type differs: {artifact['path']}")
        expected_target = expected.resolve().relative_to(roots[0].resolve()).as_posix()
        actual_target = actual.resolve().relative_to(roots[1].resolve()).as_posix()
        if expected_target != actual_target:
            raise ComparisonError(f"symlink target differs: {artifact['path']}")
        return
    if classification == "byte_exact":
        if expected.read_bytes() != actual.read_bytes():
            raise ComparisonError(f"bytes differ: {artifact['path']}")
    elif classification == "parsed_exact":
        if expected.suffix == ".bam":
            left = _bam_records(expected, artifact)
            right = _bam_records(actual, artifact)
        else:
            left = _parse_exact(expected)
            right = _parse_exact(actual)
        if left != right:
            raise ComparisonError(f"parsed content differs: {artifact['path']}")
    elif classification == "numeric_tolerance":
        _compare_numeric(
            expected,
            actual,
            float(tolerance["absolute"]),
            float(tolerance["relative"]),
        )
    elif classification == "canonicalized_html_log":
        approved = artifact.get("canonicalize", [])
        if _canonical_text(expected, roots, approved) != _canonical_text(actual, roots, approved):
            raise ComparisonError(f"canonical content differs: {artifact['path']}")
    else:
        raise ComparisonError(f"unsupported classification: {classification}")


def compare_tree(expected_root: Path, actual_root: Path, manifest: dict, scenario: str | None = None):
    expected_root = Path(expected_root)
    actual_root = Path(actual_root)
    artifacts = [
        artifact
        for artifact in manifest.get("artifacts", [])
        if scenario is None or artifact.get("scenario") == scenario
    ]
    if not artifacts:
        raise ComparisonError("manifest has no artifacts for comparison")
    declared = set()
    oracle_declared = set()
    expected_present = set()
    actual_present = set()
    for artifact in artifacts:
        if set(artifact) - {
            "canonicalize", "classification", "oracle_path", "path", "required_tags", "require_index", "samtools", "scenario", "sha256"
        }:
            raise ComparisonError(f"unknown artifact keys for {artifact.get('path')!r}")
        relative = _relative_path(artifact.get("path"))
        path_text = relative.as_posix()
        if path_text in declared:
            raise ComparisonError(f"duplicate artifact path: {path_text}")
        declared.add(path_text)
        oracle_relative = _relative_path(artifact.get("oracle_path", artifact.get("path")))
        oracle_text = oracle_relative.as_posix()
        if oracle_text in oracle_declared:
            raise ComparisonError(f"duplicate oracle artifact path: {oracle_text}")
        oracle_declared.add(oracle_text)
        classification = artifact.get("classification")
        if classification not in CLASSIFICATIONS:
            raise ComparisonError(f"invalid classification for {path_text}")
        expected = _validate_path(expected_root, oracle_relative)
        actual = _validate_path(actual_root, relative)
        if classification == "intentionally_absent":
            if expected.exists() or expected.is_symlink() or actual.exists() or actual.is_symlink():
                raise ComparisonError(f"intentionally absent path exists: {path_text}")
            continue
        expected_present.add(oracle_text)
        actual_present.add(path_text)
        if not expected.exists() and not expected.is_symlink():
            raise ComparisonError(f"missing oracle path: {path_text}")
        if not actual.exists() and not actual.is_symlink():
            raise ComparisonError(f"missing actual path: {path_text}")
        expected_hash = artifact.get("sha256", "")
        if sha256_path(expected) != expected_hash:
            raise ComparisonError(f"frozen oracle checksum differs: {path_text}")
        _compare_artifact(
            expected,
            actual,
            artifact,
            (expected_root, actual_root),
            manifest["numeric_tolerance"],
        )
    expected_paths = _listed_paths(expected_root)
    actual_paths = _listed_paths(actual_root)
    if expected_paths != expected_present:
        missing = sorted(expected_present - expected_paths)
        extra = sorted(expected_paths - expected_present)
        raise ComparisonError(f"oracle paths differ: missing={missing}, extra={extra}")
    if actual_paths != actual_present:
        missing = sorted(actual_present - actual_paths)
        extra = sorted(actual_paths - actual_present)
        raise ComparisonError(f"actual paths differ: missing={missing}, extra={extra}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--expected", required=True, type=Path)
    parser.add_argument("--actual", required=True, type=Path)
    parser.add_argument("--scenario")
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    try:
        compare_tree(args.expected, args.actual, manifest, args.scenario)
    except ComparisonError as error:
        parser.exit(1, f"workflow comparison failed: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
