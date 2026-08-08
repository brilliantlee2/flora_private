#!/usr/bin/env python3
"""Compile Flora runtime scripts without embedding private source paths."""

from __future__ import annotations

import argparse
import importlib.util
import marshal
import py_compile
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--public-prefix", default="flora/scripts")
    parser.add_argument("sources", nargs="+", type=Path)
    return parser.parse_args()


def code_filenames(code: object) -> set[str]:
    filenames: set[str] = set()
    if hasattr(code, "co_filename"):
        filenames.add(code.co_filename)
        for value in code.co_consts:
            filenames.update(code_filenames(value))
    return filenames


def load_code(path: Path) -> object:
    with path.open("rb") as handle:
        handle.read(16)
        return marshal.load(handle)


def main() -> int:
    args = parse_args()
    if sys.implementation.name != "cpython" or sys.version_info[:2] != (3, 11):
        raise SystemExit(
            "Flora release bytecode must be compiled with CPython 3.11; "
            f"found {sys.implementation.name} {sys.version.split()[0]}"
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for source in args.sources:
        source = source.resolve(strict=True)
        output = args.output_dir / f"{source.stem}.pyc"
        public_name = f"{args.public_prefix.rstrip('/')}/{source.name}"
        py_compile.compile(
            str(source),
            cfile=str(output),
            dfile=public_name,
            doraise=True,
            invalidation_mode=py_compile.PycInvalidationMode.CHECKED_HASH,
        )
        filenames = code_filenames(load_code(output))
        if filenames != {public_name}:
            raise SystemExit(
                f"Unexpected code filenames in {output}: {sorted(filenames)}"
            )

    cache_tag = sys.implementation.cache_tag or "unknown"
    (args.output_dir.parent / "PYTHON_ABI.txt").write_text(
        f"implementation=cpython\nmajor_minor=3.11\ncache_tag={cache_tag}\n",
        encoding="ascii",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
