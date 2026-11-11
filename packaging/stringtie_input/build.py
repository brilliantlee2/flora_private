#!/usr/bin/env python3
"""Build optional, separately named stock/accelerated StringTie 2.2.2 binaries.

No installation, PATH manipulation, workflow changes or implicit downloads.
System zlib/bzip2/lzma development libraries are required. libdeflate is built
from the exact commit selected by the upstream StringTie build script.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import tarfile

HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parents[1]
STRINGTIE = {
    "url": "https://codeload.github.com/gpertea/stringtie/tar.gz/refs/tags/v2.2.2",
    "sha256": "17382efe7240c0efb2ac5efd6493882d58797062b3cc713344fb1a0211cb8ce1",
    "root": "stringtie-2.2.2",
}
LIBDEFLATE = {
    "url": "https://codeload.github.com/ebiggers/libdeflate/tar.gz/9b565afd996d8b798fc7b94cddcc7cfa49293050",
    "sha256": "8bc23999e454da31b8117db932b859107f1023d144ce89b4165cb819bfacdf2d",
    "root": "libdeflate-9b565afd996d8b798fc7b94cddcc7cfa49293050",
}
SOURCE_HASHES = {
    "rlink.cpp": "c0623a6f02f130ecfe914f60fcfbf05bf049ca734be220d9f8ede435307c0430",
    "rlink.h": "32478a9ca71ef7b5232a0ad772e75fc81575a4474c7ff8e7e3804c5b53c75d43",
}
OLD_LOOP = """\tif(!mergeMode) while(n>-1 && readlist[n]->start==brec.start) {
\t\tif(strand==readlist[n]->strand && (readlist[n]->longread==longr) && (!isunitig || (unitig_cov>0) == readlist[n]->unitig)) {
\t\t\tmatch=exonmatch(readlist[n]->segs,brec.exons);
\t\t\tif(match && (longreads || mixedMode)) match=deljuncmatch(readlist[n],brec.juncsdel); //DEL AWARE
\t\t}
\t\t//if(strand==readlist[n]->strand) match=exonmatch(readlist[n]->segs,brec.exons);
\t\tif(match) break; // this way I make sure that I keep the n of the matching readlist
\t\tn--;
\t}
"""
NEW_LOOP = """\tif(!mergeMode) {
\t\tif (!StringTieSameStartIndex::enabled()) {
""" + OLD_LOOP.replace("\tif(!mergeMode) while(", "\t\twhile(", 1) + """\t\t}
\t\telse {
\t\tconst std::vector<int>& candidates=bdata.same_start_candidates.find(
\t\t\treadlist, brec.exons, strand, longr, brec.start);
\t\tfor (auto candidate=candidates.rbegin(); candidate!=candidates.rend(); ++candidate) {
\t\t\tn=*candidate;
\t\t\tif(strand==readlist[n]->strand && (readlist[n]->longread==longr) && (!isunitig || (unitig_cov>0) == readlist[n]->unitig)) {
\t\t\t\tmatch=exonmatch(readlist[n]->segs,brec.exons);
\t\t\t\tif(match && (longreads || mixedMode)) match=deljuncmatch(readlist[n],brec.juncsdel); //DEL AWARE
\t\t\t}
\t\t\tif(match) break;
\t\t}
\t\t}
\t}
"""


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def verify(path, expected):
    actual = sha256(path)
    if actual != expected:
        raise ValueError(f"SHA256 mismatch for {path}: expected {expected}, got {actual}")


def extract_verified(archive, destination, spec):
    verify(archive, spec["sha256"])
    if destination.exists():
        raise ValueError(f"refusing existing extraction directory: {destination}")
    with tarfile.open(archive, "r:gz") as tar:
        members = tar.getmembers()
        for member in members:
            path = Path(member.name)
            if (path.is_absolute() or ".." in path.parts or not path.parts
                    or path.parts[0] != spec["root"] or not (member.isfile() or member.isdir())):
                raise ValueError(f"unsafe archive member: {member.name}")
        destination.mkdir(parents=True)
        tar.extractall(destination, members=members)
    return destination / spec["root"]


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError("upstream patch anchor is not unique; refusing patch")
    return text.replace(old, new, 1)


def apply_accelerator(source):
    for name, expected in SOURCE_HASHES.items():
        verify(source / name, expected)
    cpp = replace_once((source / "rlink.cpp").read_text(), OLD_LOOP, NEW_LOOP)
    header = (source / "rlink.h").read_text()
    header = replace_once(header, '#include "GArgs.h"', '#include "same_start_index.h"\n#include "GArgs.h"')
    header = replace_once(header, " GList<CReadAln> readlist;", " GList<CReadAln> readlist;\n StringTieSameStartIndex same_start_candidates;")
    header = replace_once(header, " void getReady(int currentstart, int currentend) {",
                          " void getReady(int currentstart, int currentend) {\n\t same_start_candidates.clear();")
    header = replace_once(header, "\treadlist.Clear();", "\tsame_start_candidates.clear();\n\treadlist.Clear();")
    (source / "rlink.cpp").write_text(cpp)
    (source / "rlink.h").write_text(header)
    shutil.copyfile(HERE / "same_start_index.h", source / "same_start_index.h")


def obtain(path, spec, fetch):
    if not path.exists():
        if not fetch:
            raise ValueError(f"missing {path}; supply archive or explicitly use --fetch")
        path.parent.mkdir(parents=True, exist_ok=True)
        partial = path.with_suffix(path.suffix + ".part")
        subprocess.run(["curl", "--fail", "--location", "--proto", "=https", "--retry", "2",
                        "--output", str(partial), spec["url"]], check=True)
        verify(partial, spec["sha256"])
        partial.replace(path)
    verify(path, spec["sha256"])
    return path


def run(command, cwd, log):
    print("Building:", " ".join(command), flush=True)
    with log.open("ab") as stream:
        subprocess.run(command, cwd=cwd, stdout=stream, stderr=subprocess.STDOUT, check=True)


def build_flags(environment):
    return {
        "cppflags": environment.get("CPPFLAGS", ""),
        "cflags": "-g -O3 " + environment.get("CFLAGS", ""),
        "cxxflags": "-g -O3 " + environment.get("CXXFLAGS", ""),
        "ldflags": environment.get("LDFLAGS", "-g"),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-root", type=Path, required=True, help="new directory outside the source repository")
    parser.add_argument("--archive", type=Path, help="verified StringTie source archive")
    parser.add_argument("--libdeflate-archive", type=Path)
    parser.add_argument("--fetch", action="store_true", help="allow downloads of the two hash-pinned archives")
    parser.add_argument("--prepare-only", action="store_true", help="extract and patch without compiling")
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--cc", default=os.environ.get("CC", "cc"))
    parser.add_argument("--cxx", default=os.environ.get("CXX", "c++"))
    args = parser.parse_args(argv)
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    root = args.build_root.resolve()
    if root.exists() or root == REPOSITORY or REPOSITORY in root.parents or root in REPOSITORY.parents:
        parser.error("--build-root must be a new directory outside this repository")
    if any(char.isspace() for char in str(root)):
        parser.error("upstream makefiles require a --build-root without whitespace")
    flags = build_flags(os.environ)
    # Make paths survive cwd changes; no upstream build scripts are invoked.
    root.mkdir(parents=True)
    archive = obtain((args.archive or root / "stringtie.tar.gz").resolve(), STRINGTIE, args.fetch)
    stock = extract_verified(archive, root / "stock", STRINGTIE)
    accelerated = extract_verified(archive, root / "accelerated", STRINGTIE)
    apply_accelerator(accelerated)
    metadata = {"stringtie": STRINGTIE, "libdeflate": LIBDEFLATE,
                "patch": "same-start-candidate-index-v1", "platform": platform.platform(),
                "lookup_control": {"variable": "STRINGTIE_INPUT_LOOKUP", "default": "indexed",
                                   "linear": "upstream reverse loop in the same binary",
                                   "unknown": "linear"},
                "index_sha256": sha256(HERE / "same_start_index.h"),
                "builder_sha256": sha256(Path(__file__)),
                "upstream_file_sha256": SOURCE_HASHES,
                "patched_file_sha256": {name: sha256(accelerated / name) for name in SOURCE_HASHES},
                "system_compression_libraries": ["zlib", "bzip2", "lzma"],
                "flags": flags,
                "commands": [], "binaries": {}}
    if not args.prepare_only:
        dependency_archive = obtain((args.libdeflate_archive or root / "libdeflate.tar.gz").resolve(), LIBDEFLATE, args.fetch)
        dependency = extract_verified(dependency_archive, root / "dependency", LIBDEFLATE)
        command = ["make", f"-j{args.jobs}", "libdeflate.a", f"CC={args.cc}",
                   f"CPPFLAGS={flags['cppflags']}", f"CFLAGS={flags['cflags']}"]
        metadata["commands"].append(command)
        run(command, dependency, root / "dependency.log")
        for name, source in [("stock", stock), ("same-start", accelerated)]:
            # This explicit target bypasses upstream build_lib.sh (which fetches dependencies).
            command = ["make", f"-j{args.jobs}", "lib-static", f"CC={args.cc}",
                       f"CPPFLAGS=-I{dependency} {flags['cppflags']}",
                       f"CFLAGS={flags['cflags']} -Wall -fvisibility=hidden -DHAVE_LIBDEFLATE -DHAVE_LIBLZMA",
                       f"LDFLAGS={flags['ldflags']}"]
            metadata["commands"].append(command)
            run(command, source / "htslib", root / f"{name}.log")
            command = ["make", f"-j{args.jobs}", "release", f"CXX={args.cxx}", f"LINKER={args.cxx}",
                       "LIBBZ2=-lbz2", "LIBLZMA=-llzma", f"LIBDEFLATE={dependency / 'libdeflate.a'}"]
            # Upstream appends its required C++ flags to environment CXXFLAGS,
            # but ignores CPPFLAGS. Fold both into its environment, not a make
            # command-line CXXFLAGS override (which would suppress BASEFLAGS).
            metadata["commands"].append(command)
            release_env = dict(os.environ, CXXFLAGS=f"{flags['cxxflags']} {flags['cppflags']}",
                               LDFLAGS=flags["ldflags"])
            with (root / f"{name}.log").open("ab") as stream:
                print("Building:", " ".join(command), flush=True)
                subprocess.run(command, cwd=source, env=release_env, stdout=stream,
                               stderr=subprocess.STDOUT, check=True)
            output = root / f"stringtie-2.2.2-{name}"
            shutil.copy2(source / "stringtie", output)
            version = subprocess.check_output([str(output), "--version"], text=True).strip()
            if version != "2.2.2":
                raise ValueError(f"unexpected binary version: {version}")
            metadata["binaries"][name] = {"path": str(output), "sha256": sha256(output), "version": version,
                                          "patched": name != "stock"}
        metadata["compiler"] = subprocess.check_output(shlex.split(args.cxx) + ["--version"], text=True).splitlines()[0]
    (root / "build-provenance.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata["binaries"], indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, subprocess.CalledProcessError) as error:
        raise SystemExit(f"StringTie optional build failed: {error}; inspect build-root logs")
