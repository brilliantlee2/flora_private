#!/usr/bin/env python3
"""Provision the pinned CPU Bonefish and Flora's indexed StringTie."""

import argparse
import hashlib
from http.client import IncompleteRead
import json
import os
from pathlib import Path
import platform
import shlex
import shutil
import ssl
import subprocess
import sys
import tarfile
import tempfile
import time
from urllib.error import URLError
from urllib.parse import urlsplit
from urllib.request import HTTPRedirectHandler, build_opener


ASSET_ROOT = Path(__file__).resolve().parents[1]
BONEFISH = "bonefish"
STRINGTIE = "stringtie-2.2.2-same-start"
RECEIPT = "installed-tools.json"


def resolve_tools_dir(explicit):
    if explicit is not None:
        return explicit.expanduser().resolve()
    if override := os.environ.get("FLORA_TOOLS_DIR"):
        return Path(override).expanduser().resolve()
    data = Path(os.environ.get("XDG_DATA_HOME") or Path.home() / ".local" / "share")
    return (data / "flora" / "tools").expanduser().resolve()


def platform_key():
    machine = platform.machine().lower()
    machine = {"amd64": "x86_64", "aarch64": "arm64"}.get(machine, machine)
    system = platform.system().lower()
    return f"{'macos' if system == 'darwin' else system}-{machine}"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class HTTPSRedirectHandler(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        if urlsplit(newurl).scheme != "https":
            raise ValueError(f"Refusing non-HTTPS redirect for {req.full_url}")
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def download_archive(url, path):
    if urlsplit(url).scheme != "https":
        raise ValueError(f"Refusing non-HTTPS archive URL for {path.name}")
    opener = build_opener(HTTPSRedirectHandler())
    print(f"Downloading {path.name}", flush=True)
    for attempt in range(3):
        try:
            with opener.open(url, timeout=60) as source, path.open("xb") as output:
                shutil.copyfileobj(source, output)
            return
        except (URLError, IncompleteRead, ConnectionError, TimeoutError, ssl.SSLEOFError) as error:
            path.unlink(missing_ok=True)
            if attempt == 2:
                raise ValueError(
                    f"Could not download {path.name} from {url}: {error}. "
                    "Use --archive-dir if the archives are available locally."
                ) from error
            print(f"Retrying {path.name}: {error}", file=sys.stderr, flush=True)
            time.sleep(2 ** attempt)


def obtain_archive(spec, archive_dir, work):
    name = spec["archive"]
    path = work / name
    if archive_dir is None:
        download_archive(spec["url"], path)
    else:
        cached = archive_dir / name
        if not cached.is_file():
            raise ValueError(f"Missing offline archive: {cached}")
        shutil.copyfile(cached, path)
    if sha256(path) != spec["sha256"]:
        raise ValueError(f"SHA-256 mismatch for {path}; expected {spec['sha256']}")
    return path.resolve()


def extract_bonefish(archive, name, output):
    with tarfile.open(archive, "r:gz") as tar:
        matches = [member for member in tar.getmembers() if member.name == name]
        if len(matches) != 1 or not matches[0].isfile():
            raise ValueError(f"Expected one regular Bonefish executable member at {name}")
        with tar.extractfile(matches[0]) as source, output.open("xb") as dest:
            shutil.copyfileobj(source, dest)
    output.chmod(0o755)


def build_environment():
    env = dict(os.environ)
    prefixes = []
    if conda_prefix := env.get("CONDA_PREFIX"):
        prefixes.append(Path(conda_prefix))
    if platform.system() == "Darwin":
        homebrew_xz = Path(env.get("HOMEBREW_PREFIX") or "/opt/homebrew") / "opt" / "xz"
        if homebrew_xz.is_dir():
            prefixes.append(homebrew_xz)
    includes, libraries = [], []
    for prefix in prefixes:
        if (prefix / "include").is_dir():
            includes.append(shlex.quote(f"-I{prefix / 'include'}"))
        if (prefix / "lib").is_dir():
            libraries.extend(shlex.quote(flag) for flag in (f"-L{prefix / 'lib'}", f"-Wl,-rpath,{prefix / 'lib'}"))
    for name, args in (("CPPFLAGS", includes), ("LDFLAGS", libraries)):
        if args:
            env[name] = " ".join([env.get(name, ""), *args]).strip()
    return env


def build_stringtie(asset_root, source, libdeflate, build, jobs, output):
    builder_dir = asset_root / "packaging" / "stringtie_input"
    builder = next((p for p in (builder_dir / "build.py", builder_dir / "build.pyc") if p.is_file()), None)
    if builder is None or not (builder_dir / "same_start_index.h").is_file():
        raise ValueError(f"Missing Flora StringTie build assets in {builder_dir}")
    env = build_environment()
    for tool in ("make", env.get("CC") or "cc", env.get("CXX") or "c++"):
        if shutil.which(shlex.split(tool)[0]) is None:
            raise ValueError(f"Missing build tool: {tool}. Install make and a C/C++ compiler, then rerun setup-tools.")
    command = [sys.executable, str(builder), "--build-root", str(build),
               "--archive", str(source), "--libdeflate-archive", str(libdeflate), "--jobs", str(jobs)]
    result = subprocess.run(command, env=env, check=False)
    if result.returncode:
        for path in sorted(build.glob("*.log")):
            print(f"{path.name} tail:", file=sys.stderr)
            print("\n".join(path.read_text(errors="replace").splitlines()[-12:]), file=sys.stderr)
        raise ValueError("StringTie build failed. Install the zlib, bzip2, and xz (liblzma) development libraries "
                         "or activate the Flora Conda environment and rerun setup-tools.")
    provenance = json.loads((build / "build-provenance.json").read_text())
    shutil.copy2(build / STRINGTIE, output)
    return {key: provenance[key] for key in ("patch", "compiler", "index_sha256", "patched_file_sha256", "flags")}


def executable_version(path):
    result = subprocess.run([str(path), "--version"], capture_output=True, text=True, check=False)
    output = (result.stdout + result.stderr).strip()
    if result.returncode:
        raise ValueError(f"{path.name} could not run (exit {result.returncode}): {output}")
    return output


def exists(path):
    return path.exists() or path.is_symlink()


def verify_existing(receipt, bin_dir, installation, names):
    if not receipt.is_file() or receipt.is_symlink():
        raise ValueError(f"Invalid setup-tools receipt: {receipt}")
    installed = json.loads(receipt.read_text())
    if installed.get("installation") != installation:
        raise ValueError(f"Existing setup-tools receipt does not match these pins: {receipt}. Use a new --tools-dir.")
    for name in names:
        tool = bin_dir / name
        recorded = installed.get("binaries", {}).get(name, {})
        if (tool.is_symlink() or not tool.is_file() or not os.access(tool, os.X_OK)
                or sha256(tool) != recorded.get("sha256")):
            raise ValueError(f"Installed executable integrity check failed: {tool}. Use a new --tools-dir.")


def publish(stage, bin_dir, receipt, names):
    bin_dir.mkdir(exist_ok=True)
    paths = [(stage / name, bin_dir / name) for name in names] + [(stage / RECEIPT, receipt)]
    try:
        for source, destination in paths:
            # Linking fails atomically if the user created a destination during the build.
            os.link(source, destination)
    except BaseException:
        for source, destination in reversed(paths):
            if not destination.is_symlink() and destination.exists() and source.samefile(destination):
                destination.unlink()
        raise


def install(tools, archive_dir=None, jobs=2, asset_root=ASSET_ROOT, target=None):
    config = json.loads((asset_root / "runtime-tools.json").read_text())
    target = target or platform_key()
    bonefish = config["bonefish"]["releases"].get(target)
    if bonefish is None:
        supported = ", ".join(sorted(config["bonefish"]["releases"]))
        raise ValueError(f"No pinned Bonefish release for {target}. Supported platforms: {supported}.")
    spec = config["stringtie"]
    names = (BONEFISH, STRINGTIE)
    installation = {"platform": target, "bonefish_archive_sha256": bonefish["sha256"],
                    "stringtie_source_sha256": spec["source"]["sha256"],
                    "libdeflate_source_sha256": spec["libdeflate"]["sha256"], "stringtie_patch": spec["patch"],
                    "stringtie_index_sha256": spec["index_sha256"]}
    tools.mkdir(parents=True, exist_ok=True)
    lock = tools / ".setup-tools.lock"
    try:
        lock.mkdir()
    except FileExistsError as error:
        raise ValueError(f"Setup lock exists at {lock}. Wait for setup to finish, or remove the empty lock if it stopped.") from error

    bin_dir, receipt = tools / "bin", tools / RECEIPT
    try:
        if bin_dir.is_symlink() or (bin_dir.exists() and not bin_dir.is_dir()):
            raise ValueError(f"Expected a regular installation bin directory: {bin_dir}")
        if exists(receipt):
            verify_existing(receipt, bin_dir, installation, names)
            return False
        for name in names:
            if exists(bin_dir / name):
                raise ValueError(f"{bin_dir / name} already exists outside setup-tools management. Use a new --tools-dir.")
        # The StringTie makefiles require a whitespace-free build path.
        with tempfile.TemporaryDirectory(prefix="flora-build-tools-", dir="/tmp") as scratch, \
                tempfile.TemporaryDirectory(prefix=".stage-", dir=tools) as staged:
            work, stage = Path(scratch), Path(staged)
            archives = [obtain_archive(archive, archive_dir, work)
                        for archive in (bonefish, spec["source"], spec["libdeflate"])]
            extract_bonefish(archives[0], bonefish["member"], stage / BONEFISH)
            build = build_stringtie(asset_root, archives[1], archives[2], work / "stringtie", jobs, stage / STRINGTIE)
            if (build["patch"] != spec["patch"] or build["index_sha256"] != spec["index_sha256"]
                    or b"STRINGTIE_INPUT_LOOKUP" not in (stage / STRINGTIE).read_bytes()):
                raise ValueError("Built StringTie is missing the pinned same-start patch")
            versions = dict(zip(names, (config["bonefish"]["version"], spec["version"])))
            for name, version in versions.items():
                if executable_version(stage / name) != version:
                    raise ValueError(f"Unexpected executable version for {name}; expected {version}")
            data = {"installation": installation, "stringtie_build": build,
                    "binaries": {name: {"version": versions[name], "sha256": sha256(stage / name)} for name in names}}
            (stage / RECEIPT).write_text(json.dumps(data, sort_keys=True, indent=2) + "\n")
            publish(stage, bin_dir, receipt, names)
            return True
    finally:
        lock.rmdir()


def main(argv=None):
    prog = "flora setup-tools" if Path(__file__).suffix == ".pyc" else None
    parser = argparse.ArgumentParser(prog=prog, description=__doc__)
    parser.add_argument("--tools-dir", type=Path, help="Override FLORA_TOOLS_DIR and the per-user installation path")
    parser.add_argument("--archive-dir", type=Path, help="Offline directory with the three pinned source and binary archives")
    parser.add_argument("--jobs", type=int, default=2, help="Parallel StringTie build jobs (default: 2)")
    args = parser.parse_args(argv)
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    try:
        tools = resolve_tools_dir(args.tools_dir)
        changed = install(tools, args.archive_dir, args.jobs)
    except (OSError, ValueError, tarfile.TarError, subprocess.SubprocessError) as error:
        parser.exit(1, f"setup-tools: {error}\n")
    print(f"{'Installed' if changed else 'Verified existing'} Flora tools in {tools / 'bin'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
