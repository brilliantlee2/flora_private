[![English](https://img.shields.io/badge/Language-English-2563eb)](README.md)
[![中文](https://img.shields.io/badge/Language-中文-0f766e)](README_zh-CN.md)

# Flora Private Development Repository

This repository contains the private source, tests, packaging tools, and release documentation for Flora. Do not mirror or push this repository to the public `brilliantlee2/Flora` repository.

The public repository contains user documentation only. Compiled archives are distributed through GitHub Releases.

## Source layout

```text
Cargo.toml / Cargo.lock       Rust package and locked dependencies
src/                          Flora and integrated Glycine Rust source
vendor/                       Patched rust-htslib and edlib_rs
scripts/                      Python analysis, QC, plotting, and report source
run_all.sh                    Private single-species regression baseline
run_all_mixed_species.sh      Private mixed-species regression baseline
environment.yml               Private build and test environment
environment.runtime.yml       Public binary runtime environment
packaging/                    Release packaging tools
tests/                        Rust/Python/release tests
docs/repository-templates/    Public/private GitHub documentation
```

## Create the private build environment

```bash
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export PYTHONUTF8=1

conda env create -f environment.yml
conda activate flora
```

The private environment includes Python 3.11, Rust/Cargo, GCC/G++, Clang/libclang, CMake, samtools, minimap2, bedtools, and all Python dependencies.

Verify it before building:

```bash
python --version
rustc --version
cargo --version
cmake --version
gcc --version
clang --version
find "$CONDA_PREFIX" -name 'libclang.so*' | head
```

If Bindgen cannot locate libclang:

```bash
export LIBCLANG_PATH="$(dirname "$(find "$CONDA_PREFIX" -name 'libclang.so*' -print -quit)")"
```

## Development build and tests

```bash
cargo build --release --locked
cargo test --locked
bash -n run_all.sh run_all_mixed_species.sh packaging/build_binary_release.sh
python -m unittest discover -s tests -v
```

Verify the main CLI:

```bash
./target/release/flora --version
./target/release/flora glycine --help
./target/release/flora analyze --help
bash run_all.sh -h
```

Before producing a single-binary public release, `flora run --help` and
`flora run-mixed --help` must also succeed. The packaging script deliberately
stops when either command is unavailable.

## Version update checklist

1. Update the package version in `Cargo.toml`.
2. Run `cargo check` or `cargo build --release` to update `Cargo.lock`.
3. Update release URLs and version text in the public README templates.
4. Add `docs/repository-templates/public/RELEASE_NOTES_vX.Y.Z.md`.
5. Run all tests before packaging.

Flora follows semantic versioning. Use versions such as `0.1.0`, `0.1.1`, and `0.2.0`.

## Build the Linux x86_64 release

Build on Linux x86_64, preferably the target HPC or an older compatible Linux system. Do not build the public Linux archive on macOS.

```bash
uname -s
uname -m
ldd --version | head -n 1
python --version
```

Expected architecture:

```text
Linux
x86_64
Python 3.11.x
```

Build and package:

```bash
conda activate flora
export LIBCLANG_PATH="$(dirname "$(find "$CONDA_PREFIX" -name 'libclang.so*' -print -quit)")"

cargo test --locked
PYTHON_BIN="$CONDA_PREFIX/bin/python" \
  bash packaging/build_binary_release.sh
```

The archive is generated as:

```text
dist/Flora-<version>-linux-x86_64.tar.gz
```

The packaging script builds only the `flora` binary, copies
`environment.runtime.yml` into the archive as `environment.yml`, uses the
public README templates, compiles an allowlist of Python files to deterministic
Python 3.11 bytecode, and excludes shell workflow sources, standalone stage
binaries, and Rust/Python source.

## Validate the release archive

```bash
ARCHIVE="dist/Flora-0.1.0-linux-x86_64.tar.gz"

tar -tzf "$ARCHIVE" | head -50
tar -tzf "$ARCHIVE" | \
  grep -E '(^|/)(src|vendor|tests)/|Cargo\.(toml|lock)$|\.py$' && {
    echo "ERROR: source file found in public archive" >&2
    exit 1
  } || true
```

Extract and test independently:

```bash
rm -rf /tmp/flora-release-test
mkdir -p /tmp/flora-release-test
tar -xzf "$ARCHIVE" -C /tmp/flora-release-test
cd /tmp/flora-release-test/Flora-0.1.0-linux-x86_64

file flora
ldd flora
./flora --version
./flora glycine --help
./flora run --help
./flora run-mixed --help
cat PYTHON_ABI.txt
cat BUILD_INFO.txt
```

`file flora` should report an ELF x86-64 executable and `stripped`. Do not
publish an archive whose main executable is reported as `not stripped`.

## Generate the checksum

On Linux:

```bash
cd /path/to/Flora
sha256sum dist/Flora-0.1.0-linux-x86_64.tar.gz \
  > dist/Flora-0.1.0-linux-x86_64.tar.gz.sha256
```

Verify it:

```bash
cd dist
sha256sum -c Flora-0.1.0-linux-x86_64.tar.gz.sha256
```

## Publish to the public GitHub repository

Update the public repository README from this private repository:

```bash
cp docs/repository-templates/public/README.md /path/to/Flora-public/README.md
cp docs/repository-templates/public/README_zh-CN.md /path/to/Flora-public/README_zh-CN.md
```

Commit and push only the public documentation. Upload the `.tar.gz` and `.sha256` files as GitHub Release assets; do not commit them to Git history.

With GitHub CLI:

```bash
gh release create v0.1.0 \
  dist/Flora-0.1.0-linux-x86_64.tar.gz \
  dist/Flora-0.1.0-linux-x86_64.tar.gz.sha256 \
  --repo brilliantlee2/Flora \
  --title "Flora v0.1.0" \
  --notes-file docs/repository-templates/public/RELEASE_NOTES_v0.1.0.md
```

## Confidentiality rules

- Keep this repository private.
- Never push `src/`, `scripts/`, `vendor/`, Cargo manifests, tests, or packaging code to the public repository.
- Inspect every archive before upload.
- Python bytecode is obfuscation, not encryption; determined users may still reverse engineer it.
- Build artifacts are platform-specific. Publish separate assets for each supported operating system and architecture.
