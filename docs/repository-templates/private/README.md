[![English](https://img.shields.io/badge/Language-English-2563eb)](README.md)
[![中文](https://img.shields.io/badge/Language-中文-0f766e)](README_zh-CN.md)

# Flora Private Development Repository

This repository contains the private source, tests, packaging tools, and release documentation for Flora. Do not mirror or push this repository to the public `brilliantlee2/Flora` repository.

The public repository contains user documentation only. Compiled archives are distributed through GitHub Releases.

## Clone the private repository

You must have access to the private `brilliantlee2/flora_private` repository.
Clone it with HTTPS:

```bash
git clone https://github.com/brilliantlee2/flora_private.git
cd flora_private
```

Or use SSH after adding an SSH key to your GitHub account:

```bash
git clone git@github.com:brilliantlee2/flora_private.git
cd flora_private
```

Because this is a private repository, GitHub may request browser authorization,
a personal access token, or an authorized SSH key. The repository includes the
patched crates under `vendor/`; no additional Git submodule command is needed.

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

The full workflow accepts multiple paths after one `--fastq`, or a
non-recursive `--fastq-dir`. Parallel Glycine defaults are 10 concurrent jobs
sharing 64 total threads; verify these options in `flora --help` after changes.

Before producing a single-binary public release, `flora run --help` and
`flora run-mixed --help` must also succeed. The packaging script deliberately
stops when either command is unavailable.

## Workflow output policy

Default light-output mode retains the final tagged BAM and index, expression
matrices, cell/barcode mapping, clustering results, QC tables,
`qc/metrics_summary.xlsx`, logs, and the HTML report. After report generation,
Flora removes reproducible read-level tables, the pre-tagging aligned BAM,
assignment TSVs, and other large intermediates. Use `--save-intermediate` for
debugging files or `--full-output` for all upstream FASTQ outputs and
intermediates. These options do not alter analytical results.

Use `--remove-final-bam` to remove every alignment/matrix BAM and index after
the report succeeds, including the final tagged BAM. Expression matrices, QC,
`metrics_summary.xlsx`, and the HTML report remain available.

## Version update checklist

1. Update the package version in `Cargo.toml`.
2. Run `cargo check` or `cargo build --release` to update `Cargo.lock`.
3. Update release URLs and version text in the public README templates.
4. Add `docs/repository-templates/public/RELEASE_NOTES_vX.Y.Z.md`.
5. Run all tests before packaging.

Flora follows semantic versioning. Use versions such as `0.1.0`, `0.1.1`, and `0.2.0`.

## Build the Linux x86_64 release

Build on Linux x86_64, preferably the target HPC. Do not build the public Linux
archive on macOS. The following sequence starts from a fresh clone.

### 1. Enter the repository and check the host

```bash
cd /path/to/flora_private

uname -s
uname -m
ldd --version | head -n 1
```

Expected architecture:

```text
Linux
x86_64
```

### 2. Create the complete build environment

```bash
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export PYTHONUTF8=1

conda env create -f environment.yml
conda activate flora
```

If the environment already exists, update it from the checked-in definition:

```bash
conda env update -n flora -f environment.yml --prune
conda activate flora
```

Verify the required toolchain:

```bash
python --version
rustc --version
cargo --version
cmake --version
gcc --version
clang --version
samtools --version | head -n 2
minimap2 --version
bedtools --version

export LIBCLANG_PATH="$(dirname "$(find "$CONDA_PREFIX" -name 'libclang.so*' -print -quit)")"
test -n "$LIBCLANG_PATH"
```

Python must report `3.11.x`.

### 3. Validate the source tree

```bash
cargo metadata --locked --no-deps --format-version 1 >/dev/null
cargo test --locked
bash -n run_all.sh
bash -n run_all_mixed_species.sh
bash -n packaging/build_binary_release.sh
bash -n packaging/refresh_binary_release_metadata.sh
python -m unittest discover -s tests -v
```

For a development-only build:

```bash
cargo build --release --locked
./target/release/flora --version
./target/release/flora glycine --help
./target/release/flora analyze --help
```

### 4. Build and package the public binary release

The packaging script performs its own locked release build, so a separate
`cargo build` command is not required here:

```bash
PYTHON_BIN="$CONDA_PREFIX/bin/python" \
  bash packaging/build_binary_release.sh
```

To enforce a specific maximum glibc requirement, set it explicitly:

```bash
FLORA_MAX_GLIBC=2.34 \
PYTHON_BIN="$CONDA_PREFIX/bin/python" \
  bash packaging/build_binary_release.sh
```

Do not set a glibc ceiling lower than the symbols provided by the build host and
linked libraries. The detected requirement is recorded in `BUILD_INFO.txt`.

### 5. Locate the generated files

```text
dist/Flora-<version>-linux-x86_64.tar.gz
dist/Flora-<version>-linux-x86_64.tar.gz.sha256
```

The script builds only `flora`, stages it at the archive root, compiles an
allowlist of Python files as deterministic CPython 3.11 bytecode, copies the
runtime environment and public documentation, removes extended attributes, and
generates and verifies SHA256 automatically.

The script deliberately refuses to package until both `flora run --help` and
`flora run-mixed --help` succeed. During the migration period this failure is
expected and prevents publication of an incomplete binary.

## Build a Singularity image

When the target HPC host provides an older glibc than the Flora binary requires,
package the release and its Python/bioinformatics dependencies in a Singularity
SIF. The validated base is Ubuntu 22.04 with glibc 2.35. The build host requires
Linux x86_64, Singularity, and working `--fakeroot`; Docker is needed only to
prepare the base image. The Flora Conda environment does not need to be active.

Create an isolated build directory and copy the inputs:

```bash
mkdir -p singularity_build
cp dist/Flora-0.1.1-linux-x86_64.tar.gz singularity_build/
cp packaging/singularity/Flora.def singularity_build/
cd singularity_build
```

Validate and export the Ubuntu 22.04 base image:

```bash
docker run --rm quay.io/nf-core/ubuntu:22.04 \
  bash -c 'cat /etc/os-release; getconf GNU_LIBC_VERSION'

docker save quay.io/nf-core/ubuntu:22.04 \
  -o flora-base-ubuntu22.04.tar

singularity build --fakeroot \
  flora-base-ubuntu22.04.sif \
  docker-archive://flora-base-ubuntu22.04.tar
```

If Singularity cannot access `/var/run/docker.sock`, do not make the socket
world-writable. Use `docker save` plus `docker-archive://`, or ask the system
administrator to configure Docker group access.

Download the Miniforge installer once from a nearby mirror:

```bash
curl -L --fail --retry 5 --retry-delay 5 \
  -o Miniforge3-Linux-x86_64.sh \
  https://mirror.nju.edu.cn/github-release/conda-forge/miniforge/LatestRelease/Miniforge3-Linux-x86_64.sh
```

The four build inputs must be in the same directory:

```text
Flora.def
Flora-0.1.1-linux-x86_64.tar.gz
Miniforge3-Linux-x86_64.sh
flora-base-ubuntu22.04.sif
```

Build and validate the SIF:

```bash
env -u LD_LIBRARY_PATH \
singularity build --fakeroot \
  Flora-0.1.1-linux-x86_64.sif \
  Flora.def

singularity run --cleanenv \
  Flora-0.1.1-linux-x86_64.sif \
  --version

singularity run --cleanenv \
  Flora-0.1.1-linux-x86_64.sif \
  --help

sha256sum Flora-0.1.1-linux-x86_64.sif \
  > Flora-0.1.1-linux-x86_64.sif.sha256
```

The tested image contains Flora, Python 3.11, samtools, minimap2, bedtools, and
all Python dependencies, and is approximately 820 MiB. `Flora.def` uses
per-user writable cache directories under `/tmp` for Fontconfig, Matplotlib,
and Numba.

## Validate the release archive

```bash
ARCHIVE="dist/Flora-0.1.1-linux-x86_64.tar.gz"

tar -tzf "$ARCHIVE" | head -50
tar -tzf "$ARCHIVE" | \
  grep -E '(^|/)(src|vendor|tests)/|Cargo\.(toml|lock)$|run_all|\.(rs|py|sh)$' && {
    echo "ERROR: source file found in public archive" >&2
    exit 1
  } || true
```

Extract and test independently:

```bash
rm -rf /tmp/flora-release-test
mkdir -p /tmp/flora-release-test
tar -xzf "$ARCHIVE" -C /tmp/flora-release-test
cd /tmp/flora-release-test/Flora-0.1.1-linux-x86_64

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

The executable count must be exactly one:

```bash
find . -type f -perm /111 -print
```

The only result should be `./flora`.

## Generate the checksum

On Linux:

```bash
cd /path/to/Flora
sha256sum dist/Flora-0.1.1-linux-x86_64.tar.gz \
  > dist/Flora-0.1.1-linux-x86_64.tar.gz.sha256
```

Verify it:

```bash
cd dist
sha256sum -c Flora-0.1.1-linux-x86_64.tar.gz.sha256
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
gh release create v0.1.1 \
  dist/Flora-0.1.1-linux-x86_64.tar.gz \
  dist/Flora-0.1.1-linux-x86_64.tar.gz.sha256 \
  --repo brilliantlee2/Flora \
  --title "Flora v0.1.1" \
  --generate-notes
```

Commit only README files, `Flora.def`, licenses, and release notes to Git
history. Upload the binary archive and `.sha256` as GitHub Release assets. The
SIF and its checksum are optional Release assets; omit them when users will
build the image themselves. Never commit the base SIF, Docker tar, Miniforge
installer, `dist/`, or private source code to the public repository.

## Confidentiality rules

- Keep this repository private.
- Never push `src/`, `scripts/`, `vendor/`, Cargo manifests, tests, or packaging code to the public repository.
- Inspect every archive before upload.
- Python bytecode is obfuscation, not encryption; determined users may still reverse engineer it.
- Build artifacts are platform-specific. Publish separate assets for each supported operating system and architecture.
