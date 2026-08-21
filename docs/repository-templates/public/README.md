[![English](https://img.shields.io/badge/Language-English-2563eb)](README.md)
[![中文](https://img.shields.io/badge/Language-中文-0f766e)](README_zh-CN.md)

# Flora

Flora is an end-to-end workflow for full-length single-cell RNA sequencing. It is optimized for MGI C4 single-cell libraries sequenced with the Cyclone long-read platform.

Flora performs:

- full-length cDNA identification, orientation, trimming, and chimeric-read rescue;
- dual-end barcode extraction, correction, merging, and cell assignment;
- UMI extraction and directional clustering;
- gene and isoform expression matrix generation;
- RNA QC and saturation analysis;
- Scanpy RNA clustering and UMAP visualization;
- generation of a self-contained HTML report.

[Glycine](https://github.com/CycloneSEQ-Bioinformatics/Glycine) is integrated into Flora. Users do not need to install Glycine separately or provide `--glycine-bin-dir`.

## Platform

For native execution, the current binary release requires:

- Linux x86_64;
- glibc 2.35 or newer (the release is built on Ubuntu 22.04/glibc 2.35);
- Conda, Miniforge, Mambaforge, or Micromamba;
- Python 3.11, installed by the bundled environment file.

The release is not compatible with macOS, ARM Linux, or Windows.

Check the system architecture and glibc version before downloading:

```bash
uname -m
ldd --version | head -n 1
```

The expected architecture is `x86_64`. Older HPC operating systems, including
systems with glibc 2.17, cannot run this binary release natively; use the
Singularity procedure below instead.

## Download

Open [GitHub Releases](https://github.com/brilliantlee2/Flora/releases) and
download these two files from the latest release:

```text
Flora-<version>-linux-x86_64.tar.gz
Flora-<version>-linux-x86_64.tar.gz.sha256
```

Use the accompanying `.sha256` file to verify the download, then extract the
archive. Replace `<version>` with the version that was downloaded:

```bash
sha256sum -c Flora-<version>-linux-x86_64.tar.gz.sha256
tar -xzf Flora-<version>-linux-x86_64.tar.gz
cd Flora-<version>-linux-x86_64
```

## Install the runtime environment

```bash
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export PYTHONUTF8=1

conda env create -f environment.yml
conda activate flora
```

Verify the installation:

```bash
python --version
samtools --version
minimap2 --version
bedtools --version

./flora --version
./flora --help
./flora mixed --help
./flora glycine --help
```

The embedded Python analysis modules require Python 3.11. Do not replace the environment with Python 3.10, 3.12, 3.13, or 3.14.

For reproducible runs, always provide the 10 bp barcode whitelist explicitly with `--barcode-list-10bp /path/to/BC_1536.txt`. Automatic discovery of a whitelist placed next to the executable is not guaranteed in the current release.

## Use Singularity

For older HPC systems such as glibc 2.17 hosts, build a Singularity SIF on a
separate Linux x86_64 build machine and copy the resulting single SIF to the
cluster. The container uses Ubuntu 22.04 with glibc 2.35, so Flora no longer
depends directly on the host's older glibc. The host must still be able to run
Singularity, and the build machine must support `--fakeroot`; running the final
SIF normally does not require root.

The SIF includes Flora, Python 3.11, samtools, minimap2, bedtools, and the Python
dependencies. Do not run `conda env create` or `conda activate flora` when using
the SIF.

Create a build directory and place the Flora archive downloaded from the GitHub
Release in it:

```bash
mkdir -p flora_singularity_build
cd flora_singularity_build

wget https://raw.githubusercontent.com/brilliantlee2/Flora/main/Flora.def
```

Prepare the validated Ubuntu 22.04/glibc 2.35 base image:

```bash
docker run --rm quay.io/nf-core/ubuntu:22.04 \
  bash -c 'cat /etc/os-release; getconf GNU_LIBC_VERSION'

docker save quay.io/nf-core/ubuntu:22.04 \
  -o flora-base-ubuntu22.04.tar

singularity build --fakeroot \
  flora-base-ubuntu22.04.sif \
  docker-archive://flora-base-ubuntu22.04.tar
```

If the build host cannot access the Docker socket, run `docker save` on a
machine with Docker access and copy the tar. Do not make
`/var/run/docker.sock` world-writable.

Download the Miniforge installer:

```bash
curl -L --fail --retry 5 --retry-delay 5 \
  -o Miniforge3-Linux-x86_64.sh \
  https://mirror.nju.edu.cn/github-release/conda-forge/miniforge/LatestRelease/Miniforge3-Linux-x86_64.sh
```

The build directory must contain these four files. The Flora archive name must
match the version referenced by `Flora.def`:

```text
Flora.def
Flora-0.1.1-linux-x86_64.tar.gz
Miniforge3-Linux-x86_64.sh
flora-base-ubuntu22.04.sif
```

Build and validate the image:

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
```

Bind the directories containing inputs, references, and outputs for a real run:

```bash
singularity run --cleanenv \
  --bind /data:/data \
  Flora-0.1.1-linux-x86_64.sif \
  --fastq /data/sample.fastq.gz \
  --barcode-list-10bp /data/BC_1536.txt \
  --ref-dir /data/GRCh38_flora \
  --out-dir /data/sample_output \
  --sample-id sample
```

The tested SIF is approximately 820 MiB. If the cluster provides Apptainer,
replace `singularity` with `apptainer` in the commands above.

## Prepare a reference

`--ref-dir` must contain:

```text
reference/
├── genome.fa
├── genes.gtf
├── genes.bed
└── chrom_sizes.tsv
```

Generate the auxiliary files with:

```bash
samtools faidx genome.fa
paftools.js gff2bed -j genes.gtf > genes.bed
cut -f1,2 genome.fa.fai | sort -V > chrom_sizes.tsv
```

`paftools.js` is distributed with minimap2.

Chromosome or contig names must match across `genome.fa`, `genes.gtf`, `genes.bed`, and `chrom_sizes.tsv`. For example, do not mix `chr1` with `1` between files. An explicit `--isoform-gtf` can be supplied when the isoform annotation differs from `genes.gtf`; otherwise Flora reuses `genes.gtf`.

## Analyze raw FASTQ

The integrated Glycine stage runs automatically:

```bash
./flora \
  --fastq /data/sample.fastq.gz \
  --barcode-list-10bp /data/BC_1536.txt \
  --ref-dir /data/GRCh38_flora \
  --out-dir ./sample_output \
  --sample-id sample
```

## Analyze an existing full-length FASTQ

```bash
./flora \
  --skip-glycine \
  --full-length-fastq /data/sample.full-length-plus-rescued.fq.gz \
  --barcode-list-10bp /data/BC_1536.txt \
  --ref-dir /data/GRCh38_flora \
  --out-dir ./sample_output \
  --sample-id sample
```

## Mixed-species analysis

```bash
./flora mixed \
  --skip-glycine \
  --full-length-fastq /data/mixed.full-length-plus-rescued.fq.gz \
  --barcode-list-10bp /data/BC_1536.txt \
  --ref-dir /data/merged_reference \
  --out-dir ./mixed_output \
  --sample-id mixed_sample
```

Mixed-species runs additionally generate `qc/barnyard_qc/barnyard_summary.tsv`, `barnyard_per_cell.tsv`, and a Barnyard QC section in the HTML report.

## Measured resource usage

A monitored memory-optimized `flora mixed --skip-glycine` run used a 149 GB compressed full-length FASTQ containing 138,615,368 reads. With 32 workflow threads and 16 cluster threads, the process tree completed in approximately 884 minutes (14 h 44 min). The observed peak PSS was 94.964 GiB and the observed peak summed RSS was 94.973 GiB.

## Outputs

The output directory contains:

```text
upstream/    Barcode correction, cell assignment, and knee plots
alignment/   Aligned and tagged BAM files
matrix/      Gene/isoform matrices and RNA clustering coordinates
qc/          RNA QC, saturation, and report inputs
logs/        Per-stage logs
```

Key outputs include `read_assigned_cell.csv`, `barcode_to_cell.csv`, tagged BAM files, gene and isoform expression matrices, Scanpy UMAP coordinates, and `<sample>.single_cell_report.html`. Mixed-species analysis also reports the per-cell human/mouse UMI classification under `qc/barnyard_qc/`.

## Support

Please report reproducible problems through [GitHub Issues](https://github.com/brilliantlee2/Flora/issues). Include the Flora version, command, operating system, input size, requested resources, and the relevant stage log. Do not upload private FASTQ/BAM data.

## License notices

No project-wide license has been declared for Flora. The integrated Glycine component retains its MIT license and attribution; third-party notices are included in every release archive.

## GitHub release files

Keep only the README files, `Flora.def`, release notes, and license files in the
Git repository. Upload at least these assets to each formal GitHub Release:

```text
Flora-<version>-linux-x86_64.tar.gz
Flora-<version>-linux-x86_64.tar.gz.sha256
```

A prebuilt SIF and its `.sha256` may be added as optional Release assets. Do not
upload the base SIF, Docker tar, Miniforge installer, Conda environment, or any
private source code to the public repository.
