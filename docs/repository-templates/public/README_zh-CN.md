[![English](https://img.shields.io/badge/Language-English-2563eb)](README.md)
[![中文](https://img.shields.io/badge/语言-中文-0f766e)](README_zh-CN.md)

# Flora

Flora 是一个面向三代全长单细胞 RNA 测序的端到端分析流程，针对 MGI C4 单细胞建库和 Cyclone 长读长测序进行了优化。

Flora 包括：

- 全长 cDNA 识别、定向、修剪和嵌合 read 拯救；
- 双端 barcode 提取、矫正、合并和 cell assignment；
- UMI 提取与 directional clustering；
- 基因和转录本表达矩阵生成；
- RNA QC 和饱和度分析；
- Scanpy RNA 聚类和 UMAP 可视化；
- 自包含 HTML 报告。

[Glycine](https://github.com/CycloneSEQ-Bioinformatics/Glycine) 已集成到 Flora 中，用户不需要单独安装 Glycine，也不需要传入 `--glycine-bin-dir`。

## 适用平台

当前二进制发行包在原生运行方式下需要：

- Linux x86_64；
- glibc 2.35 或更高版本（发行包在 Ubuntu 22.04/glibc 2.35 上构建）；
- Conda、Miniforge、Mambaforge 或 Micromamba；
- 由环境文件安装的 Python 3.11。

该发行包不支持 macOS、ARM Linux 或 Windows。

下载前请检查系统架构和 glibc 版本：

```bash
uname -m
ldd --version | head -n 1
```

预期架构为 `x86_64`。使用 glibc 2.17 等较旧系统的 HPC 节点无法直接原生运行
当前二进制发行包，可以使用下文的 Singularity 方案。

## 下载

请前往 [GitHub Releases](https://github.com/brilliantlee2/Flora/releases) 下载
最新版本的以下两个文件：

```text
Flora-<version>-linux-x86_64.tar.gz
Flora-<version>-linux-x86_64.tar.gz.sha256
```

下载后使用对应的 `.sha256` 文件验证完整性，再解压压缩包。命令中的
`<version>` 替换为实际下载的版本号：

```bash
sha256sum -c Flora-<version>-linux-x86_64.tar.gz.sha256
tar -xzf Flora-<version>-linux-x86_64.tar.gz
cd Flora-<version>-linux-x86_64
```

## 安装运行环境

```bash
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export PYTHONUTF8=1

conda env create -f environment.yml
conda activate flora
```

验证：

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

内嵌的 Python 分析模块需要 Python 3.11，不要替换为 Python 3.10、3.12、3.13 或 3.14。

为保证运行可复现，建议始终使用 `--barcode-list-10bp /path/to/BC_1536.txt` 显式传入 10 bp barcode whitelist。当前版本不能保证自动发现放在可执行文件旁边的 whitelist。

## 使用 Singularity

对于 glibc 2.17 等较旧的 HPC，可以在另一台 Linux x86_64 构建机上将 Flora
封装成 Singularity SIF，再把单个 SIF 复制到计算集群。镜像内部使用 Ubuntu
22.04/glibc 2.35，因此 Flora 不再直接依赖计算节点的旧 glibc。宿主机仍需能够
运行 Singularity，并支持 `--fakeroot` 构建；SIF 运行本身通常不需要 root。

SIF 已包含 Flora、Python 3.11、samtools、minimap2、bedtools 和 Python 依赖，
使用 SIF 时不需要再执行 `conda env create` 或 `conda activate flora`。

在构建机上准备目录，并放入从 GitHub Release 下载的 Flora 压缩包：

```bash
mkdir -p flora_singularity_build
cd flora_singularity_build

wget https://raw.githubusercontent.com/brilliantlee2/Flora/main/Flora.def
```

准备经验证的 Ubuntu 22.04/glibc 2.35 基础镜像：

```bash
docker run --rm quay.io/nf-core/ubuntu:22.04 \
  bash -c 'cat /etc/os-release; getconf GNU_LIBC_VERSION'

docker save quay.io/nf-core/ubuntu:22.04 \
  -o flora-base-ubuntu22.04.tar

singularity build --fakeroot \
  flora-base-ubuntu22.04.sif \
  docker-archive://flora-base-ubuntu22.04.tar
```

如果构建机不能访问 Docker socket，可以在有 Docker 权限的机器执行
`docker save`，再复制 tar；不要把 `/var/run/docker.sock` 修改为全局可写。

下载 Miniforge 安装器：

```bash
curl -L --fail --retry 5 --retry-delay 5 \
  -o Miniforge3-Linux-x86_64.sh \
  https://mirror.nju.edu.cn/github-release/conda-forge/miniforge/LatestRelease/Miniforge3-Linux-x86_64.sh
```

构建目录中需要以下四个文件，其中 Flora 压缩包名称必须与 `Flora.def` 中的版本
一致：

```text
Flora.def
Flora-0.1.1-linux-x86_64.tar.gz
Miniforge3-Linux-x86_64.sh
flora-base-ubuntu22.04.sif
```

构建并验证：

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

真实数据运行时显式绑定输入、参考和输出所在目录，例如：

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

构建成功的实测 SIF 大约为 820 MiB。若集群管理员已提供 Apptainer，可将上述
`singularity` 命令替换为 `apptainer`。

## 准备参考基因组

`--ref-dir` 需要包含：

```text
reference/
├── genome.fa
├── genes.gtf
├── genes.bed
└── chrom_sizes.tsv
```

生成辅助文件：

```bash
samtools faidx genome.fa
paftools.js gff2bed -j genes.gtf > genes.bed
cut -f1,2 genome.fa.fai | sort -V > chrom_sizes.tsv
```

`paftools.js` 由 minimap2 提供。

`genome.fa`、`genes.gtf`、`genes.bed` 和 `chrom_sizes.tsv` 中的染色体或 contig 名称必须一致，例如不能在不同文件中混用 `chr1` 和 `1`。如果转录本注释与 `genes.gtf` 不同，可以显式传入 `--isoform-gtf`；否则 Flora 会复用 `genes.gtf`。

## 分析原始 FASTQ

内置 Glycine 会自动运行：

```bash
./flora \
  --fastq /data/sample.fastq.gz \
  --barcode-list-10bp /data/BC_1536.txt \
  --ref-dir /data/GRCh38_flora \
  --out-dir ./sample_output \
  --sample-id sample
```

## 分析已有全长 FASTQ

```bash
./flora \
  --skip-glycine \
  --full-length-fastq /data/sample.full-length-plus-rescued.fq.gz \
  --barcode-list-10bp /data/BC_1536.txt \
  --ref-dir /data/GRCh38_flora \
  --out-dir ./sample_output \
  --sample-id sample
```

## Mixed-species 分析

```bash
./flora mixed \
  --skip-glycine \
  --full-length-fastq /data/mixed.full-length-plus-rescued.fq.gz \
  --barcode-list-10bp /data/BC_1536.txt \
  --ref-dir /data/merged_reference \
  --out-dir ./mixed_output \
  --sample-id mixed_sample
```

Mixed-species 模式还会生成 `qc/barnyard_qc/barnyard_summary.tsv`、`barnyard_per_cell.tsv`，并在 HTML 报告中增加 Barnyard QC 部分。

## 实测资源消耗

一次受监控的内存优化版 `flora mixed --skip-glycine` 运行使用了 149 GB 的压缩全长 FASTQ，共包含 138,615,368 条 reads。在流程线程为 32、cluster 线程为 16 时，进程树约运行 884 分钟（14 小时 44 分）。实际观测到的 PSS 峰值为 94.964 GiB，RSS 求和峰值为 94.973 GiB。

## 主要输出

```text
upstream/    barcode 矫正、cell assignment 和 knee plot
alignment/   比对与标签 BAM
matrix/      基因/转录本矩阵和 RNA 聚类坐标
qc/          RNA QC、饱和度和报告输入
logs/        各步骤日志
```

关键结果包括 `read_assigned_cell.csv`、`barcode_to_cell.csv`、带标签 BAM、基因和转录本表达矩阵、Scanpy UMAP 坐标和 `<sample>.single_cell_report.html`。Mixed-species 分析还会在 `qc/barnyard_qc/` 中输出逐细胞的人/鼠 UMI 分类结果。

## 问题反馈

请通过 [GitHub Issues](https://github.com/brilliantlee2/Flora/issues) 提交可复现的问题，并附上 Flora 版本、运行命令、操作系统、输入大小、资源申请和相关日志。请勿上传私有 FASTQ/BAM 数据。

## 许可声明

Flora 尚未声明项目整体许可证。内置 Glycine 保留 MIT 许可和作者声明，每个发行包中均包含第三方许可说明。

## GitHub 发布文件

Git 仓库中仅保存 README、`Flora.def`、release notes 和许可文件。每个正式
GitHub Release 至少上传：

```text
Flora-<version>-linux-x86_64.tar.gz
Flora-<version>-linux-x86_64.tar.gz.sha256
```

预构建 SIF 和对应 `.sha256` 可以作为可选 Release 附件。基础 SIF、Docker tar、
Miniforge 安装器、Conda 环境目录及任何私有源码不应上传到公开仓库。
