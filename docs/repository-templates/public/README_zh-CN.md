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
- glibc 达到发行包 `BUILD_INFO.txt` 中的 `required_glibc` 要求；
- Conda、Miniforge、Mambaforge 或 Micromamba；
- 由环境文件安装的 Python 3.11。

该发行包不支持 macOS、ARM Linux 或 Windows。

下载前请检查系统架构和 glibc 版本：

```bash
uname -m
ldd --version | head -n 1
```

预期架构为 `x86_64`。如果宿主机的 glibc 低于发行包要求，可用下文的 Singularity 方案。

## 下载

请前往 [GitHub Releases](https://github.com/brilliantlee2/Flora/releases) 下载
最新版本的以下两个文件：

```text
Flora-<version>-linux-x86_64.tar.gz
Flora-<version>-linux-x86_64.tar.gz.sha256
```

将 `version` 设为实际下载的版本号，校验压缩包后再解压。

```bash
version=0.1.2
archive="Flora-${version}-linux-x86_64.tar.gz"
sha256sum -c "${archive}.sha256" &&
  tar -xzf "$archive" &&
  cd "Flora-${version}-linux-x86_64"
```

## 安装运行环境

```bash
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export PYTHONUTF8=1

conda env create -f environment.yml
conda activate flora
./flora setup-tools
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

`./flora setup-tools` 自动下载 Bonefish CPU 0.5.0 并构建 StringTie 2.2.2
same-start 补丁版。安装时按随包提供的 `runtime-tools.json` 校验各归档的
SHA-256。Flora 会自动在 `${XDG_DATA_HOME:-$HOME/.local/share}/flora/tools`
下找到这两个程序。在安装和运行前设置 `FLORA_TOOLS_DIR` 可选用其他目录。

Flora 默认用 Bonefish 做基因组比对，可用 `--genome-aligner minimap2` 切换为
minimap2。ONT 特征指派的转录组重比对始终使用 minimap2。已有二进制程序可以用
`--bonefish-bin` 或 `--stringtie-bin` 指定。

`--bonefish-seam-retry` 是默认关闭的实验选项，启用时须用 `--bonefish-bin`
显式指定打过补丁的 Bonefish 0.5.0。`setup-tools` 安装的是标准版。

流程的大型临时文件默认放在 `<out-dir>/.flora-scratch`。用
`--ont-temp-dir /path/to/scratch` 可指定单独的工作盘。ONT 特征指派的磁盘剩余
空间和存储错误可在 `logs/04_gene.log` 查看。

内嵌的 Python 分析模块需要 Python 3.11。如果 Flora 报
`Bad magic number in .pyc file`，重新激活 Flora 环境，并确认
`python3 --version` 显示 3.11.x。

为保证运行可复现，建议始终使用 `--barcode-list-10bp /path/to/BC_1536.txt` 显式传入 10 bp barcode whitelist。当前版本不能保证自动发现放在可执行文件旁边的 whitelist。

## 使用 Singularity

对于 glibc 2.17 等较旧的 HPC，可以在另一台 Linux x86_64 构建机上将 Flora
封装成 Singularity SIF，再把单个 SIF 复制到计算集群。镜像内部使用 Ubuntu
22.04/glibc 2.35，因此 Flora 不再直接依赖计算节点的旧 glibc。宿主机仍需能够
运行 Singularity，并支持 `--fakeroot` 构建；SIF 运行本身通常不需要 root。

构建机需要联网下载 Ubuntu 软件包、Conda 依赖及工具压缩包。配套定义文件调用
`flora setup-tools`，将 Bonefish CPU 0.5.0 和 StringTie 2.2.2 same-start 版
安装到镜像内的 `/opt/flora/tools`。SIF 内置流程使用的 Conda 环境和工具，默认
使用镜像内的 Bonefish 做基因组比对。传入 `--genome-aligner minimap2` 可选用
镜像内的 minimap2。

### 1. 创建构建目录

```bash
mkdir -p flora_singularity_build
cd flora_singularity_build
```

后续命令都在这个目录中执行。

### 2. 下载 Flora 发行包

打开 [GitHub Releases](https://github.com/brilliantlee2/Flora/releases)，下载与
`Flora.def` 版本一致的压缩包和校验文件，并保存到 `flora_singularity_build`
目录。以下命令使用与配套 `Flora.def` 一致的 0.1.2 版本。

```text
Flora-0.1.2-linux-x86_64.tar.gz
Flora-0.1.2-linux-x86_64.tar.gz.sha256
```

进入构建目录后验证压缩包：

```bash
sha256sum -c Flora-0.1.2-linux-x86_64.tar.gz.sha256
```

构建 SIF 时不需要解压这个压缩包，`Flora.def` 会自动将其复制并解压到镜像中。

### 3. 下载 Flora.def

仍然在 `flora_singularity_build` 目录执行：

```bash
wget https://raw.githubusercontent.com/brilliantlee2/Flora/main/Flora.def
```

构建前核对 `Flora.def` 的 `%files` 段，确认其中的压缩包文件名与本地文件一致。
使用旧版发行包时，请使用该版本配套的定义文件。

### 4. 准备 Ubuntu 22.04 基础 SIF

以下命令从 Docker 镜像生成经验证的 Ubuntu 22.04/glibc 2.35 基础 SIF：

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
`docker save`，再把 `flora-base-ubuntu22.04.tar` 复制到当前构建目录；不要把
`/var/run/docker.sock` 修改为全局可写。如果已经取得
`flora-base-ubuntu22.04.sif`，直接放入构建目录即可跳过本步骤。

### 5. 下载 Miniforge 安装器

```bash
curl -L --fail --retry 5 --retry-delay 5 \
  -o Miniforge3-Linux-x86_64.sh \
  https://mirror.nju.edu.cn/github-release/conda-forge/miniforge/LatestRelease/Miniforge3-Linux-x86_64.sh
```

### 6. 检查构建输入

此时 `flora_singularity_build` 目录至少应包含以下文件：

```text
Flora.def
Flora-0.1.2-linux-x86_64.tar.gz
Flora-0.1.2-linux-x86_64.tar.gz.sha256
Miniforge3-Linux-x86_64.sh
flora-base-ubuntu22.04.sif
```

其中 `.sha256` 仅用于校验；真正构建 SIF 使用另外四个文件。

### 7. 构建并验证最终 SIF

```bash
env -u LD_LIBRARY_PATH \
singularity build --fakeroot \
  Flora-0.1.2-linux-x86_64.sif \
  Flora.def

singularity test --cleanenv Flora-0.1.2-linux-x86_64.sif
singularity run --cleanenv Flora-0.1.2-linux-x86_64.sif --version
```

真实数据运行时显式绑定输入、参考和输出所在目录，例如：

```bash
singularity run --cleanenv \
  --bind /data:/data \
  Flora-0.1.2-linux-x86_64.sif \
  --fastq /data/sample.fastq.gz \
  --barcode-list-10bp /data/BC_1536.txt \
  --ref-dir /data/GRCh38_flora \
  --out-dir /data/sample_output \
  --sample-id sample
```

若集群管理员已提供 Apptainer，可将上述 `singularity` 命令替换为 `apptainer`。

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

`paftools.js` 随 minimap2 分发，通过 `k8` JavaScript 解释器运行。

`genome.fa`、`genes.gtf`、`genes.bed` 和 `chrom_sizes.tsv` 中的染色体或 contig 名称必须一致，例如不能在不同文件中混用 `chr1` 和 `1`。如果转录本注释与 `genes.gtf` 不同，可以显式传入 `--isoform-gtf`；否则 Flora 会复用 `genes.gtf`。

`--genome-mmi` 指定的预构建剪接索引须与参考 FASTA 匹配。Bonefish 处理分片索引时
需要普通 cDNA FASTQ。遇到这种情况请添加 `--full-output`，让 Flora 先写出供比对的
cDNA FASTQ。

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

对于拆分为多个文件的下机数据，可以在一个 `--fastq` 后用空格传入多个文件：

```bash
./flora \
  --fastq /data/sample_1.fastq.gz /data/sample_2.fastq.gz /data/sample_10.fastq.gz \
  --barcode-list-10bp /data/BC_1536.txt \
  --ref-dir /data/GRCh38_flora \
  --out-dir ./sample_output \
  --sample-id sample
```

也可以直接传入目录。Flora 搜索该目录当前层的 `.fastq.gz` 和 `.fq.gz` 文件，
文件名判断忽略大小写。

```bash
./flora --fastq-dir /data/sample_chunks \
  --barcode-list-10bp /data/BC_1536.txt \
  --ref-dir /data/GRCh38_flora \
  --out-dir ./sample_output \
  --sample-id sample
```

Glycine 默认最多同时运行10个任务，并在这些任务之间共享64个线程。可使用
`--glycine-jobs` 和 `--glycine-threads` 分别调整。最终只生成一个合并后的
`*.full-length-plus-rescued.fq.gz` 和一个
`*.identifying_statistic.txt`；计数逐项求和，比例按合并后的总数重新计算。

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

混合物种模式会在 HTML 报告中增加 Barnyard QC 部分，并写出
`qc/barnyard_qc/barnyard_summary.tsv` 和 `qc/barnyard_qc/barnyard_per_cell.tsv`。

## 实测资源消耗

一次受监控的内存优化版 `flora mixed --skip-glycine` 运行使用了 149 GB 的压缩全长 FASTQ，共包含 138,615,368 条 reads。在流程线程为 32、cluster 线程为 16 时，进程树约运行 884 分钟（14 小时 44 分）。实际观测到的 PSS 峰值为 94.964 GiB，RSS 求和峰值为 94.973 GiB。

## 结果目录

Flora 将结果写入 `--out-dir`，文件名以 `--sample-id` 为前缀。例如使用
`--out-dir ./sample_output --sample-id sample` 时，打开
`./sample_output/qc/sample.single_cell_report.html` 查看报告。下表列出了默认
轻量模式成功完成后的主要结果，路径均相对于 `--out-dir`。

| 路径 | 内容 |
| --- | --- |
| `qc/<sample>.single_cell_report.html`、`qc/metrics_summary.xlsx` | HTML 报告和指标工作簿。 |
| `qc/<sample>.rna_qc_metrics.tsv`、`qc/<sample>.per_cell_qc.tsv`、`qc/<sample>.saturation.tsv` | RNA QC、逐细胞指标和饱和度数据。 |
| `qc/<sample>.parameters.tsv` | 此次运行实际使用的参数。 |
| `matrix/<sample>.gene_expression.tsv`、`matrix/<sample>.isoform_expression.tsv` | 基因和转录本表达量。 |
| `matrix/<sample>.rna_cluster.tsv` | RNA 聚类结果及可用的 UMAP 坐标。 |
| `matrix/<sample>.filtered.tagged.sorted.bam`、`matrix/<sample>.filtered.tagged.sorted.bam.bai` | 最终带标签的 BAM 及索引。 |
| `upstream/barcode_to_cell.csv` | Barcode 到细胞的对应关系。 |
| `glycine/<sample>.full-length-plus-rescued.fq.gz` | Glycine 运行且使用默认输出路径时生成的全长 reads。 |
| `logs/` | 各阶段日志。基因组比对看 `01_alignment.log`，特征指派看 `04_gene.log`。 |

默认的 ONT 特征指派还会写出组装转录本
`transcriptome/<sample>.transcriptome.fa.gz`、
`transcriptome/<sample>.transcriptome.gtf.gz`，以及
`transcriptome/tool_versions.json`。各 MEX 目录位于 `matrix/` 下。

| 矩阵 | 原始 UMI 计数 | 过滤、归一化后的表达量 |
| --- | --- | --- |
| 基因 | `<sample>.gene_raw_feature_bc_matrix/` | `<sample>.gene_processed_feature_bc_matrix/` |
| 转录本 | `<sample>.transcript_raw_feature_bc_matrix/` | `<sample>.transcript_processed_feature_bc_matrix/` |

每个 MEX 目录包含 `matrix.mtx.gz`、`features.tsv.gz` 和保存 Flora 细胞 ID 的
`barcodes.tsv.gz`。`--skip-isoform` 跳过转录本表达量和转录本 MEX 矩阵。
`--upstream-only` 在细胞分配完成后结束，并列出 `upstream/` 中可用的文件。

轻量模式成功完成后，Flora 会清理初始比对 BAM 和逐 read 指派表等中间文件，
以及清理后为空的 `alignment/` 目录。

- `--save-intermediate` 保留修剪后的 cDNA FASTQ、初始比对 BAM 等调试文件。
- `--full-output` 还会在 `upstream/` 写出匹配、未匹配和已分到细胞的 FASTQ。
- `--remove-final-bam` 在分析成功后删除 `alignment/` 和 `matrix/` 下流程生成的 BAM 和索引，包括最终带标签 BAM。对这些文件，删除优先于两个保留选项。`transcriptome/intermediates/` 中的 BAM 按中间文件保留设置处理。

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
