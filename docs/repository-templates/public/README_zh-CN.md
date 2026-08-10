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

Glycine 已集成到 Flora 中，用户不需要单独安装 Glycine，也不需要传入 `--glycine-bin-dir`。

## 适用平台

当前二进制发行包需要：

- Linux x86_64；
- 与 glibc 兼容的 Linux 系统（发行包在 Ubuntu 22.04/glibc 2.35 上构建）；
- Conda、Miniforge、Mambaforge 或 Micromamba；
- 由环境文件安装的 Python 3.11。

该发行包不支持 macOS、ARM Linux 或 Windows。

## 下载

请从 [GitHub Releases](https://github.com/brilliantlee2/Flora/releases) 下载最新版本。

Flora v0.1.0：

```bash
wget https://github.com/brilliantlee2/Flora/releases/download/v0.1.0/Flora-0.1.0-linux-x86_64.tar.gz
wget https://github.com/brilliantlee2/Flora/releases/download/v0.1.0/Flora-0.1.0-linux-x86_64.tar.gz.sha256

sha256sum -c Flora-0.1.0-linux-x86_64.tar.gz.sha256
tar -xzf Flora-0.1.0-linux-x86_64.tar.gz
cd Flora-0.1.0-linux-x86_64
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

## 分析原始 FASTQ

内置 Glycine 会自动运行：

```bash
./flora \
  --fastq /data/sample.fastq.gz \
  --barcode-list-10bp /data/BC_1536.txt \
  --ref-dir /data/GRCh38_flora \
  --out-dir ./sample_output \
  --sample-id sample \
  --threads 32 \
  --cluster-threads 8 \
  --top1-alpha 0.1 \
  --max-ed 2
```

## 分析已有全长 FASTQ

```bash
./flora \
  --skip-glycine \
  --full-length-fastq /data/sample.full-length-plus-rescued.fq.gz \
  --barcode-list-10bp /data/BC_1536.txt \
  --ref-dir /data/GRCh38_flora \
  --out-dir ./sample_output \
  --sample-id sample \
  --threads 32 \
  --cluster-threads 8 \
  --top1-alpha 0.1 \
  --max-ed 2
```

## Mixed-species 分析

```bash
./flora mixed \
  --skip-glycine \
  --full-length-fastq /data/mixed.full-length-plus-rescued.fq.gz \
  --barcode-list-10bp /data/BC_1536.txt \
  --ref-dir /data/merged_reference \
  --out-dir ./mixed_output \
  --sample-id mixed_sample \
  --threads 32 \
  --cluster-threads 8
```

## 初始资源建议

资源消耗受 read 数、read 长度、参考基因组、barcode 多样性、存储速度和线程数影响。下表是保守的首次投递建议，不是保证上限。

| 压缩 FASTQ | 线程 | 建议内存 | 首次任务时间 |
|---:|---:|---:|---:|
| 不超过 5 GB | 16-24 | 96-128 GB | 4-8 h |
| 5-20 GB | 24-32 | 192-256 GB | 12-24 h |
| 20-50 GB | 32 | 384-512 GB | 24-48 h |
| 50-100 GB | 32-48 | 768 GB-1 TB | 48-96 h |

增加线程可能提高内存峰值，且不一定带来线性加速。降低任务内存前，请先使用代表性样本做基准测试。

## 主要输出

```text
upstream/    barcode 矫正、cell assignment 和 knee plot
alignment/   比对与标签 BAM
matrix/      基因/转录本矩阵和 RNA 聚类坐标
qc/          RNA QC、饱和度和报告输入
logs/        各步骤日志
```

关键结果包括 `read_assigned_cell.csv`、`barcode_to_cell.csv`、带标签 BAM、基因和转录本表达矩阵、Scanpy UMAP 坐标和 `<sample>.single_cell_report.html`。

## 问题反馈

请通过 [GitHub Issues](https://github.com/brilliantlee2/Flora/issues) 提交可复现的问题，并附上 Flora 版本、运行命令、操作系统、输入大小、资源申请和相关日志。请勿上传私有 FASTQ/BAM 数据。

## 许可声明

Flora 尚未声明项目整体许可证。内置 Glycine 保留 MIT 许可和作者声明，每个发行包中均包含第三方许可说明。
