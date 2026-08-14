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

下面的命令会自动识别最新的正式 Release，后续发布新版本时无需再修改
README 中的版本号：

```bash
REPO="brilliantlee2/Flora"
LATEST_URL="$(curl -fsSL -o /dev/null -w '%{url_effective}' "https://github.com/${REPO}/releases/latest")"
LATEST_TAG="${LATEST_URL##*/}"

case "${LATEST_TAG}" in
  v[0-9]*) ;;
  *) echo "Unable to resolve the latest Flora release" >&2; exit 1 ;;
esac

VERSION="${LATEST_TAG#v}"
ARCHIVE="Flora-${VERSION}-linux-x86_64.tar.gz"

wget "https://github.com/${REPO}/releases/download/${LATEST_TAG}/${ARCHIVE}"
wget "https://github.com/${REPO}/releases/download/${LATEST_TAG}/${ARCHIVE}.sha256"

sha256sum -c "${ARCHIVE}.sha256"
tar -xzf "${ARCHIVE}"
cd "Flora-${VERSION}-linux-x86_64"
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
  --cluster-threads 16 \
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
  --cluster-threads 16 \
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
  --cluster-threads 16
```

## 初始资源建议

资源消耗受 read 数、read 长度、参考基因组、barcode 多样性、存储速度、流程模式和线程数影响。下表以当前内存优化版的实测结果为基准，属于保守的首次投递建议，不是保证上限。除约 150 GB 的实测范围外，其他范围均为外推建议。

| 压缩 FASTQ | CPU 线程 | 建议申请内存 | 首次任务时间 |
|---:|---:|---:|---:|
| 不超过 20 GB | 16-24 | 64-96 GB | 不超过 12 h |
| 20-75 GB | 24-32 | 96-128 GB | 12-20 h |
| 75-160 GB | 32 | 160 GB | 24-30 h |
| 超过 160 GB | 首次使用 32 | 192-256 GB | 36-48 h，建议先做基准测试 |

### 大数据实测基准

一次受监控的内存优化版 `flora mixed --skip-glycine` 运行使用了 149 GB 的压缩全长 FASTQ，共包含 138,615,368 条 reads。在流程线程为 32、cluster 线程为 16 时，受监控的进程树约运行 884 分钟（14 小时 44 分）。PSS 峰值为 94.964 GiB，RSS 求和峰值为 94.973 GiB，PSS 的第 95 百分位为 66.846 GiB。

对类似的约 150 GB `--skip-glycine` 任务，首次正式投递建议申请 **32 个 CPU slot、160 GB 内存和 24-30 小时**。如果希望保留更大余量、共享存储性能不稳定，或者提高了线程数，建议申请 192 GB。该基准只来自一个数据集和一套运行环境，不应视为硬性上限。

上述基准不包含 Glycine。对使用内置 Glycine 处理大型原始 FASTQ 的首次任务，建议在上述时间基础上至少再预留 8-12 小时。各阶段是串行执行的，内存不会简单相加；但当前优化后的未跳过 Glycine 全流程尚未完成端到端实测，在降低内存申请前应先监控一次正式任务。

增加线程可能提高内存峰值，且不一定带来线性加速。降低任务内存前，请先使用代表性样本做基准测试。

PSS 更适合表示进程树实际占用的物理内存，因为共享页会按比例分摊；RSS 直接求和可能重复计算共享页。本次监控的采样间隔为 1 分钟，因此报告峰值是已观测到的下限，可能遗漏持续时间短于采样间隔的尖峰。

使用默认 light-output 模式时，建议初次投递至少准备压缩 FASTQ 大小 3-5 倍的可写临时空间。对于 149 GB 输入，约需要 450-750 GB；条件允许时建议准备 1 TB。开启 full output 或保留中间文件时需要更多空间。首次正式运行期间可使用 `df -h` 和 `df -i` 同时监控磁盘容量与 inode。

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
