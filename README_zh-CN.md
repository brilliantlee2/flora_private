[![English](https://img.shields.io/badge/Language-English-2563eb)](README.md)
[![中文](https://img.shields.io/badge/语言-中文-0f766e)](README_zh-CN.md)

# Flora 私有开发仓库

本仓库保存 Flora 私有源码、测试和打包工具，必须保持 private。
公开仓库 `brilliantlee2/Flora` 保存用户文档，二进制压缩包通过 GitHub Releases 发布。

## 克隆私有仓库

```bash
git clone https://github.com/brilliantlee2/flora_private.git
cd flora_private
```

需要私有仓库访问权限，使用 GitHub Token 或已授权的 SSH key。
SSH 地址为 `git@github.com:brilliantlee2/flora_private.git`。
修改过的依赖已包含在 `vendor/` 中，不需要配置 submodule。

## 构建环境

在 Linux x86_64 的仓库根目录执行：

```bash
conda env create -f environment.yml
conda activate flora

export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export PYTHONUTF8=1
```

已有环境时，用 `conda env update -n flora -f environment.yml` 代替创建命令，
然后重新激活。必须使用 Python 3.11.x。
如果 Flora 报 `Bad magic number in .pyc file`，重新激活 Flora 环境，并确认
`python3 --version` 显示 3.11.x。
`environment.yml` 包含构建工具和运行依赖；`environment.runtime.yml` 用于
已打包的可执行文件。`requirements.txt` 只列出 Python 包。
[运行工具清单](runtime-tools.json) 固定 Bonefish 和补丁版 StringTie 的来源与校验值。
Conda 环境包含安装所需的依赖，以及注释工具 `gffread` 和 `gffcompare`。

编译前为 `rust-htslib`/`bindgen` 设置 libclang 路径。
这仅是构建依赖，运行已打包的 Flora 不需要单独设置。

```bash
LIBCLANG_SO="$(find -L "$CONDA_PREFIX/lib" -maxdepth 1 -type f \
  \( -name 'libclang.so' -o -name 'libclang.so.*' \) -print -quit)"
if [[ -n "$LIBCLANG_SO" ]]; then
  export LIBCLANG_PATH="$(dirname "$LIBCLANG_SO")"
  echo "LIBCLANG_PATH=$LIBCLANG_PATH"
else
  unset LIBCLANG_PATH
  echo "未找到 libclang，请更新 Conda 环境后再编译。" >&2
fi
```

找不到时，执行 `conda install -c conda-forge libclang`，然后重新查找；成功后再编译。
上述示例不会退出当前 SSH 会话。新开终端后需重新设置，也可在查找成功后执行
`conda env config vars set LIBCLANG_PATH="$LIBCLANG_PATH"`，再重新激活 `flora`。

## 开发构建与测试

每个准备发版的源码版本执行一次以下检查即可：

```bash
cargo test --locked
bash -n run_all.sh
bash -n run_all_mixed_species.sh
bash -n packaging/build_binary_release.sh
bash -n packaging/refresh_binary_release_metadata.sh
python -m unittest discover -s tests -v
```

每条命令均成功后再继续。`bash -n` 需要逐个检查脚本，不能在后面列出多个路径
来检查所有文件。这些测试不能代替代表性数据的全流程验证。

只生成开发用可执行文件、不打包时：

```bash
cargo build --release --locked
./target/release/flora --version
./target/release/flora --help
./target/release/flora mixed --help
./target/release/flora glycine --help
./target/release/flora analyze --help
```

`flora` 运行单物种全流程，`flora mixed` 运行混合物种全流程；
`flora glycine` 只做全长识别，`flora analyze` 只做 barcode/UMI 提取与细胞分配。
全流程支持 `--fastq` 后输入多个文件，或用 `--fastq-dir` 非递归读取目录。
Glycine 默认 10 个并发任务，共享 64 个总线程。

运行默认流程前，在已激活的环境中安装一次这两个工具。

```bash
./target/release/flora setup-tools
```

安装程序会校验并安装 Bonefish CPU 0.5.0，并构建性能对照中使用的
StringTie 2.2.2 same-start 补丁版。支持 Linux x86_64 和 Apple Silicon macOS。
也可在编译 Flora 前从源码目录执行 `python3 packaging/setup_tools.py`。

基因组比对默认使用 Bonefish 0.5.0。Flora 会自动在
`${XDG_DATA_HOME:-$HOME/.local/share}/flora/tools` 下查找工具。
需要用共享目录或自定义目录时，在安装和运行前都设置 `FLORA_TOOLS_DIR`。
可通过 `--bonefish-bin` 或 `--stringtie-bin` 指定已有的二进制。
`--genome-aligner minimap2` 可选择 minimap2 做基因组比对。
ONT 转录组重比对使用 minimap2。全量对照命令见 [PBMC 性能验证](docs/pbmc-perf.md)。

默认轻量模式通过管道传递剪切后的 cDNA。使用 `--genome-mmi` 指定多分片索引时，
Bonefish 需要重新读取 cDNA，请加 `--save-intermediate` 或 `--full-output`，
让流程在比对前先写出 FASTQ 文件。请使用由相同参考 FASTA 构建的索引，
并为保留的 FASTQ 预留磁盘空间。

超长窗口的接缝重试默认关闭。要在基因组比对中做实验，可在运行命令中同时加入
`--bonefish-seam-retry` 和 `--bonefish-bin /path/to/patched/bonefish`，
选择接缝重试 A/B 中使用的、支持 `BONEFISH_SEAM_RETRY` 的 Bonefish 0.5.0
实验版。该重试缩短了 9011 子集的比对时间，同时改变了部分剪接链和独立组装后的
逐 read 指派。Flora 会在参数 TSV 中记录是否启用。

流程的大型临时文件默认放在 `<out-dir>/.flora-scratch`。可用
`--ont-temp-dir /path/to/scratch` 选择单独的工作盘。ONT 指派日志记录启动时
的临时盘余量，排序或读写失败时也会记录临时盘和结果盘的剩余空间。

## 构建 Linux x86_64 发行包

完成环境准备和上述测试后执行：

```bash
PYTHON_BIN="$CONDA_PREFIX/bin/python" \
  bash packaging/build_binary_release.sh
```

打包脚本会自行执行锁定依赖的 release 编译，无需先单独执行 `cargo build`。
它会检查 Linux x86_64 可执行文件和 CLI，复制运行环境与公开文档，
并自动生成、验证 SHA256。Python 资源嵌入在 `flora` 内，不以独立源码或字节码文件分发。

输出如下，版本号自动读取 `Cargo.toml`：

```text
dist/Flora-<version>-linux-x86_64/
dist/Flora-<version>-linux-x86_64.tar.gz
dist/Flora-<version>-linux-x86_64.tar.gz.sha256
```

发行包内的 `BUILD_INFO.txt` 记录所需 glibc 版本。请在兼容的 Linux 主机上构建，
此脚本不能在 macOS 上生成 Linux 发行包。可选前缀 `FLORA_MAX_GLIBC=2.34`
会拒绝依赖更高版本符号的构建，但不会降低二进制的 glibc 要求。
较旧系统需使用兼容的构建环境或 [Singularity](docs/singularity-build_zh-CN.md)。

发布前，将压缩包解压到新目录，用包内的 `environment.yml` 准备并激活
Python 3.11 运行环境，独立检查：

```bash
./flora --version
./flora --help
./flora mixed --help
./flora glycine --help
./flora analyze --help
./flora setup-tools
cat BUILD_INFO.txt
```

发版前还应运行代表性数据；帮助命令只能验证 CLI，不能证明分析结果正确。
压缩包应当只有 `flora` 一个可执行文件，不包含私有源码。

## 结果目录与文件保留

建议每次运行用独立的 `--out-dir`，方便区分结果和日志。流程直接写入该目录。
例如设置 `--out-dir /data/results/sample --sample-id sample`，终端显示
`Pipeline completed successfully.` 后，打开
`/data/results/sample/qc/sample.single_cell_report.html` 查看报告。

下表路径均相对于 `--out-dir`。`<sample>` 是 `--sample-id` 的值，默认是 `sample`。

| 路径 | 内容 |
| --- | --- |
| `qc/<sample>.single_cell_report.html` | 交互式报告。 |
| `qc/metrics_summary.xlsx` | 汇总指标（含全深度饱和度）。 |
| `qc/<sample>.parameters.tsv`、`logs/` | 本次运行参数和各阶段日志。基因组比对诊断见 `logs/01_alignment.log`。 |
| `qc/<sample>.rna_qc_metrics.tsv`、`qc/<sample>.per_cell_qc.tsv`、`qc/<sample>.saturation.tsv` | 样本 QC、单细胞 QC 和饱和度。 |
| `upstream/barcode_to_cell.csv` | Barcode 与细胞 ID 的对应关系。 |
| `matrix/<sample>.gene_expression.tsv`、`matrix/<sample>.isoform_expression.tsv` | 基因和异构体表达计数。 |
| `matrix/<sample>.rna_cluster.tsv` | RNA 聚类结果。 |
| `matrix/<sample>.filtered.tagged.sorted.bam` 和 `.bam.bai` | 最终带标签的比对结果及索引。 |
| `glycine/<sample>.full-length-plus-rescued.fq.gz` | 运行 Glycine 时，全长和挽救 reads 的默认保存路径。 |
| `qc/barnyard_qc/` | `flora mixed` 的物种混合质控。 |

默认 ONT 特征识别流程还会输出以下文件和目录。

| 路径 | 内容 |
| --- | --- |
| `transcriptome/<sample>.transcriptome.fa.gz` | 重构的转录本序列。 |
| `transcriptome/<sample>.transcriptome.gtf.gz` | 重构的转录本注释。 |
| `transcriptome/tool_versions.json` | 实际工具版本。 |
| `transcriptome/assignment_summary.json` | 基因、转录本分配状态计数。 |
| `matrix/<sample>.gene_raw_feature_bc_matrix/` | 基因层面的整数 UMI 计数。 |
| `matrix/<sample>.gene_processed_feature_bc_matrix/` | 过滤、归一化后的基因表达矩阵。 |
| `matrix/<sample>.transcript_raw_feature_bc_matrix/` | 转录本层面的整数 UMI 计数。 |
| `matrix/<sample>.transcript_processed_feature_bc_matrix/` | 过滤、归一化后的转录本表达矩阵。 |
| `matrix/<sample>.matrix_processing.json` | 矩阵过滤参数及维度。 |

每个 MEX 目录包含 `matrix.mtx.gz`、`features.tsv.gz` 和
`barcodes.tsv.gz`（细胞 ID）。Processed 矩阵先按每细胞至少 200 个特征、
每特征至少 3 个细胞、线粒体比例不超过 20% 过滤，再按每细胞 10,000 计数归一化并执行
`log10(1+x)`。其中的细胞集合可能与报告的原始细胞判定不同。
`--skip-isoform` 会跳过异构体表达 TSV 和转录本 MEX 矩阵。ONT 流程仍会构建转录组，
供联合分配使用。

默认轻量模式在流程成功完成后清理可重建的逐 read 表、分配 TSV 以及加标签前的
比对 BAM 和索引。`alignment/` 变空时也会移除。下游分析请用 `matrix/` 中的最终 BAM 和索引。

| 选项 | 作用 |
| --- | --- |
| `--save-intermediate` | 保留 `alignment/<sample>.cdna.fastq.gz`、加标签前的 BAM 和索引及其他中间文件。 |
| `--full-output` | 保留比对用 FASTQ 和中间文件，并在 `upstream/` 下写出 `matched_reads.fastq.gz`、`unmatched_reads.fastq.gz` 和 `cell_reads.fastq.gz`。 |
| `--remove-final-bam` | 成功分析后删除 `alignment/` 与 `matrix/` 下由流程生成的 BAM 和索引，包括最终带标签的 BAM。 |

这些参数控制文件保留。BAM 删除选项与保留选项同时使用时，仍会清理 `alignment/`
与 `matrix/` 下的 BAM 和索引。保留的 ONT 调试文件位于 `transcriptome/intermediates/`。
`--upstream-only` 在 barcode/UMI 提取与细胞分配后结束。
若需要 Glycine 其他分类的 FASTQ，请单独运行 `flora glycine` 并添加 `--keep-all-outputs`。

## 发布新版本

1. 修改 `Cargo.toml`，然后执行不带 `--locked` 的 `cargo check`，
   更新 `Cargo.lock` 中 Flora 自身的版本号。
2. 同步更新公开文档中的版本、release notes 和 `packaging/singularity/Flora.def`。
3. 按上述步骤测试、打包，并验证发行包。
4. 将 `docs/repository-templates/public/` 中审核后的文件同步到公开仓库。
   `.tar.gz` 和 `.sha256` 作为 Release 附件上传；SIF 和校验文件为可选附件。

Tag 与标题使用同一版本，例如 `v0.1.2` 和 `Flora v0.1.2`。
不要将压缩包、SIF、Docker tar、Miniforge 安装器或 `dist/` 提交到 Git。
不要向公开仓库复制私有的 `src/`、`scripts/`、`vendor/`、Cargo 文件、测试和打包代码。
编译后的可执行文件并非加密，仍可能被逆向。

## 源码目录

| 路径 | 用途 |
| --- | --- |
| `src/`、`vendor/`、`Cargo.toml`、`Cargo.lock` | Rust 源码与锁定依赖。 |
| `scripts/` | Python 分析、QC 与报告资源。 |
| `run_all.sh`、`run_all_mixed_species.sh` | 单物种与混合物种流程脚本。 |
| `tests/` | 回归与发行检查。 |
| `packaging/` | 二进制与容器打包工具。 |
| `docs/repository-templates/` | 公开/私有仓库文档模板。 |

容器详细步骤：[构建 Singularity 镜像](docs/singularity-build_zh-CN.md)。
