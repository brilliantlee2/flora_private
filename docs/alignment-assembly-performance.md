# 人类参考比对与 StringTie 热点性能验证

日期：2026-11-11。基线为本地 main `291de65`；测试分支为 `codex/1111-alignment-assembly-performance`。本分支不是已验证的新版发行包，也未修改 main。

## 证据与修改边界

采样结果显示，chr2 的 927,114 条比对中，634,090 条集中于 `chr2:88538720-89330679`。StringTie 2.2.2 / 2.2.3 总耗时 338.53 / 529.34 秒，其中开始组装前的空档为 322 / 513 秒；空档中主线程约 99.7% CPU，四个消费者基本等待。因此更应检查输入结构归并，不应只提高消费者线程或升级版本。chr22 也有开始前长空档，但尚没有与 chr2 同等完整的线程/函数级证据。

同起点逐项比较在官方源码中确实存在；目前的证据尚未证明它是服务器真实数据上唯一的耗时函数。可选加速器的意义是用保留原比较条件的索引化实验验证这个机制，不是宣称作者已发布了 chr2/chr22 专用补丁。

| 项目 | main | 本测试分支 | 是否改变检出规则 |
|---|---|---|---|
| minimap2 query 批次 | 固定 `-K 10M` | 默认 `--minimap2-batch 500M`，可恢复 10M | 目标不改变；必须对照输出 |
| mapper I/O | 一个 I/O 线程 | 默认 `--minimap2-io-threads 2`，传入 `-2` | 目标不改变；峰值内存可能增加 |
| 保存的 gzip cDNA | mapper 内部解压 | 默认外部 `gzip -dc` 供应 stdin | 保留输入；可恢复 direct |
| 默认 cDNA FIFO | 未压缩流 | 仍为未压缩流，不再压缩/解压 | 不改变 |
| BAM 过滤管道 | 两次 samtools view | 一次 `view -uh --no-PG -F 256` | 相同 secondary 过滤，PG 头可不同 |
| StringTie | PATH 上的 2.2.2 | 默认不替换；可用 `--stringtie-bin` 指定隔离构建 | 不升级模型或改变阈值 |
| StringTie 同起点查找 | 逆序遍历 | 可选索引筛候选后保留原比较顺序/条件 | 实验目标为等价；需真实结果验证 |

未更改 `-ax splice -uf --MD`、`--junc-bed`、`-I 16G`、k15/w5、`--cap-kalloc 100m`、基因 MAPQ、转录本候选规则、`-L -c 2` 或 UMI 逻辑。不会删除免疫球蛋白区域、低 MAPQ、supplementary，或任意按窗口切开 bundle。`--ont-jobs` 仍默认 1，可在资源允许时使用既有染色体并行，但它不解决单个热点内部的串行计算。

## 官方资料如何支持这次修改

- minimap2 作者在 [#771 回复](https://github.com/lh3/minimap2/issues/771#issuecomment-868865187) 建议内存足够时提高批次。[v2.31 手册](https://github.com/lh3/minimap2/blob/3c28777e7e2dcc90f825de1b9f17a89cca7d4452/minimap2.1#L591) 的默认批次是 500M **query 碱基**，不是 5 亿条 reads；`-2` 将输入/输出分到两个 I/O 线程，可能增加内存。本分支采用该默认批次，而不是盲目加大线程数。
- [#1021 作者回复](https://github.com/lh3/minimap2/issues/1021#issuecomment-1446595850) 解释 zlib 单线程解压及外部解压流。对默认 FIFO，这项 gzip 改动不会凭空提速；主要用于 `--save-intermediate` 的落盘 gzip 输入。
- 当前稳定版 ONT wf-single-cell 为 [v3.3.4](https://github.com/epi2me-labs/wf-single-cell/releases/tag/v3.3.4)。[实际调用](https://github.com/epi2me-labs/wf-single-cell/blob/cb39932e3dfca08ad6dae1e6d08579c957874fcf/subworkflows/preprocess.nf#L36) 也用 `-K 10M`，所以不能把小批次描述为“Flora 没有对齐 ONT”。ONT 可调度多个 FASTQ chunk mapper，合并后再按染色体组装；本轮没有引入新的 chunk 调度器。
- ONT [无 preset 预建索引](https://github.com/epi2me-labs/wf-single-cell/blob/cb39932e3dfca08ad6dae1e6d08579c957874fcf/modules/local/common.nf#L24) 导致 stock minimap2 使用 w10，而 Flora 的 splice FASTA 使用 w5。该差异可能影响速度，也影响 seeding/敏感性，因此本轮不悄悄改成 w10；这不是免费的等价提速。
- StringTie 作者在 [#464](https://github.com/gpertea/stringtie/issues/464#issuecomment-2721034880) 确认单输入生产者架构；[后续回复](https://github.com/gpertea/stringtie/issues/464#issuecomment-2721051151) 支持按染色体拆分，但警告任意切分会破坏 bundle。Flora 已按染色体运行，不能把这个既有行为重新算成修复。
- [#221 作者回复](https://github.com/gpertea/stringtie/issues/221#issuecomment-488778093) 提到 `-j 3`，同时明确低表达检出可能下降。[ONT #168](https://github.com/epi2me-labs/wf-single-cell/issues/168#issuecomment-3017923731) 是用户个案，不是通用作者修复。本轮不默认调高 `-j`。若单独做阈值实验，须用 `-c 2 -j 3`，避免用 `-j 3` 替换整串后意外恢复默认 `-c 1`。
- ONT 的 [2.2.3 升级记录](https://github.com/epi2me-labs/wf-single-cell/blob/cb39932e3dfca08ad6dae1e6d08579c957874fcf/CHANGELOG.md#L62) 不等于本数据已提速；该补丁没有替换同起点扫描。本轮保留 main 的 2.2.2，隔离构建同版本 stock 和实验版本，避免混淆版本算法与数据结构变化。

## 环境与构建

可复用 main 的构建环境或现有 `flora-ont330`（StringTie 2.2.2）。不需要修改正在运行的环境。先从本分支按 README 设置 libclang，构建新的 Flora；旧打包二进制不认识新增选项。

```bash
conda activate flora-ont330
PYTHON_BIN="$CONDA_PREFIX/bin/python" cargo build --release --locked
./target/release/flora --help
```

可选 StringTie 构建器会生成两个同编译环境的独立程序，使用哈希固定的官方源码，不覆盖 Conda 中的 stringtie。需要 C/C++ 编译器、make 和压缩库开发文件；详见 [构建器说明](stringtie-input-acceleration.md)。

```bash
BUILD="$HOME/0.software/stringtie-input-1111"
python packaging/stringtie_input/build.py --build-root "$BUILD" --fetch --jobs 4
STOCK="$BUILD/stringtie-2.2.2-stock"
FAST="$BUILD/stringtie-2.2.2-same-start"
```

`BUILD` 必须为新目录。两个程序仍报告上游版本 2.2.2，来源区别记录在 `build-provenance.json` 和 Flora 输出的 StringTie 路径/SHA256；实验构建不是官方发布版本，不能冒充系统安装的原版。

## 指定 sample_dir 的对照

以下使用同一目录，不重新抽样。`OUTBASE` 可调整，三个输出目录都必须是新目录。首次比较顺序运行，避免多个基准互相争 CPU/磁盘。记录同一实际 minimap2 版本和参考文件；不要在运行中升级工具。

```bash
FLORA="$(pwd)/target/release/flora"
OUTBASE="$HOME/2.project/C4_V3/performance_1111"
mkdir -p "$OUTBASE"
COMMON=(
  --fastq-dir /home/liyy/1.data/C4_Cyclone/gufen_1103/260F800579011/sample_dir
  --barcode-list-10bp /home/liyy/2.project/C4_V3/Strint2.6/BC_1536.txt
  --ref-dir /home/liyy/1.data/REF/GRCH38
  --sample-id 260F800579011 --exp-cells 7500
  --threads 32 --cluster-threads 16 --glycine-jobs 15
  --ont-jobs 1 --feature-backend ont --save-intermediate
  --stringtie-opts='-c 2 -v'
)

```

COMMON 数组在同一 Bash 终端中保留，再依次运行：

```bash
# A: previous mapping settings, isolated stock StringTie
"$FLORA" "${COMMON[@]}" --out-dir "$OUTBASE/A_previous" \
  --minimap2-batch 10M --minimap2-io-threads 1 --alignment-gzip-input direct \
  --stringtie-bin "$STOCK" > "$OUTBASE/A.log" 2>&1

# B: mapper changes only, same stock StringTie
"$FLORA" "${COMMON[@]}" --out-dir "$OUTBASE/B_mapping" \
  --stringtie-bin "$STOCK" > "$OUTBASE/B.log" 2>&1

# C: mapper changes plus isolated same-start accelerator
"$FLORA" "${COMMON[@]}" --out-dir "$OUTBASE/C_mapping_assembly" \
  --stringtie-bin "$FAST" > "$OUTBASE/C.log" 2>&1
```

这三次是完整的分层对照，可能耗时。若不想重复上游，StringTie A/B 可直接复用已完成测试中相同的完整 chr2/chr22 BAM 和 GTF，只替换两个独立程序；不可任意剪碎热点边界。只有同一编译器/库下的 stock 与 same-start 比较，才能分离编译环境差异。

同一个新 Flora 可恢复旧 mapper 行为，无需反复重新打包：`--minimap2-batch 10M --minimap2-io-threads 1 --alignment-gzip-input direct`。单独检验 batch 时保持其他选项相同，再检验 I/O / gzip；A/B/C 组合验证不等于各参数的独立因果实验。

## 验收与仍需验证的内容

1. 记录日志 Stage elapsed、minimap2 自报 Real time/CPU/Peak RSS、每染色体 StringTie elapsed，不能把父子计时重复相加。`--cap-kalloc` 不是整个进程的硬内存上限。
2. 核对 cDNA 解压序列一致，BAM 的 read ID、flag、位置、CIGAR、MAPQ、AS/NM 等内容一致；BAM 压缩字节和 PG/CL 头不是等价标准。
3. 核对 chr2/chr22 重构模型的 exon/junction 结构、覆盖、参考关联和 novel 模型；多线程可能改变输出顺序和生成的临时模型 ID，不能只比文件大小。
4. 比较 gene/transcript assignment、gene/isoform UMI 矩阵和 QC。任何差异须解释，不能只看到提速就 merge。
5. 本地已运行 170 项 Python 测试（3 项因完整工具链/全流程设置不可用而跳过）、109 项 Rust 测试，以及两个 shell 脚本各自的语法检查。真实 minimap2 2.31 合成剪接数据的新旧管道解码一致；StringTie 已做原版/实验版线性/实验版索引的逐条 processRead 状态、组装、丰度和 merge 对照。服务器 sample_dir 和完整 chr2/chr22 输入不在本机，尚未证明真实数据的加速倍数或全流程完全等价。

本轮不改 saturation 图或统计口径，避免与性能对照混在一起。暂不将可选 StringTie 加速器作为正式发行默认值。
