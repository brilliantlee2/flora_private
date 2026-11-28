# PBMC 性能分支

`codex/perf` 以 Flora `main` 为基线。下面的实测收益来自各自的 A/B 对照。全量 PBMC 的总耗时和组装差异由集成对照衡量。

| 环节 | 分支实现 | 之前的对照 |
| --- | --- | --- |
| Step 2 基因组比对 | 默认 Bonefish CPU 0.5.0，使用 `-a -x cdna -u f --MD -I 16G -K 10M` 和原 junction BED。不指定 MMI 时根据输入 FASTA 建索引。 | 120 万条 cDNA 第一轮到 BAM 和 BAI 就绪为 412.545 → 202.447 秒。第二轮 minimap2 耗时异常。chr19 重新组装后，0.51% 的共同 reads 改变了基因指派。 |
| Step 2 BAM 索引 | 保留 `-F 256`。坐标排序时用 `samtools sort --write-index` 生成对应的 `.bam.bai`。 | 120 万条 cDNA 的输出子阶段平均 15.841 → 11.428 秒，解码 BAM 和索引查询一致。 |
| Step 4 StringTie2 | 自动使用 `flora setup-tools` 构建的 2.2.2 same-start，或通过 `--stringtie-bin` 指定，默认 `STRINGTIE_INPUT_LOOKUP=indexed`。 | chr19 两轮组装 945 → 627 秒，归一化组装内容和第一轮逐 read 指派一致。 |
| Step 4 GN/TX 回填 | 对 `(chrom, read_id)` 主键做有界的 `IN` 批量查询，分块遵守 SQLite 参数上限。 | chr19 回填至 BAM 和 BAI 就绪的耗时分别减少约 63% 和 40%，六项输出核对一致。 |
| Step 4 UMI 聚类 | 原样使用 `edf000a` 的 15 bp、编辑距离 3 的流式精确实现，其他输入仍走原有算法。 | 全量从 7,840 秒降到 1,879.58 秒，解码 BAM 和基因矩阵一致，新版峰值 RSS 3.26 GiB。 |

`--genome-aligner minimap2` 选择基因组 minimap2。它仍使用 `-ax splice -uf --MD` 和原内存限制。Step 4 转录组重比对使用 minimap2。旧版 `--secondary=no`、`-K 100M`、`--ont-jobs 2` 和 Glycine 的 Edlib 实验均属于单独对照。

## 全量验证结果

- PBMC 两套全流程分别用时 25 小时 34 分 07 秒和 18 小时 33 分 46 秒。两套结果识别到相同的 6,894 个细胞。以逐项 UMI 计数差的绝对值之和除以两套总 UMI，基因矩阵差异为 1.71%，转录本矩阵差异为 4.69%。两次运行使用了不同路径的全长 FASTQ，新流程同时改变了比对器、StringTie 构建和其他环节。上述结果用于比较两套全流程。
- 9011 的 75,086,206 条 reads 使用同一份剪切后 cDNA、同一个单片 MMI 和 32 线程。实验版 Bonefish 开启接缝重试后，从开始比对到 BAM 和 BAI 就绪由 19,661 秒降到 6,314 秒，耗时减少 67.89%。全量有 396,000 条 reads 的主比对剪接链改变。在 chr2、chr14、chr22 独立组装后的共同 reads 中，12.54% 改变了基因指派，13.01% 改变了参考转录本指派。`--bonefish-seam-retry` 默认为关闭。

## 运行依赖

运行使用包含 samtools 1.24 和 minimap2 2.31 的 Flora 环境。构建环境与安装步骤见 [README](../README_zh-CN.md)。
下面的 `flora setup-tools` 会按 [运行工具清单](../runtime-tools.json) 下载
Bonefish CPU 0.5.0，并用原脚本构建 StringTie 2.2.2 same-start。
Flora 自动使用安装好的程序。需要复用 A/B 测试中的二进制时，可用
`--bonefish-bin` 和 `--stringtie-bin` 指定。可以通过
`transcriptome/tool_versions.json` 中的路径、SHA-256 和模式区分两种
`--version` 均显示 2.2.2 的 StringTie 构建。

## cycloneserver 上的全量 PBMC 集成对照

在服务器的 `codex/perf` 工作树执行。命令复用 11 月 23 日基线的参数表，保留调试中间文件。若基线记录的全长 FASTQ 已移走，在启动前将 `PBMC_FULL_LENGTH_FASTQ` 设置为当前保存的全长 FASTQ 路径。命令先用基线 FASTA 准备单片 k15/w5 MMI，单独记录索引准备时间。

```bash
conda activate flora-ont330
export PYTHON_BIN="$CONDA_PREFIX/bin/python"
if [[ -z ${LIBCLANG_PATH:-} ]]; then
  clang_lib=$(find -L "$CONDA_PREFIX/lib" -maxdepth 1 -type f -name 'libclang.so*' -print -quit)
  if [[ -n $clang_lib ]]; then export LIBCLANG_PATH="$(dirname "$clang_lib")"; fi
fi
cargo build --release --locked --bin flora
./target/release/flora setup-tools --jobs 8
```

```bash
(
  set -euo pipefail
  S=260F700591011_barcode6
  DIST=/home/liyy/0.software/test_space/flora_private_saturation_re/dist/Flora-0.1.2-linux-x86_64
  BASE="$DIST/${S}_perf_baseline_20261123_123741"
  PARAMS="$BASE/results/qc/$S.parameters.tsv"
  FLORA="$PWD/target/release/flora"
  [[ -s "$PARAMS" && -x "$FLORA" ]]

  get() { awk -F '\t' -v k="$1" '$1==k {sub(/^[^\t]*\t/,""); print; exit}' "$PARAMS"; }
  [[ "$(get skip_glycine)" == 1 && "$(get feature_backend)" == ont ]]
  args=(--skip-glycine --save-intermediate --sample-id "$S")
  for key in full_length_fastq barcode_list_10bp gene_fasta junction_bed chrom_sizes gene_gtf isoform_gtf; do
    value="$(get "$key")"
    if [[ "$key" == full_length_fastq ]]; then value="${PBMC_FULL_LENGTH_FASTQ:-$value}"; fi
    [[ -r "$value" && -s "$value" ]] || { printf '缺少 %s %s\n' "$key" "$value" >&2; exit 1; }
    args+=("--${key//_/-}" "$value")
  done
  for key in threads cluster_threads exp_cells min_q max_ed barcode_extract_mode pair_min \
      auto_pair_min_floor auto_pair_min_quantile top1_alpha top1_alpha_umi dominance_min \
      drop_umiA_ratio_gt gene_assign_mapq transcript_assign_mapq ref_interval cell_gene_max_reads ont_jobs; do
    value="$(get "$key")"
    if [[ -n "$value" ]]; then args+=("--${key//_/-}" "$value"); fi
  done
  if [[ "$(get require_pass_both_ends)" == 1 ]]; then args+=(--require-pass-both-ends); fi
  args+=(--feature-backend ont "--stringtie-opts=$(get stringtie_opts)")
  RUN=$(mktemp -d "$BASE/full_pbmc_perf.XXXXXX")
  MMI="$RUN/genome.splice.mmi"
  args+=(--out-dir "$RUN/results" --genome-mmi "$MMI")
  cmd=(env STRINGTIE_INPUT_LOOKUP=indexed "$FLORA" "${args[@]}")
  printf '%q ' "${cmd[@]}" >"$RUN/command.txt"
  printf '\n' >>"$RUN/command.txt"
  cat >"$RUN/run.sh" <<'BASH'
set -euo pipefail
run=$1
threads=$2
genome=$3
mmi=$4
shift 4
trap 'status=$?; printf "%s\n" "$status" >"$run/exit.status"' EXIT
/usr/bin/time -v -o "$run/genome-index.time" \
  minimap2 -x splice -I 16G -t "$threads" -d "$mmi" "$genome" >"$run/genome-index.log" 2>&1
"$@"
BASH
  nohup bash "$RUN/run.sh" "$RUN" "$(get threads)" "$(get gene_fasta)" "$MMI" "${cmd[@]}" \
    </dev/null >"$RUN/nohup.log" 2>&1 &
  printf '%s\n' "$!" >"$RUN/pid"
  printf 'PID %s\n结果目录 %s\n日志 %s\n' "$!" "$RUN/results" "$RUN/nohup.log"
)
```

`exit.status` 为 0 时表示流程完成。首次运行的耗时对照应计入 `genome-index.time`，后续复用索引时单列这项一次性成本。输出目录内的 `qc/<sample>.parameters.tsv`、`logs/01_alignment.log`、`logs/04_gene.log`、`transcriptome/tool_versions.json` 和 nohup 日志可用于计时及版本对照。重新组装的基因和 isoform 原始 UMI 矩阵用于判断全量结果变化。
