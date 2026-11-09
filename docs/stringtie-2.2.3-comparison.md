# StringTie 2.2.3 对照测试

日期：2026-11-09（Asia/Shanghai）。
分支：`codex/1109-stringtie-2.2.3`。
基线：本地 `main`，`291de65fa6b38a3abadd012095076124bc8d9e25`。

## 修改范围

- 构建/运行环境固定 `stringtie=2.2.3`，单物种、混合物种与 Python 后端均检查该版本。
- `transcriptome/tool_versions.json` 保留实际工具版本，并记录 `required_stringtie_version`；
  `qc/*.parameters.tsv` 在使用组装后端时也记录该要求，legacy 不记录。
- 不改 minimap2、samtools、gffread、gffcompare 的版本范围或调用参数；
  不改 Glycine、barcode、UMI、MAPQ、`-L -c 2`、线程/jobs、候选选择、矩阵或报告逻辑。
- 不合入独立 `codex/1106-ont-hotspot-fix` 的 UMI/调度优化，避免混入第二个实验变量。
- 历史 ONT commit/container 信息仍表示 feature-selection 参考来源，不表示正在运行旧容器。

StringTie 升级可能改变重构模型，进而改变 gene/transcript assignment 和表达矩阵。
这里保证的是 Flora 其他逻辑没有修改，不承诺生物学输出逐项相同或异常耗时必然消失。

ONT 在 [v3.3.1 更新日志](https://github.com/epi2me-labs/wf-single-cell/blob/master/CHANGELOG.md)
记录了升级到 2.2.3 以修复组装停滞；是否解决本样本仍需实测。

## 获取与环境隔离

分支推送后，在 Linux 的新目录克隆，不覆盖 `main` 的源码或发行包：

```bash
git clone --branch codex/1109-stringtie-2.2.3 --single-branch \
  https://github.com/brilliantlee2/flora_private.git flora_private_stringtie223
cd flora_private_stringtie223
git rev-parse HEAD
```

做严格对照时，优先克隆原来成功运行 main 的环境。下面假设原环境名为 `flora`，
若实际使用 `flora-ont330`，只把 `BASE_ENV` 改成该名称：

```bash
BASE_ENV=flora
conda create -n flora-stringtie223 --clone "$BASE_ENV"
conda list -n "$BASE_ENV" --explicit > main_environment_explicit.txt
conda install -n flora-stringtie223 --freeze-installed \
  --override-channels -c conda-forge -c bioconda 'stringtie=2.2.3'
conda list -n flora-stringtie223 --explicit > stringtie223_environment_explicit.txt
conda activate flora-stringtie223
stringtie --version
```

应显示 `2.2.3`。检查安装计划和两份环境清单，确认没有额外改变关键工具/依赖。
`--freeze-installed` 限制更新已安装的依赖，但仍可能安装新增依赖；若无法求解，
先保留原环境，不用 `--no-deps` 强行安装。[Conda 官方说明](https://docs.conda.io/projects/conda/en/stable/commands/install.html)

没有可克隆环境时，可直接创建：

```bash
conda env create -n flora-stringtie223 -f environment.yml
conda activate flora-stringtie223
```

这种方式会重新求解其他包的版本范围，不自动构成“只有 StringTie 不同”的环境。
原环境不更新；主流程仍在原环境使用 StringTie 2.2.2。

## 编译与运行

按 README 设置 `LIBCLANG_PATH` 后，在新仓库、新环境内执行：

```bash
PYTHON_BIN="$CONDA_PREFIX/bin/python" cargo build --release --locked
./target/release/flora --help
./target/release/flora mixed --help
```

如需原来的 Linux 发行包布局：

```bash
PYTHON_BIN="$CONDA_PREFIX/bin/python" bash packaging/build_binary_release.sh
```

Cargo 软件版本保持 0.1.2，两组发行包文件名和 `flora --version` 因此相同。
必须放在不同目录，并用明确的可执行文件路径运行；用源码 commit、日志及
`tool_versions.json` 区分组别，不能覆盖原来的发行包。旧二进制嵌入了旧检查，
不能只换 Conda 环境而继续运行未重新编译的 main 二进制。

保持相同输入（采样只生成一次）、参考、sample ID、线程、jobs 和其他参数，
只改输出目录和执行文件/环境。两组都建议增加：

```bash
--save-intermediate \
--ont-jobs 1 \
--stringtie-opts='-c 2 -v'
```

不要加 `--remove-final-bam`。不要同时测试 `-j 3` 或改变 MAPQ。
两组分别使用独立的输出目录和 scratch 子目录。若并行执行，使用独立 shell
激活各自环境，或用 `conda run -n <环境名> --no-capture-output <可执行文件> ...`。
功能比较可并行；可靠 wall-time 比较应串行运行或分配独立资源，避免 CPU、内存、
磁盘争用。两组都保留中间文件也会增加 I/O，不能与以前 FIFO 模式总耗时直接比较。
新 Conda 环境不会修改旧 SIF 镜像里的工具；这轮应使用本分支重新编译的原生可执行文件。

重点比较 `logs/04_gene.log` 的 chr2/chr22 StringTie 计时、转录组模型数量、
assignment 状态和最终矩阵。采样没有卡顿不证明全量问题已解决；升级 StringTie
也不负责解决 minimap2 或 UMI 阶段的其他瓶颈。

## 验证记录

版本检查覆盖：两个 Bash 入口接受 2.2.3，拒绝 2.2.2/3.0.3；Python 后端在组装前
拒绝旧版本并保留实际版本记录。候选选择、BAM 标签/顺序、矩阵和嵌入运行时沿用既有回归测试。
本分支没有用户真实异常 FASTQ/BAM，实际 GRCh38 大样本速度与结果差异由服务器对照验证。
