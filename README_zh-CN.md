[![English](https://img.shields.io/badge/Language-English-2563eb)](README.md)
[![中文](https://img.shields.io/badge/语言-中文-0f766e)](README_zh-CN.md)

# Flora 私有开发仓库

本仓库保存 Flora 的私有源码、测试、打包工具和发行文档。不要将本仓库镜像或推送到公开的 `brilliantlee2/Flora` 仓库。

公开仓库只保留用户文档，编译后的压缩包通过 GitHub Releases 发布。

## 克隆私有仓库

你必须具有私有仓库 `brilliantlee2/flora_private` 的访问权限。使用 HTTPS 克隆：

```bash
git clone https://github.com/brilliantlee2/flora_private.git
cd flora_private
```

或者在 GitHub 账户中添加 SSH key 后使用 SSH：

```bash
git clone git@github.com:brilliantlee2/flora_private.git
cd flora_private
```

由于这是私有仓库，GitHub 可能要求浏览器授权、personal access token 或已授权的
SSH key。修改过的 crate 已保存在 `vendor/` 中，不需要再执行 Git submodule 命令。

## 源码目录

```text
Cargo.toml / Cargo.lock       Rust 包和锁定依赖
src/                          Flora 与内置 Glycine Rust 源码
vendor/                       修改过的 rust-htslib 和 edlib_rs
scripts/                      Python 分析、QC、作图和报告源码
run_all.sh                    私有单物种回归基准
run_all_mixed_species.sh      私有 mixed-species 回归基准
environment.yml               私有构建与测试环境
environment.runtime.yml       公开二进制运行环境
packaging/                    发行打包工具
tests/                        Rust/Python/发行测试
docs/repository-templates/    公开/私有 GitHub 文档
```

## 创建私有构建环境

```bash
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export PYTHONUTF8=1

conda env create -f environment.yml
conda activate flora
```

私有环境包含 Python 3.11、Rust/Cargo、GCC/G++、Clang/libclang、CMake、samtools、minimap2、bedtools 和全部 Python 依赖。

构建前验证：

```bash
python --version
rustc --version
cargo --version
cmake --version
gcc --version
clang --version
find "$CONDA_PREFIX" -name 'libclang.so*' | head
```

如果 Bindgen 找不到 libclang：

```bash
export LIBCLANG_PATH="$(dirname "$(find "$CONDA_PREFIX" -name 'libclang.so*' -print -quit)")"
```

## 开发构建与测试

```bash
cargo build --release --locked
cargo test --locked
bash -n run_all.sh run_all_mixed_species.sh packaging/build_binary_release.sh
python -m unittest discover -s tests -v
```

验证 CLI：

```bash
./target/release/flora --version
./target/release/flora glycine --help
./target/release/flora analyze --help
bash run_all.sh -h
```

生成单二进制公开发行包前，`flora run --help` 和
`flora run-mixed --help` 也必须成功。如果任意命令不可用，打包脚本会主动停止。

## 版本更新清单

1. 修改 `Cargo.toml` 中的版本号。
2. 执行 `cargo check` 或 `cargo build --release` 更新 `Cargo.lock`。
3. 更新公开 README 模板中的版本号和下载链接。
4. 添加 `docs/repository-templates/public/RELEASE_NOTES_vX.Y.Z.md`。
5. 打包前运行全部测试。

Flora 使用语义化版本，例如 `0.1.0`、`0.1.1` 和 `0.2.0`。

## 构建 Linux x86_64 发行包

请在 Linux x86_64 上构建，最好使用目标 HPC。不要在 macOS 上构建 Linux
压缩包。下面是从新克隆仓库开始的完整流程。

### 1. 进入仓库并检查主机

```bash
cd /path/to/flora_private

uname -s
uname -m
ldd --version | head -n 1
```

预期：

```text
Linux
x86_64
```

### 2. 创建完整构建环境

```bash
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export PYTHONUTF8=1

conda env create -f environment.yml
conda activate flora
```

如果环境已经存在，根据仓库中的环境文件更新：

```bash
conda env update -n flora -f environment.yml --prune
conda activate flora
```

检查必要工具：

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

Python 必须显示为 `3.11.x`。

### 3. 检查源码与测试

```bash
cargo metadata --locked --no-deps --format-version 1 >/dev/null
cargo test --locked
bash -n run_all.sh
bash -n run_all_mixed_species.sh
bash -n packaging/build_binary_release.sh
bash -n packaging/refresh_binary_release_metadata.sh
python -m unittest discover -s tests -v
```

如果只需要开发构建：

```bash
cargo build --release --locked
./target/release/flora --version
./target/release/flora glycine --help
./target/release/flora analyze --help
```

### 4. 构建并打包公开二进制发行包

打包脚本会自行执行锁定版本的 release 构建，因此这里不需要先单独执行
`cargo build`：

```bash
PYTHON_BIN="$CONDA_PREFIX/bin/python" \
  bash packaging/build_binary_release.sh
```

如果需要强制指定最高 glibc 需求：

```bash
FLORA_MAX_GLIBC=2.34 \
PYTHON_BIN="$CONDA_PREFIX/bin/python" \
  bash packaging/build_binary_release.sh
```

不要将 glibc 上限设置得低于构建主机和链接库所需的版本。检测到的需求会记录在
`BUILD_INFO.txt` 中。

### 5. 查看输出文件

```text
dist/Flora-<version>-linux-x86_64.tar.gz
dist/Flora-<version>-linux-x86_64.tar.gz.sha256
```

脚本只构建 `flora`，将其放在压缩包根目录，将白名单中的 Python 脚本编译为确定性
CPython 3.11 字节码，复制运行环境和公开文档，清理扩展属性，并自动生成与验证 SHA256。

只有 `flora run --help` 和 `flora run-mixed --help` 都成功时才允许打包。在迁移尚未完成的
阶段，这个失败是预期行为，用于防止发布不完整的二进制包。

## 检查发行包

```bash
ARCHIVE="dist/Flora-0.1.0-linux-x86_64.tar.gz"

tar -tzf "$ARCHIVE" | head -50
tar -tzf "$ARCHIVE" | \
  grep -E '(^|/)(src|vendor|tests)/|Cargo\.(toml|lock)$|run_all|\.(rs|py|sh)$' && {
    echo "ERROR: source file found in public archive" >&2
    exit 1
  } || true
```

独立解压测试：

```bash
rm -rf /tmp/flora-release-test
mkdir -p /tmp/flora-release-test
tar -xzf "$ARCHIVE" -C /tmp/flora-release-test
cd /tmp/flora-release-test/Flora-0.1.0-linux-x86_64

file flora
ldd flora
./flora --version
./flora glycine --help
./flora run --help
./flora run-mixed --help
cat PYTHON_ABI.txt
cat BUILD_INFO.txt
```

`file flora` 应显示 ELF x86-64 可执行文件且包含 `stripped`。如果显示
`not stripped`，不要发布该压缩包。

可执行文件数量必须正好为1：

```bash
find . -type f -perm /111 -print
```

唯一结果应当是 `./flora`。

## 生成 SHA256

Linux：

```bash
cd /path/to/Flora
sha256sum dist/Flora-0.1.0-linux-x86_64.tar.gz \
  > dist/Flora-0.1.0-linux-x86_64.tar.gz.sha256
```

验证：

```bash
cd dist
sha256sum -c Flora-0.1.0-linux-x86_64.tar.gz.sha256
```

## 发布到公开 GitHub 仓库

从私有仓库更新公开 README：

```bash
cp docs/repository-templates/public/README.md /path/to/Flora-public/README.md
cp docs/repository-templates/public/README_zh-CN.md /path/to/Flora-public/README_zh-CN.md
```

只提交公开文档。`.tar.gz` 和 `.sha256` 通过 GitHub Release 附件上传，不要提交进 Git 历史。

使用 GitHub CLI：

```bash
gh release create v0.1.0 \
  dist/Flora-0.1.0-linux-x86_64.tar.gz \
  dist/Flora-0.1.0-linux-x86_64.tar.gz.sha256 \
  --repo brilliantlee2/Flora \
  --title "Flora v0.1.0" \
  --notes-file docs/repository-templates/public/RELEASE_NOTES_v0.1.0.md
```

## 保密规则

- 私有开发仓库必须保持 private。
- 不要向公开仓库推送 `src/`、`scripts/`、`vendor/`、Cargo manifest、测试或打包代码。
- 每次上传前都要检查压缩包。
- Python 字节码是混淆，不是加密，仍可能被逆向。
- 二进制与平台相关，不同操作系统和架构需要分别发布。
