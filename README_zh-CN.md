[![English](https://img.shields.io/badge/Language-English-2563eb)](README.md)
[![中文](https://img.shields.io/badge/语言-中文-0f766e)](README_zh-CN.md)

# Flora 私有开发仓库

本仓库保存 Flora 的私有源码、测试、打包工具和发行文档。不要将本仓库镜像或推送到公开的 `brilliantlee2/Flora` 仓库。

公开仓库只保留用户文档，编译后的压缩包通过 GitHub Releases 发布。

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

请在 Linux x86_64 上构建，最好使用目标 HPC 或版本较老的兼容 Linux。不要在 macOS 上构建 Linux 压缩包。

```bash
uname -s
uname -m
ldd --version | head -n 1
python --version
```

预期：

```text
Linux
x86_64
Python 3.11.x
```

构建和打包：

```bash
conda activate flora
export LIBCLANG_PATH="$(dirname "$(find "$CONDA_PREFIX" -name 'libclang.so*' -print -quit)")"

cargo test --locked
PYTHON_BIN="$CONDA_PREFIX/bin/python" \
  bash packaging/build_binary_release.sh
```

输出：

```text
dist/Flora-<version>-linux-x86_64.tar.gz
```

打包脚本只构建 `flora` 一个程序，将 `environment.runtime.yml` 作为
发行包的 `environment.yml`，使用公开 README 模板，将白名单中的 Python
脚本编译为确定性 Python 3.11 字节码，并排除流程 Shell、独立阶段程序和
Rust/Python 源码。

## 检查发行包

```bash
ARCHIVE="dist/Flora-0.1.0-linux-x86_64.tar.gz"

tar -tzf "$ARCHIVE" | head -50
tar -tzf "$ARCHIVE" | \
  grep -E '(^|/)(src|vendor|tests)/|Cargo\.(toml|lock)$|\.py$' && {
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
