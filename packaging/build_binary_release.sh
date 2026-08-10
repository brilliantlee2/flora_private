#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VERSION="$(sed -n 's/^version = "\([^"]*\)"/\1/p' "${ROOT_DIR}/Cargo.toml" | head -n 1)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
RUST_TARGET="x86_64-unknown-linux-gnu"
RELEASE_PLATFORM="linux-x86_64"
RELEASE_NAME="Flora-${VERSION}-${RELEASE_PLATFORM}"
DIST_DIR="${ROOT_DIR}/dist"
STAGE_DIR="${DIST_DIR}/${RELEASE_NAME}"
ARCHIVE_PATH="${DIST_DIR}/${RELEASE_NAME}.tar.gz"
BUILD_BIN="${ROOT_DIR}/target/${RUST_TARGET}/release/flora"
MAX_GLIBC_VERSION="${FLORA_MAX_GLIBC:-}"

PYTHON_RUNTIME_ASSETS=(
  Saturation.py
  barnyard_qc.py
  build_report.py
  generate_knee_plots.py
  rna_cluster_analysis.py
  rna_violin_plot.py
)

die() {
  echo "[ERROR] $*" >&2
  exit 1
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

version_gt() {
  [[ "$(printf '%s\n%s\n' "$1" "$2" | sort -V | tail -n 1)" == "$1" && "$1" != "$2" ]]
}

[[ "$(uname -s)" == "Linux" ]] || die "Linux release packages must be built on Linux"
[[ "$(uname -m)" == "x86_64" ]] || die "Expected x86_64 host, found $(uname -m)"
[[ -n "${VERSION}" ]] || die "Unable to read the package version from Cargo.toml"

for cmd in cargo rustc file tar strings readelf sort find install; do
  require_cmd "${cmd}"
done
require_cmd "${PYTHON_BIN}"

"${PYTHON_BIN}" -c '
import sys
if sys.implementation.name != "cpython" or sys.version_info[:2] != (3, 11):
    raise SystemExit("Flora releases require CPython 3.11")
'

cd "${ROOT_DIR}"
export RUSTFLAGS="${RUSTFLAGS:+${RUSTFLAGS} }-C target-cpu=x86-64 --remap-path-prefix=${ROOT_DIR}=flora-src"
cargo build --release --locked --target x86_64-unknown-linux-gnu --bin flora

[[ -x "${BUILD_BIN}" ]] || die "Missing release executable: ${BUILD_BIN}"
file "${BUILD_BIN}" | grep -q 'ELF 64-bit.*x86-64' || die "flora is not a Linux x86-64 ELF executable"
if file "${BUILD_BIN}" | grep -q 'not stripped'; then
  die "flora still contains symbols; check [profile.release] strip configuration"
fi

PRIVATE_PATHS=("${ROOT_DIR}" "${HOME:-}")
for private_path in "${PRIVATE_PATHS[@]}"; do
  [[ -n "${private_path}" ]] || continue
  if strings "${BUILD_BIN}" | grep -Fq "${private_path}"; then
    die "Private build path is embedded in flora: ${private_path}"
  fi
done

mapfile -t GLIBC_VERSIONS < <(
  readelf --version-info "${BUILD_BIN}" \
    | grep -oE 'GLIBC_[0-9]+\.[0-9]+' \
    | sed 's/^GLIBC_//' \
    | sort -Vu
)
for glibc_version in "${GLIBC_VERSIONS[@]}"; do
  if [[ -n "${MAX_GLIBC_VERSION}" ]] && version_gt "${glibc_version}" "${MAX_GLIBC_VERSION}"; then
    die "flora requires GLIBC_${glibc_version}; release ceiling is GLIBC_${MAX_GLIBC_VERSION}"
  fi
done
REQUIRED_GLIBC_VERSION="$(printf '%s\n' "${GLIBC_VERSIONS[@]}" | tail -n 1)"

rm -rf "${STAGE_DIR}"
mkdir -p "${STAGE_DIR}/scripts" "${STAGE_DIR}/licenses"
install -m 0755 "${BUILD_BIN}" "${STAGE_DIR}/flora"
{
  printf 'target=%s\n' "${RUST_TARGET}"
  printf 'rustc=%s\n' "$(rustc --version)"
  printf 'required_glibc=%s\n' "${REQUIRED_GLIBC_VERSION:-unknown}"
  printf 'python=%s\n' "$("${PYTHON_BIN}" --version 2>&1)"
} > "${STAGE_DIR}/BUILD_INFO.txt"

# Refuse to publish a partial binary while the shell-to-Rust migration is incomplete.
# Checking only the exit status is insufficient: clap may treat an unknown positional
# command followed by --help as valid help for the analyze compatibility path.
TOP_LEVEL_HELP="$("${STAGE_DIR}/flora" --help 2>&1)" \
  || die "failed to inspect the staged flora CLI"
printf '%s\n' "${TOP_LEVEL_HELP}" | grep -Eq '^  flora run([[:space:]]|$)' \
  || die "flora run is not advertised by the CLI; complete the Rust workflow migration before packaging"
printf '%s\n' "${TOP_LEVEL_HELP}" | grep -Eq '^  flora run-mixed([[:space:]]|$)' \
  || die "flora run-mixed is not advertised by the CLI; complete the Rust workflow migration before packaging"

RUN_HELP="$("${STAGE_DIR}/flora" run --help 2>&1)" \
  || die "flora run help failed"
printf '%s\n' "${RUN_HELP}" | grep -Eq '^Usage: flora run([[:space:]]|$)' \
  || die "flora run resolved to a different command; refusing to package an incomplete workflow"

RUN_MIXED_HELP="$("${STAGE_DIR}/flora" run-mixed --help 2>&1)" \
  || die "flora run-mixed help failed"
printf '%s\n' "${RUN_MIXED_HELP}" | grep -Eq '^Usage: flora run-mixed([[:space:]]|$)' \
  || die "flora run-mixed resolved to a different command; refusing to package an incomplete workflow"

PYTHON_SOURCES=()
for asset in "${PYTHON_RUNTIME_ASSETS[@]}"; do
  source_path="${ROOT_DIR}/scripts/${asset}"
  [[ -f "${source_path}" ]] || die "Missing Python runtime asset: ${source_path}"
  PYTHON_SOURCES+=("${source_path}")
done
"${PYTHON_BIN}" "${ROOT_DIR}/packaging/compile_python_assets.py" \
  --output-dir "${STAGE_DIR}/scripts" \
  "${PYTHON_SOURCES[@]}"

install -m 0644 "${ROOT_DIR}/runtime_manifest.json" "${STAGE_DIR}/runtime_manifest.json"
install -m 0644 "${ROOT_DIR}/scripts/report_template.html" "${STAGE_DIR}/scripts/report_template.html"
install -m 0644 "${ROOT_DIR}/scripts/plotly-2.26.0.min.js" "${STAGE_DIR}/scripts/plotly-2.26.0.min.js"
install -m 0644 "${ROOT_DIR}/docs/repository-templates/public/README.md" "${STAGE_DIR}/README.md"
install -m 0644 "${ROOT_DIR}/docs/repository-templates/public/README_zh-CN.md" "${STAGE_DIR}/README_zh-CN.md"
install -m 0644 "${ROOT_DIR}/environment.runtime.yml" "${STAGE_DIR}/environment.yml"
install -m 0644 "${ROOT_DIR}/requirements.txt" "${STAGE_DIR}/requirements.txt"
install -m 0644 "${ROOT_DIR}/THIRD_PARTY_NOTICES.md" "${STAGE_DIR}/THIRD_PARTY_NOTICES.md"
install -m 0644 "${ROOT_DIR}/licenses/Glycine-MIT.txt" "${STAGE_DIR}/licenses/Glycine-MIT.txt"
if [[ -f "${ROOT_DIR}/BC_1536.txt" ]]; then
  install -m 0644 "${ROOT_DIR}/BC_1536.txt" "${STAGE_DIR}/BC_1536.txt"
fi

if find "${STAGE_DIR}" -type l -print -quit | grep -q .; then
  die "Release staging directory contains a symbolic link"
fi
if find "${STAGE_DIR}" -type f -perm /111 ! -name flora -print -quit | grep -q .; then
  die "Release staging directory contains an unexpected executable"
fi
if find "${STAGE_DIR}" -type f \( -name '*.rs' -o -name '*.py' -o -name '*.sh' \
  -o -name Cargo.toml -o -name Cargo.lock \) -print -quit | grep -q .; then
  die "Release staging directory contains source or build files"
fi

if command -v xattr >/dev/null 2>&1; then
  xattr -cr "${STAGE_DIR}" 2>/dev/null || true
fi

rm -f "${ARCHIVE_PATH}" "${ARCHIVE_PATH}.sha256"
TAR_OPTIONS=()
for option in --no-xattrs --no-acls --no-selinux; do
  if tar --help 2>&1 | grep -q -- "${option}"; then
    TAR_OPTIONS+=("${option}")
  fi
done

# The staged tree is already allowlisted and link-free; archive that exact tree.
COPYFILE_DISABLE=1 tar "${TAR_OPTIONS[@]}" \
  -C "${DIST_DIR}" -czf "${ARCHIVE_PATH}" "${RELEASE_NAME}"

mapfile -t ARCHIVE_ENTRIES < <(tar -tzf "${ARCHIVE_PATH}")
for entry in "${ARCHIVE_ENTRIES[@]}"; do
  [[ "${entry}" == "${RELEASE_NAME}" || "${entry}" == "${RELEASE_NAME}/"* ]] \
    || die "Archive entry escapes release root: ${entry}"
  [[ "${entry}" != *'/../'* && "${entry}" != '../'* ]] \
    || die "Archive contains path traversal: ${entry}"
done

(cd "${DIST_DIR}" && sha256sum "${RELEASE_NAME}.tar.gz" > "${RELEASE_NAME}.tar.gz.sha256")
(cd "${DIST_DIR}" && sha256sum -c "${RELEASE_NAME}.tar.gz.sha256")

echo "Created ${ARCHIVE_PATH}"
echo "Checksum ${ARCHIVE_PATH}.sha256"
