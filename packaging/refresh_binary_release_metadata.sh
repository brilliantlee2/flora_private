#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ARCHIVE="${1:-}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

if [[ -z "${ARCHIVE}" || ! -f "${ARCHIVE}" ]]; then
  echo "Usage: bash packaging/refresh_binary_release_metadata.sh /path/to/Flora-<version>-linux-x86_64.tar.gz" >&2
  exit 1
fi
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || {
  echo "[ERROR] Python is required to validate release archives" >&2
  exit 1
}

ARCHIVE="$(cd "$(dirname "${ARCHIVE}")" && pwd)/$(basename "${ARCHIVE}")"

validate_archive() {
  local archive_path="$1"
  "${PYTHON_BIN}" - "${archive_path}" <<'PY'
import pathlib
import sys
import tarfile

archive = pathlib.Path(sys.argv[1])
forbidden_names = {"Cargo.toml", "Cargo.lock", "run_all.sh", "run_all_mixed_species.sh"}
forbidden_suffixes = {".rs", ".py", ".sh"}
with tarfile.open(archive, "r:gz") as handle:
    members = handle.getmembers()
    if not members:
        raise SystemExit("[ERROR] Empty release archive")
    roots = {pathlib.PurePosixPath(member.name).parts[0] for member in members}
    if len(roots) != 1:
        raise SystemExit(f"[ERROR] Archive must contain one root: {sorted(roots)}")
    root = next(iter(roots))
    expected = f"{root}/flora"
    executable_count = 0
    for member in members:
        path = pathlib.PurePosixPath(member.name)
        if path.is_absolute() or ".." in path.parts:
            raise SystemExit(f"[ERROR] Unsafe archive path: {member.name}")
        if member.issym() or member.islnk() or member.isdev():
            raise SystemExit(f"[ERROR] Unsupported archive member type: {member.name}")
        if path.name in forbidden_names or path.suffix in forbidden_suffixes:
            raise SystemExit(f"[ERROR] Source/build file in release: {member.name}")
        if member.isfile() and member.mode & 0o111:
            executable_count += 1
            if member.name != expected:
                raise SystemExit(f"[ERROR] Unexpected executable: {member.name}")
    if executable_count != 1:
        raise SystemExit(f"[ERROR] Expected one executable, found {executable_count}")
PY
}

validate_archive "${ARCHIVE}"
TOP_DIR="$(tar -tzf "${ARCHIVE}" | sed -n '1{s:/$::;p;q;}')"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/flora-release-refresh.XXXXXX")"
trap 'rm -rf "${TMP_DIR}"' EXIT
tar -xzf "${ARCHIVE}" -C "${TMP_DIR}"
STAGE_DIR="${TMP_DIR}/${TOP_DIR}"
[[ -d "${STAGE_DIR}" ]] || { echo "[ERROR] Missing release root: ${STAGE_DIR}" >&2; exit 1; }

install -m 0644 "${ROOT_DIR}/docs/repository-templates/public/README.md" "${STAGE_DIR}/README.md"
install -m 0644 "${ROOT_DIR}/docs/repository-templates/public/README_zh-CN.md" "${STAGE_DIR}/README_zh-CN.md"
install -m 0644 "${ROOT_DIR}/environment.runtime.yml" "${STAGE_DIR}/environment.yml"
install -m 0644 "${ROOT_DIR}/requirements.txt" "${STAGE_DIR}/requirements.txt"
install -m 0644 "${ROOT_DIR}/runtime_manifest.json" "${STAGE_DIR}/runtime_manifest.json"

if command -v xattr >/dev/null 2>&1; then
  xattr -cr "${STAGE_DIR}" 2>/dev/null || true
fi

NEW_ARCHIVE="${ARCHIVE}.new"
rm -f "${NEW_ARCHIVE}"
TAR_OPTIONS=()
for option in --no-xattrs --no-acls --no-selinux; do
  if tar --help 2>&1 | grep -q -- "${option}"; then
    TAR_OPTIONS+=("${option}")
  fi
done
COPYFILE_DISABLE=1 tar "${TAR_OPTIONS[@]}" -C "${TMP_DIR}" -czf "${NEW_ARCHIVE}" "${TOP_DIR}"
validate_archive "${NEW_ARCHIVE}"
mv "${NEW_ARCHIVE}" "${ARCHIVE}"

ARCHIVE_DIR="$(dirname "${ARCHIVE}")"
ARCHIVE_NAME="$(basename "${ARCHIVE}")"
if command -v sha256sum >/dev/null 2>&1; then
  (cd "${ARCHIVE_DIR}" && sha256sum "${ARCHIVE_NAME}" > "${ARCHIVE_NAME}.sha256")
else
  (cd "${ARCHIVE_DIR}" && shasum -a 256 "${ARCHIVE_NAME}" > "${ARCHIVE_NAME}.sha256")
fi

echo "Refreshed ${ARCHIVE}"
echo "Checksum ${ARCHIVE}.sha256"
