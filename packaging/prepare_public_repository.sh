#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST_DIR="${1:-}"

if [[ -z "${DEST_DIR}" ]]; then
  echo "Usage: bash packaging/prepare_public_repository.sh /path/to/Flora-public" >&2
  exit 1
fi

mkdir -p "${DEST_DIR}/licenses"
install -m 0644 "${ROOT_DIR}/docs/repository-templates/public/README.md" "${DEST_DIR}/README.md"
install -m 0644 "${ROOT_DIR}/docs/repository-templates/public/README_zh-CN.md" "${DEST_DIR}/README_zh-CN.md"
install -m 0644 "${ROOT_DIR}/docs/repository-templates/public/RELEASE_NOTES_v0.1.0.md" "${DEST_DIR}/RELEASE_NOTES_v0.1.0.md"
install -m 0644 "${ROOT_DIR}/docs/repository-templates/public/.gitignore" "${DEST_DIR}/.gitignore"
install -m 0644 "${ROOT_DIR}/THIRD_PARTY_NOTICES.md" "${DEST_DIR}/THIRD_PARTY_NOTICES.md"
install -m 0644 "${ROOT_DIR}/licenses/Glycine-MIT.txt" "${DEST_DIR}/licenses/Glycine-MIT.txt"

echo "Prepared public repository files in ${DEST_DIR}"
