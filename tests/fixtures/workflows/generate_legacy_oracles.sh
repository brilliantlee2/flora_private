#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
ENV_DIR="${ROOT}/.conda-env"
VERSIONS_FILE="${SCRIPT_DIR}/tool_versions.tsv"
REFRESH_ENVIRONMENT=0
LEGACY_COMMIT="9f41ae3"

usage() {
  cat <<'EOF'
Usage: bash tests/fixtures/workflows/generate_legacy_oracles.sh [--refresh-environment]

Regenerate frozen legacy workflow outputs using only this worktree's .conda-env.
The command refuses any tool-version mismatch. Use --refresh-environment only
when intentionally reviewing and accepting a new pinned environment; this
rewrites tool_versions.tsv as part of the reviewable oracle diff.
EOF
}

if [[ "${1:-}" == "--refresh-environment" ]]; then
  REFRESH_ENVIRONMENT=1
  shift
fi
if [[ $# -ne 0 ]]; then
  usage >&2
  exit 2
fi

for tool in python python3 cargo rustc minimap2 samtools bedtools; do
  if [[ ! -x "${ENV_DIR}/bin/${tool}" ]]; then
    echo "[oracle] missing pinned tool: ${ENV_DIR}/bin/${tool}" >&2
    exit 2
  fi
done
export PATH="${ENV_DIR}/bin:/usr/bin:/bin"
export MPLCONFIGDIR="${ROOT}/target/oracle-matplotlib"
export LC_ALL=C
export LANG=C
export PYTHONHASHSEED=0
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

if ! git -C "${ROOT}" diff --quiet "${LEGACY_COMMIT}" -- \
  Cargo.toml Cargo.lock src scripts main.py args_parser.py utils.py \
  run_all.sh run_all_mixed_species.sh; then
  echo "[oracle] legacy workflow sources differ from ${LEGACY_COMMIT}; refusing regeneration" >&2
  exit 2
fi

TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/flora-legacy-oracles.XXXXXX")"
trap 'chmod -R u+w "${TMP_ROOT}" 2>/dev/null || true; rm -rf "${TMP_ROOT}"' EXIT
CURRENT_VERSIONS="${TMP_ROOT}/tool_versions.tsv"

python_packages="$(${ENV_DIR}/bin/python - <<'PY'
import importlib.metadata

packages = {
    f"{distribution.metadata['Name']}=={distribution.version}"
    for distribution in importlib.metadata.distributions()
    if distribution.metadata.get("Name")
}
print(";".join(sorted(packages, key=str.casefold)))
PY
)"
flora_version="$(${ROOT}/target/release/flora --version 2>/dev/null || printf 'flora 0.1.0 (pre-build)')"
cargo_lock_sha="$(shasum -a 256 "${ROOT}/Cargo.lock" | awk '{print $1}')"
architecture="$(uname -m)"
os_version="$(sw_vers -productName) $(sw_vers -productVersion) ($(sw_vers -buildVersion))"
darwin_version="$(uname -r)"
libc_version="$(otool -L /usr/bin/env | awk '/libSystem/{sub(/^[[:space:]]*/, ""); print; exit}')"
{
  printf 'tool\tversion\n'
  printf 'python\t%s\n' "$(${ENV_DIR}/bin/python --version 2>&1)"
  printf 'python_packages\t%s\n' "${python_packages}"
  printf 'flora_build\t%s;legacy_commit=%s;Cargo.lock.sha256=%s;%s;%s\n' \
    "${flora_version}" "${LEGACY_COMMIT}" "${cargo_lock_sha}" \
    "$(${ENV_DIR}/bin/rustc --version)" "$(${ENV_DIR}/bin/cargo --version)"
  printf 'minimap2\t%s\n' "$(${ENV_DIR}/bin/minimap2 --version)"
  printf 'samtools\t%s\n' "$(${ENV_DIR}/bin/samtools --version | head -1)"
  printf 'bedtools\t%s\n' "$(${ENV_DIR}/bin/bedtools --version)"
  printf 'libc\t%s;Darwin=%s\n' "${libc_version}" "${darwin_version}"
  printf 'architecture\t%s\n' "${architecture}"
  printf 'os\t%s;Darwin=%s\n' "${os_version}" "${darwin_version}"
} > "${CURRENT_VERSIONS}"

if [[ ! -f "${VERSIONS_FILE}" ]] || ! cmp -s "${VERSIONS_FILE}" "${CURRENT_VERSIONS}"; then
  if [[ "${REFRESH_ENVIRONMENT}" -ne 1 ]]; then
    echo "[oracle] pinned environment mismatch; refusing regeneration" >&2
    diff -u "${VERSIONS_FILE}" "${CURRENT_VERSIONS}" >&2 || true
    echo "[oracle] rerun with --refresh-environment only for an intentional environment refresh" >&2
    exit 2
  fi
  chmod u+w "${VERSIONS_FILE}" 2>/dev/null || true
  cp "${CURRENT_VERSIONS}" "${VERSIONS_FILE}"
fi

cd "${ROOT}"
cargo build --release --locked --bins

RUST_STAGES=(
  flora generate_26bp_whitelists prepare_read_tags add_cb_ur_tags assign_genes
  add_gene_tags cluster_umis_allbam cell_umi_gene_table gene_expression
  assign_transcripts isoform_expression rna_qc_metrics read_qc_summary
)
for stage in "${RUST_STAGES[@]}"; do
  if [[ ! -x "${ROOT}/target/release/${stage}" ]]; then
    echo "[oracle] missing expected Rust release stage: target/release/${stage}" >&2
    exit 2
  fi
done

GUARD_BIN="${TMP_ROOT}/guard-bin"
mkdir -p "${GUARD_BIN}"
cat > "${GUARD_BIN}/python3" <<'EOF'
#!/usr/bin/env bash
case "${1:-}" in
  */main.py|*/generate_26bp_whitelists.py|*/prepare_read_tags.py|*/add_cb_ur_tags.py|\
  */assign_genes.py|*/add_gene_tags.py|*/cluster_umis_allbam.py|\
  */cell_umi_gene_table.py|*/gene_expression.py|*/assign_transcripts.py|\
  */isoform_expression.py|*/rna_qc_metrics.py|*/rna_qc_metrics_mixed.py|\
  */read_qc_summary.py)
    echo "[oracle] forbidden Python fallback selected: $1" >&2
    exit 97
    ;;
esac
exec "${FLORA_ORACLE_PYTHON}" "$@"
EOF
cat > "${GUARD_BIN}/bedtools" <<'EOF'
#!/usr/bin/env bash
if [[ "${FLORA_FORCE_BEDTOOLS_FAILURE:-0}" == "1" ]]; then
  echo "[oracle] forced bedtools failure" >&2
  exit 86
fi
exec "${FLORA_ORACLE_BEDTOOLS}" "$@"
EOF
chmod +x "${GUARD_BIN}/python3" "${GUARD_BIN}/bedtools"
export FLORA_ORACLE_PYTHON="${ENV_DIR}/bin/python"
export FLORA_ORACLE_BEDTOOLS="${ENV_DIR}/bin/bedtools"
export PATH="${GUARD_BIN}:${ENV_DIR}/bin:/usr/bin:/bin"

STAGING="${TMP_ROOT}/staging"
mkdir -p "${STAGING}"

run_case() {
  local workflow="$1"
  local scenario="$2"
  local expected_status="$3"
  shift 3
  local case_dir="${STAGING}/${workflow}/${scenario}"
  local output_dir="${case_dir}/output"
  local entrypoint="${ROOT}/run_all.sh"
  [[ "${workflow}" == "mixed" ]] && entrypoint="${ROOT}/run_all_mixed_species.sh"
  mkdir -p "${output_dir}"
  {
    printf 'stage\trelease_binary\tsha256\n'
    for stage in "${RUST_STAGES[@]}"; do
      printf '%s\t%s\t%s\n' \
        "${stage}" \
        "${ROOT}/target/release/${stage}" \
        "$(shasum -a 256 "${ROOT}/target/release/${stage}" | awk '{print $1}')"
    done
  } > "${case_dir}/rust_stage_selection.log"
  set +e
  bash "${entrypoint}" "$@" --out-dir "${output_dir}" > "${case_dir}/workflow.log" 2>&1
  local status=$?
  set -e
  printf 'status\t%d\n' "${status}" > "${case_dir}/exit_status.tsv"
  if [[ "${status}" -ne "${expected_status}" ]]; then
    echo "[oracle] ${workflow}/${scenario}: expected exit ${expected_status}, got ${status}" >&2
    tail -80 "${case_dir}/workflow.log" >&2
    exit 1
  fi
  if grep -Fq "forbidden Python fallback" "${case_dir}/workflow.log"; then
    echo "[oracle] ${workflow}/${scenario}: a Rust release stage used Python fallback" >&2
    exit 1
  fi
  if [[ "$(wc -l < "${case_dir}/rust_stage_selection.log")" -ne $(( ${#RUST_STAGES[@]} + 1 )) ]]; then
    echo "[oracle] ${workflow}/${scenario}: incomplete Rust release stage audit" >&2
    exit 1
  fi
}

for workflow in single mixed; do
  fixture="${SCRIPT_DIR}/${workflow}"
  common=(
    --barcode-list-10bp "${fixture}/barcodes_10bp.txt"
    --ref-dir "${fixture}/ref"
    --isoform-gtf "${fixture}/ref/isoforms.gtf"
    --sample-id "fixture_${workflow}"
    --threads 1 --cluster-threads 1 --exp-cells 1 --pair-min 1 --min-q 0
    --gene-assign-mapq 0 --transcript-assign-mapq 0
  )
  run_case "${workflow}" light 0 \
    --fastq "${fixture}/raw_reads.fastq" --min-len 100 --umi-len 0 \
    --glycine-err 0,0 --glycine-shift 20,20 "${common[@]}"
  run_case "${workflow}" skip_glycine 0 \
    --skip-glycine --full-length-fastq "${fixture}/reads.fastq" "${common[@]}"
  run_case "${workflow}" skip_isoform 0 \
    --skip-glycine --skip-isoform --full-length-fastq "${fixture}/reads.fastq" "${common[@]}"
  run_case "${workflow}" upstream_only 0 \
    --skip-glycine --upstream-only --full-length-fastq "${fixture}/reads.fastq" "${common[@]}"
  run_case "${workflow}" full 1 \
    --skip-glycine --full-output --full-length-fastq "${fixture}/reads.fastq" "${common[@]}"

  stale_dir="${STAGING}/${workflow}/stale_output/output"
  mkdir -p "${stale_dir}"
  printf 'legacy stale marker\n' > "${stale_dir}/stale.marker"
  run_case "${workflow}" stale_output 0 \
    --skip-glycine --full-length-fastq "${fixture}/reads.fastq" "${common[@]}"

  malformed_dir="${STAGING}/${workflow}/malformed_input"
  mkdir -p "${malformed_dir}"
  printf '@malformed\nACGT\n+\nIII\n' > "${malformed_dir}/malformed.fastq"
  run_case "${workflow}" malformed_input 1 \
    --skip-glycine --full-length-fastq "${malformed_dir}/malformed.fastq" "${common[@]}"

  export FLORA_FORCE_BEDTOOLS_FAILURE=1
  run_case "${workflow}" forced_failure 86 \
    --skip-glycine --full-length-fastq "${fixture}/reads.fastq" "${common[@]}"
  unset FLORA_FORCE_BEDTOOLS_FAILURE
done

for workflow in single mixed; do
  oracle_dir="${SCRIPT_DIR}/${workflow}/oracles"
  chmod -R u+w "${oracle_dir}" 2>/dev/null || true
  mkdir -p "${oracle_dir}"
  rsync -a --delete "${STAGING}/${workflow}/" "${oracle_dir}/"
  "${ENV_DIR}/bin/python" - "${ROOT}" "${STAGING}/${workflow}" "${oracle_dir}" <<'PY'
import os
import sys
from pathlib import Path

repository_root = Path(sys.argv[1])
staging_root = Path(sys.argv[2])
oracle_root = Path(sys.argv[3])
for scenario_root in oracle_root.iterdir():
    if not scenario_root.is_dir():
        continue
    staging_scenario = staging_root / scenario_root.name
    for path in scenario_root.rglob("*"):
        if path.is_symlink():
            target = Path(os.readlink(path))
            if not target.is_absolute():
                continue
            try:
                suffix = target.relative_to(staging_scenario)
            except ValueError as error:
                raise SystemExit(f"oracle symlink target escapes scenario: {path} -> {target}") from error
            destination = scenario_root / suffix
            relative_target = os.path.relpath(destination, path.parent)
            path.unlink()
            path.symlink_to(relative_target)
            continue
        if not path.is_file():
            continue
        if path.suffix not in {".html", ".log"} and not path.name.endswith("parameters.tsv"):
            continue
        text = path.read_text(encoding="utf-8")
        text = text.replace(str(repository_root), "<REPOSITORY_ROOT>")
        text = text.replace(str(staging_scenario), "<WORKFLOW_ROOT>")
        path.write_text(text, encoding="utf-8")
PY
  while IFS= read -r -d '' report; do
    gzip -n -9 "${report}"
  done < <(find "${oracle_dir}" -type f -name '*.html' -print0)
done

"${ENV_DIR}/bin/python" - "${ROOT}" "${SCRIPT_DIR}" <<'PY'
import hashlib
import json
import os
import sys
from pathlib import Path

root = Path(sys.argv[1])
fixtures = Path(sys.argv[2])
scenarios = [
    "forced_failure", "full", "light", "malformed_input", "skip_glycine",
    "skip_isoform", "stale_output", "upstream_only",
]
known_outputs = [
    "output/logs/00_glycine.log",
    "output/upstream/matched_reads.fastq.gz",
    "output/upstream/unmatched_reads.fastq.gz",
    "output/upstream/cell_reads.fastq.gz",
    "output/matrix/fixture_{workflow}.isoform_expression.tsv",
    "output/matrix/fixture_{workflow}.read_transcript_assigns.tsv",
    "output/qc/barnyard_qc/per_cell_barnyard.tsv",
    "output/qc/fixture_{workflow}.single_cell_report.html",
]

def digest(path):
    payload = os.readlink(path).encode() if path.is_symlink() else path.read_bytes()
    return hashlib.sha256(payload).hexdigest()

def classify(oracle_path, logical_path):
    name = logical_path.name
    if oracle_path.is_symlink():
        return "byte_exact", {}
    if logical_path.suffix == ".bam":
        tags = []
        if "cb_ur" in name or "tagged" in name:
            tags = ["CB", "CR", "UR", "C5", "C3", "U5", "U3"]
        if ".gn." in name or "tagged" in name:
            tags.append("GN")
        if "tagged" in name:
            tags.append("UB")
        return "parsed_exact", {
            "required_tags": tags,
            "require_index": Path(str(oracle_path) + ".bai").exists(),
            "samtools": ".conda-env/bin/samtools",
        }
    if logical_path.suffix in {".html", ".log"} or name.endswith("parameters.tsv"):
        return "canonicalized_html_log", {
            "canonicalize": ["workspace_root", "elapsed_seconds"]
        }
    if any(token in name for token in ("saturation", "rna_cluster", "per_cell_qc", "rna_qc_metrics")) and logical_path.suffix == ".tsv":
        return "numeric_tolerance", {}
    if logical_path.suffix in {".csv", ".json", ".tsv", ".txt", ".bed", ".gtf", ".gz"}:
        return "parsed_exact", {}
    return "byte_exact", {}

for workflow in ("single", "mixed"):
    oracle_root = fixtures / workflow / "oracles"
    artifacts = []
    for scenario in scenarios:
        scenario_root = oracle_root / scenario
        actual = set()
        for path in sorted(scenario_root.rglob("*")):
            if not (path.is_file() or path.is_symlink()):
                continue
            oracle_relative = path.relative_to(scenario_root).as_posix()
            relative = oracle_relative[:-3] if oracle_relative.endswith(".html.gz") else oracle_relative
            actual.add(relative)
            classification, metadata = classify(path, Path(relative))
            artifact = {
                "scenario": scenario,
                "path": relative,
                "classification": classification,
                "sha256": digest(path),
                **metadata,
            }
            if oracle_relative != relative:
                artifact["oracle_path"] = oracle_relative
            artifacts.append(artifact)
        for template in known_outputs:
            relative = template.format(workflow=workflow)
            if relative not in actual:
                artifacts.append({
                    "scenario": scenario,
                    "path": relative,
                    "classification": "intentionally_absent",
                })
    manifest = {
        "schema_version": 1,
        "workflow": workflow,
        "oracle_source_commit": "9f41ae3",
        "known_deviations": {
            "legacy_full_output_orchestration": {
                "approved": True,
                "legacy_scenario": "full",
                "legacy_expected_exit": 1,
                "legacy_diagnostic_contains": "Missing cell_reads.fastq.gz",
                "future_command": (
                    ["flora", "--full-output"]
                    if workflow == "single"
                    else ["flora", "mixed", "--full-output"]
                ),
                "future_expected_exit": 0,
                "future_required_outputs": [
                    "upstream/matched_reads.fastq.gz",
                    "upstream/unmatched_reads.fastq.gz",
                    "upstream/cell_reads.fastq.gz",
                ],
                "scope": "Legacy run_all orchestration only; do not reproduce in the new Rust workflow.",
            }
        },
        "numeric_tolerance": {"absolute": 1e-9, "relative": 1e-7},
        "canonicalization": {
            "workspace_root": "comparison roots and this repository root only",
            "elapsed_seconds": "legacy step and total elapsed fields only",
        },
        "scenarios": {scenario: {"oracle_root": f"oracles/{scenario}"} for scenario in scenarios},
        "artifacts": artifacts,
    }
    target = root / "tests" / f"artifact_manifest_{workflow}.json"
    target.chmod(0o644) if target.exists() else None
    target.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

for workflow in ("single", "mixed"):
    for path in (fixtures / workflow / "oracles").rglob("*"):
        if path.is_file() and not path.is_symlink():
            path.chmod(0o444)
for target in root.glob("tests/artifact_manifest_*.json"):
    target.chmod(0o444)
(fixtures / "tool_versions.tsv").chmod(0o444)
PY

echo "[oracle] regeneration complete under the pinned environment"
echo "[oracle] review with: git diff --stat && git diff -- tests/artifact_manifest_*.json tests/fixtures/workflows/tool_versions.tsv"
