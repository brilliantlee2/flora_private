# Flora Single-Binary Workflow Design

## Goal

Replace the public `run_all.sh` and `run_all_mixed_species.sh` entrypoints with
`flora run` and `flora run-mixed`, while preserving the current analysis
algorithms, command-line defaults, output files, logs, and result semantics.
The public archive must contain one Rust executable plus required Python
bytecode and report assets. The legacy shell scripts remain only in the private
source repository as regression references.

## Scope

### Included

- Add complete single-species and mixed-species workflow subcommands to `flora`.
- Move argument parsing, path resolution, validation, stage orchestration,
  logging, and external command pipelines from Bash into Rust.
- Refactor existing Rust stage binaries so the workflow can invoke their shared
  library functions directly without changing their algorithms.
- Retain Python-only stages as packaged `.pyc` assets invoked by the Rust
  workflow.
- Remove `run_all.sh`, `run_all_mixed_species.sh`, and independent Rust stage
  executables from the public binary archive.
- Preserve legacy scripts and optional standalone Rust binaries in the private
  source tree for parity testing and developer troubleshooting.

### Excluded

- Rewriting Scanpy, plotting, saturation, or HTML report generation in Rust.
- Bundling Python, minimap2, samtools, bedtools, or reference files.
- Changing barcode correction, UMI clustering, cell assignment, gene
  assignment, transcript assignment, matrix generation, QC, or clustering
  algorithms.
- Guaranteeing that a client-side binary cannot be reverse engineered.

## User Interface

The public interface is:

```text
flora glycine ...
flora analyze ...
flora run ...
flora run-mixed ...
```

`flora run` accepts the same options, defaults, aliases, and validation rules as
the current `run_all.sh`. `flora run-mixed` does the same for
`run_all_mixed_species.sh`, including `--singlet-threshold` and mixed-species
QC. Existing `glycine` and `analyze` behavior remains backward compatible.

Compatibility means that every invocation accepted by a legacy script is
accepted with the same normalized value and every legacy-invalid value remains
invalid. Rust/Clap may additionally accept `--option=value`; this is the only
intentional syntax expansion. Exact help layout and exact parser error wording
may change, but errors remain on stderr with a nonzero exit status. A generated
CLI contract fixture records every option, alias, type, default, conflict,
required condition, and legacy accepted/rejected example for both workflows.

The executable resolves its installation root from `current_exe()`. A packaged
layout with `flora` at the archive root is preferred. `FLORA_HOME` is supported
as an explicit override for development and unusual installation layouts.

## Architecture

### Top-Level CLI

`src/main.rs` dispatches typed Clap subcommands. The current `analyze` parser is
moved behind a reusable command module rather than duplicated. Top-level help
lists all four public modes.

### Workflow Modules

```text
src/workflow/
  mod.rs             public workflow entrypoints
  config.rs          shared options, defaults, aliases, and validation
  paths.rs           input, reference, output, and runtime asset resolution
  command.rs         child processes, pipelines, logs, and exit-status checks
  stages.rs          common single/mixed stage sequence
  single_species.rs  single-species configuration
  mixed_species.rs   mixed QC and barnyard additions
```

Single and mixed workflows use one shared stage engine. Species-specific
behavior is represented by configuration and explicit hooks rather than copied
or condition-heavy orchestration.

### Internal Rust Stages

Logic currently implemented directly in `src/bin/*.rs` is extracted into
library modules with typed configuration structs and `run(...) -> Result<()>`
entrypoints. Standalone binary files become thin private wrappers around those
same functions. `flora run` calls the library functions directly, so the public
archive does not need separate Rust executables and does not duplicate statically
linked dependencies.

The refactor is mechanical: parsing and `main()` move to wrappers; computation
and file-format behavior remain unchanged.

Stage APIs do not call `process::exit`, mutate process-wide environment or the
working directory, install global thread pools, or write through unmanaged
global stdout/stderr. Each accepts explicit paths, thread settings, and log
writers. Per-stage thread pools and large resources are owned by stage-local
scopes and dropped before the next stage. The workflow catches a stage panic at
the boundary, records it as a failed stage, and exits nonzero. Wrapper-versus-
library tests run each Rust stage both ways and compare outputs and diagnostics.

### Python Stages and Assets

Python-only stages continue to run through CPython 3.11. Interpreter selection
is `--python` when supplied, then `FLORA_PYTHON`, then `python3` on `PATH`.
The workflow selects hash-based `.pyc` assets from `<FLORA_HOME>/scripts` and
verifies implementation, major/minor ABI, and imports/version constraints from
the packaged `runtime_manifest.json` before starting expensive analysis. Its
minimum schema contains `python.implementation`, `python.major_minor`, and a
`packages` object mapping each required import name to a PEP 440 version
specifier. Bytecode is
compiled with deterministic public `co_filename` values such as
`flora/scripts/build_report.py`; package-relative imports and `__file__` asset
lookups are tested after archive extraction. HTML templates and the vendored
Plotly JavaScript remain external runtime assets.

The Rust workflow invokes Python with an argument vector, never a shell command
string. Paths containing spaces therefore remain valid and shell injection is
avoided.

### External Tools

minimap2, samtools, and bedtools remain external dependencies. Preflight is
derived from the resolved stage plan, so a tool or Python package used only by a
skipped stage is not required. It validates required inputs before the first
expensive stage. Resolved executable paths and versions are written to a new,
approved `logs/runtime_provenance.tsv` artifact rather than changing the legacy
parameter TSV. Parity canonicalizes installation-dependent executable paths in
this file while comparing tool names and versions exactly. Tests use one pinned
tool environment.

The alignment pipeline is constructed with `std::process::Command` and piped
standard streams:

```text
minimap2 -> samtools view -> samtools sort
```

All child statuses are checked. Failure of any pipeline member fails the stage,
even if a downstream process happens to exit successfully.

The process runner closes unused parent pipe handles immediately and drains all
stderr streams concurrently to bounded in-memory buffers plus stage log files.
Children run in a dedicated Unix process group. SIGINT and SIGTERM are forwarded
to that group; on failure or cancellation, remaining descendants receive TERM,
then KILL after a bounded grace period, and every child is reaped. If an
intentional downstream shutdown causes upstream SIGPIPE, the original downstream
failure is reported as primary. Tests cover stderr larger than pipe capacity,
early downstream exit, a hung child, cancellation, and surviving-process checks.

## Data Flow and Compatibility

The stage order, working directories, filenames, output-directory structure,
and arguments passed to each algorithm match the legacy scripts. Existing skip
switches retain their current meaning. The legacy scripts do not implement a
general resume protocol, so the Rust workflow does not infer completion merely
from an output file's existence and does not introduce a new generic resume
mode. Outputs normally overwritten by the legacy flow remain overwritten. The
one legacy conditional plot fallback is preserved and covered explicitly.

For each run, Rust writes the same parameter TSV and per-stage logs. User-facing
progress messages retain stable step names so existing monitoring remains
useful. The workflow does not silently reuse a truncated BAM or incomplete
output from a failed stage.

Raw input files are never copied into memory by the orchestration layer. Child
process pipes are streamed, so migration itself does not introduce input-size
dependent buffering.

## Error Handling

- Missing tools, inputs, references, Python assets, or incompatible Python ABI
  fail during preflight with a concrete path or command name.
- Every stage error includes the stage name and log path.
- Pipeline children are terminated and reaped when another member fails.
- A failed stage returns a nonzero `flora` exit code immediately.
- Output directories are created deterministically; unsafe implicit cleanup is
  not performed.
- Existing user files are never deleted merely because a stage failed.
- Where a stage implementation supports a temporary output, it writes beside
  the destination and atomically renames only after validation. Stages whose
  legacy tools write directly to fixed outputs continue to overwrite those
  outputs, but downstream stages never run after a nonzero exit. Interrupted
  runs are restarted with the same legacy overwrite behavior; malformed stale
  files are not treated as completion markers.

## Packaging

The public archive contains:

```text
Flora-<version>-linux-x86_64/
  flora
  scripts/*.pyc
  scripts/report_template.html
  scripts/plotly-2.26.0.min.js
  PYTHON_ABI.txt
  runtime_manifest.json
  environment.yml
  requirements.txt
  README.md
  README_zh-CN.md
  THIRD_PARTY_NOTICES.md
  licenses/
```

It does not contain shell workflow sources, Rust sources, Cargo metadata,
standalone Rust stage binaries, tests, vendor sources, or private build scripts.
Release binaries are stripped and archives exclude macOS extended attributes.

Packaging uses an allowlisted manifest rather than copying directories. It
rejects symlinks and archive paths that escape the release root. Rust builds use
path-prefix remapping, and bytecode compilation supplies deterministic public
source filenames. Release validation scans archive names, ELF sections/strings,
and Python code-object filenames for private absolute paths and source-tree
names. The Linux x86-64 build targets the baseline x86-64 instruction set and
glibc 2.17 or newer; dynamic dependencies are inspected and tested on a clean
compatible host.

The present archive is 28 MB compressed and 75 MB unpacked, including about
71 MB of separate Rust executables. Consolidation is expected to reduce the
archive rather than enlarge it. Exact size is an acceptance measurement, not a
hard compatibility requirement.

## Verification

### Unit Tests

- CLI defaults and aliases match each legacy script.
- Reference and runtime asset resolution covers packaged and development paths.
- Command construction preserves every legacy argument.
- Pipeline status handling detects failure in each child position.
- Single/mixed stage plans contain the expected stages in the expected order.

### Integration Tests

- Run legacy `run_all.sh` and `flora run` on the same small single-species input.
- Run legacy `run_all_mixed_species.sh` and `flora run-mixed` on the same small
  mixed-species input.
- Compare every artifact according to the manifest below. Nondeterministic
  metadata is allowed only where listed.
- Compare exit behavior for missing tools, malformed references, and a forced
  failure in each external pipeline member.
- Interrupt at every stage, then restart into the same output directory and
  compare the final artifact set with a clean run. Include stale zero-byte,
  truncated, malformed, and mismatched-parameter files.

### Release Tests

- Archive contains exactly one ELF executable named `flora`.
- Archive contains no `.rs`, `.py`, `.sh`, `Cargo.toml`, `Cargo.lock`, `src`,
  `tests`, or `vendor` content.
- `flora --help`, `flora run --help`, and `flora run-mixed --help` work after
  extraction.
- `runtime_manifest.json` is present, schema-valid, and agrees with the bytecode
  ABI and documented dependency constraints.
- The executable is Linux x86-64 and stripped.
- The archive contains no extended-attribute headers and its SHA256 verifies.
- No private source path appears in ELF strings, Python code-object filenames,
  HTML/JavaScript assets, documentation, or archive metadata.

## Artifact Parity Manifest

The implementation records the legacy single and mixed output trees from fixed
fixtures and classifies every generated path. A test fails when either workflow
produces an unclassified extra or missing path.

- FASTQ, BAM, BAI, BED, PNG, JavaScript, and compressed outputs: byte-for-byte
  comparison when the producing tool is deterministic in the pinned environment;
  otherwise compare decompressed records or image dimensions/pixel hashes with
  only the documented producer metadata removed.
- TSV/CSV matrices, read assignments, barcode mappings, QC tables, parameter
  tables, and cluster tables: headers and row order must match where order is
  part of the legacy output; explicitly order-insensitive tables are parsed and
  compared by their documented key. Integers and strings are exact. Floating
  values use an absolute/relative tolerance of `1e-9` unless the legacy file
  prints fewer decimal places, in which case printed text is exact.
- JSON: recursively compare parsed values with the same numeric rule; object key
  order is ignored and array order is preserved.
- HTML reports: canonicalize only generated timestamp, elapsed-time text, and
  nondeterministic HTML identifier values, then compare embedded payloads,
  section order, labels, and referenced assets exactly.
- Logs: compare stage presence, command arguments, summaries, warnings, and
  errors after removing timestamps, elapsed durations, PIDs, temporary paths,
  tool banner ordering caused by concurrent stderr, and private installation
  roots. Log filenames and stdout/stderr destination are exact.
- `logs/runtime_provenance.tsv` is an approved new artifact. Executable paths
  are canonicalized to command basenames; resolved tool and package versions
  are compared exactly in the pinned test environment.
- Empty optional outputs and absent outputs are part of the manifest for every
  skip/light/full-output combination. Symlink names and targets used by QC are
  compared exactly.
- Failure fixtures compare exit code, completed-stage logs, absence of later
  stage outputs, and presence/validity of any partial output permitted by the
  underlying legacy external tool.

The manifest covers, at minimum, barcode whitelists and knee data/plots,
Glycine statistics/FASTQ outputs, aligned and tagged BAM/BAI/BED artifacts,
read-tag and gene/transcript assignment tables, barcode-to-cell mappings,
gene/isoform matrices, RNA cluster output, QC TSV/JSON/PNG files, saturation,
barnyard outputs for mixed species, parameter/provenance tables, logs, symlinks,
and the final HTML report. The exact filenames are generated from the two legacy
fixtures and checked into the private test suite.

## Acceptance Criteria

- New full-workflow commands finish successfully on representative single and
  mixed test data.
- Scientifically meaningful outputs match the legacy workflows according to the
  Artifact Parity Manifest.
- Existing `glycine` and `analyze` commands remain compatible.
- Public release layout contains one Rust executable and no workflow shell
  source.
- The public archive is no larger than the current 28 MB archive unless a
  measured dependency change is documented and approved.
