# Optional StringTie 2.2.2 Input Lookup Experiment

This is an opt-in, locally built performance patch, not an upstream release or a
global StringTie upgrade. The workflow's default remains the stock PATH binary.
Select the separate patched executable only through `--stringtie-bin`. Both
executables intentionally report `2.2.2`; use their file names, SHA256 hashes and
`build-provenance.json` to distinguish them.

## Evidence and Scope

The reported chr2 observation was 322/513 seconds before assembly, with the
producer near 99.7% CPU and the consumer mostly idle. This identifies a producer
bottleneck, **not a function-level profile proving `processRead` is responsible**.
The same-start reverse scan is a source-derived hypothesis. No server-input
speedup is claimed here.

[Official issue 464](https://github.com/gpertea/stringtie/issues/464) describes
single-producer limitations and deep-coverage bundle loading, but does not prove
this particular scan is the bottleneck. [Official issue 221](https://github.com/gpertea/stringtie/issues/221)
reports a historical dense-region stall, without establishing a general fix.
Neither issue is a justification for filtering reads, changing `-j`/`-c`, or
splitting a locus. This patch does none of those things.

## Pinned Sources

The complete StringTie source archive includes its bundled gclib and HTSlib:

- [Official v2.2.2 archive](https://codeload.github.com/gpertea/stringtie/tar.gz/refs/tags/v2.2.2)
- SHA256: `17382efe7240c0efb2ac5efd6493882d58797062b3cc713344fb1a0211cb8ce1`
- [Official rlink.cpp](https://github.com/gpertea/stringtie/blob/v2.2.2/rlink.cpp)
- `rlink.cpp` SHA256: `c0623a6f02f130ecfe914f60fcfbf05bf049ca734be220d9f8ede435307c0430`
- [Official rlink.h](https://github.com/gpertea/stringtie/blob/v2.2.2/rlink.h)
- `rlink.h` SHA256: `32478a9ca71ef7b5232a0ad772e75fc81575a4474c7ff8e7e3804c5b53c75d43`

The existing research snapshot's `rlink.cpp` hash matches the official archive.
Archive and file hashes are checked before patching; changed tags, drifted files,
and reapplying the patch fail closed. There is no fuzzy patch application.

The builder also pins [libdeflate commit 9b565afd](https://codeload.github.com/ebiggers/libdeflate/tar.gz/9b565afd996d8b798fc7b94cddcc7cfa49293050),
the commit selected by [upstream's dependency builder](https://github.com/gpertea/stringtie/blob/v2.2.2/htslib/build_lib.sh).
Its archive SHA256 is `8bc23999e454da31b8117db932b859107f1023d144ce89b4165cb819bfacdf2d`.
Upstream's network-fetching build script is deliberately bypassed. zlib, bzip2,
and lzma use the active toolchain's development libraries, not newly installed
packages. Thus source pinning does **not** imply byte-reproducible builds across
different compilers or compression-library versions.

## Exact Candidate Semantics

Only the candidate enumeration in non-merge `processRead` changes. The original
reverse loop remains available in the patched binary as a control. The final
strand, long-read, unitig, `exonmatch` and `deljuncmatch` predicates remain intact.
All read creation, NH minimization, YC/YK weighting, mismatch accounting, junction
updates, pair hash/count mutations and assembly code remain upstream code.

The cache key contains the exon count and **every exon start/end**, strand, and
the computed `longr` value. The cache covers only the contiguous trailing
same-start run that upstream would scan. Each bucket stores original readlist
indices in insertion order; candidates are visited in reverse order. Therefore
the last matching readlist entry wins, including when several entries have the
same key and a later entry fails the deletion/unitig predicate. Hash collisions
cannot cause a match: vector equality distinguishes keys and the original
predicates are rechecked.

Deletion signatures are deliberately **not** keys. Guide-matched junctions can
retain unadjusted coordinates, so the stored deletion signature need not equal
the incoming signature. Several separately created reads can consequently have
the same stored signature. Unitig is also deliberately not a key: upstream can
set it after insertion, including for nonzero negative YK. Both predicates are
checked against the live stored object on every lookup.

Assumptions are limited to pinned upstream input loading: readlist is append-only
and stored exon coordinates, strand and long-read flags do not change before
`getReady`. `evalReadAln` touches the incoming record/guide counts, not those stored
key fields. Assembly-time mutations occur only after producer handoff.

The index belongs to each `BundleData`, not a process-global cache. It indexes new
readlist entries lazily, releases all allocations on a start change, and clears
both at `getReady` (before consumer handoff) and `Clear` (including destruction
and reuse). It retains no read pointers. A shrinking list also invalidates it.
Merge-mode thresholds and enumeration are untouched.

Space is proportional to exon keys and read indices for one same-start run per
loading bundle. This adds memory; it is not a memory-reduction patch. Identical
exons with many distinct deletion/unitig outcomes still require a reverse scan
within that bucket and can remain quadratic. Small runs may be slightly slower.
The patch does not address BAM I/O, guide evaluation, or assembly graph costs.

## Build and Select

Run from this dedicated worktree, with Python 3.9+, make, a C/C++11 compiler, and
zlib/bzip2/lzma development libraries already available. No environment mutation,
package install, system install, or PATH replacement is performed.

```bash
python3 packaging/stringtie_input/build.py \
  --build-root /tmp/stringtie-input-222 --fetch --jobs 4
```

The build root must not exist, must be outside the repository, and must not
contain whitespace (upstream makefile limitation). Downloads require explicit
`--fetch`; offline builds accept `--archive /path/stringtie.tar.gz` and
`--libdeflate-archive /path/libdeflate.tar.gz`. `--prepare-only` verifies, extracts,
and patches without compiling. `--cc`/`--cxx` override activated `CC`/`CXX`.

The builder preserves environment `CPPFLAGS`, `CFLAGS`, `CXXFLAGS`, and `LDFLAGS`,
including Conda include/library prefixes and rpaths. It adds release `-O3` defaults
before user flags and upstream's required flags. Its bundled libdeflate include
directory takes precedence. If a toolchain does not automatically set prefixes,
supply e.g. `CPPFLAGS=-I$CONDA_PREFIX/include` and
`LDFLAGS="-L$CONDA_PREFIX/lib -Wl,-rpath,$CONDA_PREFIX/lib"` explicitly. Never compare
binaries built with different optimization/math flags as evidence of model parity.
Linux/Conda execution is not tested by the local macOS validation.

Outputs under the chosen build root:

- `stringtie-2.2.2-stock`: unmodified upstream model and input implementation.
- `stringtie-2.2.2-same-start`: patched input implementation, unchanged model.
- `build-provenance.json`: source/patch/binary hashes, patched markers, compiler,
  flags, commands, system-library requirements, and the lookup control contract.
- `stock.log`, `same-start.log`, `dependency.log`, and both extracted source trees.

The optional binary accepts a process-lifetime environment control:

```bash
STRINGTIE_INPUT_LOOKUP=linear /tmp/stringtie-input-222/stringtie-2.2.2-same-start -L -o linear.gtf reads.bam
STRINGTIE_INPUT_LOOKUP=indexed /tmp/stringtie-input-222/stringtie-2.2.2-same-start -L -o indexed.gtf reads.bam
```

Unset means indexed **in the patched binary only**. Unknown values fall back to
linear. The stock binary ignores this variable. Set it before process startup;
changing it during a process is unsupported. Use identical inputs and arguments
for all three controls: upstream stock, patched-linear, patched-indexed. Actual
workflow selection is `--stringtie-bin /absolute/path/stringtie-2.2.2-same-start`.
Keep all other assembly options unchanged. No standalone server benchmark helper
is added here; the parent task owns full workflow/chr2/chr22 benchmarking.

## Validation

```bash
STRINGTIE_INPUT_ARCHIVE=/path/stringtie.tar.gz \
STRINGTIE_INPUT_BUILD=/tmp/stringtie-input-222 \
python3 -m unittest discover -s tests -p 'test_stringtie_input*.py' -v
```

The integration tests additionally require `pysam`; they skip explicitly unless
`STRINGTIE_INPUT_BUILD` names a completed build. Ordinary offline tests require no
download and exercise hash rejection, patch drift, build-root constraints, flag
preservation, same-binary control, and 40,000 compiled differential candidate
queries. Initial regressions were run before implementation and failed for the
missing builder/index/control/flag handling.

Integration validation is stronger than a lookup-only comparison:

- A harness links the **real upstream `processRead` and release objects**, renaming
  only CLI `main`. It compares exact full state after every one of 1,080 calls in
  each of short, long, mixed, unitig and no-multimap-correction modes across all
  three controls. Snapshots include exon order, strand/longread/unitig fields,
  read count/NH/length, junction coordinates/guide flags/support/mismatch values,
  pair indices/counts, coverage arrays and bundle counters/status.
- Adversarial inputs include guide-adjusted deletion duplicates followed by a
  query that must update the newest matching read, negative-to-positive YK
  transitions, different strands and mixed read types, secondary/supplementary
  flags, soft clipping, weighted NH/YC, actual mate hash/count updates, start
  changes, separated loci, and chromosome/bundle reuse with repeated positions.
- Actual synthetic BAMs contain 1,080 alignments on chr2/chr22 in six loci.
  Seven assembly scenarios cover short, `-L`, guides, `-e`, `-u`, two worker
  threads, and `--mix`. GTF data rows and `-A` abundance output must match exactly;
  only GTF command-line comment headers are omitted from comparisons.
- A dense BAM has 24,000 records with 12,000 distinct same-start single-exon
  geometries and legal `SEQ=*`. It checks nonempty GTF parity, and then checks
  `--merge` parity. This is a focused geometry stress case, not a representative
  ONT runtime benchmark. Test timings have no hard speed threshold.

Local macOS release builds succeeded with Apple Clang. The final local run
measured stock 0.284 s, same-binary linear 0.277 s, indexed 0.013 s. Small spliced
fixtures were roughly 5-12 ms per run, with no meaningful speed conclusion. These
are local synthetic observations, not the chr2/chr22 server-input result.

The completed local build is `/tmp/flora-stringtie-input-reviewed-1111`; its
`differential-tests.log` records **10 tests passed, zero skipped**, including the
real pair-count mutation checks. The test interpreter was
`/tmp/flora-stringtie223-env/bin/python3` (Python 3.11, pysam 0.24.1). Build logs
also verify propagation of a CPPFLAGS marker and C/C++ `-fPIC` flags. This local
build/log path is evidence from this workstation, not a distributable runtime.

Before adoption, the parent must run the actual server BAMs with unchanged
arguments in all controls, verify output/quantification parity, collect repeated
wall/CPU/RSS measurements, and profile the producer to test the same-start
hypothesis. Retain stock as the default and rollback path.
