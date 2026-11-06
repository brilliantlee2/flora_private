# ONT human-sample hotspots: isolated fix and official evidence

Date: 2026-11-06 (Asia/Shanghai).
Branch: `codex/1106-ont-hotspot-fix`.
Base: `codex/1028-ont-v330-performance`, `9025ff6`.

## What changes, and what does not

- Large UMI groups use exact Levenshtein trie searches during directional graph
  traversal. No sample-sized list of all candidate pairs or graph edges is built.
  Already assigned nodes and impossible abundance matches are pruned. Small
  groups (at most 64 unique UMIs) retain exhaustive pairwise clustering.
- The edit threshold remains 3; substitutions, insertions, deletions, mixed UMI
  lengths, abundance rules and descending-count/lexicographic roots are retained.
  Fixed UMI prefixes are not stripped or reinterpreted.
- UMI logs list the eight largest cell/gene groups per chromosome, including read
  counts, distinct UMIs and length range. Groups with at least 5,000 distinct
  UMIs also log start/end and elapsed time. These are diagnostic thresholds,
  not biological filters or read caps.
- With `--ont-jobs >1`, workers refill when tasks complete. Outputs remain in
  annotation order. At most `2 * jobs` tasks/results remain unexported, so a
  slow early chromosome cannot create an unbounded completed-result backlog.
  Backpressure can still pause workers when this bounded window is full.
- StringTie logs identify chromosome, options and scratch directory at start.
  Its version, default options and input data are unchanged.

No change to dependencies, MAPQ, gene/transcript selection, read caps, UMI
sequences, matrices, report definitions or default job/thread budgets is made.
The same `flora-ont330` environment can be used; rebuild/repackage the executable
to include the new Rust code and embedded Python script.

The supplied human log's 4.88-hour chr14 clustering is consistent with the old
candidate-enumeration failure mode, but the triggering cell/gene is not proven
without group data. Exact neighbor search can still be expensive for adversarial
inputs: removing the quadratic pair allocation is not a universal linear-time
guarantee. This update does not claim to fix the separate 12.19-hour genome
alignment or the 18.78-hour StringTie assembly observed in that log.

## What upstream projects report

1. ONT [issue #168](https://github.com/epi2me-labs/wf-single-cell/issues/168)
   reports a single chromosome stuck in StringTie for a week, including with
   wf-single-cell v3.3.0. Thus native ONT can encounter this failure pattern too.
2. ONT [issue #187](https://github.com/epi2me-labs/wf-single-cell/issues/187)
   quotes a user's successful `-j 3` workaround; `-j 2` did not complete for that
   dataset. This is an empirical workaround, not a universally validated fix.
   StringTie 2.2.2 defines `-j` as minimum junction coverage (default 1), and `-v`
   as bundle-processing diagnostics: see its
   [version-pinned source/help](https://github.com/gpertea/stringtie/blob/v2.2.2/stringtie.cpp).
   Raising `-j` changes transcript assembly evidence and can change gene and
   transcript assignments. It is deliberately not enabled by this branch.
3. StringTie [v3.0.2 release notes](https://github.com/gpertea/stringtie/releases/tag/v3.0.2)
   describe a fix for a rare regression involving complex long-read splice
   graphs, alongside other assembly changes. The notes do not prove that the
   older 2.2.2 runtime in these Flora logs has that exact regression. An upgrade
   requires a separate assembly/assignment accuracy and timing comparison.
4. In StringTie [issue #477](https://github.com/gpertea/stringtie/issues/477),
   maintainers ask the reporter to check a recent version and share input data
   for a severe sample-specific slowdown. This supports isolating the failing
   input, not presuming that more threads alone will solve it.

The evidence above establishes similar reported symptoms, not a diagnosis of
the actual human chr2/chr22 bundles. Immunoglobulin regions remain a hypothesis
to inspect, not a reason to remove those reads or split genes arbitrarily.

## Server validation and useful logs

Use the normal command in a **new output directory**, with the same data and
parameters as the old run. No extra flag is needed for the UMI fix. To obtain
StringTie bundle details, add this existing option:

```bash
--stringtie-opts='-c 2 -v'
```

Retain any other StringTie options used in the original run; only append `-v`.
It increases diagnostic output, not assembly thresholds. For initial diagnosis,
keep `--ont-jobs 1` to avoid interleaved bundle messages. The scheduler change
only helps concurrent runs; the supplied human/mouse logs both used one job.
`--save-intermediate` is optional for completed intermediate files, but retains
large BAM/FASTQ files and can increase I/O. Prefer disk-backed local scratch with
`--ont-temp-dir`; do not place large scratch files in a small RAM-backed `/tmp`.

Extract the main diagnostics after or during the run:

```bash
grep -E 'Stage stringtie|cluster_umis\.|\[umi_cluster\]' logs/04_gene.log
```

The new group messages distinguish a large UMI group from BAM collection/writing.
StringTie verbose output helps identify slow bundles for an isolated replay
using the same BAM, chromosome annotation, options and version. Compare decoded
BAM records/tags, assignment TSVs, matrices and counts, not compressed BAM hashes.
Do not apply `-j 3`, a read cap or a version upgrade in a behavior-preserving
timing comparison; test those independently if the default StringTie still stalls.

## Local verification

The test suite includes an exhaustive graph oracle, randomized and dense UMI
groups, mixed lengths, indels, zero counts and active-node removal. Scheduler
tests cover completion-driven refill, ordering, boundedness and error cleanup.
A separate BAM/matrix integration test verifies UB tags against an exhaustive
Python oracle with 64/65-UMI boundary groups and large constant-prefix groups.

Reproduce the optional optimized synthetic timing test:

```bash
cargo test --release --locked --lib constant_prefix_benchmark -- --ignored --nocapture
```

One local macOS run (not real FASTQ) gave:

| Distinct 15bp UMIs | Count per UMI | Old candidate pairs | Old seconds | New seconds |
| ---: | ---: | ---: | ---: | ---: |
| 1,000 | 1 | 499,500 | 0.507456 | 0.000752 |
| 1,000 | 2 | 499,500 | 0.479823 | 0.000359 |
| 2,000 | 1 | 1,999,000 | 2.827254 | 0.002563 |
| 2,000 | 2 | 1,999,000 | 2.743118 | 0.000771 |

All four UMI maps matched the original indexed implementation. These deliberately
pathological shared-prefix inputs demonstrate the targeted improvement only;
they do not predict full-run speedup, real-data memory peaks or StringTie speed.
