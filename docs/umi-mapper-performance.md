# UMI performance candidate and mapper follow-up

Branch: `codex/1113-umi-mapper-performance`.
Base: `codex/1111-alignment-assembly-performance`, `0b704e2`.

## Scope

Only `src/umi_cluster.rs` changes normal execution. StringTie sources/options,
mapper defaults, MAPQ thresholds, barcode selection, UMI extraction, output
formats, expression counting and dependencies remain unchanged.

The fast path applies only when every UR in a cell/gene group is exactly 15
bytes and the existing edit-distance threshold is 3. Other lengths, mixed
lengths and other thresholds use the original indexed-graph implementation.
This restriction avoids regressing sparse long-UMI groups whose legacy
segment index is highly selective. The legacy fallback can still materialize
large candidate sets; it is not a universal memory fix.

## Unchanged assignment rules

- Global byte Levenshtein distance <= 3, including insertions/deletions.
- Directed edge from a to b when count(a) >= 2*count(b)-1.
- Roots ordered by descending abundance, then ascending lexical sequence.
- Directed reachability, including transitive chains, not merely neighbors
  within distance 3 of the root.
- No removal of singleton UMIs, fixed bases, low-MAPQ reads or large groups.

## Implementation

The legacy code enumerates all segment candidates, stores a de-duplication set
and pair vector, verifies distances, stores the graph, and finally traverses it.
On the measured large immunoglobulin groups the candidate index still emits
roughly 80% of all possible pairs.

The new path visits roots in the same order and explores outgoing edges on
demand, testing only globally unvisited targets. A queued target is marked as
visited immediately and will be expanded. A completed component is forward-
closed: otherwise its outgoing traversal would already have discovered the
target. Thus skipping edges into visited nodes cannot lose a new assignment.
The current root is the first-ranked node in its component, exactly as in the
legacy component sort.

Short-sequence distance checks use a precomputed single-word Myers bit-vector
kernel without allocating DP rows for each comparison. It computes global
byte distance, not Hamming distance. Arithmetic wraps explicitly and length 64
uses the highest valid bit. Intermediate score can decrease; early rejection
uses a remaining-text lower bound, not score alone. The independent legacy
DP remains the test oracle.

Auxiliary pair/graph storage on the optimized path is replaced with linear
node, position, queue and result structures. Worst-case distance scanning is
still quadratic. Real-group performance and assignment equivalence must be
measured; no blanket speedup or final-matrix equivalence is implied.

## Validation

The regression test first failed on the original implementation because normal
clustering invoked pair materialization. It passes on the optimized dispatch.
Tests compare full UMI-to-root maps against the original graph/DP implementation.
They cover fixed motifs, abundance direction, lexical ties, transitive paths,
indels, insertion-order changes and mixed/empty inputs. Distance tests cover
all binary strings of length 0..6, seeded random edits, UTF-8 bytes and
63/64/65-byte boundaries. Dispatch tests retain the legacy path for other
lengths and thresholds.

Local tests are not a full Linux FASTQ-to-matrix acceptance run. Use the supplied
independent UMI replay package before a full rerun; its report distinguishes
full groups, deterministic subsets, skipped baselines and mismatches.

## Mapper experiment

The accompanying standalone mapper diagnostic reuses the previous 200,000-read
`sample.fastq.gz`, without re-sampling or re-running upstream stages.

1. Memory mode: cap-kalloc 100m vs 500m; K100M, mapper/sort 32 threads, one I/O
   thread, direct gzip input. Two forward/reverse rounds.
2. Threads mode: mapper 8/16/32 threads; K100M, cap500m, sort fixed at 32 threads,
   one I/O thread, direct gzip input. Run after memory mode, not concurrently.

It records versions, commands, wall/user/system CPU, faults, RSS, thread/wait
and host PSI snapshots, and decoded-record equivalence. No production mapper
parameter has been changed based on this untested hypothesis.

Official parameter definition:
https://github.com/lh3/minimap2/blob/v2.31/minimap2.1

## Build

Use the existing Flora environment. No new Conda package or Rust dependency is
required. The Linux binary must be built on Linux; a local macOS build is not a
Linux deployment artifact.

```bash
conda activate flora-ont330
PYTHON_BIN="$CONDA_PREFIX/bin/python" cargo build --release --locked
```

Existing run arguments and the optional isolated StringTie candidate remain
valid. The new UMI path is automatic for qualifying groups and does not need a
new biological or workflow option. Keep `--save-intermediate` for initial
end-to-end equivalence comparison.
