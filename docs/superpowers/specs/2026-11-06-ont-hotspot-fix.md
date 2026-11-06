# Behavior-preserving ONT hotspot fixes

## Scope

Base: `codex/1028-ont-v330-performance`, commit `9025ff6`.
Work only on `codex/1106-ont-hotspot-fix`; do not move or edit the original branch.

The human/mouse log audit identified pathological UMI candidate enumeration,
head-of-line blocking in parallel chromosome scheduling, and insufficient
StringTie/UMI hotspot visibility. Runtime causes in the supplied server BAMs
cannot all be established from logs alone.

## Design

1. Replace the materialized UMI candidate-pair list and adjacency graph with
   exact Levenshtein trie queries during directional component traversal.
   Query only unassigned UMIs satisfying the unchanged abundance inequality.
   Preserve descending-count/lexicographic component roots, indels, mixed lengths,
   threshold 3 in the caller, raw UMI sequences, and output mappings.
2. Log the largest cell/gene UMI groups and start/end of large group clustering.
3. Refill chromosome workers on completion, not input-order completion. Keep
   exports ordered and completed results bounded to avoid unbounded disk use.
4. Add per-chromosome StringTie start/progress diagnostics without changing its
   options or version. Document official issue reports, optional verbose mode,
   and why changing `-j` or upgrading StringTie requires a separate experiment.

## Non-goals

No changes to MAPQ, gene/transcript assignment, UMI edit/count rules, read caps,
StringTie junction thresholds, tool dependencies, or QC counting semantics.
No automatic retries with altered scientific parameters, downsampling, or
removal of immunoglobulin reads. No promise that a StringTie algorithmic stall
is fixed by this branch.

## Verification

Use exhaustive all-pair directional clustering as an independent test oracle,
including substitutions, indels, constant prefixes, zero counts, mixed lengths,
and deterministic random DNA. Test scheduler ordering, bounded consumption and
failure cleanup. Run Rust and Python regression suites. Record synthetic timing
separately from real-data performance, which requires server validation.
