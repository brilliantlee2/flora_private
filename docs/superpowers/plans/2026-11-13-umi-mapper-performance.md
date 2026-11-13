# UMI and Mapper Performance Implementation Plan

**Goal:** Remove the measured UMI candidate-materialization bottleneck without changing any UMI-to-root assignments; prepare controlled mapper memory-pool/thread experiments without changing production mapping defaults.

**Base:** codex/1111-alignment-assembly-performance at 0b704e2. Work only in codex/1113-umi-mapper-performance; main, the base branch and the StringTie patch remain unchanged.

**Approved scope:** The user approved continuing the previous report's lossless UMI optimization and single-variable mapper diagnosis. No MAPQ, read filtering, UMI trimming, distance threshold, StringTie or environment changes.

## Design

The legacy implementation first builds all segment candidates, then an adjacency graph, then deterministic directed components. For the measured 15-byte UR / threshold-3 workload ONLY, the new implementation computes the same directed reachability on demand: visit roots in the existing abundance/lexical order, expand a queue, and test only globally unvisited targets. Already visited nodes cannot expose new reachability: they were, or will be, expanded in the same traversal, and all prior components have been fully expanded. This avoids materializing the quadratic candidate set and the graph on the optimized path. All other lengths/thresholds retain the original graph/index implementation: sparse long sequences can have far fewer legacy candidates than exhaustive streaming scans. The fallback retains the legacy memory risks.

Use a precomputed single-word Myers byte-Levenshtein kernel for source sequences of 1..64 bytes; verify exact distances exhaustively/randomly against the existing DP. Retain DP fallback for other lengths. Preserve indels, arbitrary bytes, empty strings, directional abundance rules and component root identity. The graph implementation is retained as a test/replay oracle, not called by normal runs.

Alternatives considered: another short-segment index remains poorly selective on shared fixed motifs; dropping fixed sequences, UMI singletons, or capping groups changes behavior and is rejected. On-demand traversal is preferred because it needs no new dependency and avoids constructing edges that cannot change an already assigned component.

## Tasks

- [x] Establish baseline Rust UMI tests in the new worktree; fix Python interpreter selection without code/environment changes.
- [x] Add regression tests first: original-graph assignment oracle, dense fixed-motif groups, directional branches, ties, indels, mixed lengths and short/long boundaries; observe the materialization assertion failing before implementation.
- [x] Implement the short exact-distance kernel and streamed directed traversal in src/umi_cluster.rs; retain the legacy graph/index for non-qualifying inputs and as independent DP-based test oracle.
- [x] Add a bounded standalone replay harness using the actual modules, exporting per-UR assignments and hashes. All 40 groups / 411,073 group-UR entries match a full independent DP reference; actual old/new subsets (512 per group) also match. Full actual-legacy large-group cases are explicitly skipped rather than risking unbounded memory.
- [x] Extend the independent mapper benchmark with cap-kalloc and thread modes, reusing the existing 200k FASTQ instead of re-sampling the original cDNA. Sequential reverse-order plans, decoded-record equality and input-change protections tested. Linux performance execution is deliberately left to the user.
- [x] Run Rust tests (97 library plus 20 stage tests), release build, workflow/Python tests (170, 8 external-tool skips), mapper tests (25, 2 Linux-only skips with actual mapper smoke enabled), shell syntax checks; independent source correctness review approved the scoped dispatch. Python fixture BAM decoding used an explicitly test-only pysam dispatcher, not stock samtools.
- [x] Record measured local replay results, limitations and Linux commands. Prepare portable read-only UMI replay and mapper diagnostic ZIPs. Local branch commit created; no automatic push/merge. Full production Linux acceptance remains future work, not an achieved claim.

## Acceptance

- Every tested original UR maps to the same root, not merely the same UMI count.
- No graph/candidate-pair materialization on the 15-byte/threshold-3 optimized path; no biological-rule changes. Other inputs retain the original path.
- Same short-kernel exact byte distance as DP, including indels and 64-bit boundary cases.
- Mapper plans vary only the stated factor, hold sort settings fixed during mapper-thread tests, and compare all decoded records and tags.
- Failed, changed-input, timeout or mismatch trials must not emit a successful completion marker.
- The existing StringTie candidate files and production mapper defaults remain byte-identical to base.
