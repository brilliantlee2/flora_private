# Pre-Merge Barcode Whitelist Threshold Experiment

Base: main 291de65. This branch does not include the separate alignment/assembly
performance experiment.

Current rule: `main_window_min_log_log_slope_v2`. This supersedes the strict
NGS-derived guards in commit 6782643. Cross-window instability and broad cutoff
ties now warn without rejecting or moving the main-window cutoff.

Only the observed 3p/5p whitelist cutoff changes when explicitly requested.
Extraction, minQ, theoretical whitelist orientation, correction distance/tie
rules, pair support, UMI pruning, graph merging and downstream analysis remain
unchanged. A barcode rank is not a final cell rank.

## Full Workflow

Default remains `--whitelist-method legacy`. To test the new mode, add:

```bash
--whitelist-method steepest \
--whitelist-slope-window 100 \
--whitelist-slope-min-count 100 \
--whitelist-slope-max-rank-spread 0.10
```

Supported by `flora`, `flora mixed` and `flora analyze`. Use a NEW output
directory for every steepest run. No new Conda packages are needed; rebuild
Flora because Rust logic and embedded scripts changed. The legacy Python
`main.py` fallback rejects steepest rather than maintaining a second algorithm.

## Offline Preview

This reads existing high-quality barcode count tables, not FASTQ/BAM. Use the
same theoretical whitelist files and orientation as the original run:

```bash
RESULT=/path/to/previous_sample_result
PREVIEW=/path/to/new_barcode_preview
./flora barcode-thresholds \
  --counts-3p "$RESULT/upstream/barcode_counts_3p.tsv" \
  --counts-5p "$RESULT/upstream/barcode_counts_5p.tsv" \
  --full-whitelist-3p "$RESULT/barcodes/barcode_3prime_26bp.txt" \
  --full-whitelist-5p "$RESULT/barcodes/barcode_5prime_26bp.txt" \
  --exp-cells 7500 \
  --whitelist-method steepest \
  --out-dir "$PREVIEW"
```

The preview directory must not exist. `--exp-cells` affects only the legacy
comparison, not the steepest search. Add `--no-revcomp-whitelist` only if the
original run used it. Do not use observed `whitelist_3p.csv` or
`whitelist_5p.csv` as theoretical inputs: that truncates the search curve.
The production pipeline and preview share the Rust selector.

Plot the authoritative preview from the source checkout:

```bash
python scripts/generate_knee_plots.py \
  --counts-3p "$RESULT/upstream/barcode_counts_3p.tsv" \
  --counts-5p "$RESULT/upstream/barcode_counts_5p.tsv" \
  --full-whitelist-3p "$RESULT/barcodes/barcode_3prime_26bp.txt" \
  --full-whitelist-5p "$RESULT/barcodes/barcode_5prime_26bp.txt" \
  --threshold-report "$PREVIEW/barcode_thresholds.json" \
  --out-3p "$PREVIEW/knee_plot_3p.png" \
  --out-5p "$PREVIEW/knee_plot_5p.png"
```

Use the same orientation option for plotting. Full workflows generate plots
automatically from this run's report, not a recomputed cutoff. Red is the
applied threshold; dotted lines show available alternatives. Preview does not
update prior expression matrices. Applying a new cutoff requires rerunning
correction/merge and downstream stages; an existing Glycine FASTQ can be reused.

## Algorithm and Failure Policy

- Legacy: Q95 among the highest max(exp_cells,1) observed theoretical barcode
  counts, divided by 20; retain strictly `count > threshold`.
- Steepest: sort the same complete eligible counts independently for 3p/5p.
  Calculate centered log(count)/log(rank) secants. Choose the most negative
  finite main-window slope; equal slopes choose the earliest rank.
  Half/double/fivefold windows provide advisory location diagnostics only. All
  endpoints must strictly exceed the minimum-count guard.
- The former 3x peak/median descent and 1.5x endpoint-drop requirements are
  removed. Their raw metrics remain descriptive, not calibrated evidence of
  cell/background separation. They were local NGS engineering heuristics, not
  validated biological cutoffs for these pre-merge counts.
- `--whitelist-slope-max-rank-spread` is now advisory: a spread above the
  limit generates `cross_window_peak_spread` warnings only. A valid main window
  is sufficient even if auxiliary windows are unavailable; this generates
  `insufficient_auxiliary_windows` instead of rejection.
- Retain all `count >= selected center count` ties. If tie expansion exceeds
  max(main half-window, ceil(spread limit * main rank)), emit
  `broad_tied_plateau` but keep the cutoff and ALL tied barcodes unchanged.
  This may retain a large background plateau and must not be treated as a
  validated set of cells.
- The main window must exist and have an actual finite negative slope.
  Too few eligible points, a completely flat eligible main-window curve or
  invalid parameters still stop the run. Sorted counts cannot legitimately
  produce a positive slope. Do not invent a threshold when no decrease exists.
- Both sides must yield valid main-window thresholds. Otherwise diagnostics
  are saved and neither whitelist is published. There is no silent legacy
  fallback; the new mode still does not change correction/pair/graph rules.
- Legacy keeps its original empty/asymmetric-side behavior. Steepest refuses
  nonempty standalone output directories and existing full-workflow output
  directories to prevent stale successful results masking a failed rerun.
  Directory creation/reservation is atomic; a hidden reservation file remains
  in standalone upstream outputs. Do not share an output directory between jobs.

This retains the main-window secant calculation from the local NGS package's
`engines/cdna/call_cells.py:steepest`, but not its strict rejection guards.
Inputs are **pre-merge single-end barcode read counts**, not post-merge cDNA
read-pair support. It is not
EmptyDrops, an FDR-controlled classifier, or final cell calling. The floor is
a search-domain guard, not proof that lower-count barcodes are empty cells.
Excluding low-support real barcode members or retaining extra noisy barcodes
can change graph groups. Higher median genes may reflect cell selection alone.

## Outputs and Validation

`upstream/barcode_thresholds.json` and `.tsv` (or the preview directory)
record old/new/applied thresholds, retained barcode counts, comparison rules,
the selection-rule version, peak ranks, cross-window stability, tie expansion
and advisory warnings. Warnings are also printed in pipeline/preview stderr.
Light output retains the reports. Observed whitelist CSVs remain sorted as before.

After updating and rebuilding, reuse the failed run's existing count tables
with the offline preview command above; no FASTQ rescan is needed. Use a new
preview directory. A successful preview changes no prior merge or expression
results. Rerun upstream only when the preview cutoff has been inspected.

Compare old/new whitelist overlap, barcode-valid and cell-assigned reads,
group membership, final cell count and expression. Run legacy as the control;
larger cell counts or higher median genes alone do not establish improvement.
Local algorithm/CLI/synthetic parity tests do not replace Linux real-data
validation; no sensitivity gain is claimed before that validation.
