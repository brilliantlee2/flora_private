# Pre-Merge Barcode Whitelist Threshold Experiment

Base: main 291de65. This branch does not include the separate alignment/assembly
performance experiment.

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
  Calculate centered log(count)/log(rank) secants. The main half-window gives
  the cutoff; half/double/fivefold windows check location stability. All
  endpoints must strictly exceed the minimum-count guard.
- Require at least two computed windows including the main window, relative
  peak-rank spread <= configured limit, peak descent >=3x median descent, and
  endpoint count drop >=1.5x. Equal slopes choose the earliest rank.
- Retain all `count >= selected center count` ties. Reject if tie expansion
  exceeds max(main half-window, ceil(spread limit * main rank)).
- Both sides must pass in steepest mode. Failure saves diagnostics and stops
  before publishing either whitelist. No silent fallback. Review the curve
  before explicitly choosing legacy or adjusting search parameters.
- Legacy keeps its original empty/asymmetric-side behavior. Steepest refuses
  nonempty standalone output directories and existing full-workflow output
  directories to prevent stale successful results masking a failed rerun.
  Directory creation/reservation is atomic; a hidden reservation file remains
  in standalone upstream outputs. Do not share an output directory between jobs.

This ports the guards from the local NGS package's
`engines/cdna/call_cells.py:steepest`, but the inputs are **pre-merge single-end
barcode read counts**, not post-merge cDNA read-pair support. It is not
EmptyDrops, an FDR-controlled classifier, or final cell calling. The floor is
a search-domain guard, not proof that lower-count barcodes are empty cells.
Excluding low-support real barcode members or retaining extra noisy barcodes
can change graph groups. Higher median genes may reflect cell selection alone.

## Outputs and Validation

`upstream/barcode_thresholds.json` and `.tsv` (or the preview directory)
record old/new/applied thresholds, retained barcode counts, comparison rules,
peak ranks, cross-window stability and failure reasons. Light output retains
them. Observed whitelist CSVs remain sorted as before.

Compare old/new whitelist overlap, barcode-valid and cell-assigned reads,
group membership, final cell count and expression. Run legacy as the control;
larger cell counts or higher median genes alone do not establish improvement.
Local algorithm/CLI/synthetic parity tests do not replace Linux real-data
validation; no sensitivity gain is claimed before that validation.
