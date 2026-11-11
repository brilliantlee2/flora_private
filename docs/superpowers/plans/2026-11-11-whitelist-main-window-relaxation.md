# Main-Window Steepest Threshold Relaxation Plan

**Approved change:** The user selected direct main-window steepest selection, with cross-window instability and broad cutoff ties reported as warnings, not rejection gates. Continue on codex/1111-barcode-steepest-threshold, baseline 6782643. Main and the local NGS project remain unchanged.

**Selection rule:** Keep the oriented theoretical-whitelist intersection, configured count floor, main half-window, log-log secant calculation, earliest-rank tie break and inclusive count cutoff. Use the main-window minimum finite negative slope only. Remove 3x prominence and 1.5x endpoint-drop gates entirely. Insufficient main-window points or no actual decline remain explicit failures, not arbitrary thresholds. Auxiliary window availability/stability and large tied plateaus are advisory only.

**Diagnostics:** Preserve raw prominence/drop metrics without biological significance claims. Add an algorithm identifier, advisory warnings and tie expansion to the JSON/TSV, and print warnings in both production and offline preview. The max-rank-spread option controls warnings only, never the selected rank. Keep legacy as the top-level default; no new flags or dependencies.

- [x] Add failing tests for low prominence, shallow decline, unstable windows, broad tied plateaus, missing auxiliary windows and warning-limit invariance.
- [x] Replace the decision gates with main-window selection and explicit advisory diagnostics; preserve invalid-input/flat-curve rejection.
- [x] Wire warnings into pipeline/preview logs and reports; document changed semantics and false inclusion/exclusion risks.
- [x] Run focused/full Rust/Python tests, legacy parity, preview/production parity and independent review.

**Delivery:** Commit on the existing branch without pushing; real-data validation remains pending.

The former 3x/1.5x guards were local NGS engineering heuristics, not published or calibrated biological thresholds. This revision is an exploratory numerical cutoff, not a validated cell/background separator. Counts tied at the cutoff are all retained, even if this includes a large background plateau; the warning must disclose that risk.

Verification: 115 Rust tests passed across all targets. Python ran 179 tests with no failures and 3 existing external/opt-in skips. All 26 threshold tests ran, including 20 seeded curves checked against an independent main-window reference. A synthetic example with unstable auxiliary windows produced byte-identical production/preview reports and whitelists, with warnings in both logs. Default legacy key outputs still match the original-main executable byte-for-byte. Both shell runners passed syntax checks. Independent focused review found no actionable defects in the diff. No real Linux FASTQ/count table from the failing run was available locally; validating the actual sample's selected point remains pending.
