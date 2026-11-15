# Original controller-aware calibration

This runner uses the completed DeepSeek V4 Pro R=10, T=0.7 outputs and the frozen
semantic 70/30 split. It makes no model calls. The loss is exactly
`g_k * 1[Phi_k does not imply reference]`, where
`g_k = product(tilde_g_j for j=1..k)`.

The synthesis backend and plant are the existing all-admissible-initial-state
Büchi-first / Spot-parity pipeline. Issuance requires a nonempty projected
admissible winning region, infinite-horizon verification, and controller export.
Per-formula time limit: 120 seconds. Representative simulation horizon: 1,200
steps; a finite simulation is not the proof of satisfaction.

The sampling unit is a semantic class: average the original binary losses over
its rewritten NL prompts, then apply `(sum(class_losses)+1)/(n_classes+1)` to 788
calibration classes. NL-pair-weighted risks are also reported as diagnostics.
Levels are evaluated in increasing order; the first passing level is frozen
before any test synthesis. The 338 test classes use only that fixed threshold.

Synthesis is skipped where implication already proves loss zero, or where an
earlier failed synthesis forces cumulative issuance to zero. Unknown issuance
flags remain null; they are not fabricated. These shortcuts compute the exact
original loss but do not measure overall issuance coverage for skipped cases.

Basic Spot normalization deduplicates equivalent spellings. Verified old
controllers or known empty winning regions are reused only when plant and code
identities match; old computational failures are not imported. New per-formula
results, artifacts, and computational failures are checkpointed and reused on
resume. Source experiments and the R=7/R=10 semantic reports remain unchanged.

From the Dionysos root:

```sh
python3 ControllerCalibration/run.py --root . --alpha 0.20 --workers 2
```

Repeat the identical command to resume. Outputs go to
`Results/deepseek-r10-crc-original-alpha0.20/`. Inspect `progress.json`,
`summary.json`, `selected-threshold.json`, worker logs, per-level case JSONs,
job results, winning regions, and controller artifacts. No incomplete calibration
level is used to select a threshold. Startup failures stop the experiment rather
than being turned into zero losses for every pending case.

This is an empirical assessment on the existing benchmark, including previously
evaluated development examples. The population targeted by the class-average
loss is a uniform semantic class followed by a uniform NL instance in that class.
