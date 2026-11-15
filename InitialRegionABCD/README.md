# All-admissible-initial-state development ABCD

This version preserves the earlier fixed-start experiment unchanged. It uses
all states of the frozen free-space abstraction (excluding any `obs` labels)
as possible initial states. The original cached TS initial state is ignored.

Both backends compute the label-consistent region

`X_win = {x admissible : (x, delta(q_init, label(x))) is product-winning}`.

Spot is invoked with global solving enabled. Its Boolean return value for the
superinitial state is not the issuance criterion. Each winning candidate start
is projected and all projected starts are independently verified. Büchi uses
the full product winning set but likewise tests only correctly initialized
product states. Issuance requires a nonempty projected region, universal
verification, and controller export. There is no check of `any(raw_product_W)`.

The source models, candidate ranking, SAT filtering and saturation are fixed.
The 120-second unique-formula deadline, parity arena caps and 1,200-step finite
simulation horizon are retained. Unknown regions from computational abstentions
remain null; they are never treated as mathematically empty.

```sh
python3 InitialRegionABCD/run.py --dionysos-root . \
  --source Results/task-specific-temperature-prompt-v2 \
  --plant-cache Results/deepseek-development-abcd-v1/plant.jls \
  --self-test
```

Outputs default to `Results/deepseek-development-abcd-all-initial-v3/`:

- `regions/*.json`: exact abstract state IDs and sizes, one per unique formula.
- `controllers/*.jls`: verified policy, region and bounded representative trace.
- `merged-cases.json` and `.csv`: all 200 prompt/budget cases, implication,
  issuance, loss, abstention reason and region size (null if unknown).
- `joint-risk.md` and `.csv`: requested risk table. Conditional diagnostic is
  losses divided by issuance; joint risk divides by 40; CRC correction is
  `(losses + 1) / 41`.
- `monotonicity.json`: all 160 adjacent prompt pairs, loss and exact-region
  inclusion checks, saturation identity checks and any violations.

The correction is a development statistic, not a held-out CRC guarantee.
Computational abstention may cause operational loss monotonicity violations;
the audit reports these rather than imposing monotonicity after computation.
No model API calls are made.
