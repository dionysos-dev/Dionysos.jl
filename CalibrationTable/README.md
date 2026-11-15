# Multi-alpha controller-aware calibration and issuance

Uses the completed DeepSeek V4 Pro T=0.7 R=10 collection. No model/API calls.
The existing random split has 300 calibration semantic classes and 826 test
classes. All NL variants and equivalent references remain grouped. This is a
resplit of previously evaluated data, not a pristine held-out test experiment.

Parse, allowed APs and individual SAT precede semantic equivalence clustering.
Clusters are ranked by frequency in the ten queries, breaking ties by the
earliest query. At k, conjoin representatives of the first min(k,m) clusters.
No clusters means abstention, not an empty true conjunction. The saved ranking
and every saved conjunction are checked before running.

For alpha=0.05,0.10,0.15,0.20,0.25, select the first k=1,...,10 for which
`(sum(class_mean_losses)+1)/301 <= alpha`, using calibration only. A semantic
class loss is the mean original binary loss over that class's rewritten prompts.
The original prompt loss is `g_k * !implies_reference`, with cumulative
`g_k=product(tilde_g_j,j<=k)`. Computational failures cause abstention.
If no finite k passes, infinity always abstains, giving zero test risk and
zero issuance. The example table's infinity and historical 16.32% are not inputs.

Once all thresholds are frozen, evaluate complete test issuance through the
largest selected finite k. Calibration may skip synthesis when implication
proves loss zero; test issuance never uses that shortcut. All main table rates
are means over the 826 test classes, with equal weighting of prompts inside
each class. Pair-weighted rates are saved separately as diagnostics.
Thus test risk + guaranteed issuance = issuance, exactly before rounding.

The production scan-based nested Buchi fixed-point solver remains unchanged;
if Spot cannot produce a deterministic Buchi automaton, use the unchanged
native Spot parity game solver. Issuance requires a nonempty projected
admissible initial-state winning region, independent infinite-run verification,
and controller export. Simulation is illustrative and bounded to 1200 steps.

New cases run sequentially, with 600 seconds for translation, process startup,
game solving, independent verification and export. Parity limits: 4M nodes,
60M edges. Process-group RSS guard 3 GiB; free-system-memory floor 15%.
These are polled guards, not temperature sensors or OS allocation quotas.
Buchi uses four Julia threads; separate parity verification uses one.
Processes exit after each new case, followed by 60 seconds idle. After around
30 minutes of active case processing, take a 300-second break instead.

Reuse is proof-based and recorded per formula: compatible resolved reference
decisions, verified old positive controllers (their independent proof remains
valid after the race fix), and resolved decisions from the fixed-solver retry.
No negative Buchi decision from the old race-affected run is imported. Old
120-second or 2M-node computational failures are not imported. The latest
600-second / 4M-node reference failures may be reused as resource abstentions.
Each source hash, plant identity, code hash and new resource policy is frozen.
This is an explicitly checkpointed, proof-reusing evaluation; resource failures
are algorithmic abstentions, not proofs of plant unrealizability.

From the Dionysos root:

```sh
python3 CalibrationTable/run.py --root .
```

Repeat the same command to resume. SIGINT safely stops the active child group;
interrupted cases are not scored as failures. `progress.json` reports activity
or rest, `cases/*/result.json` checkpoints decisions, and `table.md` / `table.tex`
fill as complete levels become available. `summary.json` complete=true means
the whole table is done. Original source results are kept unchanged.

Run tests with `python3 -m unittest discover -s CalibrationTable -p 'test_*.py'`.

To resume an existing frozen experiment with shorter ordinary cooldowns:

```sh
python3 CalibrationTable/resume_with_cooldown.py --root . --cooldown 30
```

Stop the existing runner safely with SIGINT first. This launcher preserves the
original runner and frozen configuration, and records the effective rest policy
in `runtime-policy.json` and the run log. Five-minute rests after 30 minutes of
active processing, memory waits, solver budgets, and checkpoint reuse remain
unchanged. Use the same launcher command for subsequent resumes to retain the
30-second setting; the original command uses the frozen cooldown.

## Single-query baseline

```sh
python3 CalibrationTable/single_query_baseline.py --root . --cooldown 30
```

This evaluates the original `query_index == 1` for each of the same 3,480 test
prompts across 826 semantic classes. It never selects a later valid response
or the most frequent cluster. Invalid, disallowed-AP, and unsatisfiable first
responses abstain. No LLM calls are made.

Both reported risks are joint probabilities, not errors conditional on issuance:
`issued AND NOT equivalent_to_reference` and `issued AND NOT implies_reference`.
Each rate averages prompts within each semantic class and then averages classes
uniformly. The issuance rate is shared between the two baseline rows.

The runner validates cached controller certificates and artifact hashes. It
reuses normalized formulas, proven reference equivalence, and equal Spot level-3
simplifications (which establish language equivalence). Source eligibility and
new-synthesis resource limits match the completed calibration table. The
five-minute rests and memory guards are retained. Missing cases are synthesized
with the same backend and saved individually for safe interruption/resumption.

Outputs are separate under `Results/deepseek-single-query-baseline-cal300-v1`:
`inventory.json`, `test-rows.json`, `summary.json`, `progress.json`, and per-case
checkpoints. Until complete, the summary reports exact lower and upper bounds
from unresolved issuance decisions. Once complete, `table.tex` and `table.md`
contain the combined baseline and NL2RCtrl table. Repeat the same command to
resume; `--prepare-only` validates caches and reports the missing work without
starting new synthesis.

For shorter baseline rests while preserving the frozen experiment and checkpoints:

```sh
python3 CalibrationTable/resume_baseline.py --root . --cooldown 5 --long-rest 30
```

Stop the existing baseline process safely first. This launcher changes only
scheduled post-synthesis waits: five seconds normally, and 30 seconds after
30 minutes of active processing. It records the effective settings in
`runtime-policy.json` and the log. Memory waits and solver budgets remain
unchanged. Use this launcher for subsequent resumes to keep these settings.
