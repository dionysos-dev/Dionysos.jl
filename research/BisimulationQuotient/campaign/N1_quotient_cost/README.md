# N1 — what a bisimulation quotient costs

Two ways to get a bisimulation quotient out of a path-complete Lyapunov function handed over by an
oracle:

| | |
| :-- | :-- |
| **route 1** | determinise the PCLF into a common Lyapunov function (`build_common_lyapunov`), then run the single-node construction on it |
| **route 2** | run the construction directly on the PCLF, one partition per node |

Both certify the same rate — determinisation preserves it, which is what makes a *cost* comparison
between them meaningful in the first place. This campaign measures what each costs.

## The finished experiments are in `research/HSCC2027/`

Three self-contained scripts with their figures, and the reported counts in that folder's README.
They are the same three examples as here, re-measured under the construction's own stopping rule:

| example here | experiment there | system |
| :-- | :-- | :-- |
| `C_gol_lazar_belta` | `experiment1_gol_lazar_belta.jl` | Example 3.1 of Gol, Ding, Lazar & Belta |
| `A_identical_pieces` | `experiment2_primal_and_dual.jl` | the two-mode system, on both De Bruijn graphs |
| `B_diverse_pieces` | `experiment3_diverse_pieces.jl` | the same working set, with pieces that separate |

The dynamics, working sets and regions are identical in both places. **Quote numbers from there, not
from here** — the rows in this folder that were measured under a forced rung count were void, and
have been dropped rather than carried.

## What is left here

| | |
| :-- | :-- |
| `plan.md` | why the two De Bruijn graphs behave oppositely (§2, §7), the per-node threshold that predicts the sign of the gain (§4a), the check that both routes certify the same region (§4b), and what is *not* claimed (§6) |
| `comparison.jl` | the three examples in one matrix, at two ladder depths |
| `paper_experiments.md` | the paper's narrative chain, and which experiment serves which link |

`comparison.jl` reports three things `research/HSCC2027/` does not, which is why it survives:

* the **transition count** — a co-safe solve is a fixed point over transitions, not over cells, so
  this is what predicts downstream cost;
* the **covered volume** of each arm — the two routes do not tile the same region, and volume is the
  only check that they nevertheless certify the same one;
* a **second ladder depth** — which is how the cell ratio is shown not to move with problem size.

```
julia --project=test research/BisimulationQuotient/campaign/N1_quotient_cost/comparison.jl
```

`SAMPLES=n` sets the number of paired rounds, `--figures` draws instead of timing. Figures are
written next to the script, one folder per example and depth, and are not in git.

## Two traps this folder has already fallen into

**Their working set is not the ±5.9 box.** `gol_lazar_belta_problem` carries a box of area 139 that
is a strict subset of their `{‖Lx‖_∞ ≤ 10}` (area 308). A quotient whose outer level is taken from
that box tiles under half their domain and leaves two of the three observation regions largely
outside it. The outer level has to come from the regions instead — which is what
`region_level` does in `research/HSCC2027/experiment1_gol_lazar_belta.jl`.

**The terminal set must clear the observation regions.** The construction descends until the terminal
level no longer meets R1, R2 or R3, and that is how their `Γ_D = 5.063` is fixed. Forcing a rung
count instead — whether a fixed number, or their band ratio `⌈log(Γ_X/Γ_D)/log(1/γ)⌉` evaluated at
our rate — produces a `D` that overlaps the regions, mislabels every point of the overlap, and makes
the controllable sets incomparable to theirs. Both wrong rules were used here before being caught.

## Measurement

Timings come from paired rounds: both arms run back-to-back inside one round, the ratio taken per
round, the median over rounds. A machine that slows down mid-run then moves both arms of a round
together rather than one arm of a pool. Figures are drawn strictly after the last measurement.

Treat a row whose per-round ratios straddle 1 as undecided rather than as a win, and note that a
single round (`SAMPLES=1`) gives no spread at all, so it is a probe and not a measurement.
