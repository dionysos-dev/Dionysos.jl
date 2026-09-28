# N1 — determinise first, or build on the PCLF directly?

Two ways to get a bisimulation quotient out of a path-complete Lyapunov function handed over by an
oracle:

| | |
| :-- | :-- |
| **route 1** | determinise the PCLF into a common Lyapunov function (`build_common_lyapunov`), then run the single-node construction on it |
| **route 2** | run the construction directly on the PCLF, one partition per node |

Both certify the same rate — determinisation preserves it, which is what makes a *cost* comparison
between them meaningful in the first place. The campaign measures what each costs.

## Layout

One folder per example, each holding **its own scripts and its own figures**. There is no shared
`figures/` tree: a figure lives next to the example it belongs to.

```
├── plan.md                           the argument, the results, and what is NOT claimed
├── README.md
├── 1_comparison.jl                   the full matrix: every example × graph × case × size
├── 2_piece_diversity.jl              design screen: which piece geometries separate the routes
├── demo_primal_vs_dual_quotients.jl  [glb|A] — the 2×3 grid + paired timing
│
├── A_identical_pieces/               two_mode_problem, rotated 2n-face template
│   ├── demo_route1_determinise_first.jl
│   ├── demo_route2_direct_on_pclf.jl
│   ├── demo_route{1,2}_{primal,dual}.png      from the demo pair
│   ├── primal_vs_dual_quotients.png           from the grid script
│   └── large_D/, mid_D/                       from 1_comparison.jl
├── B_diverse_pieces/                 a second system chosen so the node pieces differ
│   └── large_D/, mid_D/
├── C_gol_lazar_belta/                Example 3.1 of arXiv:1208.5471, conic order-2 template
│   ├── demo_route1_determinise_first.jl
│   ├── demo_route2_direct_on_pclf.jl           both also solve the co-safe LTL, ∃ and ∀
│   ├── demo_route{1,2}_{quotient,synthesis,verification}.png
│   └── primal_vs_dual_quotients.png
│
├── logs/
└── retired/                          superseded scripts, kept for provenance
```

`A_identical_pieces` and `C_gol_lazar_belta` are the same system *and* the same certificate that the
demo pairs in those folders use — the folders previously called `canonical/` and `glb/` were
duplicates of them and have been merged in.

Run anything from the repository root:

```
julia --project=test research/BisimulationQuotient/campaign/N1_incomplete_speedup/C_gol_lazar_belta/demo_route1_determinise_first.jl
SAMPLES=1 julia --project=test research/.../demo_primal_vs_dual_quotients.jl glb
```

The `demo_*` pairs take the same arguments and must be run with the same ones, or they tile
different regions and the comparison measures coverage instead of efficiency.

## What the three examples are for

**`A_identical_pieces`** is the synthetic two-mode system. It is where route 2 wins on *both*
orientations of the De Bruijn graph, and it is the example that exhibits the dual graph's mechanism
cleanly. Its name is about the certificate, not the system: the node pieces are indistinguishable on
the dual graph (support-function gap below 1e-4) and 2.8 % apart on the primal.

**`B_diverse_pieces`** keeps A's working set and regions and changes only the dynamics, so that the
node pieces come out genuinely different. One difference between A and B, by construction — giving B
its own geometry too would make every cross-example reading ambiguous.

**`C_gol_lazar_belta`** is Example 3.1 of Gol, Ding, Lazar & Belta (arXiv:1208.5471) — their
dynamics, their working set `{‖Lx‖_∞ ≤ 10}`, their three observation regions, their co-safe formula
and their initial point. What is *not* theirs is the certificate: instead of their common Lyapunov
function they are handed a path-complete one from an oracle, and the two routes are compared on it.
Using a benchmark a reader already knows is the point; it is also where the claim is hardest,
because route 2 loses on the dual orientation here.

Its pieces are the most diverse of the three — gap 0.178 on the primal, 0.107 on the dual, against
A's 0.028 and below-1e-4 — which is why its induced common on the complete graph is a union of 17
parts where A's is 3.

## Two traps this folder has already fallen into

**The working set is not the ±5.9 box.** `gol_lazar_belta_problem` carries a box of area 139 that is
a strict subset of their `{‖Lx‖_∞ ≤ 10}` (area 308). A quotient built on the box tiles under half
their domain and leaves two of the three observation regions largely outside it.

**The terminal set must clear the observation regions.** The construction descends until the
terminal level no longer meets R1, R2 or R3 — that is how their `Γ_D = 5.063` is fixed. Forcing a
rung count instead (whether a fixed number, or their band ratio `⌈log(Γ_X/Γ_D)/log(1/γ)⌉` evaluated
at our rate) produces a `D` that overlaps the regions, mislabels every point of the overlap, and
makes the controllable sets incomparable to theirs. Both wrong rules were used in this folder before
being caught; `rungs_clearing_regions` in the grid script is the right one.

## State of the evidence

`plan.md` carries a banner at the top marking which of its rows are superseded; read it before
quoting anything. In short:

A row is void when the terminal set `D` still meets an observation region. The property is **monotone
in the rung count** — more rungs, smaller `D`, so once it clears it stays clear — hence forcing *more*
rungs than needed is harmless and only forcing *fewer* breaks.

| | rungs used / needed | status |
| :-- | :-- | :-- |
| example A, both graphs | 7 / 5 | valid |
| any "no regions" row | any | valid — nothing for `D` to overlap |
| `large_D`, any case | natural stop | valid |
| GLB, primal | 5 / 7 | **void** — superseded by the 7-rung run |
| GLB, primal, 7 rungs | 7 / 7 | valid; one round only, warm-up anomaly unresolved |
| GLB, dual | 9 / 13 | **void**; the 13-rung run is estimated 4–6 h |

No figure in this folder was produced under the corrected rule: the ones that had been were deleted
in the clean-up rather than left to be mistaken for current. Re-run the scripts to regenerate them.

## Measurement

Timings come from paired rounds: both arms run back-to-back inside one round, the ratio taken per
round, the median over rounds. A machine that slows down mid-run then moves both arms of a round
together rather than one arm of a pool. Figures are drawn strictly after the last measurement.
Treat a row whose per-round ratios straddle 1 as undecided rather than as a win — and note that a
single round (`SAMPLES=1`) gives no spread at all, so it is a probe, not a measurement.
