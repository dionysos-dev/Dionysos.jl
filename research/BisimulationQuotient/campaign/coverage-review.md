# Coverage review: do the experiments cover the cases of interest?

Checked against the repository, not only against the notes. This review is what added N3, N4 and N7
to the campaign, and **three of its seven gaps have since closed**. The four that remain are all
cheap, and none of them is a cost experiment.

## What is covered

| claim | where | grade |
| :--- | :--- | :--- |
| the graph does not change the certified answer | E1 (4 graphs, 101 probes, 0 disagreements) | proved + measured |
| $\exists$ / $\forall$ from one quotient, with a $\forall$ counterexample | E5 (lasso `1^ω`; script retired) | measured |
| bounded per-node budget forces memory | N5 / `memory_vs_geometry.jl` | analytic |
| the same, for a **single** linear system | N6 | measured |
| validation against the published numbers | E4 (rate 0.940008, 11 slices) | measured |
| complexity redistribution | E2, E8 | measured |
| when memory *cannot* help | N12 (s.m.p. predictor), E8 as negative control | measured |
| **the method against what a practitioner would reach for** | N3 | measured |
| **the object really is a bisimulation** | N4 | measured |
| API-level correctness of the quotient statistics | `test/optim/PCLFBisimulationQuotient/` (1 223 lines) | tested |

## The three gaps that closed

**G1 — no baseline outside the PCLF family.** Every comparison used to be route 1 against route 2, or
against Gol–Lazar–Belta's published certificate, and nothing was ever compared against the uniform
grid abstraction Dionysos ships in
[`src/optim/continuous_systems/UniformGridAbstraction/`](../../../src/optim/continuous_systems/UniformGridAbstraction/).
That is the first question any referee asks — *why not just grid it?* — and the answer was never
given. **Closed by [N3](N3_bisimulation_vs_simulation/plan.md):** the grid is fine on reachability and
collapses on co-safe LTL, which is the qualitative statement, not a cost ratio. A grid yields an
alternating simulation, conservative and one-directional; this construction yields a bisimulation,
exact in both semantics.

**G2 — the conservatism separation had no home.** The one measured certificate-quality result was
written up only in a design note, with no folder and no script. **Closed by
[N7](N7_conservatism_separation/plan.md)** — and the move shrank the claim, because a fair 41-draw
baseline search cut the gap from 0.0447 to 0.0046. That correction is the reason the folder exists.

**G3 — solver soundness was unaddressed.** Three anomalies sat in [`../README.md`](../README.md) and
were never chased: `max_slices` not capping the slice count, deadend states in a bisimulation of a
total system, and a non-monotone JSR bound over nested conic orders — the bug that forced a
retraction. A referee who finds deadend states in an object the paper calls a bisimulation stops
reading. **Closed by [N4](N4_soundness_regressions/plan.md):** all three resolved, 0 violations in
8 367 samples. **One thing remains**: N4's monotonicity table is meant to become a regression test in
[`test/runtests.jl`](../../../test/runtests.jl) and has not been wired in, so nothing stops the
retraction bug returning.

## The gaps still open

### G4 — $r$, the number of observation regions, is never swept

It appears **explicitly** in the complexity bound, $M_{i+1} \le (r+1)^{\sum_j \Delta^j}$, and every
experiment uses whatever the problem happens to define — one region, two, or three. The theorem's own
parameter has never been varied. One sweep on a fixed system would test the bound's shape and cost an
afternoon.

### G5 — specification complexity is never varied

Synthesis is a fixed point over (abstract states × monitor states), so the scLTL formula multiplies
the cost. Every experiment uses a single fixed formula. This matters specifically for the cost story:
route 2's cell-count penalty is multiplied by the monitor size, so E8's 4× becomes 4× *at every
monitor state*. A sweep over formulas of growing monitor size belongs as columns in N1's and N13's
tables rather than as standalone work.

### G6 — `atol` sensitivity is documented but is not an experiment

The README records that `atol = 1e-3` silently drops ~3.7 % of the domain and `1e-4` about 0.4 %.
**Every number in the folder depends on it**, and the graph-invariance result explicitly relies on an
exactness that `atol > 0` breaks. One sensitivity figure — a headline quantity against `atol` over
$10^{-3} \ldots 10^{-5}$ — would immunize the whole paper for an afternoon's work. This is the
cheapest insurance left.

### G7 — no comparison to $l$-complete abstractions

The paper positions itself against them in `rem:l-complete-contrast` and never measures anything.
Lowest priority, and arguably correct: that remark argues the two are incomparable *in mechanism* —
memory tightens an over-approximation there, whereas here the quotient already attains the exact set.
That is a defensible reason not to benchmark, but say it explicitly rather than leaving it silent.

## Where the balance stands now

The original complaint was structural: every planned experiment answered "PCLF versus the common
function it induces", the internal question, and neither external question — *is this better than
gridding?*, *is the object actually a bisimulation?* — had an experiment at all. Both now do, and both
are done. The imbalance that remains is milder and runs the other way: what is measured is
concentrated on two-mode, order-1, planar systems, which is what [N10](N10_dimension_sweep/plan.md)
and [N11](N11_graph_order_sweep/plan.md) exist to widen.

## What to do next

1. **Wire N4's monotonicity table into `test/runtests.jl`** — the one piece of G3 still missing, and
   the cheapest thing on this list.
2. **Write up what is already measured**: N5, N6, N12, N7, and N2's qualitative separation on their
   own matrices.
3. **G6**, `atol` sensitivity — an afternoon, and it underwrites every other number.
4. **N8**, then N1 with a real harness.
5. **N10**, and **N11** — only once N8 says the facet term is live for N10. N11 is independent of it.
6. **G4 and G5** as columns in the tables above rather than as separate experiments.
