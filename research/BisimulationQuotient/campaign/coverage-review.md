# Coverage review: do the experiments cover the cases of interest?

Checked against the repository, not only against the notes. Short answer: **the *claims* are now well
covered, the *baselines* and the *soundness* are not**, and there is a structural imbalance — twelve
planned experiments all compare a PCLF against the common function it induces, and none asks whether
the method beats the thing a practitioner would actually reach for instead.

## What is well covered

| claim | where | grade |
| :--- | :--- | :--- |
| the graph does not change the certified answer | E1 (4 graphs, 101 probes, 0 disagreements) | proved + measured |
| $\exists$ / $\forall$ from one quotient, with a $\forall$ counterexample | `augmented_showcase.jl` (lasso `1^ω`) | measured |
| bounded per-node budget forces memory | N2 / `memory_vs_geometry.jl` | analytic |
| the same, for a **single** linear system | N12 | measured |
| validation against the published numbers | E4 (rate 0.940008, 11 slices) | measured |
| complexity redistribution | E2, E8 | measured |
| when memory *cannot* help | N3 (s.m.p. predictor), E8 as negative control | measured |
| API-level correctness of the quotient statistics | `test/optim/PCLFBisimulationQuotient/` (1 223 lines) | tested |

## The gaps, ranked by how much they would hurt

### G1 — No baseline outside the PCLF family *(the biggest gap)*

Every comparison in the folder is route 1 versus route 2, or against Gol–Lazar–Belta's published
certificate. **Nothing is ever compared against uniform grid abstraction**, which Dionysos ships
([`src/optim/continuous_systems/UniformGridAbstraction/`](../../../src/optim/continuous_systems/UniformGridAbstraction/))
and which is the first question any referee or user asks: *why not just grid it?*

The answer is strong and it is never given. A grid yields an **alternating simulation** — an
over-approximation, conservative, one direction of transfer — while this construction yields a
**bisimulation**: exact, both semantics, no conservatism. On the same problem and specification, the
comparison writes itself:

- grid cells vs quotient cells at a grid resolution fine enough to certify anything at all;
- the certified set of each — **ours should be strictly larger**, since the grid is conservative and
  we are exact, and that is a *qualitative* advantage, not a cost ratio;
- build time, and how the grid degrades as the resolution is refined to close the gap.

This is in-repo, cheap, and by far the most persuasive missing experiment for a control audience. See
[N14](N14_bisimulation_vs_simulation/plan.md).

### G2 — The conservatism separation has no home

The one measured certificate-quality result — on E2, four nodes reach 0.86896 (0.02 % above the JSR
lower bound) where a single node at the same template cannot beat 0.873586 (0.55 % above) — is written
up in [`../notes/paper-advantages.md`](../notes/paper-advantages.md) §1 and has **no experiment folder
and no script**. It is one of only three positive results in hand. See
[N13](N13_conservatism_separation/plan.md).

### G3 — Solver soundness is unaddressed, and it undermines every measurement

Three anomalies are recorded in [`../README.md`](../README.md) and never chased:

1. **`max_slices` does not cap the number of slices** — 12 requested, 20 built.
2. **Quotients contain deadend states** (minimum outgoing degree 0) where a bisimulation of a *total*
   system must have none.
3. **The JSR bound was non-monotone over nested conic orders** — the bug that forced a retraction.

The test suite covers the statistics API thoroughly, including a `deadend_states` unit test, but
**nothing asserts these semantic properties on a real quotient**.
[`../notes/complexity-vs-memory.md`](../notes/complexity-vs-memory.md) §8 already concluded that the
first thing to add "is not a benchmark but a regression test asserting that the bound is non-increasing
over nested conic-partition orders" — still not done.

A referee who finds deadend states in an object the paper calls a bisimulation does not read the rest
of the paper. This outranks every cost experiment. See [N15](N15_soundness_regressions/plan.md).

### G4 — $r$, the number of observation regions, is never swept

It appears **explicitly** in the complexity bound, $M_{i+1} \le (r+1)^{\sum_j \Delta^j}$, and every
experiment uses whatever the problem happens to define — one region, two, or three. The theorem's own
parameter has never been varied. One sweep on a fixed system would test the bound's shape and cost an
afternoon.

### G5 — Specification complexity is never varied

Synthesis is a fixed point over (abstract states × monitor states), so the scLTL formula multiplies
the cost. Every experiment uses a single fixed formula. This matters specifically for the cost story:
route 2's cell-count penalty is multiplied by the monitor size, so E8's 4× becomes 4× *at every
monitor state*. A sweep over formulas of growing monitor size belongs in N1's and N6's tables.

### G6 — `atol` sensitivity is documented but is not an experiment

The README records that `atol = 1e-3` silently drops ~3.7 % of the domain and `1e-4` about 0.4 %.
**Every number in the folder depends on it**, and the graph-invariance result explicitly relies on
exactness that `atol > 0` breaks. One sensitivity figure — a headline quantity against `atol` over
$10^{-3} \ldots 10^{-5}$ — would immunize the whole paper for an afternoon's work.

### G7 — No comparison to $l$-complete abstractions

The paper positions itself against them in `rem:l-complete-contrast` and never measures anything.
Lowest priority: that remark argues the two are incomparable *in mechanism* — memory tightens an
over-approximation there, whereas here the quotient already attains the exact set — which is a
defensible reason not to benchmark. But say so explicitly rather than leaving it silent.

## The structural imbalance

Twelve campaign experiments, and eleven of them answer "PCLF versus the common function it induces".
That is the internal question. The external questions — *is this better than gridding?* and *is the
object actually a bisimulation?* — have zero experiments between them. G1 and G3 are worth more than
N4, N5, N9 combined.

## Revised priority

1. **N15** (soundness) — nothing else matters if the object is not what the paper says it is.
2. **N14** (bisimulation vs alternating simulation) — the missing external comparison, and the most persuasive one. **Done:** the grid is fine on reachability and collapses on co-safe LTL.
3. **N2, N12, N3, N13** — the results already in hand; write them up.
4. **N10** — the qualitative separation on their matrices.
5. **G6** (`atol` sensitivity) — cheap insurance.
6. **N7**, then N1 with a real harness.
7. **N4, N5, N9, N11** — only if N7 says the facet term is live.
8. **G4, G5** — sweeps that belong as columns in the above rather than as standalone work.
