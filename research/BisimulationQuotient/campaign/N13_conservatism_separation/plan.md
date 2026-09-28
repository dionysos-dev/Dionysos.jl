# N13 — The conservatism separation: a tighter certificate at equal per-node budget

Already measured and currently homeless. It is one of only three positive results in hand, and it is
the direct evidence for "less conservative at the same template" — the claim the introduction makes
and no experiment currently backs.

## What we want to show

> At an identical per-node template budget, a path-complete family certifies a strictly better
> contraction rate than any single common function. Both certificates exist; one is tighter.

**This is a *conservatism* separation and must not be confused with the *feasibility* separation of
[N2](../N2_feasibility_separation/plan.md)**, where the single node certifies nothing at all. Both are
worth claiming; conflating them is an overclaim a referee will catch, because they are different
quantities.

## The measurement — and the correction the fair search forced

**An earlier version of this file quoted a gap of 0.0447 from a single-orientation baseline. That was
wrong by a factor of ten, and it is the third time in this work that an under-searched baseline has
inflated a separation.** The honest numbers, from
[`conservatism.jl`](conservatism.jl) with 41 orientations:

| | rate | slack over the JSR lower bound |
| :--- | --: | --: |
| JSR lower bound (periodic, $L \le 10$), word $211$ | 0.868807 | — |
| **4-node observer graph** | **0.868964** | **0.02 %** |
| 1 node, best of 41 orientations | **0.873586** | **0.55 %** |
| **gap** | **0.004622** | |

Best-so-far against search budget, which is why the script reports it:

| draws | 1 | 2 | 5 | 10 | 20 | 41 |
| :-- | --: | --: | --: | --: | --: | --: |
| single-node rate | 0.913680 | 0.913680 | 0.880081 | 0.880081 | 0.875496 | **0.873586** |

The identity orientation alone gives 0.9137; ten times the search closes 90 % of the apparent gap.

**The separation is real but small.** The right claim is that the PCLF is *essentially tight* against
the JSR (0.02 %) where a properly searched single node is 0.55 % off — not that there is a large gap.
The four pieces are genuinely distinct (pairwise support-function distance 0.025–0.470) and the graph
is incomplete.

**Whether this is worth a paragraph in the paper is now a judgement call.** A 0.0046 rate difference
is defensible but unexciting; the feasibility separations of [N2](../N2_feasibility_separation/plan.md)
and [N12](../N12_single_mode_clock/plan.md), where the single node certifies *nothing*, are far
stronger and carry no search-budget risk at all. Consider reporting this as supporting evidence for
the s.m.p. predictor of [N3](../N3_smp_predictor/plan.md) rather than as a headline.

**And the mechanism is understood**: the s.m.p. is the length-3 word $211$, 8.6 % worse than the best
single mode — so there is word structure for memory to remember. Contrast Gol–Lazar–Belta, s.m.p.
length $\approx 1$, gap exactly 0. See [N3](../N3_smp_predictor/plan.md).

## What we will try

1. **Give the single node a fair search.** The number above is one template orientation. Search over
   many, as `memory_vs_geometry.jl`'s `single_node_rate` does with 40 draws plus the identity, and
   report best-so-far against search budget. **Under-searching the baseline has manufactured
   separations twice in this work**; it must not happen a third time.
2. **Sweep the node count**: 1, 2, 4 nodes on the same template. If the rate improves monotonically
   the story is clean; if 2 nodes already captures it, say so — that is a *better* result, since it
   bounds how much memory is needed.
3. **Sweep the template budget**: conic orders 1, 2, 3. The interesting question is whether the single
   node ever catches up, and at what cost in facets. If it catches up at order 3, the honest claim
   becomes "equal certificate quality at one third of the per-node geometry", which is still the
   thesis.
4. **Connect it to the abstraction.** A better rate means fewer slices to reach the terminal set;
   report slice counts and cells alongside, so the certificate result and the cost result are one
   story rather than two.

## What would falsify it

The single node reaching 0.86896 once searched properly. That would make E2 a second memory-free
benchmark alongside Gol–Lazar–Belta, and would leave N2 and N12 as the only positive results — which
would be worth knowing immediately, since it would decide the shape of the paper.

## Cost

Low. Only certificates are computed, no quotients, so each point is a conic program of a few seconds.
