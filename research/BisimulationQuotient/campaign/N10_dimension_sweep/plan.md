# N10 — The facet advantage is dimensional

**Targets mechanism B** (same total, simpler cells). This is the experiment the histogram story
actually needs: the per-cell facet distribution is measured and striking, but it has never produced a
clock win, and the code says why.

It absorbs what was a separate plan for timing the $\mathbb{Z}_q$ family. That family is the right
vehicle — it is *designed* so that per-node geometry stays fixed while the induced common's grows
without bound — but in the plane it cannot win, for the reason below. The two questions are one
experiment: run that family, and sweep the dimension.

## What we want to show

> The per-cell facet advantage of the lifted construction does not pay in the plane because the
> polytope primitives are too cheap there. In dimension 4 and above it should pay, and the threshold
> is predictable from the implementation rather than guessed.

On the $\mathbb{Z}_q$ construction the per-node template is four facets for every $q$, while the
induced common needs $4q$. If the cost of the primitives is superlinear in the facet count, the
lifted construction must eventually win as $q$ grows. Find the crossing point in $(q, n)$, or show
there is none.

Note this is the *cost* half of N5. The *feasibility* half — that at four facets per node route 1
does not exist at all — is already proved analytically and needs no timing.

## Why dimension 4 and not 3 — read off the code

`_screen_vertices` ([`semilinear_set.jl:185`](../../../../src/utils/sets/semilinear_set.jl#L185))
begins `LazySets.dim(P) <= 3 || return nothing`. That screen, in its own comment, "retires the large
majority of candidate pieces before any LP is spent". So:

- **$n = 2$:** the screen is active and cheap. A planar polytope with $F$ facets has $F$ vertices, so
  screening is $O(F)$ and most candidate pieces of a set difference die before an LP is allocated.
- **$n = 3$:** still active and still cheap — McMullen's upper bound gives at most $2F-4$ vertices,
  linear in $F$.
- **$n \ge 4$:** **the screen switches off entirely.** Every candidate piece of every set difference
  costs a full `clean_poly` pass, which is one LP per constraint.

Route 1's cells have far more facets than route 2's (measured 1 728 vs 146 at the tail on E2, 1 557
vs 187 on E8), so route 1 generates far more candidate pieces per difference. Losing the screen
therefore punishes route 1 disproportionately.

**Prediction:** the route-1/route-2 time ratio is flat near 1 at $n = 2, 3$ and jumps at $n = 4$.


## The setup

$A_1 = \rho R(2\pi/q)$, $A_2 = \rho R(-2\pi/q)$ on the $\mathbb{Z}_q$ cycle, $V_s = \|R^{-s}x\|_\infty$
(see [`../N5_feasibility_separation/plan.md`](../N5_feasibility_separation/plan.md)). Route 1 is the
induced common at $4q$ facets, which *does* exist; route 2 is the cycle at four facets per node.

## What the existing numbers say, and why this is borderline

Measured earlier, over $M = 2 \ldots 5$:

| $M$ | max facets/cell, PCLF | max facets/cell, CLF | cells, PCLF | cells, CLF |
| --: | --: | --: | --: | --: |
| 2 | 225 | 259 | 260 | 99 |
| 3 | 210 | 346 | 585 | 144 |
| 4 | 198 | 393 | 1 662 | 376 |
| 5 | 264 | 493 | 2 680 | 431 |

Per-cell complexity stays flat for the PCLF and grows for the CLF, exactly as designed. But route 2
carries 2.6–6.2× more cells against only ~1.9× simpler ones. With a cost model $c(F) = F^\alpha$,
route 2 wins at $M = 5$ only if $(493/264)^\alpha > 6.2$, that is $\alpha > 2.9$. Nearly cubic.

That is not absurd — a set difference emits $O(F_Q)$ pieces, each cleaned at one LP per constraint,
so $O(F^3)$ is a defensible model — but it is tight, and the cell ratio is growing faster than $M$,
which works against us.

**Also note** the facet counts above are dominated by *refinement* cuts, not by the template: 264
facets per cell against a template of 8 halfspaces. Refinement cuts are largely shared between the two
routes, which is why the facet ratio (1.9×) is so much smaller than the template ratio ($q$). This is
the structural reason the family may not deliver.


## What we will try

1. **Run N8 first.** If bookkeeping dominates the build at all dimensions, this experiment cannot
   show anything and should not be attempted.
2. **Push $q$ much further** — $q = 4, 8, 16, 32$ — and plot build time for both arms against $q$,
   in the plane. The existing table stops at 5, which is far too early to see a crossing.
3. **Sweep the dimension at fixed $q$**, $n = 2, 3, 4$, on the block-rotation generalization of the
   $\mathbb{Z}_q$ construction, which keeps the per-node template at $2n$ facets in every dimension
   while the induced common grows. N9 carries the trap: coordinate permutations do **not** work,
   since $\|Px\|_\infty = \|x\|_\infty$ makes the infinity norm already common.
4. **Plot the time ratio against $n$**, with the facet ratio and cell ratio beside it, so the reader
   can see which of the three moves.
5. **Instrument the screen.** Count how many candidate pieces `_screen_vertices` retires and how many
   reach an LP, per arm, per dimension. That turns the explanation from a plausible story into a
   measured one, and it is a few counters.
6. **Report the crossing point or its absence explicitly.** "The curves do not cross for $q \le 32$
   in the plane" is a publishable, useful sentence.

## Risks and fallbacks

- **The curse of dimensionality may make $n = 4$ intractable.** Every benchmark in the folder is
  planar and two-mode for a reason. Mitigate with a small working set, few observation regions and a
  low `max_slices`; a quotient of a few thousand cells is enough to time.
- **If $n = 4$ is out of reach**, the honest fallback is to run N8's instrumentation at $n = 3$ with
  the screen *artificially disabled*, which isolates the screen's contribution without paying the
  dimensional cost. That is a weaker experiment but it tests the same mechanism.
- **`atol` behaves differently in higher dimension** and the tolerance study in
  [`../../README.md`](../../README.md) is planar. Expect to re-tune, and report what was used.


## What would falsify it

A flat time ratio at $n = 4$ with the screen confirmed off, or the curves diverging rather than
converging as $q$ grows — the cell penalty outrunning the facet saving. Either would mean the
per-cell facet count is not what the build time depends on, and mechanism B should be dropped from
the paper's claims, leaving mechanism A (N1) as the whole cost story.

## Relation to the rest

If the crossing exists, the claim becomes "the facet advantage is dimensional" and mechanism B is
live. If it does not, the claim becomes "memory buys feasibility (N5) and smaller abstractions on
incomplete graphs (N1); it does not buy cheaper geometry" — which is narrower, fully supported, and
still a paper. Discovering that early is worth the experiment.
