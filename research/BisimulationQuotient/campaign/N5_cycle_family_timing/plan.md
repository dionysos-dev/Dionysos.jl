# N5 — Time the $\mathbb{Z}_q$ family against the common function it induces

**Targets mechanism B.** This family is the best home for it in the folder, because it is *designed*
so that per-node geometry stays fixed while the induced common's grows without bound.

## What we want to show

> On the $\mathbb{Z}_q$ construction the per-node template is four facets for every $q$, while the
> induced common function needs $4q$. If the cost of the polytope primitives is superlinear in the
> facet count, the lifted construction must eventually win on the clock as $q$ grows. Find the
> crossing point, or show there is none in the plane.

Note this is the *cost* half of N2. The *feasibility* half — that at four facets per node route 1
does not exist at all — is already proved analytically and needs no timing.

## The setup

$A_1 = \rho R(2\pi/q)$, $A_2 = \rho R(-2\pi/q)$ on the $\mathbb{Z}_q$ cycle, $V_s = \|R^{-s}x\|_\infty$
(see [`../N2_feasibility_separation/plan.md`](../N2_feasibility_separation/plan.md)). Route 1 is the
induced common at $4q$ facets, which *does* exist; route 2 is the cycle at four facets per node.

## What the existing numbers say, and why this is borderline

From [`../../notes/paper-experiments.md`](../../notes/paper-experiments.md) §E4, $M = 2 \ldots 5$:

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

1. **Run N7 first.** Same gate as N4.
2. **Push $q$ much further** — $q = 4, 8, 16, 32$ — and plot build time for both arms against $q$.
   The existing table stops at 5, which is far too early to see a crossing.
3. **Fit the exponent directly** rather than inferring it: this is N7's micro-benchmark applied to
   the actual cells produced here.
4. **Run it at $n = 4$** where `_screen_vertices` switches off (see N4). If the crossing exists
   anywhere, that is where it should appear first.
5. **Report the crossing point or its absence explicitly.** "The curves do not cross for $q \le 32$ in
   the plane" is a publishable, useful sentence.

## What would falsify it

The curves diverging rather than converging as $q$ grows — the cell penalty outrunning the facet
saving — which would mean bounded per-node geometry is not a cost advantage at all in this family,
only a feasibility one. In that case N2 stands alone and the cost story rests entirely on mechanism A
(N1). Discovering that early is worth the experiment.

## Relation to the rest

If N5 fails and N4 succeeds, the claim becomes "the facet advantage is dimensional". If both fail, the
claim becomes "memory buys feasibility (N2) and smaller abstractions on incomplete graphs (N1); it
does not buy cheaper geometry" — which is narrower, fully supported, and still a paper.
