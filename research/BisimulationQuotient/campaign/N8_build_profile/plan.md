# N8 — Where the build time actually goes

**The gate.** It decides whether N10 is worth attempting, and it is cheap. Run it first.

## What we want to show

> What fraction of a quotient build is dependent on the per-cell facet count, and what fraction is
> fixed per-cell overhead. If overhead dominates, no facet-based argument can produce a speedup at
> current benchmark sizes, however it is measured.

## Why the question is open

The structural argument says route 1's cells are 8–12× fatter at the tail and cost is superlinear in
facet count, so route 2 should be much faster. Measured, it is not: E2 is a tie (40.7 s vs 41.5 s)
and E8 is a 1.30× loss. Three explanations are live and this experiment separates them:

1. **Per-cell overhead dominates** — bookkeeping, allocation, emptiness LPs on tiny polytopes — so
   the facet-dependent term is a small fraction of the build. The facet argument would then be
   correct but not yet *active*.
2. **The vertex screen is doing the work.** `_screen_vertices` retires most candidate pieces before
   any LP, and it is active only at $n \le 3$, which is where every benchmark lives.
3. **Route 2's cell-count penalty simply exceeds its facet saving** on these problems.

## RESULT, step 3 — the exponent is measured, and `rem:complexity` is vindicated

[`1_primitive_scaling.jl`](1_primitive_scaling.jl). Random bounded polytopes at a controlled facet
count $F$, overlapping operands, nanosecond timing with adaptive batching.

| primitive | $n=2$ | $n=3$ | $n=4$ | ms at $F=64$, $n=2$ |
| :--- | --: | --: | --: | --: |
| `set_difference_decompose` | **$F^{2.45}$** (r²=.988) | **$F^{2.17}$** (r²=.987) | **$F^{2.72}$** (r²=.998) | **190** |
| `clean_poly` | $F^{1.56}$ (.981) | $F^{1.61}$ (.996) | $F^{1.64}$ (1.000) | 6.4 |
| `preimage_linear` | $F^{1.56}$ (.973) | $F^{1.64}$ (.994) | $F^{1.80}$ (.989) | 4.6 |
| `isempty` | $F^{0.68}$ (.941) | $F^{0.45}$ (.872) | $F^{0.83}$ (.923) | 0.04 |

**Three things follow, and they are not what the stalled timing numbers suggested.**

1. **The claim is true and quantified.** Every refinement primitive is superlinear, and the dominant
   one is roughly quadratic-to-cubic: at $n=2$ the set difference runs 0.50 ms at $F=4$ and
   **17.0 s at $F=256$**, a factor of 33 700 for a factor of 64 in facets. `rem:complexity` can now
   cite a number instead of an assertion, and this is `fig:complexity-comparison`.
2. **Mechanism B is live, not dormant.** Simpler cells are worth paying for. Whatever explains E2's
   tie and E8's 1.30× loss, it is *not* that facet count fails to matter.
3. **Emptiness LPs are not the bottleneck** — the one primitive that is *sub*linear, and 0.04 ms
   against the set difference's 190 ms at $F=64$, roughly 4 500× cheaper. That removes the "emptiness
   LPs on tiny polytopes" half of hypothesis 1 without needing a build profile.

**Read the dimension column with care.** The $F$ ranges differ per dimension (capped so the sweep
terminates: 4–256, 6–96, 8–64), so the exponents are **not** directly comparable across $n$ and the
apparent rise from 2.45 to 2.72 is confounded. What the column does support is the robust statement:
the set difference is strongly superlinear in every dimension tested.

**What this does *not* settle.** The exponent says cost climbs steeply *when facet counts are large*.
It says nothing about how often they are. Real cells average 4.00 facets (Euler) with a tail reaching
142–243, and at $F=4$ the set difference costs 0.50 ms — flat and cheap. So the live question is now
sharper than the original three hypotheses: **what fraction of build time is spent on the fat tail?**
That is steps 1–2, and it is the number the decision rule below actually needs.

## What we will try

1. **Profile one build per arm** on E2 and on E8, attributing time to: `set_difference_decompose`,
   `clean_poly` / `remove_redundant_constraints`, emptiness LPs, `_screen_vertices`,
   `preimage_linear`, and everything else.
2. **Count, as well as time.** Number of set differences, candidate pieces generated, pieces retired
   by the screen, LPs solved, and the facet-count distribution of the polytopes each primitive was
   handed. Counts survive machine noise in a way timings do not.
3. **Micro-benchmark the primitives separately** — this is the old X2. Random polytopes at controlled
   facet counts $F \in \{4, 8, \ldots, 512\}$; time `set_difference_decompose(P, Q)` over
   $(F_P, F_Q)$ and fit the exponent; separately fix $F$ and sweep the number of pieces $k$ of the
   subtrahend union. Expect roughly quadratic in $F$ and multiplicative in $k$. Repeat at $n = 4$,
   where the screen is off.
4. **Report a single number**: the facet-dependent fraction of the build, per arm, per benchmark.
   That number is what gates N10.

## Decision rule

- **Facet-dependent fraction above ~50 %** — mechanism B is live and N10 is worth running.
- **Below ~20 %** — mechanism B cannot show anything at these sizes. Run N10 only as the
  screen-disabled experiment, since the dimensional jump is then the only way to move the fraction.
- **In between** — run N10's dimension sweep first, since it changes the regime, and let its result
  decide whether the $\mathbb{Z}_q$ timing is worth pushing.

## Why this also improves the paper directly

The measured exponent from step 3 is what `rem:complexity`'s claim that the refinement operations are
"superlinear in the facet count of the cells involved" currently asserts without evidence, and it
supplies the missing `fig:complexity-comparison`. Even if the decision rule kills N10, this
experiment pays for itself.

## What would falsify the framing

Finding that the dominant cost is neither facet-dependent nor per-cell, but something structural such
as the slice bookkeeping or the transition construction. That would redirect the optimization effort
entirely, and it is worth knowing.
