# N4 — The facet advantage should switch on at dimension 4

**Targets mechanism B** (same total, simpler cells). This is the experiment the histogram story
actually needs: the per-cell facet distribution is measured and striking, but it has never produced a
clock win, and the code says why.

## What we want to show

> The per-cell facet advantage of the lifted construction does not pay in the plane because the
> polytope primitives are too cheap there. In dimension 4 and above it should pay, and the threshold
> is predictable from the implementation rather than guessed.

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

## What we will try

1. **Run N7 first.** If bookkeeping dominates the build at all dimensions, this experiment cannot
   show anything and should not be attempted.
2. **One PCLF and its induced common, at $n = 2, 3, 4$.** The certificate family must be comparable
   across dimensions; the natural choice is a block-rotation generalization of the $\mathbb{Z}_q$
   construction of N2, which keeps the per-node template at $2n$ facets in every dimension while the
   induced common grows.
3. **Plot the time ratio against $n$**, with the facet ratio and cell ratio beside it, so the reader
   can see which of the three moves.
4. **Instrument the screen.** Count how many candidate pieces are retired by `_screen_vertices` and
   how many reach an LP, per arm, per dimension. That converts the explanation from a plausible story
   into a measured one, and it is a few counters.

## Risks and fallbacks

- **The curse of dimensionality may make $n = 4$ intractable.** Every benchmark in the folder is
  planar and two-mode for a reason. Mitigate with a small working set, few observation regions and a
  low `max_slices`; a quotient of a few thousand cells is enough to time.
- **If $n = 4$ is out of reach**, the honest fallback is to run N7's instrumentation at $n = 3$ with
  the screen *artificially disabled*, which isolates the screen's contribution without paying the
  dimensional cost. That is a weaker experiment but it tests the same mechanism.
- **`atol` behaves differently in higher dimension** and the tolerance study in
  [`../../README.md`](../../README.md) is planar. Expect to re-tune, and report what was used.

## What would falsify it

A flat time ratio at $n = 4$ with the screen confirmed off. That would mean the per-cell facet count
is simply not what the build time depends on, and mechanism B should be dropped from the paper's
claims — leaving mechanism A (N1) as the whole cost story, which is a perfectly good outcome to
discover early.
