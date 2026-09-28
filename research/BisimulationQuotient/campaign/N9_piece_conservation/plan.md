# N9 — Why is $\Sigma$ pieces conserved on some certificates and not others?

Closes `rem:when-memory-pays` from the geometry side, where N1 approaches it from the cell-count side.

## The puzzle

Total geometric complexity is conserved to ~1 % on two benchmarks and provably inflated by $|S|$ on a
third:

| case | route 2 | route 1 | ratio |
| :--- | --: | --: | --: |
| E2, $\Sigma$ pieces | 24 780 | 24 580 | 1.008 |
| E8, $\Sigma$ pieces | 40 979 | 40 463 | 1.013 |
| complete graph, identical pieces | $|S| \cdot n$ | $n$ | $|S|$ — **proved** |
| N1 dual De Bruijn, $\Sigma$ facets | 3 519 | 5 439 | **0.65** |

So conservation is not a law: `prop:complete-redundant` refutes it exactly, and N1's dual De Bruijn
goes the *other* way. Something about the certificate decides which regime holds.

## Two things that must be said before any of this is published

1. **Conservation is not "almost by construction"**, as
   [`../../notes/complexity-vs-memory.md`](../../notes/complexity-vs-memory.md) §6 currently puts it.
   The complete-graph case is the construction, and it refutes conservation.
2. **The facet totals are not independent evidence.** Every arm measures **4.00 facets per convex
   piece** (3.995, 3.998, 3.999, 4.002), which is what Euler's formula forces for a generic planar
   subdivision. So every facet ratio equals the corresponding piece ratio to three digits, and the
   conserved quantity is $\Sigma$ *pieces*. Report pieces as primary and state the identity, rather
   than presenting two apparently independent conserved quantities. The identity is planar and will
   not survive N4.

## The hypothesis

Conservation holds when the per-node refining families are largely *disjoint* — each cut used at one
node only — and degrades toward the $|S|$ factor as they become identical. The complete-graph case
with identical pieces is one endpoint; E2 and E8, with distinct pieces, sit at the other.

## What we will try

1. **The interpolation.** A 2-node graph; rotate the second node's template by $\theta$ and sweep
   $\theta \to 0$, so the pieces go from distinct to identical. Plot
   $\Sigma\text{pieces}(\text{route 2}) / \Sigma\text{pieces}(\text{route 1})$ against $\theta$.
   *Expected:* a monotone run from ~1 to 2.
2. **Clear the coverage confound first.** Two certificates induce different working sets and terminal
   sets, and on Gol–Lazar–Belta ~86 % of every quotient's cells lie outside $\mathcal X$. Restrict
   both arms to a common reference set before comparing totals. If the 1 % conservation does not
   survive restriction, it was compensating coverage differences and the whole observation dissolves
   — better to know.
3. **Measure the overlap directly.** For each pair of nodes, count how many refining cuts are shared.
   That turns "overlap" from a story into a number and lets it be plotted against the ratio.

## What would falsify it

Non-monotone behaviour in $\theta$, or a ratio below 1 somewhere (N1's dual De Bruijn already achieves
0.65 on facets, so the range is wider than the hypothesis allows — this may already be a
counterexample and should be checked first). Either outcome means overlap is not the mechanism, which
matters because one candidate explanation for `rem:when-memory-pays` — the mode set of a node — has
already been refuted.

## Relation to N1

N1 asks *when does lifting reduce the number of cells*; N9 asks *when does it reduce total geometric
complexity*. The dual De Bruijn answers yes to both at once, which suggests the two questions have one
answer. If the interpolation confirms the overlap mechanism, that answer is in hand.
