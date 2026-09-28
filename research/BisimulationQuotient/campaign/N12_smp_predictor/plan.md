# N11 — The spectrum-maximizing product predicts whether memory can help

Costs a handful of eigenvalue computations and explains every data point in the folder. It is a
*contribution*, not a benchmark: it tells a reader when to reach for a PCLF before building anything.

## What we want to show

> Memory encodes which words have been seen. If the extremal behaviour of the system is a single
> mode repeated, there is nothing to remember, and a path-complete family cannot beat a common one at
> any template budget. The length of the spectrum-maximizing product is therefore a cheap *a priori*
> test of whether the method can help on a given system.

## The evidence, already measured — and the predictor is the *ratio*, not the word length

| | JSR lower bound | attained by | $\max_m \rho(A_m)$ | ratio | 1-node vs PCLF gap |
| :--- | --: | :-- | --: | --: | --: |
| **Gol–Lazar–Belta** | 0.85585 | $2111111$ | 0.85580 | **1.00006** | **0.0** |
| **E2 observer-graph** | 0.86881 | $211$ | 0.80016 | **1.086** | **0.0046** |
| **two_mode** | 0.70395 | $2111111111$ | 0.70000 | **1.00564** | **0.0** |

Three systems, three consistent calls — and `two_mode` is the one that sharpens the rule. It has a
**length-10** s.m.p. word and yet memory buys nothing on it: N1 measures the dual De Bruijn and the
single node certifying the identical rate 0.717572. The word is long but beats the best single mode
by only 0.56 %.

**So the instrument is $\max_w \rho(A_w)^{1/|w|} \big/ \max_m \rho(A_m)$, not the length of $w$.**
When that ratio is near 1 the extremal behaviour is essentially "stay in one mode", there is no word
structure to remember, and no path-complete family can beat a common one. On Gol–Lazar–Belta this is
confirmed independently by the budget sweep: one, two and four nodes certify at identical rates at
conic orders 1, 2 and 3 alike.

The same reasoning explains an earlier negative result: for random planar systems
a four-facet polytope is essentially an extremal norm, so the s.m.p. is short and the common
certificate is already tight.

## What we will try

1. **Tabulate across every system in the folder**: `two_mode_problem`, `observer_graph_problem`,
   `gol_lazar_belta_problem`, the $\mathbb{Z}_q$ family of N5, and the constrained system of N2.
   Columns: s.m.p. word and length, JSR lower bound, $\max_m\rho(A_m)$, the ratio, the 1-node rate,
   the PCLF rate, the gap.
2. **Plot the gap against the ratio.** If the correlation holds across five or six systems, the rule
   is publishable as a heuristic. If it does not, the counterexample is more interesting than the
   rule.
3. **Test the rule prospectively.** Generate systems with a deliberately long s.m.p. — the standard
   knob is a family whose extremal norm's facet count grows with s.m.p. length — and check that the
   gap appears where predicted. This doubles as the family N10 needs.
4. **State the caveat.** A periodic lower bound over words of length $\le L$ is a lower bound, not the
   JSR; a long s.m.p. beyond $L$ would be invisible. Report $L$ and the best word.

## Why it matters beyond this paper

Every user of the method faces the question "should I use one node or many?", and until now the only
answer was "try both". This gives an answer costing a few eigenvalues, and it turns the E5 negative
result from an embarrassment into a characterization of where the method's regime is.

## What would falsify it

A system with a short s.m.p. where memory nevertheless buys a large gap, or a long s.m.p. where it
buys nothing. Either kills the rule as stated — and would be worth chasing, because it would point at
the real mechanism.
