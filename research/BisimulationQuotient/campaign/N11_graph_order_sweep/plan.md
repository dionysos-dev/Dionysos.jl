# N11 — What happens as the graph order grows

Everything measured in this folder is De Bruijn order 1 over two modes. N1's own threats section says
so, and lists order 2 as untested. This is the other axis of the design space, and it is the one a
referee reaches for first: *you showed it on the smallest non-trivial graph — does it survive?*

## What we want to show

> As the De Bruijn order $k$ grows, the certified rate improves or saturates, and the cost advantage
> of building on the PCLF directly either holds or degrades in a way the mechanism predicts. State
> which, rather than leaving the method as an order-1 result.

## Why the answer is not obvious, and what the mechanism predicts

A De Bruijn graph of order $k$ over $M$ modes has $\lvert S\rvert = M^k$ nodes, so route 2 lifts the
abstraction over exponentially many nodes as $k$ grows. The two orientations should then separate
*harder*, not less, and in opposite directions.

**Primal, complete.** N1 §2b: if every node has an outgoing edge for every mode, lifting produces
$\lvert S\rvert$ near-isomorphic copies of one partition — the node records what was played and
constrains nothing about what follows, so it must still serve the whole alphabet. At order 1 that
costs route 2 a factor of 2 and it wins anyway, on cell simplicity. At order $k$ it costs $M^k$.
**Prediction: route 2's primal advantage shrinks with $k$ and reverses at some order.** Where, and
whether the fragmentation of the induced common grows fast enough to compensate, is the question.

**Dual, co-complete.** Every node commits to the mode played next, so $\lvert\mathrm{enabled}(s)\rvert
= 1$ at every order. The per-node partition should keep shrinking as the committed future lengthens.
But N1 §2c's break-even is $g > \lvert S\rvert$, and $\lvert S\rvert = M^k$ grows exponentially.
**Prediction: the dual holds only if $g$ grows at least as fast as $M^k$.** Measured at order 1,
$g = 4.81$ against $\lvert S\rvert = 2$ on example A — comfortable. At order 2, $\lvert S\rvert = 4$
and $g$ must exceed it. There is no reason yet to believe it does.

**The rate may not move at all.** The campaign was reset partly on this: Gol–Lazar–Belta at conic
orders 1, 2 and 3 — one, two and four nodes — certifies at *identical* rates to six digits. That
system is memory-free, and N12's s.m.p. ratio says why. So the sweep needs a system where memory
demonstrably pays, which means example A or the $\mathbb{Z}_q$ family of N5, not the benchmark.

## What we will try

1. **Rate against $k$, for $k = 1, 2, 3$**, on a system N12's predictor says has something to
   remember. Cheap: it is a certificate search, no quotient. This alone answers whether higher order
   buys anything, and it is the prerequisite for everything below.
2. **Both orientations, $k = 1, 2$**, with the full cost comparison of N1 — cells, polytopes per
   cell, facets per cell, build time — so the two predictions above are tested against each other on
   one system.
3. **Track $g$ and $\lvert S\rvert$ separately** at each order and plot $g / \lvert S\rvert$. That
   ratio is the whole prediction, and it is the number to report whether it holds or fails.
4. **Report the order at which the primal reverses**, if it does within reach.

[`../../experiments/debruijn_sweep.jl`](../../experiments/debruijn_sweep.jl) already sweeps
$k = 0, 1, 2$ and draws the quotient; it is a demo, not a measurement, but it is the starting point
and it establishes that the machinery accepts $k > 1$.

## Risks and fallbacks

- **$M^k$ partitions is a memory problem before it is a time problem.** At $k = 3$ over two modes
  route 2 builds eight partitions of the working set. Keep the ladder short and the working set
  small; a few thousand cells per node is enough to compare.
- **The certificate search gets harder with $k$**, since the LP has $M^k$ pieces to fit
  simultaneously. If order 3 does not solve reliably, order 2 against order 1 is still a result — one
  doubling is enough to see the direction.
- **A null result is publishable and cheap to state.** "The advantage is an order-1 phenomenon"
  narrows the claim honestly, and it is better found by us than by a referee.

## What would falsify the method's generality

The dual's $g / \lvert S\rvert$ falling below 1 at $k = 2$, with the primal already reversed. That
would make the whole cost story an order-1 result, and the paper would have to say so. Note this does
**not** touch N5's feasibility separation, which is analytic and carries no graph order.
