# N11 — Designing system families where memory pays

The rotation family of [N2](../N2_feasibility_separation/plan.md) is not a lucky example; it
instantiates one of **four** distinct structural mechanisms. Naming all four turns "we found a family"
into "here is how to find families", which is worth more than any single benchmark.

## The four mechanisms

### (i) Group symmetry — what the rotation family actually uses

If the mode set is closed under a finite group action $G$, a common Lyapunov function must be
invariant under $G$: the sublevel set $B$ satisfies $gB \subseteq (\gamma/\rho)B$ for every generator,
and as $\gamma \to \rho$ this forces $B$ to be $G$-invariant. A polytopic $G$-invariant body needs at
least $|G|$ facets. A path-complete family over the **Cayley graph of $G$** sidesteps it: node $s$
carries one phase of the orbit, $V_s = \|R^{-s}x\|_\infty$, at $O(1)$ facets for every $|G|$.

This is exactly `memory_vs_geometry.jl`: $A_{1,2} = \rho R(\pm2\pi/q)$ on $\mathbb{Z}_q$.

**Design rule.** Pick a finite group with a faithful linear action that does *not* preserve a simple
polytopic norm. Rotations by $2\pi/q$ qualify. **Coordinate permutations do not** — $\|Px\|_\infty =
\|x\|_\infty$ for any permutation matrix, so $\|\cdot\|_\infty$ is already a common Lyapunov function
and the construction collapses. This is the trap when generalizing to $n$ dimensions.

**$n$-dimensional version** (what [N4](../N4_dimension_sweep/plan.md) needs): block-diagonal rotations
on $\mathbb{R}^{2k}$, $R(\theta_1)\oplus\cdots\oplus R(\theta_k)$ with $\theta_i = 2\pi/q_i$. The
$\infty$-norm is not invariant, so the orbit argument survives, and the per-node template stays at
$2n$ facets.

**Its weakness.** The Cayley graph $\mathbb{Z}_q$ is **complete**, so route 2 pays $q\times$ in cells
(measured: 260 / 585 / 1662 / 2680 against 99 / 144 / 376 / 431 for $M = 2\ldots5$). The family wins
on feasibility and on facets, and loses on cell count. See mechanism (v).

### (ii) A long spectrum-maximizing product

The complexity of a polytopic extremal norm tracks the combinatorial complexity of the s.m.p.: a
system whose worst behaviour is a long periodic word needs a unit ball shaped by that whole orbit,
while a path-complete family over the cycle of that length localizes it to one node per position.

This is [N3](../N3_smp_predictor/plan.md)'s predictor read as a design rule. Measured: E2 has s.m.p.
$211$ (length 3, 8.6 % above the best single mode) and memory buys 0.0046; Gol–Lazar–Belta has
s.m.p. $\approx$ length 1 (0.006 % above the best single mode) and memory buys exactly 0.

**Distinct from (i).** In the rotation family *every* product of length $L$ has spectral radius
$\rho^L$, so every word is extremal and the s.m.p. is degenerate. That family works by symmetry, not
by s.m.p. length. The two mechanisms are genuinely different and a good benchmark suite should have
one of each.

**Design rule.** Take a family whose s.m.p. lengthens while the stability margin stays fixed, so that
hardness grows without the system merely becoming less stable. (An earlier attempt — making the modes
pull in opposite directions with growing $\mathrm{cond}(S)$ — failed for exactly this reason: it
destroyed *stability*, not certificate simplicity, and past $s = 2.5$ nobody certified at all.)

**Literature to mine:** `Jung09` (the JSR book), `ProJunBlo10` (conic programming for joint spectral
characteristics), `PhiAth19` (path-complete Lyapunov functions: geometry and comparison). The
invariant-polytope algorithm of Guglielmi and Protasov constructs polytopic extremal norms whose
vertex count tracks the s.m.p. — **not in the bib, verify the reference before citing it.**

### (iii) Determinization blow-up — the one nobody has tried

This targets the *only* genuinely exponential step in the cost model. `th:induced_common` routes the
common function through the **observer graph**, up to $2^{|S|}$ nodes, and its sublevel set becomes
$\bigcup_{Q\in S_{\mathrm{obs}}}\bigcap_{s\in Q}\mathcal P^{(s)}_\Gamma$ — a union of that many
pieces, and differencing against a union of $k$ pieces multiplies the piece count.

**Our benchmarks never trigger it.** Measured: the 4-node observer-graph certificate determinizes to
**3** observer states out of 16, and the 2-node De Bruijn to **3** out of 4. The bound is a worst case
that the certificates in this folder do not approach — which is why the exponential argument in the
paper is currently theoretical.

**Design rule.** Build the path-complete graph as an NFA whose subset construction really is
exponential, using the classical automata-theoretic constructions (`HMU2001`, already in the bib).
Then route 1's induced common has exponentially many pieces *by design* and route 2 has $|S|$ simple
ones.

**What must be checked first.** The graph has to be path-complete **and** have exponentially many
*reachable* subsets starting from the full node set — the observer here begins at $S$, not at a single
initial state. These are compatible in principle; the construction is the work.

### (iv) A constrained switching language

No common function exists at any complexity. This is [N10](../N10_constrained_switching/plan.md) and
it is qualitatively different from (i)–(iii): not a cost separation but a feasibility one, with the
graph supplied by the problem rather than chosen. Families: any pair of Schur matrices with unstable
product, and the dwell-time / average-dwell-time literature (`ZHANG201342`, `PEDJ:16`,
`Thesis:MPhilippe`).

## (v) The combination nobody has built — highest value

| family | graph | pieces | cells | facets |
| :--- | :--- | :--- | :--- | :--- |
| rotation ($\mathbb{Z}_q$) | **complete** | **distinct** | loses $q\times$ | wins |
| dual De Bruijn | **incomplete** | **identical** | **wins 2.18×** | ties |
| **the missing one** | **incomplete** | **distinct** | should win | should win |

Mechanism A (fewer cells) needs an *incomplete* graph; the facet advantage needs *distinct* pieces.
Every family in the folder has exactly one of the two. **A family with both should win on cells and
facets simultaneously**, which is the largest speedup available.

**First attempt: tried, and it fails — but the test was badly designed.**
`rotation_templates(nodes; mode = :alternating)` assigns different templates per node, so dual De
Bruijn + `:alternating` was a one-line test. Measured, with the README's loosened tolerances
(`tol_feas = tol_gap_* = 1e-6`):

| $\theta$ | $k{=}0$ identical | $k{=}1$ dual identical | $k{=}1$ dual **alt** | $k{=}2$ dual **alt** |
| --: | --: | --: | --: | --: |
| $\pi/12$ | 0.734833 | 0.734833 | **0.761890** | **0.794147** |
| $\pi/6$ | 0.717560 | 0.717566 | **INFEASIBLE** | **INFEASIBLE** |
| $\pi/4$ | 0.742542 | 0.742542 | **INFEASIBLE** | **INFEASIBLE** |
| $\pi/3$ | 0.758582 | 0.758588 | **INFEASIBLE** | **INFEASIBLE** |

Distinct templates make the certificate *worse* where they work and infeasible everywhere else.

**Why the test was wrong, and what to run instead.** `:alternating` gives one node the rotated
template and the other the **identity** — so it does not diversify the two nodes, it *degrades* one
of them. And mismatched templates are strictly harder to satisfy along an edge, since
$V_d(A_mx)\le\gamma V_s(x)$ couples the two, so an over-constrained conic program is exactly what one
should expect.

The right test diversifies without degrading: give node $s$ the template $R(s\theta)$ — distinct,
equally expressive, and phase-shifted, which is precisely the structure that makes the $\mathbb{Z}_q$
family of mechanism (i) work. That needs a new mode in `rotation_templates` (say `:phase`) and is
still a small change. **Until that is run, the combination in the table above is untested, not
refuted.**

Also worth recording from this sweep: the dual De Bruijn's rate *equals* the single node's at every
$\theta$ (0.734833 vs 0.734833, 0.717566 vs 0.717560), confirming again that on this system memory
buys nothing in certificate quality and the 174-vs-379 cell reduction is purely structural.

## What we will try, in order

1. **Dual De Bruijn with *phase-shifted* templates** $R(s\theta)$ per node — incomplete graph,
   distinct but equally expressive pieces. Needs a `:phase` mode in `rotation_templates`. The
   `:alternating` version has already been tried and fails for a reason that does not apply here (see
   above).
2. **Sweep the two axes independently**: {complete, incomplete} × {identical, distinct}, four cells of
   a table, same system, same budget. That table *is* the characterization the paper is missing, and
   it subsumes `rem:when-memory-pays`.
3. **The $n$-dimensional block-rotation family** for N4.
4. **Mine `PJ:19` and `DDJ:21`** — these papers characterize when one path-complete method dominates
   another at fixed template, which is precisely our question. Their separating examples are
   ready-made candidate families and were built for exactly this purpose.
5. **The determinization-blow-up construction** (iii) — the most novel, the most work.

## Families to avoid, recorded so they are not tried twice

- **Positive switched systems** (`MasSho07`). They admit linear copositive common Lyapunov functions,
  so the common certificate is already cheap and there is nothing to redistribute.
- **Generic random low-dimensional systems.** Searched: 24 planar, 20 in dimensions 3–4, zero
  separations. A four-facet polytope is essentially an extremal norm for a random planar pair.
- **Misalignment families** ($A_1 = cR(+\theta)S$, $A_2 = cR(-\theta)S^{-1}$ with growing
  $\mathrm{cond}(S)$). They destroy stability rather than certificate simplicity.
- **Gol–Lazar–Belta's own system**, for any cost claim: one, two and four nodes certify at identical
  rates at every template budget tried.
