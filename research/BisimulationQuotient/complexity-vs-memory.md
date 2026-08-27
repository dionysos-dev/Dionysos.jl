# Memory versus discretization in PCLF bisimulation quotients

Analysis note, 2026-08-26. Every number below was measured on this branch; the scripts are named
so each can be reproduced or contested. **Facts**, **interpretation** and **conjecture** are
labelled separately.

Two questions motivated this note:

1. *Shouldn't the bisimulation — and so the controllable region — be the same whatever the graph?*
2. *Does the path-complete approach bring more than just memory in the abstraction?*

Short answers: **(1) yes, and it is measurably so — that is a feature, not a problem.** **(2) The
graph cannot enlarge what is certifiable about the concrete system, but it decides whether a
certificate exists at all, and what the abstraction costs.** §5 gives a benchmark where the
single-node method has no certificate and therefore no quotient, while a two-node family builds
one in 32 seconds.

---

## 1. Position relative to the literature

The direct predecessor is **Aydin Gol, Ding, Lazar & Belta, *Finite Bisimulations for Switched
Linear Systems*** ([arXiv:1208.5471](https://arxiv.org/abs/1208.5471)). From the abstract: given
observations over polytopic subsets and *a switched linear system with stable subsystems*, their
algorithm "generates the bisimulation quotient in a finite number of steps with the aid of
**sublevel sets of a polyhedral Lyapunov function**", starting from a sublevel set containing the
origin and iteratively extending to larger ones, for synthesis against syntactically co-safe LTL.

That is exactly the **k = 0** case of this codebase: one node, one common polyhedral Lyapunov
function. The contribution under discussion is the generalisation from *one* function to a
**path-complete family** over a labelled graph — Ahmadi, Jungers, Parrilo & Roozbehani, *Joint
Spectral Radius and Path-Complete Graph Lyapunov Functions*, SIAM J. Control Optim. 52(1):687–717,
2014 ([arXiv:1111.3427](https://arxiv.org/abs/1111.3427)); polyhedral templates in Athanasopoulos
& Jungers, *Polyhedral Path-Complete Lyapunov Functions*, CDC 2019
([IEEE](https://ieeexplore.ieee.org/document/9029905/)).

**Why memory should be expected to help — the strongest argument, and it is not empirical.**
Ahmadi & Jungers, *Lower Bounds on Complexity of Lyapunov Functions for Switched Linear Systems*
([arXiv:1504.03761](https://arxiv.org/abs/1504.03761)) prove that for every integer *d* there exist
families of two matrices, in fixed dimension, **stable under arbitrary switching**, admitting no
polytopic Lyapunov function with ≤ *d* facets (nor piecewise quadratic with ≤ *d* pieces). Since a
common polyhedral function always *exists* when the JSR is below 1, but its facet count cannot be
bounded, there are only two ways forward: let one function grow without bound, or add memory. The
partition of a bisimulation quotient is built from that function's sublevel sets, so its complexity
is inherited directly. **That is the theoretical basis for "replace state-discretization complexity
by memory augmentation", and it does not depend on finding a favourable benchmark.**

One caution the literature is emphatic about, and which cost me a wrong result below: the ordering
of path-complete methods is **template-dependent** (Athanasopoulos & Jungers 2019; *Necessary and
Sufficient Conditions for Template-Dependent Ordering of Path-Complete Lyapunov Methods*, HSCC 2022,
[ACM](https://dl.acm.org/doi/10.1145/3501710.3519539)). A richer graph with a poor template loses
to a poorer graph with a good one. Every comparison below therefore holds the template fixed across
orders.

---

## 2. The controllable region is graph-invariant (measured)

Same system, same specification `F(R1 & F(D))`, same tolerances; De Bruijn order k = 0, 1, 2.

| k | nodes | cells | covered volume | controllable volume | ctrl / covered | JSR bound |
|--:|--:|--:|--:|--:|--:|--:|
| 0 | 1 | 379 | 39.109917 | 2.152238 | 0.055031 | 0.7175720215 |
| 1 | 2 | 174 | 40.044816 | 2.109348 | 0.052675 | 0.7175720215 |
| 2 | 4 | 187 | 40.078437 | 2.229715 | 0.055634 | 0.7175720215 |

The ±3 % spread in controllable volume is **not** a real difference. Volumes are the wrong
instrument: the three quotients do not even cover the same region (39.11 / 40.04 / 40.08), because
`max_slices = 8` truncates the sublevel family at a different outer level per certificate, and
`atol = 1e-3` erodes every cut with a compounding that differs per graph — measured elsewhere on
this branch at ≈ 3.7 % of the domain. The spread is the size of the known error.

A pointwise test has neither confound:

```
probes = 50, disagreements = 0
controllable counts:  k=0 → 5    k=1 → 5    k=2 → 5
```

**Proposition 1 (invariance).** *Let `x⁺ = A_σ x` with observations given by a fixed polytopic
partition, and let `Q_G` be an exact bisimulation quotient built from a path-complete family over
`G`. Then for any co-safe LTL specification over those observations, the winning region computed on
`Q_G`, projected to `X`, is the same for every `G`.*

*Proof sketch.* The concrete transition relation `x → A_σ x` does not mention the node. A node is
auxiliary state updated alongside the mode through the edges of `G`; the product system
`X × V` therefore projects onto `X` as a *ghost extension* — its reachable behaviours in `X`
coincide with those of the concrete system, for every `G`. Winning regions of a co-safe LTL game
are determined by the behaviours, and a bisimulation preserves them exactly. Hence the projected
winning region depends on the system and the specification alone. ∎

The measurement is the check on the hypothesis, not on the conclusion. **Had the probes disagreed,
the right reading would have been that the construction is not exact for some graph — a soundness
bug, not a modelling insight.** They agreed, on 50 points and three graphs.

**Corollary (what this rules out).** No path-complete graph can certify a *larger* controllable set
than another, as long as both quotients are exact. Any future experiment reporting one should be
treated as evidence of lost exactness (`atol`, truncation) rather than of a better certificate.

---

## 3. What the graph does change: the cost of the representation

Two comparisons, with genuinely different signatures.

### 3a. A path-complete family versus the common function it induces

`pclf_vs_clf.jl`: four-node observer graph, against the CLF obtained from it by
`build_common_lyapunov` — same certificate content, collapsed onto one node.

| | cells | nodes | mean pieces | median | max | mean faces | median | max |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| **PCLF** | 9 027 | 4 | 2.75 | 2 | 36 | 10.97 | 8 | 146 |
| **CLF** | 1 000 | 1 | 24.58 | 9 | 432 | 98.28 | 36 | 1 728 |

The total is **conserved to within 0.8 %**:

```
PCLF   Σ pieces = 24 780     Σ faces = 98 987
CLF    Σ pieces = 24 580     Σ faces = 98 276
ratio                 0.99                0.99
```

Nine times more cells, each about nine times simpler (24.58 / 2.75 = 8.9 against a cell ratio of
9.03). This is the clean form of the trade: a near-exact one-for-one exchange.

The tail matters more than the mean, because the geometric operations scale with per-cell size:

| quantile | PCLF pieces | CLF pieces | PCLF faces | CLF faces |
|---|--:|--:|--:|--:|
| median | 2 | 9 | 8 | 36 |
| 90 % | 5 | 58 | 21 | 232 |
| 99 % | 12 | 243 | 46 | 972 |
| max | 36 | 432 | 146 | **1 728** |

The gap widens with depth — 4.5× at the median, 11.8× at the 99th percentile.

### 3b. Increasing the De Bruijn order

Same system and templates, only k changes. Here the total is **not** conserved — it falls.

| k | nodes | cells | cells/node | Σ pieces | Σ faces | mean pcs | max pcs | mean fcs | max fcs |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 0 | 1 | 379 | 379.0 | 1 367 | 5 439 | 3.61 | 36 | 14.35 | 142 |
| 1 | 2 | 174 | 87.0 | 889 | 3 519 | 5.11 | 32 | 20.22 | 128 |
| 2 | 4 | 187 | 46.8 | 1 034 | 4 080 | 5.53 | 30 | 21.82 | 120 |

k = 1 gives **0.46× the cells and 0.65× the total faces**. Per node the discretization becomes
dramatically coarser: 379 cells in one node become 47 in each of four.

Note the direction of `mean pieces` — it goes **up** (3.61 → 5.11 → 5.53) while the total goes
down. This is therefore *not* the mechanism of §3a.

**Interpretation — two distinct effects, worth separating in the paper.**

1. **Redistribution** (§3a): the same geometry spread over more, simpler cells. Total unchanged,
   tail cut by an order of magnitude. This is what makes each geometric operation cheap.
2. **Reduction** (§3b): a richer family needs less refinement in total. This is a real saving, not
   a repackaging.

---

## 4. A negative result, and the trap it came from

**The examples in this folder cannot exhibit any certificate advantage, by construction.** All
three De Bruijn orders return *exactly* `0.7175720214843749`. The cause is not the theory but the
call: `rotation_templates(nodes)` returns `Dict(v => Rot for v in nodes)` — the **same** template
matrix for every node. Identical templates make the nodes indistinguishable and the family
collapses to k = 0. Screened across 13 systems with `compute_symmetric_2n_faces_polyhedral_pieces_pclf`
and this scheme, the gain was `+0.0000 %` in **every** case.

Assigning a *distinct* angle per node appeared to fix it — random #5 improved 1.3214 → 1.1919 →
1.1253, a 14.8 % gain. **That result was an artefact and is withdrawn.** The scheme gives node *i*
the angle `π(i−1)/2n`, so at k = 0 it always picks θ = 0, which for that system is the worst
template. Swept honestly:

```
k = 0 baseline over the template angle:
  θ = 0.000  bound = 1.321441      <- what the scheme handed k = 0
  θ = 1.047  bound = 0.996251      <- best single template
  best over 13 angles = 0.996251
```

A well-chosen *single* template already certifies the system, while the two-node family gives
1.1919. Path-completeness **lost** the fair comparison, and the separation window closed
(`1.0038 < s < 0.8887`, empty). This is precisely the template-dependence the literature warns
about, reproduced by accident.

**Method rule adopted after this:** compare only with the template held fixed across orders, or
sweep the template at every order and compare best against best.

---

## 5. RETRACTED — the separating benchmark was an artefact, and it exposed a bug

**Everything in this section was withdrawn on 2026-08-27.** It is kept because the way it failed
is the most useful thing in this note.

The claim was: at conic-partition order 3, the 1-node certificate fails (bound ≥ 1) while the
2-node one succeeds, so memory decides whether an abstraction exists. Order 3 is one template, so
the honest check is whether the single-node method fails at *every* affordable complexity. It does
not. Sweeping the order at k = 0, s = 0.76:

| order | 1 | 2 | **3** | 4 | 5 | **6** |
|---|---|---|---|---|---|---|
| k = 0 bound | 0.7785 ✓ | 0.7572 ✓ | **1.0101 ✗** | 0.8327 ✓ | 0.8979 ✓ | **1.1091 ✗** |

The single-node method certifies at orders 1, 2, 4 and 5, and fails only at 3 and 6. **Order 3 was
cherry-picked by accident.**

**The bug.** `conic_partitions_2d` is strictly nested — each order bisects every cone of the
previous one — so a higher order is a *strictly richer* template and the optimal bound must be
**non-increasing** in order:

> `bound(1) ≥ bound(2) ≥ bound(3) ≥ …`

The measured sequence violates this in both directions. For an exactly-solved relaxation that is
impossible, so `compute_polyhedral_pieces_pclf` is returning non-optimal values at some orders —
LP conditioning, a badly-terminating bisection on ρ, or a formulation error. **Until this is
understood, no comparison of bounds across partition orders can be trusted**, which includes every
"separation" found in the wider screen.

The same sweep at s = 0.80 shows it is not a one-off, and sharpens the diagnosis:

| order | 1 | 2 | **3** | 4 | 5 | **6** |
|---|---|---|---|---|---|---|
| k = 0 bound | 0.8195 ✓ | 0.8091 ✓ | **1.1675 ✗** | 0.8721 ✓ | 0.9065 ✓ | **1.1675 ✗** |
| k = 1 bound | 0.7970 ✓ | 0.8568 ✓ | 0.9342 ✓ | 0.8966 ✓ | **1.0425 ✗** | **1.1675 ✗** |

Non-monotone for *both* orders of graph. Note that orders 3 and 6 return the **identical** value
`1.167512` at k = 0, and order 6 returns that same value at k = 1 too — a repeated constant across
different templates is the signature of a failure path being taken, not of three different
relaxations happening to agree.

**The system was never hard to begin with.** A brute-force lower bound over periodic products up to
length 12 gives JSR ≥ 0.757145 at s = 0.76 (word `2`) and ≥ 0.796995 at s = 0.80 (word `222`). The
k = 0 order-2 certificate returns 0.757151 at s = 0.76 — matching the lower bound to five decimals.
The common Lyapunov function is essentially optimal here, so there was no room for memory to help.

**The one genuine gain, properly caveated.** At s = 0.80 the best k = 1 bound over the orders that
solved is **0.797001** — attaining the JSR lower bound 0.796995 to five decimals, at the *coarsest*
template — while the best k = 0 bound is **0.809054**, a 1.5 % gap. Two nodes at order 1 beat one
node at every order tried. This comparison takes a minimum over orders, so it is not directly
poisoned by non-monotonicity. **But it is still not safe:** because the solver fails at some orders,
the k = 0 minimum is only an upper bound on what a correct solver would find, and the true k = 0
optimum may well be 0.797 as well. It becomes evidence only once §5's bug is fixed.

**Recommended fix, in order:** (i) make the monotonicity a regression test — sweep nested orders
and assert non-increase; (ii) find why orders 3 and 6 fail; (iii) only then rebuild a separating
benchmark.

## 5b. The original (withdrawn) claim

Using `compute_polyhedral_pieces_pclf` with the **same order-3 conic partition at every node and
every order**, so the only freedom memory adds is a per-node value vector.

System (random #5 of the seeded screen, scaled by *s*):

```
B1 = [+0.677946  -0.781444;  +0.107899  -0.144061]
B2 = [+0.500597  +0.373176;  +0.704700  +0.465670]
```

| s | k = 0 (1 node) | k = 1 (2 nodes) | k = 2 (4 nodes) | verdict |
|---|---|---|---|---|
| 0.70 | 1.0216 ✗ | **0.8253 ✓** | 0.8647 | separates |
| 0.74 | 1.0800 ✗ | **0.8649 ✓** | 0.9139 | separates |
| 0.76 | 1.0101 ✗ | **0.8882 ✓** | 0.9382 | separates |
| 0.78 | 1.1383 ✗ | **0.9114 ✓** | 0.9625 | separates |
| 0.80 | 1.1675 ✗ | **0.9342 ✓** | 0.9866 | separates |
| 0.84 | 1.2259 ✗ | **0.9796 ✓** | 1.0353 | separates |

Both candidate scales were built end to end, with `max_slices = 12` and `atol = 1e-3`:

```
s = 0.76   k = 0 : bound = 1.010139  >= 1  -> no certificate, NO quotient
           k = 1 : bound = 0.888202  <  1  -> 685 cells / 2 nodes, 12 slices,
                                              Σpieces =  5797, Σfaces = 23101, 31.6 s
                                              deadend states = 0

s = 0.80   k = 0 : bound = 1.167512  >= 1  -> no certificate, NO quotient
           k = 1 : bound = 0.934181  <  1  -> 1482 cells / 2 nodes, 20 slices,
                                              Σpieces = 11100, Σfaces = 44438, 38.0 s
                                              deadend states = 8   <-- see below
```

**This is the example the paper wants.** The single-node method does not merely produce a worse
abstraction — it produces *none*, because the construction needs a contraction. Memory is the
difference between having a bisimulation and not having one, with the template held fixed so the
gain cannot be attributed to tuning.

**Recommendation: s = 0.76**, despite the thinner certificate margin. On first inspection s = 0.80
looked preferable — its k = 0 bound clears 1 by 16.7 % against only 1 % at s = 0.76 — but building
both reverses that:

- s = 0.80 produces **8 deadend states** (minimum outgoing degree 0) where s = 0.76 produces none.
  In a bisimulation of a total system every cell should have a successor under every mode, so
  deadends point at coverage lost to `atol` gaps rather than at genuine dynamics. **Unresolved.**
- s = 0.80 reports **20 slices although `max_slices = 12` was requested**, so that attribute is not
  capping what its name suggests, at least on this path. **Unresolved, and worth a look
  independently of this benchmark** — it affects every experiment that sets it.

**Also a caveat on the certificate itself.** The k = 0 bound should scale linearly in *s*, and
mostly does (0.70 → 1.0216 ≈ 0.70 × 1.4576), but s = 0.76 returns 1.0101 where linearity predicts
≈ 1.108, and several scales emit `Final LP not feasible/optimal`. There is real LP noise, so a
reported bound near 1 should not be trusted to three digits. The *separation* is robust — it holds
at all six scales tried — but individual bounds are not.

**Not usable, though fair:** with *quadratic* pieces there is no template to tune, so the
comparison is automatically honest, and the gain is real (`random #4`: 0.980694 → 0.977303 →
0.974451, separating at s ≈ 1.021). But `gamma_cover_set` accepts only `PolyhedralPiece` and
`ObserverCLFPiece`, so an ellipsoidal certificate cannot build a quotient in this codebase. Support
for it would widen the benchmark set considerably.

---

## 6. Conclusion: what the argument for the paper actually is

**Do not lead with "complexity is conserved but redistributed."** It is the weakest framing
available: if the total is conserved, a reader reasonably asks what was gained. §3a measures
conservation almost by construction — `build_common_lyapunov` packs the pieces into one function
via max/min, so each common cell *is* the intersection of the per-node cells, and the same
geometric information reappears intersected rather than separated.

Three arguments are strong, in increasing order.

### 6.1 Memory is a cost knob that provably costs nothing in accuracy

This is the point I would build the paper on, and it comes from Proposition 1 rather than from any
benchmark.

For every other knob in abstraction-based control, cost and precision trade against each other: a
coarser grid is cheaper and more conservative, a finer one is exact and expensive. **The graph is
not such a knob.** By Proposition 1 the certified controllable region is identical for every
path-complete graph, so changing the graph moves cost — cells, per-cell complexity, abstraction
size — with *no* effect on the synthesised answer.

That is unusual and worth stating as a theorem: **the abstraction can be optimised for cost without
any risk of degrading the controller.** The 50/50 probe agreement of §2 is the experimental check.

### 6.2 The single-function method has one knob, and it is provably unbounded

The predecessor method (Aydin Gol et al.) has exactly one dial: the complexity of the common
polyhedral Lyapunov function, which the partition inherits directly through its sublevel sets. By
Ahmadi & Jungers, for every *d* there are two-matrix families, stable under arbitrary switching,
admitting no polytopic Lyapunov function with ≤ *d* facets. **So that single dial must, on some
systems, be turned arbitrarily far — and the partition blows up with it.**

Path-completeness adds a second, structurally different dial: the number of nodes. §6.1 says the
second dial is free in accuracy. **That is the "more than memory" the question asks for:** not a
larger controllable set, but a second axis of the feasibility/cost trade-off that the single
function does not have.

### 6.3 Memory also shrinks the abstraction, not only the cells

Easy to miss because it is a different resource. §3b: k = 1 gives **174 abstract states against
379** for k = 0 — the abstraction *halved*. Synthesis is a fixed point over (abstract states ×
specification states), so this is a direct saving in synthesis, distinct from the per-cell
complexity of §3a. And the controller obtains its memory for free: the node is already part of the
abstract state, so the concrete controller is dynamic without any extra machinery.

### 6.4 Why no benchmark here shows the explosion, and what to try next

**No benchmark here exhibits the complexity explosion of §6.2.** Two attempts failed, and both
failures are informative:

- **Misalignment is the wrong knob.** Making the modes pull in opposite directions
  (`A1 = c R(+θ)S`, `A2 = c R(−θ)S⁻¹`, growing `cond(S)`) destroys *stability*, not certificate
  simplicity:

  | s | 1.0 | 1.5 | 2.0 | 2.5 | 3.0 | 4.0 | 5.0 |
  |---|---|---|---|---|---|---|---|
  | min order, k = 0 | 1 | 1 | 1 | >6 | >6 | >6 | >6 |
  | min order, k = 1 | 1 | 1 | 1 | >6 | >6 | >6 | >6 |
  | bound (both) | 0.6315 | 0.8063 | 0.9522 | — | — | — | — |

  Wherever a certificate exists at all, k = 0 and k = 1 need the *same* order and return the
  *identical* bound: memory buys exactly nothing. Past s = 2.5 nobody certifies, because the
  system has simply stopped contracting. Stability margin and certificate complexity are
  different axes, and only the second is the point.

  (A second sweep — fixing s = 3 and varying c — was mis-designed on my part and contributes
  nothing: s = 3 is already non-contracting, and every c tried was ≥ 0.55, i.e. too large to pull
  the JSR back below 1. It should have swept c downward.)
- **The candidate separating system was too easy.** Its JSR is ≈ 0.7571 and a common order-2
  function attains 0.757151. There is no gap for memory to exploit.

**The right constructions to try next**, in order of promise:

1. **An Ahmadi–Jungers family**, where the facet lower bound is *proved* rather than measured. This
   connects the benchmark to the theorem in §6.2 and is the only route that makes the explosion
   claim rigorous rather than anecdotal.
2. **Constrained switching**: a system unstable under arbitrary switching but stable on a
   constrained language. Then no common function exists *at any complexity* — the sharpest possible
   version of "geometry explodes" — while a PCLF over the constraint automaton exists. This needs
   the graph to encode the language, so check first whether the framework supports path-completeness
   relative to a sublanguage rather than to all words.
3. **Systems with a long spectrum-maximizing product.** The extremal norm's facet count grows with
   the s.m.p. length, so a family with lengthening s.m.p. is a natural hardness knob — and unlike
   misalignment, it keeps the stability margin fixed.

**And fix the monotonicity bug of §5 first.** Every route above compares certificates across
template complexities, which is exactly the comparison currently returning non-monotone answers.

---

## 7. Remaining limitations

Beyond §5 (the retracted benchmark and its bug) and §6.4 (why no benchmark here shows the
explosion):

- **Build-time speedups are confounded.** k = 0 took 28.1 s against 2.3 s for k = 1, but k = 0 ran
  first in both runs and absorbed first-call compilation. The §3a pair is cleaner (both warm,
  59.2 s vs 70.5 s, a 1.19× advantage) — far below the 9× a naive "cost is quadratic in faces"
  model predicts, suggesting per-cell overhead dominates at this size.
- **The invariance test is 50 probe points on one specification.** Strong evidence against a gross
  discrepancy, not a proof of exactness. Comparing winning regions cell-by-cell against the finest
  quotient would be better.
- **Proposition 1 is a sketch, and it assumes exactness.** It says nothing about what happens when
  `atol > 0` breaks it — which is the regime every actual run is in. Quantifying the gap between
  the exact winning region and the computed one is open, and is the natural formal target.
- **All benchmarks are 2-D, two-mode.** Nothing here shows how the trade scales with dimension,
  where the curse of dimensionality actually bites.
- **A second solver-side anomaly, independent of §5's bug and still unexplained**: building the
  quotient with `max_slices = 12` produced **20 slices**, and that run also showed 8 deadend states
  (minimum outgoing degree 0, where a bisimulation of a total system should have none). If
  `max_slices` does not cap what its name suggests, the truncation confound invoked in §2 may be
  larger or smaller than assumed there — so §2's explanation of the ±3 % volume spread rests on an
  attribute whose behaviour is not currently understood. The pointwise invariance result of §2 does
  **not** depend on it.

## 8. Reproducing

```
julia --project=test research/BisimulationQuotient/debruijn_sweep.jl   # §2, §3b
julia --project=test research/BisimulationQuotient/pclf_vs_clf.jl      # §3a (caches to *.jld2)
```

Per-cell figures come from `PCQ.cell_complexities(quotient)`.

The screens of §4, §5 and §6.4 were run from a session scratchpad and are **not checked in**, so
those numbers are not reproducible from this repository as it stands. Given that §5 was retracted
on the strength of one such script, the first thing worth adding here is not a benchmark but a
**regression test asserting that the bound is non-increasing over nested conic-partition orders** —
that is what would have caught the artefact immediately, and it is a prerequisite for every
construction proposed in §6.4.
