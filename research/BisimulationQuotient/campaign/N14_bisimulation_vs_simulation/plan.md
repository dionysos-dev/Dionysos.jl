# N14 — Bisimulation versus alternating simulation, as the specification gets harder

**Status: the result is established and reproducible.** Intended as the backbone of a paper extension
on why exactness matters once the specification is more than reachability.

---

## The result, in one paragraph

> A uniform grid abstraction is an **alternating simulation**: sound but over-approximating, and its
> spurious transitions give an adversary choices the real system does not have. This is survivable
> for reachability — the abstraction converges to the exact answer from below as the grid refines —
> and it is **fatal for specifications that constrain the order of events**. On Gol–Lazar–Belta's own
> system and their own co-safe LTL formula, a 26 751-cell grid certifies **0.31 %** of the domain and
> *falls* under refinement, while the 10 611-cell bisimulation certifies **83 %**. The reason
> refinement cannot rescue it is structural: for linear dynamics the branching factor of a grid
> abstraction with tight box over-approximation is
> $\prod_i \bigl(\sum_j |A_{ij}| + 1\bigr)$, **independent of the grid step**. The bisimulation has
> no spurious branching at all, and the Lyapunov slice decrease makes liveness obligations
> structural rather than something the fixed point must discover.

---

## The evidence, in the order it was built

Every run uses the **same domain** for both abstractions — the certificate's outermost sublevel set,
handed to the grid as its state set — and gives the grid the **best possible over-approximation**: for
a linear map the image of a box is a parallelogram whose tight bounding box is exactly
$\mathrm{Hyperrectangle}(Ac, |A|r)$, so `USER_DEFINED` is given that map. No growth-bound slack, no
linearization error. Targets use `MP.INNER`, which under-approximates them and so cannot inflate the
grid's certified set. The comparison is **pointwise on a probe grid**, never by volume, because
`max_slices` truncation makes covered volumes differ for reasons unrelated to the method.

### 1. Reachability on the two-mode system — the grid works, and the quotient is smaller

[`1_reachability_two_mode.jl`](1_reachability_two_mode.jl) · [figure](fig_1_reachability_two_mode.png)

| $h$ | grid cells | probes | both | quotient only | **GRID ONLY** |
| --: | --: | --: | --: | --: | --: |
| 0.40 | 209 | 992 | 4 | 49 | **0** |
| 0.20 | 919 | 992 | 27 | 26 | **0** |
| 0.10 | 3 839 | 992 | 31 | 22 | **0** |
| 0.05 | **15 701** | 992 | 43 | 10 | **0** |

Containment holds at every resolution — the falsifier stayed at zero. The grid climbs
4 → 27 → 31 → 43 toward the quotient's constant 53. **15 701 grid cells against 379, and the grid
still misses 10 of 53.**

A branching check runs first, because the two sides resolve nondeterminism differently and must: the
grid's branching is *spurious* and so is resolved adversarially, while a quotient's branching, when
it has any, is *genuine* graph structure. Measured here: **max 1 successor per (state, mode), 0 of
757 pairs branch**, so ∃ and ∀ coincide on the quotient side and comparing against the grid's worst
case is exact.

### 2. Reachability on Gol–Lazar–Belta — still works, on harder dynamics

[`2_reachability_glb.jl`](2_reachability_glb.jl) · [figure](fig_2_reachability_glb.png)

| target | $h$ | cells | certified | area | % of domain |
| :-- | --: | --: | --: | --: | --: |
| $D$ | 0.500 | 1 603 | 1 603 | 400.75 | 94.5 % |
| $D$ | 0.250 | 6 591 | 6 591 | 411.94 | 97.1 % |
| $D$ | 0.125 | 26 751 | 26 751 | 417.98 | 98.6 % |
| $R_1$ | 0.500 | 1 603 | 400 | 100.0 | 23.6 % |
| $R_1$ | 0.250 | 6 591 | 2 492 | 155.8 | 36.7 % |
| $R_1$ | 0.125 | 26 751 | 11 852 | 185.2 | **43.7 %** |

Monotone, converging from below, and the figure shows the certified region taking the *same shape* as
the bisimulation's. The grid is behaving exactly as an over-approximation should.

### 3. Their co-safe LTL formula — the grid collapses

[`3_cosafe_glb.jl`](3_cosafe_glb.jl) · [figure](fig_3_cosafe_glb.png)

$$\varphi = (\lnot R_2 \mathbin{U} D) \wedge F(R_1) \wedge \bigl((R_3 \to X \lnot R_1) \mathbin{U} D\bigr)$$

| $h$ | cells | certified | area | % of domain |
| --: | --: | --: | --: | --: |
| 1.0000 | 377 | 13 | 13.00 | 3.07 % |
| 0.5000 | 1 603 | 40 | 10.00 | 2.36 % |
| 0.2500 | 6 591 | 48 | 3.00 | 0.71 % |
| 0.1250 | 26 751 | 83 | 1.30 | 0.31 % |
| 0.0625 | **107 791** | 100 | 0.39 | **0.09 %** |
| **bisimulation** | **10 611** | **8 794 cells** | — | **83 %** |

**Identical abstraction, identical domain, identical branching — only the specification differs, and
the grid goes from 43.7 % to 0.71 %.** At 107 791 cells, ten times the bisimulation's, it manages
0.09 %.

The certified area decays roughly in proportion to $h$ (area/$h$ ≈ 13, 20, 12, 10, 6), the signature
of a **boundary layer**: the fixed point propagates about one step out from the accepting set and no
further.

### 4. Why refinement cannot help — branching is scale-invariant

[`4_branching_is_scale_invariant.jl`](4_branching_is_scale_invariant.jl)

| $h$ | cells | transitions | **mean successors** | max | cells with no transition |
| --: | --: | --: | --: | --: | --: |
| 1.000 | 377 | 3 464 | **4.59** | 6 | 0 |
| 0.500 | 1 603 | 14 802 | **4.62** | 6 | 0 |
| 0.250 | 6 591 | 60 778 | **4.61** | 6 | 0 |
| 0.125 | 26 751 | 246 770 | **4.61** | 6 | 0 |

Flat at 4.61 across a 70× increase in cell count, matching the closed form
$\prod_i(\sum_j |A_{ij}| + 1) = 1.97 \times 2.34 = 4.61$ exactly. **The image box and the cell scale
together, so refining a uniform grid buys precision in *where* states are and none at all in *how
nondeterministic* the abstraction is.** The adversary keeps ~4.6 choices per step forever.

By contrast the bisimulation quotients measure **1 successor per (state, mode)** whenever the graph is
deterministic (0 of 21 220 pairs branch on the Gol–Lazar–Belta certificate). A quotient branches only
when its *graph* does, and that branching is real structure, not error.

### 5. Confirmation that refining is not merely slow but counterproductive

[`5_refinement_does_not_help.jl`](5_refinement_does_not_help.jl) — the $h = 0.25 \ldots 0.0625$ rows
of §3, with the terminal set measured at area **141.3 of 424** (equivalent radius 6.7) against a
residual quantization uncertainty of only $h/(2(1-\rho)) = 3.77h = 0.24$ at the finest step.
Resolution is nowhere near the binding constraint; the target is a third of the domain and the grid
still cannot reach it.

---

## The mechanism, stated for the paper

Spurious transitions cannot stop the system from **reaching** a set — the contraction guarantees
arrival, and a sound over-approximation only adds arrival routes. What they destroy is **order**.
$\varphi$ is a conjunction demanding a sequence: visit $R_1$; never enter $R_2$ before $D$; never hit
$R_1$ in the step after $R_3$; discharge everything before $D$. Each of the ~4.6 spurious branches per
step is another opportunity for the adversary to violate the ordering, and no amount of refinement
removes them.

The bisimulation faces none of this. Its partition is built **from** the certificate rather than laid
over it, so it carries no spurious branching, every cell has one well-defined observation label, and
the Lyapunov slice index strictly decreases along every transition — which makes "eventually reaches
$D$" a structural property of the construction rather than something a fixed point must establish
against an adversary.

> **Over-approximating abstractions are adequate for reachability and collapse on temporal
> specifications that constrain the order of events. Exactness is not a refinement of precision; it
> is what makes ordering enforceable at all.**

---

## Reproducing

**Run from the repository root**, since `--project=test` is resolved relative to the working
directory. The scripts themselves locate `common.jl` and the caches relative to their own file, so
only the project flag cares where you are.

```
cd <repo root>
D=research/BisimulationQuotient/campaign/N14_bisimulation_vs_simulation

julia --project=test $D/1_reachability_two_mode.jl        # containment, two-mode system
julia --project=test $D/2_reachability_glb.jl             # branching constants + GLB reachability
julia --project=test $D/3_cosafe_glb.jl                   # their formula (needs the E8 cache)
julia --project=test $D/4_branching_is_scale_invariant.jl # the mechanism
julia --project=test $D/5_refinement_does_not_help.jl     # refinement sweep to h = 0.0625
julia --project=test $D/plot_1_reachability_two_mode.jl   # and plot_2_, plot_3_
```

`3_cosafe_glb.jl` and the GLB plots read `examples/gol_lazar_belta_pclf.jld2`; run
`examples/gol_lazar_belta_pclf.jl` first if it is missing. Figures are gitignored.

---

## Open, and honest about it

- **A second system for §3.** The collapse is measured on one system with one formula. Before this
  carries a paper, repeat it on the two-mode system with a comparable ordering formula, and ideally
  on a third system, to show it is a property of the specification class rather than of
  Gol–Lazar–Belta's geometry.
- **Which clause does the damage.** $\varphi$ has three conjuncts. Running each alone would say
  whether the collapse comes from the `U` obligations, from the `X` constraint, or only from their
  conjunction — a much sharper statement than "co-safe LTL".
- **The closed form for branching** is derived for tight box over-approximation of a linear map and
  confirmed numerically to three digits. State it as a proposition with that hypothesis, not as a
  general fact about grid abstractions.
- **`D` is the certificate's terminal set**, not problem data, so the grid has no natural counterpart
  and is handed the PCLF's. Say so plainly — the specification is phrased in the certificate's
  vocabulary, which is itself a point in the method's favour rather than something to hide.
- **Five explanations were offered and falsified** before the right one (polarity of the atomic
  propositions; semantics not reaching the solver; the domain shrinking under refinement; a
  probe-comparison artefact; out-of-domain handling, which turned out to be runtime quantization).
  Recorded so they are not retried.
