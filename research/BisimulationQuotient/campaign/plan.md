# Experiment campaign for the HSCC/Automatica submission

One folder per experiment, **self-contained**: a `plan.md` stating what we want to show, how, what to
try if the first attempt fails, and what would falsify it — plus the runnable script that produced
whatever numbers the plan quotes. Figures and caches stay gitignored.

Every script is standalone. Run any of them from anywhere:

```
julia --project=test research/BisimulationQuotient/campaign/N12_single_mode_clock/single_mode_clock.jl
```

They pick up the shared setup with `include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))`, so
the benchmark problems, module aliases and MOI call sequences come from
[`../common.jl`](../common.jl) rather than being duplicated. Scripts that read a `*.jld2` cache say so
and exit cleanly when it is absent, since caches are gitignored.

| folder | script | what it measures |
| :--- | :--- | :--- |
| N1 | `incomplete_speedup.jl` | cells, facets and rate for single node vs dual vs plain De Bruijn |
| N3 | `smp_predictor.jl` | s.m.p. against the memory gap, plus the Gol–Lazar–Belta budget sweep |
| N8 | `induced_common_cost.jl` | determinization cost, observer size, per-level geometry cost |
| N10 | `constrained_window.jl` | the window of $c$ where $\{cA_1, A_2\}$ needs constrained switching |
| N11 | `template_diversity.jl` | whether an incomplete graph admits distinct per-node templates |
| N12 | `single_mode_clock.jl` | the separation for a single linear system |
| N13 | `conservatism.jl` | tighter rate at equal budget, with the baseline search curve |
| N14 | `1_`..`5_`, `plot_1_`..`plot_3_` | grid vs bisimulation on reachability and on co-safe LTL |
| N15 | `deadend_audit.jl`, `monotonicity.jl`, `simulation_check.jl` | the three soundness checks |

The reasoning behind the ranking is in [`../notes/paper-advantages.md`](../notes/paper-advantages.md).
This file is the index and the run order.

## Why a new campaign

Three measurements made in September 2026 reset what the experiments should be:

1. **Gol–Lazar–Belta's system is memory-free.** At conic order 1, 2 and 3 (4, 8, 16 rows per node),
   one node, two nodes and four nodes all certify at *identical* rates to six digits. No choice of
   template budget, graph order, domain size or slice count can make it demonstrate the mechanism.
2. **A 2-D speedup already exists** and had never been timed warm: the dual De Bruijn graph of order
   1 builds in 0.72 s against 1.25 s for the single node, at an identical certified rate, with 174
   cells against 379.
3. **The facet totals carry no information beyond the piece counts** in the plane — every arm
   measures 4.00 facets per convex piece, which is Euler's formula. The per-cell *distribution* is
   the real content.

## The experiments

| # | folder | one line | status |
| :-- | :--- | :--- | :--- |
| N1 | [`N1_incomplete_speedup/`](N1_incomplete_speedup/) | incomplete lifting is faster at an identical rate | pilot measured, needs a campaign |
| N2 | [`N2_feasibility_separation/`](N2_feasibility_separation/) | at 4 facets/node the common construction does not exist | done, needs writing up |
| N3 | [`N3_smp_predictor/`](N3_smp_predictor/) | the s.m.p. predicts whether memory can help | not started, trivial |
| N4 | [`N4_dimension_sweep/`](N4_dimension_sweep/) | the facet advantage should switch on at `n = 4` | not started |
| N5 | [`N5_cycle_family_timing/`](N5_cycle_family_timing/) | bounded per-node geometry vs a `4q`-facet common, on the clock | not started |
| N6 | [`N6_online_cost/`](N6_online_cost/) | the *controller* is cheaper, not just the build | not started, cheap |
| N7 | [`N7_build_profile/`](N7_build_profile/) | where the build time actually goes | **step 3 done: set difference is $F^{2.2}$–$F^{2.7}$, so mechanism B is live**; steps 1–2 (profile a real build) next |
| N8 | [`N8_completeness_shortcut/`](N8_completeness_shortcut/) | stop handicapping route 1 | not started, a fix |
| N9 | [`N9_piece_conservation/`](N9_piece_conservation/) | why `Σ` pieces is conserved on some certificates | not started |
| N10 | [`N10_constrained_switching/`](N10_constrained_switching/) | their matrices, unstable without the constraint, route 1 cannot start | **done: ∃ 44.2→37.0 %, ∀ 18.0→26.7 %, opposite directions as predicted**; needed two library fixes; 3 proof obligations open |
| N11 | [`N11_system_families/`](N11_system_families/) | the four mechanisms that make memory pay, and how to design families | taxonomy written, one family tested |
| N12 | [`N12_single_mode_clock/`](N12_single_mode_clock/) | the separation holds for a **single** linear system; the node is a clock | **done, FEASIBILITY only**: GLB 4-facet node rate 1.134 ≥ 1 (cannot start), 5-cycle certifies at ρ=0.9, 1 694-cell quotient, ∃≡∀ exactly. **Route 1 (induced common) wins on cost 1.32× — the cycle is a complete graph, so mechanism A is off** |
| N13 | [`N13_conservatism_separation/`](N13_conservatism_separation/) | a tighter rate at equal per-node budget — 0.868964 vs 0.873586 after a fair 41-draw search | **done; the gap is 0.0046, ten times smaller than the single-orientation figure** |
| N14 | [`N14_bisimulation_vs_simulation/`](N14_bisimulation_vs_simulation/) | **bisimulation vs alternating simulation as the spec hardens** — the grid is fine on reachability and collapses on co-safe LTL | **done, result established** |
| N15 | [`N15_soundness_regressions/`](N15_soundness_regressions/) | **is the object actually a bisimulation?** three open anomalies | **done: all three closed, 0 violations in 8 367 samples** |

See [`coverage-review.md`](coverage-review.md) for why N13–N15 were added and what remains uncovered.

## Run order

0. **N15 (soundness), then N14 (bisimulation vs simulation).** Added after the coverage review. Nothing else
   matters if the object is not a bisimulation, and no PCLF-vs-common ratio answers the question a
   reader asks first, which is why not just grid it. These outrank every cost experiment below.
1. **N2, N3, N12** — cost nothing, can be written today. N2 is the safest result in the folder, and
   **N12 shows N2's separation does not need switching at all**: for the single linear system
   $x^+ = \rho R(2\pi/q)x$, one node at four facets certifies nothing (rate 1.07–1.34 over
   $q = 5\ldots8$, $\rho = 0.85\ldots0.95$) while the $q$-cycle attains $\rho$ exactly. Decide early
   whether that becomes the opening example or a remark — it changes the paper's shape.
2. **N1** — best return on effort; the pilot is already positive.
3. **N10** — the clearest statement of what memory buys, on Gol–Lazar–Belta's own matrices.
4. **N7** — gates N4 and N5. Do not invest in either before knowing what fraction of the build is
   facet-dependent.
5. **N8** — a fairness fix; do it before any timing table is published.
6. **N4, N5** — only if N7 says the facet term is live.
7. **N6, N9** — independent, cheap, do whenever.

## Two vocabulary rules for every plan here

**Route 1** = induce the common Lyapunov function from the PCLF (`build_common_lyapunov`) and run the
single-node construction. **Route 2** = run the construction on the PCLF directly. Both end with a
bisimulation certifying the same thing (`thm:graph-invariance`); they differ only in what they
compute on the way. Comparisons against a *foreign* certificate — Gol–Lazar–Belta's published `L`,
for instance — are validation runs, never cost benchmarks, because a certificate determines its own
slice family and terminal set.

**Mechanism A** = route 2 has *fewer cells*. **Mechanism B** = same total, *simpler cells*; needs cost
superlinear in per-cell facet count. Every plan below says which mechanism it targets. Mechanism B's
premise is now measured — [N7](N7_build_profile/plan.md) puts `set_difference_decompose` at
$F^{2.2}$–$F^{2.7}$ — so the open question is not whether facets cost, but how often cells are fat.

### When mechanism A is active — the sharpened rule

The old rule, "needs an incomplete graph", is the right conclusion from a vague premise, and it is not
predictive: it does not say how much, and it does not explain why a complete graph loses. The property
that actually matters is **how much of the alphabet each node still has to handle**.

Refining at node $s$ takes pre-images under the modes on $s$'s outgoing edges. So

**The full discussion — what decides a win, why the dual De Bruijn graph wins where the plain one
loses (memory of the *future* rather than of the past), the per-node decomposition, and what remains
unmeasured — is centralised in [`N1_incomplete_speedup/plan.md`](N1_incomplete_speedup/plan.md)** and
deliberately not repeated here. The short form:

$$\text{net} \;=\; \frac{\text{route 1 cells}}{\text{route 2 cells per node}} \Big/\; |S|,$$

so route 2 wins exactly when the per-node partition shrinks by more than $|S|$. Measured on
`two_mode_problem` at an identical rate, with route 1 = the induced common:

| oracle PCLF | route 1 | route 2 | route 2 per node | per-node gain | net |
| :--- | --: | --: | --: | --: | :--- |
| **dual** De Bruijn $k{=}1$ | 379 | **174** | 87.0 | 4.36× | **wins 2.18×** |
| **plain** De Bruijn $k{=}1$ | 397 | 775 | 387.5 | 1.02× | loses 1.95× |

A quantitative predictor of the form $\sum_s k_s^2 / M^2$ was proposed, fitted three retrospective
points, and **failed its first prospective test** — but that test used a *degenerate family* (modes
that are exact negatives of one another under symmetric templates), so it is not evidence either way.
Both the model's independent flaw and the broken family are documented in N1's plan. **There is no
validated quantitative predictor; do not quote one.**

### The two benefits are orthogonal, and must be claimed separately

A node can do two different things, and only the second reduces cost:

- **record the past** — different $V_s$ per node, which can certify a *better rate*. Available on a
  complete graph: [N13](N13_conservatism_separation/plan.md) measures 0.868964 against 0.873586 at
  equal per-node budget.
- **constrain the future** — the node limits which modes come next, which is what shrinks the
  partition. Requires $|\mathrm{enabled}(s)| < M$.

So a complete graph can buy a better rate *while costing more cells*, and an incomplete one buys fewer
cells *at an identical rate*. Claiming "memory helps" as a single proposition invites a referee to
test it on a complete graph and find it false. This is the same property as `PCLF.restricts_future`,
which decides **soundness** in [N10](N10_constrained_switching/plan.md)'s verification path and
**cost** here.
