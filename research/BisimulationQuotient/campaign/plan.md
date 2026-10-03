# Experiment campaign for the HSCC/Automatica submission

One folder per experiment, **self-contained**: a `plan.md` stating what we want to show, how, what to
try if the first attempt fails, and what would falsify it — plus the runnable script that produced
whatever numbers the plan quotes. Figures and caches stay gitignored.

Every script is standalone. Run any of them from anywhere:

```
julia --project=test research/BisimulationQuotient/campaign/N1_quotient_cost/comparison.jl
```

They pick up the shared setup with `include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))`, so
the benchmark problems, module aliases and MOI call sequences come from
[`../common.jl`](../common.jl) rather than being duplicated. Scripts that read a `*.jld2` cache say so
and exit cleanly when it is absent, since caches are gitignored.

| # | folder | scripts |
| :-- | :--- | :--- |
| N1 | [`N1_quotient_cost/`](N1_quotient_cost/) | `comparison.jl` |
| N2 | [`N2_constrained_switching/`](N2_constrained_switching/) | `certificate_exists.jl`, `constrained_window.jl`, `dwell_time_example.jl`, `glb_constraint_moves_both_ways.jl`, `plot_dwell_threshold.jl`, `quotient_on_dwell.jl`, `quotient_on_glb_constrained.jl` |
| N3 | [`N3_bisimulation_vs_simulation/`](N3_bisimulation_vs_simulation/) | `1_`…`5_`, `plot_1_`…`plot_3_` |
| N4 | [`N4_soundness_regressions/`](N4_soundness_regressions/) | `deadend_audit.jl`, `monotonicity.jl`, `simulation_check.jl` |
| N5 | [`N5_feasibility_separation/`](N5_feasibility_separation/) | plan only — the result is analytic, and [`../experiments/memory_vs_geometry.jl`](../experiments/memory_vs_geometry.jl) is the runnable form |
| N6 | [`N6_single_mode_clock/`](N6_single_mode_clock/) | `1_certificate_separation.jl`, `2_quotient_on_clock.jl`, `3_route1_vs_route2.jl` |
| N7 | [`N7_conservatism_separation/`](N7_conservatism_separation/) | `conservatism.jl` |
| N8 | [`N8_build_profile/`](N8_build_profile/) | `1_primitive_scaling.jl` |
| N9 | [`N9_system_families/`](N9_system_families/) | `template_diversity.jl` |
| N10 | [`N10_dimension_sweep/`](N10_dimension_sweep/) | plan only |
| N11 | [`N11_graph_order_sweep/`](N11_graph_order_sweep/) | plan only |
| N12 | [`N12_smp_predictor/`](N12_smp_predictor/) | `smp_predictor.jl` |
| N13 | [`N13_online_cost/`](N13_online_cost/) | plan only |

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
| N1 | [`N1_quotient_cost/`](N1_quotient_cost/) | building on the PCLF directly costs less than determinising first | done — the three experiments are in `research/HSCC2027/` |
| N2 | [`N2_constrained_switching/`](N2_constrained_switching/) | their matrices, unstable without the constraint, route 1 cannot start | **done: ∃ 44.2→37.0 %, ∀ 18.0→26.7 %, opposite directions as predicted**; needed two library fixes; 3 proof obligations open |
| N3 | [`N3_bisimulation_vs_simulation/`](N3_bisimulation_vs_simulation/) | **bisimulation vs alternating simulation as the spec hardens** — the grid is fine on reachability and collapses on co-safe LTL | **done, result established** |
| N4 | [`N4_soundness_regressions/`](N4_soundness_regressions/) | **is the object actually a bisimulation?** three open anomalies | **done: all three closed, 0 violations in 8 367 samples** |
| N5 | [`N5_feasibility_separation/`](N5_feasibility_separation/) | at 4 facets/node the common construction does not exist | done, needs writing up |
| N6 | [`N6_single_mode_clock/`](N6_single_mode_clock/) | the separation holds for a **single** linear system; the node is a clock | **done, FEASIBILITY only**: GLB 4-facet node rate 1.134 ≥ 1 (cannot start), 5-cycle certifies at ρ=0.9, 1 694-cell quotient, ∃≡∀ exactly. **Route 1 (induced common) wins on cost 1.32× — the cycle is a complete graph, so mechanism A is off** |
| N7 | [`N7_conservatism_separation/`](N7_conservatism_separation/) | a tighter rate at equal per-node budget — 0.868964 vs 0.873586 after a fair 41-draw search | **done; the gap is 0.0046, ten times smaller than the single-orientation figure** |
| N8 | [`N8_build_profile/`](N8_build_profile/) | where the build time actually goes | **step 3 done: set difference is $F^{2.2}$–$F^{2.7}$, so mechanism B is live**; steps 1–2 (profile a real build) next |
| N9 | [`N9_system_families/`](N9_system_families/) | the four mechanisms that make memory pay, and how to design families | taxonomy written, one family tested |
| N10 | [`N10_dimension_sweep/`](N10_dimension_sweep/) | the facet advantage is dimensional — the $\mathbb{Z}_q$ family swept over $n$ | not started; the extension to run first |
| N11 | [`N11_graph_order_sweep/`](N11_graph_order_sweep/) | what happens as the De Bruijn order grows | not started; plan only |
| N12 | [`N12_smp_predictor/`](N12_smp_predictor/) | the s.m.p. predicts whether memory can help | not started, trivial |
| N13 | [`N13_online_cost/`](N13_online_cost/) | the *controller* is cheaper, not just the build | not started, cheap |

See [`coverage-review.md`](coverage-review.md) for why N7, N3 and N4 were added and what remains uncovered.

## Run order

**The numbering is the order.** N1–N3 are the paper's experiments, N4–N7 are finished results waiting
to be written up, N8–N9 are half measured, and N10–N13 have not started. N10 and N11 are the two
directions the work extends in: higher dimension, and higher graph order.

Four things the numbering does not say:

* **N8 gates N10.** Do not invest in the dimension sweep before knowing what fraction of the build is
  facet-dependent. N8 has already measured `set_difference_decompose` at $F^{2.2}$–$F^{2.7}$, so the
  superlinearity is established; what is missing is the fraction.
* **N9 is the prerequisite for N10 and N11**, not an experiment in its own right: it says how to
  build a family where memory pays, and it carries the trap for lifting one to $n$ dimensions.
* **N6 shows N5's separation does not need switching at all**: for the single linear system
  $x^+ = \rho R(2\pi/q)x$, one node at four facets certifies nothing (rate 1.07–1.34 over
  $q = 5\ldots8$, $\rho = 0.85\ldots0.95$) while the $q$-cycle attains $\rho$ exactly. Whether that
  becomes the opening example or a remark changes the paper's shape.
* **N4 outranks everything below it.** Nothing else matters if the object is not a bisimulation, and
  no PCLF-against-common ratio answers the question a reader asks first, which is why not just grid
  it.

## Two vocabulary rules for every plan here

**Route 1** = induce the common Lyapunov function from the PCLF (`build_common_lyapunov`) and run the
single-node construction. **Route 2** = run the construction on the PCLF directly. Both end with a
bisimulation certifying the same thing (`thm:graph-invariance`); they differ only in what they
compute on the way. Comparisons against a *foreign* certificate — Gol–Lazar–Belta's published `L`,
for instance — are validation runs, never cost benchmarks, because a certificate determines its own
slice family and terminal set.

**Mechanism A** = route 2 has *fewer cells*. **Mechanism B** = same total, *simpler cells*; needs cost
superlinear in per-cell facet count. Every plan below says which mechanism it targets. Mechanism B's
premise is now measured — [N8](N8_build_profile/plan.md) puts `set_difference_decompose` at
$F^{2.2}$–$F^{2.7}$ — so the open question is not whether facets cost, but how often cells are fat.

### When mechanism A is active — the sharpened rule

The old rule, "needs an incomplete graph", is the right conclusion from a vague premise, and it is not
predictive: it does not say how much, and it does not explain why a complete graph loses. The property
that actually matters is **how much of the alphabet each node still has to handle**.

Refining at node $s$ takes pre-images under the modes on $s$'s outgoing edges, so what governs the
growth is how many modes a node must still serve — not whether the graph is complete.

**The full discussion — what decides a win, why the dual De Bruijn graph wins where the plain one
loses (memory of the *future* rather than of the past), the per-node decomposition, and what remains
unmeasured — is centralised in [`N1_quotient_cost/plan.md`](N1_quotient_cost/plan.md)** and
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
  complete graph: [N7](N7_conservatism_separation/plan.md) measures 0.868964 against 0.873586 at
  equal per-node budget.
- **constrain the future** — the node limits which modes come next, which is what shrinks the
  partition. Requires $|\mathrm{enabled}(s)| < M$.

So a complete graph can buy a better rate *while costing more cells*, and an incomplete one buys fewer
cells *at an identical rate*. Claiming "memory helps" as a single proposition invites a referee to
test it on a complete graph and find it false. This is the same property as `PCLF.restricts_future`,
which decides **soundness** in [N2](N2_constrained_switching/plan.md)'s verification path and
**cost** here.
