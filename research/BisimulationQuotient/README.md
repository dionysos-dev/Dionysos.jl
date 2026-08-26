# Bisimulation quotient (HSCC 2027)

Scripts reproducing the experiments of our **HSCC 2027** submission on path-complete Lyapunov
function (PCLF) based bisimulation quotients for switched systems.

> **Citation:** _add the full reference once the paper is finalized._

## Layout

[`common.jl`](common.jl) holds everything the experiments share — the imports and module
aliases, the two benchmark problems, the MOI call sequences for building a quotient
(`build_quotient`) and synthesising a co-safe LTL controller on it (`synthesize_cosafe_ltl`),
and the figures. Each script `include`s it and then contains only what makes that experiment
different. Run any of them directly:

```
julia --project=test research/BisimulationQuotient/paper_example_3_1.jl
```

## Experiments

| Script | What it varies |
| :--- | :--- |
| [`paper_example_3_1.jl`](paper_example_3_1.jl) | Example 3.1 of the paper: the three-region co-safe LTL specification. Reference run: 10 611 states, controllable volume 186.728. |
| [`debruijn_graph.jl`](debruijn_graph.jl) | A De Bruijn path-complete graph of order 1, with a rotated template per node. Includes the lifted 3D view. |
| [`custom_graph.jl`](custom_graph.jl) | The same system over a hand-written path-complete graph, with a different rotation per node. |
| [`debruijn_sweep.jl`](debruijn_sweep.jl) | The De Bruijn order swept over k = 0, 1, 2 — 1, 2 and 4 nodes. |
| [`pclf_vs_clf.jl`](pclf_vs_clf.jl) | A path-complete certificate against the common Lyapunov function induced by it, histogrammed by per-cell complexity. |
| [`rotated_state_space.jl`](rotated_state_space.jl) | A state-space polytope rotated off the axes. |

## Notes

**Caching.** Building a quotient costs minutes, so the expensive scripts save the solved
optimizer to a `*.jld2` beside themselves (gitignored) and reload it. `pclf_vs_clf.jl` rebuilds
only when its caches are missing — pass `force = true` to `compute_and_cache` after changing the
problem or the tolerances.

**`atol`.** The inset applied when one polytope is cut from another, and the tolerance that most
affects the result. It is bounded on both sides: at `1e-3` roughly 3.7 % of the domain is
silently dropped, at `1e-6` degenerate flat pieces survive and break plotting, and `0` diverges.
`1e-4` loses about 0.4 % and is the default.

**3D figures.** The lifted views use the `DionysosMakieExt` recipes
(`Dionysos.plot_lifted_bisimulation!` / `Dionysos.plot_lifted_trajectory!`); load a Makie backend
(`using GLMakie` for an interactive window, `using CairoMakie` for a static figure) to enable
them.
