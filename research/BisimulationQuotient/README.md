# Bisimulation quotient (HSCC 2027)

Scripts reproducing the experiments of our **HSCC 2027** submission on path-complete Lyapunov
function (PCLF) based bisimulation quotients for switched systems.

> **Citation:** _add the full reference once the paper is finalized._

## Layout

[`common.jl`](common.jl) holds everything the experiments share — the imports and module aliases,
the two benchmark problems, the MOI call sequences for building a quotient (`build_quotient`) and
synthesising a co-safe LTL controller on it (`synthesize_cosafe_ltl`), and the figures. Each script
`include`s it and then contains only what makes that experiment different. Run any of them
directly:

```
julia --project=test research/BisimulationQuotient/gol_lazar_belta_pclf.jl
```

## What is claimed, and what backs it

Four claims, four scripts. The evidence is of different kinds and the table says which is which,
because they carry different weight: a construction settles a claim, a benchmark supports one.

| # | Script | Claim it addresses | Evidence | Headline numbers |
| :-- | :--- | :--- | :--- | :--- |
| E1 | [`graph_invariance.jl`](graph_invariance.jl) | The certified answer does not depend on the graph | proved + measured | 4 graphs, 101 probes, **0 disagreements**; cells span 9.3×, facets 6.7× |
| E2 | [`complexity_distribution.jl`](complexity_distribution.jl) | Geometric complexity is redistributed, not removed | statistical | Σ facets conserved to **0.7 %**; worst cell **146 vs 1 728** (11.8×) |
| E3 | [`memory_vs_geometry.jl`](memory_vs_geometry.jl) | With a bounded per-node budget, memory decides whether an abstraction exists | formal construction | at 4 facets/node: one node **rate ≥ 1**, no quotient; ℤ_q cycle **rate = ρ**, 1 797 cells |
| E4 | [`gol_lazar_belta_example.jl`](gol_lazar_belta_example.jl) | The construction contains the predecessor as `\|S\| = 1`, and answers both their Problem 5.1 (synthesis, ∃) and 5.2 (verification, ∀) from the one quotient — the two runs differ only in who the system says owns the switching signal | validation against published numbers | rate **0.940008** (exact LP) vs their 0.94; **11** slices vs their 11; 6 332 cells; winning **4 379** cells (vol 91.88) ⊋ verified **1 867** cells (vol 51.79) |
| E5 | [`lifted_showcase.jl`](lifted_showcase.jl) | What the lifting looks like when the memory is actually used | figure for the paper | 10 steps, **9 node changes**, all 5 nodes visited (`1→2→3→2→3→4→0→1→2→1`); writes a planar and a 3-D figure |

Run order: E2 reads the caches written by [`pclf_vs_clf.jl`](pclf_vs_clf.jl), so run that first. The
others are independent.

## What the experiments do *not* show

Recorded so the claims are not read as wider than they are.

**No separation on generic systems.** A separation — the common certificate failing where a
path-complete one succeeds at equal per-node budget — was searched for over 24 random planar
systems, 20 random systems in dimensions 3 and 4, and the predecessor's own Example 3.1 at two
template budgets. **None was found.** For generic systems the common certificate already attains
the joint spectral radius to several digits, so the facet budget `𝔉(Σ)` is small and there is
nothing to redistribute. The separation of `memory_vs_geometry.jl` is a *construction*; it is not
a property of typical systems.

**Cell count does not improve for complete graphs — it degrades, provably.** Lifting along a
complete graph whose nodes carry the same piece replicates the partition exactly, so the quotient
has `|S|` times the states of the single-node one and, by node obliviousness, the copies are
pairwise bisimilar. Measured on the De Bruijn family: 379 / 775 / 1 613 cells for `|S|` = 1 / 2 / 4,
that is 1.00 / 2.05 / 4.26. The benefit of a complete lifting is confined to the geometry of the
individual cells.

**One reduction remains unexplained.** The dual De Bruijn graph of order 1 — incomplete, identical
pieces — yields 174 cells against 379 for the single node, a genuine reduction. A mode-committed
graph sharing those properties yields an *increase* (305 against 93) on the same problem. The
characterisation of when lifting reduces the quotient is open; an explanation based on the mode set
of a node was tried and refuted.

**The reproduction is partial.** E4 recovers the predecessor's contraction rate (0.9400 against
their 0.94) and slice count (11 against their 11) exactly, but not their reported quotient size:
they state **9 677** states for this example and we obtain **6 332**, some 33 % fewer. Sweeping
`atol` over `1e-3 … 5e-5` moves the count only from 6 332 to 6 468, so lost coverage does not
explain it. The cause is unidentified; candidates are the transcription of the three observation
polytopes, which the paper gives only as a figure, and a difference in how states are counted.
The certificate and the slice structure reproduce; the cell count does not.

**Two solver-side anomalies are open.** `max_slices` does not cap the number of slices built, and
some quotients contain deadend states (minimum outgoing degree 0) where a bisimulation of a total
system should have none. Both look like `atol` coverage loss rather than dynamics.

## Worked examples

| Script | What it varies |
| :--- | :--- |
| [`gol_lazar_belta_pclf.jl`](gol_lazar_belta_pclf.jl) | **Same problem as E4**, with *our* certificate (2 nodes) instead of theirs, plus the full pipeline: quotient, three-region co-safe LTL synthesis, controller, trajectory. The paper's reference run: 10 611 states, controllable volume 186.728. E4 validates against the predecessor; this shows the method working. Their cell counts are **not** comparable — different certificates induce different slice families and terminal sets. |
| [`debruijn_graph.jl`](debruijn_graph.jl) | A De Bruijn path-complete graph of order 1, with a rotated template per node. Includes the lifted 3D view. |
| [`custom_graph.jl`](custom_graph.jl) | The same system over a hand-written path-complete graph, with a different rotation per node. |
| [`debruijn_sweep.jl`](debruijn_sweep.jl) | The De Bruijn order swept over k = 0, 1, 2 — 1, 2 and 4 nodes. |
| [`pclf_vs_clf.jl`](pclf_vs_clf.jl) | A path-complete certificate against the common Lyapunov function induced by it, histogrammed by per-cell complexity. |
| [`rotated_state_space.jl`](rotated_state_space.jl) | A state-space polytope rotated off the axes. |

## Notes

**Caching.** Building a quotient costs minutes, so the expensive scripts save the solved optimizer
to a `*.jld2` beside themselves (gitignored) and reload it. `pclf_vs_clf.jl` rebuilds only when its
caches are missing — pass `force = true` to `compute_and_cache` after changing the problem or the
tolerances.

**`atol`.** The inset applied when one polytope is cut from another, and the tolerance that most
affects the result. It is bounded on both sides: at `1e-3` roughly 3.7 % of the domain is silently
dropped, at `1e-6` degenerate flat pieces survive and break plotting, and `0` diverges. `1e-4`
loses about 0.4 % and is the default.

**Solver tolerance.** Clarabel's default feasibility tolerance (`1e-8`) does not always converge on
the larger conic partitions, and a certificate that cannot be found is now reported as
`JSRapprox = Inf` rather than as a bound. Comparisons across template complexities should set
`tol_feas`/`tol_gap_*` to `1e-6`; with that setting the bound is monotone in the partition order,
as it must be for a nested family.

**Comparing graphs fairly.** When a template carries a free parameter — the orientation `Gmats`, for
instance — every graph order must be given the same search over it, and the best-so-far should be
reported against the search budget. Under-searching the single-node case manufactures separations
that vanish once it is searched properly; this happened twice during development.

**Comparing across certificates.** Two certificates induce different working sets, since the slice
family is built from the certificate's own sublevel sets. Cell counts and volumes taken over the
whole quotient are therefore not directly comparable; restrict both to a common reference set, or
report the covered volume alongside. How badly this bites: on the predecessor's own example, a
single-node conic-order-1 certificate produces a quotient of **one cell**, because its terminal set
absorbs almost the whole domain. A cell count against that is meaningless, which is why
`gol_lazar_belta_example.jl` reproduces their certificate rather than benchmarking against it.

**Figures.** Scripts do not roll their own plots: `plot_quotient` draws the structural figure (the
cells, or the slice family, over the observation regions) and `plot_synthesis_result` the synthesis
figure (winning cells over losing ones, plus the regions and the closed-loop trajectory). One
colour scheme holds across the folder — **green** is what the controller wins, **red** is what it
does not — so that two figures placed side by side in the paper mean the same thing. Three
conventions follow from it: the observation regions are drawn as outlines in colours disjoint from
green and red, so that a region neither reads as a verdict nor tints the cells beneath it; the
working set is not outlined, since the cells already tile it; and slices are drawn without contours,
because a slice is a union of polytopes and stroking it shows the internal cuts of that union.

**3D figures.** The lifted views use the `DionysosMakieExt` recipes
(`Dionysos.plot_lifted_bisimulation!` / `Dionysos.plot_lifted_trajectory!`); load a Makie backend
(`using GLMakie` for an interactive window, `using CairoMakie` for a static figure) to enable them.
