# N14 — The controller is cheaper, not just the build

Cheap, uncontroversial, and it measures a resource nobody else in this folder measures. For a control
audience it may matter more than build time, because the build happens once offline and the
controller runs forever.

## What we want to show

> The lifted construction produces a controller whose cells are an order of magnitude simpler, so
> every online membership test is an order of magnitude cheaper and the serialized controller is
> correspondingly smaller. This advantage is independent of build time and holds even on benchmarks
> where the build is a tie.

## Why it holds regardless of the timing question

A static controller must, at each step, locate the current state in the partition: given $x$, find the
cell containing it. For a cell represented as a union of $p$ polytopes with $F$ facets in total, a
membership test is $O(F)$ dot products. Measured per-cell facet counts:

| benchmark | route 2 (PCLF) | route 1 (induced CLF) | ratio |
| :--- | --: | --: | --: |
| E2, mean facets/cell | 10.97 | 98.28 | **9.0×** |
| E2, max facets/cell | 146 | 1 728 | **11.8×** |
| E8, mean facets/cell | 15.44 | 60.42 | **3.9×** |
| E8, max facets/cell | 187 | 1 557 | **8.3×** |

So the worst-case online test is 8–12× cheaper under route 2, and the average 4–9× cheaper. This is
true on E8 — where route 2 *loses* the build by 1.30× — which is exactly what makes it a useful,
independent claim.

## What we will try

1. **Measure it, do not infer it.** Time `is_defined` / point location over a grid of $10^4$ sampled
   states for both arms on E2 and E8. Report mean and worst case.
2. **Measure the serialized size.** Both controllers already round-trip through JLD2
   (`export_optimizer_jld2`); report bytes per controller and bytes per cell.
3. **Add the memory footprint of the quotient itself** — the number of halfspaces stored is
   $\Sigma$ facets, measured at 98 987 vs 98 276 (E2) and 163 854 vs 161 937 (E8), so here the two
   routes tie on storage and differ only on per-query cost. Say so; it sharpens the claim rather than
   weakening it.
4. **State the caveat.** Point location can be accelerated with spatial indexing, which would narrow
   the gap. Report the naive cost and note the caveat rather than claiming an asymptotic advantage.

## Why it is worth a column even if small

Dionysos ships controllers that run online through `control_server/`, so "the controller is 10×
cheaper to evaluate" is a statement about a deliverable, not about a research prototype. It also costs
almost nothing to produce, since both controllers already exist in the caches.

## What would falsify it

Point location being dominated by the *number* of cells rather than their facet count — route 2 has
4–9× more cells, so a naive linear scan over cells would cancel the advantage exactly. **Check which
regime the implementation is in before claiming anything**: if `is_defined` scans all cells, the
correct claim is about the per-cell test only, not the end-to-end query.
