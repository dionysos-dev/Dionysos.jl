# Path-Complete Abstractions of Constrained Switched Systems

Code and data for the experiments of the journal version of the HSCC 2027 paper. Work in progress:
the Code Ocean entry points (`run`, `environment/`) come once the Dionysos commit they pin is merged.

```
julia --project=. -e 'using Pkg; Pkg.develop(path="../.."); Pkg.instantiate()'
julia --project=. run_all.jl
```

`run_all.jl` runs every script in its own process. `CACHE=0` rebuilds the quotients instead of
reading them from `cache/`. `TIMING=0` skips the build times, and `SAMPLES=n` reports the median of
`n` timed rebuilds instead of three, each after an untimed warm-up and once every quotient of the
experiment exists, as in the conference experiments. Each script writes its numbers to
`results/numbers_<script>.tex`, which the paper includes, so that the text cannot disagree with the
code, and its figure to `results/fig_<script>.png`.

| script | paper | what it shows |
| :--- | :--- | :--- |
| `experiments/E2_graph.jl` | Section 10.2 | the primal and the dual graph with almost identical pieces; runs the conference experiment `research/HSCC2027/experiment2_primal_and_dual.jl` unchanged, with `SAMPLES=3` |
| `experiments/E3_constrained_plant.jl` | Section 10.3 | the running example when the plant never applies mode 1 twice in a row: the exact answer on the constraint graph against the arbitrary-switching quotient, whose ∀-set is sound and whose ∃-set is not |
| `experiments/E4_unstable_plant.jl` | Section 10.4 | the running example with mode 1 scaled by 1.2: no common Lyapunov function exists, and a controller that never plays mode 1 twice in a row obtains a certificate, a quotient and a sound synthesis answer |
| `experiments/cells_per_slice.jl` | Section 10.5 | the cells of each slice of the quotients of E3 and E4, read from `cache/`: geometric growth, slower under the constraint |
| `experiments/E5_dwell_time.jl` | Section 10.6 | two contracting modes whose alternation diverges: the certified rate against the minimum dwell time, and the cycle lower bound that makes the threshold a property of the system |
| `experiments/E6_coset_family.jl` | Sections 9, 10.7 | the rotation family in closed form: the quotients without regions match Proposition 9.3, and with regions the two constructions are compared in dimensions 2 and 4 |
| `experiments/fig_coset_family.jl` | Section 10.7 | the planar case of the family: square and diamond pieces, and the octagon they induce |

`common.jl` holds what the experiments share: the running example of Gol, Ding, Lazar and Belta,
the certificate search with the template of the conference experiments (`research/HSCC2027`), the
construction and timing of a quotient, the solution of the specification, the projection of
certified cells onto the plant at the points of a probe grid, and the figure helpers.
