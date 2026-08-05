# Bisimulation quotient (HSCC 2027)

Scripts reproducing the experiments of our **HSCC 2027** submission on path-complete Lyapunov
function (PCLF) based bisimulation quotients for switched systems.

> **Citation:** _add the full reference once the paper is finalized._

## Scripts

- `example1.jl`, `example2.jl`, `example4.jl`, `example_3_1.jl` — worked examples of the quotient
  construction and closed-loop control.
- `example_PCLF_vs_CLF.jl` — comparison of a path-complete versus a common (single) Lyapunov
  certificate.
- `example_state_space.jl` — state-space visualization of the quotient.

The 3D "lifted quotient" figures use the `DionysosMakieExt` recipes
(`Dionysos.plot_lifted_bisimulation!` / `Dionysos.plot_lifted_trajectory!`); load a Makie backend
(e.g. `using GLMakie`) to enable them.
