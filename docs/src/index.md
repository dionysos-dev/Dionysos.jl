```@raw html
<img class="display-light-only" src="assets/logo.png" height="240" alt="Dionysos Logo">
<img class="display-dark-only" src="assets/logo-dark.png" height="240" alt="Dionysos Logo">
```

# Dionysos

**Describe your system, state what it must do, and get back a controller — together with a proof
that it works.**

[Dionysos](https://github.com/dionysos-dev/Dionysos.jl) is a [Julia](https://julialang.org/)
framework for correct-by-construction controller synthesis through symbolic (abstraction-based)
control. It is the software of the ERC project
[Learning to Control](https://perso.uclouvain.be/raphael.jungers/content/erc-consolidator-grant) (L2C).

Designing a controller for a complex system traditionally means a team of experts hand-crafting one
over months. Dionysos pursues **Control as a Service**: an automated pipeline usable without a
dedicated control department.

> **describe the system → state the specification → pick a solver → obtain a controller and a
> certificate**

## A complete example

A pendulum, hanging down, that has to be swung up and caught at the top. This is the whole model:

```julia
using Dionysos, JuMP, StaticArrays

model = Model(Dionysos.Optimizer)
@variable(model, -π <= x1 <= 2π)                       # state bounds are the constraint set
@variable(model, -10 <= x2 <= 10)
@variable(model, -3 <= u <= 3)

@constraint(model, ∂(x1) == x2)                        # the dynamics
@constraint(model, ∂(x2) == -9.81 * sin(x1) + u)

@constraint(model, start(x1) in MOI.Interval(-0.09, 0.09))   # where it begins
@constraint(model, final(x1) in MOI.Interval(π - 0.09, π + 0.09))   # where it must end up

set_attribute(model, "state_grid", Dionysos.Mapping.GridFree(SVector(0.0, 0.0), SVector(0.05, 0.05)))
set_attribute(model, "input_grid", Dionysos.Mapping.GridFree(SVector(0.0), SVector(0.3)))
set_attribute(model, "time_step", 0.1)

optimize!(model)
trajectory = Dionysos.simulate(model, SVector(0.0, 0.0))
```

Nothing in it says *how* to swing the pendulum up. That behaviour — rock back and forth to pump in
energy, then catch it — comes out of the synthesis. [**Getting started**](generated/getting_started.md)
walks through the same model line by line and animates the result.

## How it works

The technique is **symbolic control**: the continuous system is *abstracted* into a finite-state
automaton by discretising its variables, a controller is synthesised on that finite object with
graph algorithms, and it is *concretised* back to the original system with a formal guarantee. The
price is the curse of dimensionality, and fighting it — through **smart, lazy abstractions** that
compute only the part of the abstraction actually needed — is a core research direction of the
toolbox. See [Abstraction-based control](@ref).

Dionysos is an **ecosystem, not a single algorithm**. Every solver is a
[MathOptInterface](https://jump.dev/MathOptInterface.jl) optimizer driven through
[JuMP](https://jump.dev/JuMP.jl), so a control task can be re-solved, compared and benchmarked by
*swapping the solver* rather than rewriting the model.

## What it can do

**Specifications** — reach-avoid optimal control, safety, reach-and-stay, co-safe LTL, plus the
abstraction-only problems (alternating simulation, bisimulation quotient).
See the [`Problem`](@ref Problem) reference.

**Solvers** — uniform grid abstraction (SCOTS-style), uniform and lazy ellipsoidal abstractions,
hybrid-system abstraction, a PCLF bisimulation quotient, discrete-automaton synthesis, and
optimisation-based solvers (Bemporad–Morari, branch and bound).
See [Overview](@ref) and the [`Optim`](@ref Optim) reference.

**Models** — continuous- and discrete-time systems, hybrid systems with modes, guards, reset maps
and clocks, obstacles, arbitrary bounded `LazySet` targets, and dynamics supplied as a Julia
function. See the [`Wrapper`](@ref Wrapper) reference.

## Installation

```julia
julia> import Pkg; Pkg.add("Dionysos")
```

## How to cite

If you use Dionysos in your research, please cite the paper describing the toolbox and the software
release you actually ran. The BibTeX entry is in
[How to cite](https://github.com/dionysos-dev/Dionysos.jl#how-to-cite).

## Need help?

Open an [issue](https://github.com/dionysos-dev/Dionysos.jl/issues) on GitHub.

## ERC sponsor

This project has received funding from the European Research Council (ERC) under the European
Union's Horizon 2020 research and innovation programme under grant agreement No 864017 - L2C.

```@raw html
<img class="display-light-only" src="assets/logo_erc_white.jpg" alt="ERC logo"/>
<img class="display-dark-only" src="assets/logo_erc_black.jpg" alt="ERC logo"/>
```
