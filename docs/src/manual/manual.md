# Overview

Dionysos delivers **Control as a Service**: you describe a system and a specification, and it returns a
certified controller — no bespoke, expert-crafted design required. Concretely, Dionysos designs a
controller for a system ``\mathcal{S}`` so that the closed loop satisfies a specification ``\Sigma``,
where:

- the system ``\mathcal{S}`` is a
  [`MathematicalSystems`](https://juliareach.github.io/MathematicalSystems.jl/latest/lib/types/#MathematicalSystems.AbstractSystem)
  or [`HybridSystems`](https://blegat.github.io/HybridSystems.jl/stable/lib/types/#HybridSystems.AbstractHybridSystem)
  object;
- the specification ``\Sigma`` is a [`ProblemType`](@ref Dionysos.Problem.ProblemType) object;
- the solver ``\mathcal{O}`` implements the
  [`AbstractOptimizer`](https://jump.dev/MathOptInterface.jl/stable/reference/models/#MathOptInterface.AbstractOptimizer)
  interface of [`MathOptInterface`](https://github.com/jump-dev/MathOptInterface.jl).

A control problem ``(\mathcal{S}, \Sigma)`` is therefore solved by ``\mathcal{O}`` through the
[`JuMP`](https://github.com/jump-dev/JuMP.jl) interface, so Dionysos inherits JuMP's optimization
framework. For the conceptual background, see [Abstraction-based control](@ref).

## Code structure

The core of the package lives in the [`src`](https://github.com/dionysos-dev/Dionysos.jl/tree/master/src)
folder, split into six library modules loaded in dependency order:

| Module | Description |
| :--- | :--- |
| [`Utils`](@ref Utils)       | Foundational helpers on top of LazySets: sets, cost functions, data structures, search, scalar optimization. |
| [`System`](@ref System)     | Concrete dynamical systems, their approximations, controllers, trajectories, and the simulation engine. |
| [`Problem`](@ref Problem)    | Solver-independent control-task specifications. |
| [`Mapping`](@ref Mapping)    | Concrete ↔ abstract discretization: grids, cells, mappings. |
| [`Symbolic`](@ref Symbolic)  | The finite automaton abstraction built from a system and a mapping. |
| [`Optim`](@ref Optim)       | The solver catalog. |

On top of them sits [`Wrapper`](@ref Wrapper), the JuMP front-end reached through
`Model(Dionysos.Optimizer)`. It compiles a JuMP model into a system plus a
[`ProblemType`](@ref Dionysos.Problem.ProblemType) and hands both to a solver, so it owns no
control semantics of its own.

## Systems

Dionysos supports the system types provided by:

- [`MathematicalSystems`](https://juliareach.github.io/MathematicalSystems.jl/latest/lib/types/#MathematicalSystems.AbstractSystem)
  — generic, flexible system definitions (discrete-/continuous-time, constrained, noisy). For
  instance a
  [`NoisyConstrainedAffineControlDiscreteSystem`](https://juliareach.github.io/MathematicalSystems.jl/latest/lib/types/#MathematicalSystems.NoisyConstrainedAffineControlDiscreteSystem)
  of the form
  ```math
  x(k+1) = A x(k) + B u(k) + c + D w(k), \quad x(k)\in\mathcal{X},\ u(k)\in\mathcal{U},\ w(k)\in\mathcal{W},
  ```
  where ``\mathcal{X}``, ``\mathcal{U}`` and ``\mathcal{W}`` are the state, input and noise
  constraints.
- [`HybridSystems`](https://blegat.github.io/HybridSystems.jl/stable/lib/types/#HybridSystems.AbstractHybridSystem)
  — extends the above to hybrid systems.

## Problems

The specifications currently supported are (see the [`Problem`](@ref Problem) reference for the full
definitions):

| Specification | Description |
| :--- | :--- |
| [Reach-avoid optimal control](@ref Dionysos.Problem.OptimalControlProblem) | ``(\mathcal{S},\mathcal{I},\mathcal{T},\mathcal{V},\mathcal{C},T)``: reach a target set ``\mathcal{T}`` from an initial set ``\mathcal{I}`` within horizon ``T`` while avoiding obstacles, minimizing a state cost ``\mathcal{V}`` and a transition cost ``\mathcal{C}``. |
| [Safety](@ref Dionysos.Problem.SafetyProblem) | ``(\mathcal{S},\mathcal{I},\mathcal{S}_{\text{safe}},T)``: remain inside a safe set for the whole horizon ``T``. |
| [Reach-and-stay](@ref Dionysos.Problem.ReachAndStayProblem) | Eventually reach a target set and then remain in it. |
| [Co-safe LTL](@ref Dionysos.Problem.CoSafeLTLProblem) | Satisfy a co-safe LTL formula, i.e. reach an accepting condition in finite time. |

Two abstraction-only problems parametrize the construction of a reusable abstraction without a control
objective: [`AlternatingSimulationProblem`](@ref Dionysos.Problem.AlternatingSimulationProblem) and
[`BisimulationQuotientProblem`](@ref Dionysos.Problem.BisimulationQuotientProblem).

## Solvers

**Abstraction-based solvers:**

| Solver | Discretization | Partition/Cover | Cell shape | Local controller | Reference |
| :--- | :--- | :--- | :--- | :--- | :--- |
| [Uniform grid abstraction](@ref Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer) | Full | Partition | Hyperrectangle | Piecewise constant | [SCOTS](https://dl.acm.org/doi/abs/10.1145/2883817.2883834) |
| [Uniform ellipsoid abstraction](@ref Dionysos.Optim.Abstraction.UniformEllipsoidAbstraction.Optimizer) | Full | Cover | Ellipsoid | Piecewise affine | [State-feedback abstractions](https://arxiv.org/abs/2204.00315) |
| [Lazy ellipsoids abstraction](@ref Dionysos.Optim.Abstraction.LazyEllipsoidsAbstraction.Optimizer) | Partial | Cover | Ellipsoid | Piecewise affine | — |
| [Hybrid system abstraction](@ref Dionysos.Optim.Abstraction.HybridSystemAbstraction.Optimizer) | Full | Partition | Hyperrectangle | Piecewise constant | — |
| [PCLF bisimulation quotient](@ref "PCLF bisimulation quotient") | — | Partition | Polyhedral | — | [Branch and bound Q-learning](https://proceedings.mlr.press/v144/legat21a.html) |

**Non abstraction-based solvers:**

| Solver | Description | Reference |
| :--- | :--- | :--- |
| [Bemporad–Morari](@ref Dionysos.Optim.BemporadMorari.Optimizer) | Optimal control of hybrid systems via a mixed-integer quadratic program (MIQP). | [Bemporad & Morari (1999)](https://www.sciencedirect.com/science/article/abs/pii/S0005109898001782) |
| [Branch and bound](@ref Dionysos.Optim.BranchAndBound.Optimizer) | Optimal control of hybrid systems combining branch and bound with Q-functions refined by Lagrangian duality. | [Legat et al. (2021)](https://proceedings.mlr.press/v144/legat21a.html) |

Controller synthesis on an already-built automaton is available directly through the
[discrete-system solvers](@ref "Discrete-system solvers").

### The solver interface

Every solver is a submodule exposing an `Optimizer <: MOI.AbstractOptimizer`. It is configured with
raw attributes and run through `MOI.optimize!`:

1. instantiate — `optimizer = MOI.instantiate(SomeFamily.Optimizer)`;
2. configure — `MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), problem)` and the
   solver-specific parameters (`"state_grid"`, `"input_grid"`, `"time_step"`, …);
3. run — `MOI.optimize!(optimizer)`;
4. query — `MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))`.

For abstraction-based solvers, `optimize!` follows the same conceptual pipeline:

```
concrete problem → abstraction (symbolic model) → abstract problem
                 → abstract controller → concrete controller
```

The abstraction is cached, so switching the specification on the same system (e.g. safety →
reachability) does not recompute it. Solvers compose: a high-level optimizer holds an
`abstraction_solver` and a `control_solver` and forwards attribute `set`/`get` to them.

## Running an example

The canonical entry point is a JuMP model with `Dionysos.Optimizer`, which writes the dynamics with
the `∂` operator and the target/initial constraints with `final`/`start`:

```julia
using Dionysos, JuMP, StaticArrays

model = Model(Dionysos.Optimizer)
@variable(model, x_low[i] <= x[i = 1:3] <= x_upp[i], start = x0[i])
@variable(model, -1 <= u[1:2] <= 1)
@constraint(model, ∂(x[1]) == u[1] * cos(α + x[3]) * sec(α))   # dynamics
@constraint(model, final(x[1]) in MOI.Interval(3.0, 3.6))       # target set
set_attribute(model, "time_step", 0.3)
set_attribute(model, "state_grid", Dionysos.Mapping.GridFree(x0, hx))
set_attribute(model, "input_grid", Dionysos.Mapping.GridFree(u0, hu))
optimize!(model)
concrete_controller = get_attribute(model, "concrete_controller")
```

The same problem can be driven directly through MathOptInterface by instantiating a specific family
optimizer, setting the `concrete_problem`, `state_grid` and `input_grid` attributes, and reading back
`abstract_system`, `abstract_problem`, `abstract_controller`, and `concrete_controller`.

All structures worth plotting (trajectories, discretizations, specifications, obstacles, …) carry a
[`@recipe`](https://github.com/JuliaPlots/RecipesBase.jl), so results are displayed with the single
[`plot`](https://docs.juliaplots.org/latest/) function of [`Plots.jl`](https://github.com/JuliaPlots/Plots.jl):

```julia
using Plots

plot!(concrete_system.X; color = :yellow, opacity = 0.5)
plot!(abstract_system; color = :blue, opacity = 0.5)
plot!(concrete_problem.initial_set; color = :green, opacity = 0.2)
plot!(concrete_problem.target_set; dims = [1, 2], color = :red, opacity = 0.2)
plot!(trajectory; ms = 0.5)
```

For a complete, executable walkthrough see the
[Path planning example](@ref "Example: Path planning problem solved by uniform grid abstraction") and
[Getting Started](@ref "Getting Started").
