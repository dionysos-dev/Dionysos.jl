# Overview

This page is the map of the toolbox: how the code is organised, which specifications and solvers
exist, and the one interface that ties them together. For the idea behind the method, read
[Abstraction-based control](@ref) first; to see it used, read
[Getting started](../generated/getting_started.md).

## Code structure

The library lives in [`src`](https://github.com/dionysos-dev/Dionysos.jl/tree/master/src), split into
six modules loaded in dependency order:

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
[`ProblemType`](@ref Dionysos.Problem.ProblemType) and hands both to a solver, so it owns no control
semantics of its own.

## Systems

Dionysos uses the system types of:

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
  — extends the above to hybrid systems: several modes, guarded transitions, reset maps.

## Problems

The specifications currently supported. Each is a system ``\mathcal{S}`` plus an initial set
``\mathcal{I}``, a horizon ``T``, and whatever sets the property itself needs: a target
``\mathcal{T}``, a safe set ``\mathcal{S}_{\text{safe}}``, a state cost ``\mathcal{V}`` and a
transition cost ``\mathcal{C}``. See the [`Problem`](@ref Problem) reference for the full
definitions.

| Specification | Data | Property | Description |
| :--- | :--- | :--- | :--- |
| [Reach-avoid optimal control](@ref Dionysos.Problem.OptimalControlProblem) | ``(\mathcal{S},\mathcal{I},\mathcal{T},\mathcal{V},\mathcal{C},T,\mathcal{S}_{\text{safe}})`` | ``\Box\,\mathcal{S}_{\text{safe}} \wedge \Diamond\,\mathcal{T}`` | Reach the target ``\mathcal{T}`` within horizon ``T`` without leaving ``\mathcal{S}_{\text{safe}}``, minimizing the accumulated cost. |
| [Safety](@ref Dionysos.Problem.SafetyProblem) | ``(\mathcal{S},\mathcal{I},\mathcal{S}_{\text{safe}},T)`` | ``\Box\,\mathcal{S}_{\text{safe}}`` | Never leave ``\mathcal{S}_{\text{safe}}``, for the whole horizon ``T``. |
| [Reach-and-stay](@ref Dionysos.Problem.ReachAndStayProblem) | ``(\mathcal{S},\mathcal{I},\mathcal{T},\mathcal{S}_{\text{safe}},T)`` | ``\Box\,\mathcal{S}_{\text{safe}} \wedge \Diamond\Box\,\mathcal{T}`` | Reach the target ``\mathcal{T}`` and remain in it from then on, without leaving ``\mathcal{S}_{\text{safe}}`` on the way. |
| [Co-safe LTL](@ref Dionysos.Problem.CoSafeLTLProblem) | ``(\mathcal{S},\mathcal{I},\varphi,L)`` | ``\varphi`` | Satisfy a co-safe LTL formula ``\varphi`` over the regions named by the labelling ``L``, i.e. reach an accepting condition in finite time. |

The safe set of a reach-avoid problem is optional (`nothing` means the whole state space), and it is
**not** the same as carving the unsafe region out of the state set: a region removed from the state
space is never abstracted, so the synthesis cannot reason about it, whereas a safe set keeps it
representable and lets the controller actively avoid it. The front-end writes the first as `∉` and
the second as `Always`.

Two abstraction-only problems parametrize the construction of a reusable abstraction without a
control objective: [`AlternatingSimulationProblem`](@ref Dionysos.Problem.AlternatingSimulationProblem)
and [`BisimulationQuotientProblem`](@ref Dionysos.Problem.BisimulationQuotientProblem).

## Solvers

**Abstraction-based solvers:**

| Solver | Discretization | Partition/Cover | Cell shape | Local controller | Reference |
| :--- | :--- | :--- | :--- | :--- | :--- |
| [Uniform grid abstraction](@ref Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer) | Full | Partition | Hyperrectangle | Piecewise constant | [rungger2016scots](@cite) |
| [Uniform ellipsoid abstraction](@ref Dionysos.Optim.Abstraction.UniformEllipsoidAbstraction.Optimizer) | Full | Cover | Ellipsoid | Piecewise affine | [egidio2022state](@cite) |
| [Lazy ellipsoids abstraction](@ref Dionysos.Optim.Abstraction.LazyEllipsoidsAbstraction.Optimizer) | Partial | Cover | Ellipsoid | Piecewise affine | [calbert2024smart](@cite) |
| [Hybrid system abstraction](@ref Dionysos.Optim.Abstraction.HybridSystemAbstraction.Optimizer) | Full | Partition | Hyperrectangle | Piecewise constant | — |
| [PCLF bisimulation quotient](@ref "PCLF bisimulation quotient") | — | Partition | [Semi-linear set](@ref Dionysos.Utils.SemiLinearSet) | — | — |

**Non abstraction-based solvers:**

| Solver | Description | Reference |
| :--- | :--- | :--- |
| [Bemporad–Morari](@ref Dionysos.Optim.BemporadMorari.Optimizer) | Optimal control of hybrid systems via a mixed-integer quadratic program (MIQP). | [bemporad1999control](@cite) |
| [Branch and bound](@ref Dionysos.Optim.BranchAndBound.Optimizer) | Optimal control of hybrid systems combining branch and bound with Q-functions refined by Lagrangian duality. | [legat2021abstraction](@cite) |

Controller synthesis on an already-built automaton is available directly through the
[discrete-system solvers](@ref "Discrete-system solvers").

## The solver interface

This is the architectural keystone. Every solver is a submodule exposing an
`Optimizer <: MOI.AbstractOptimizer`, configured with raw attributes and run through
`MOI.optimize!`:

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

### Why it matters: swapping the solver

Because the specification is a solver-independent object, the *same* control problem can be handed
to a different algorithm without rewriting it. Only the optimizer and its parameters change:

```julia
# One problem …
concrete_problem = PR.OptimalControlProblem(system, initial_set, target_set, nothing, nothing)

# … solved by a grid abstraction …
grid = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
MOI.set(grid, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(grid, MOI.RawOptimizerAttribute("state_grid"), MP.GridFree(x0, hx))
MOI.optimize!(grid)

# … or by an ellipsoidal one, on the very same object.
ellips = MOI.instantiate(AB.UniformEllipsoidAbstraction.Optimizer)
MOI.set(ellips, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(ellips, MOI.RawOptimizerAttribute("sdp_solver"), sdp_solver)
MOI.optimize!(ellips)
```

Both then answer `MOI.TerminationStatus` and hand back a `concrete_controller`, so results are
directly comparable. This is what makes benchmarking algorithms — rather than re-implementing
problems — the normal way to work in Dionysos. The [Solver families](../generated/dcdc_converter.md)
examples each drive one optimizer this way.

## The two entry styles

**The JuMP front-end** is the canonical one: it writes the system and the specification as
constraints and picks the solver for you. It reaches the uniform grid abstraction and the hybrid
abstraction. Start from [Getting started](../generated/getting_started.md); the full vocabulary is in
the [`Wrapper`](@ref Wrapper) reference.

**Direct MathOptInterface** builds the `ProblemType` by hand and configures one family optimizer.
It is the way to reach the solvers whose inputs the front-end cannot express — PWA systems,
ellipsoid templates, observation regions — and the way new solvers are exercised first.

Both report solution status identically, because the status is answered by the solvers themselves.

## Plotting

Every structure worth looking at — trajectories, discretizations, specifications, obstacles,
abstractions — carries a [`@recipe`](https://github.com/JuliaPlots/RecipesBase.jl), so results are
displayed with the single [`plot`](https://docs.juliaplots.org/latest/) function of
[`Plots.jl`](https://github.com/JuliaPlots/Plots.jl). Closed-loop runs can additionally be animated
as a multi-panel dashboard with
[`animate_trajectory_dashboard`](@ref Dionysos.animate_trajectory_dashboard); every example ends with
one.
