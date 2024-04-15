# Overview

Dionysos aims to design a controller for a system $\mathcal{S}$ so that the closed-loop system satisfies the specification $\Sigma$ where:
* the system $\mathcal{S}$ is specified by [`MathematicalSystems`](https://juliareach.github.io/MathematicalSystems.jl/latest/lib/types/#MathematicalSystems.AbstractSystem) or [`HybridSystems`](https://blegat.github.io/HybridSystems.jl/stable/lib/types/#HybridSystems.AbstractHybridSystem) objects;
* the specification $\Sigma$ is specified by [`ProblemType`](https://dionysos-dev.github.io/Dionysos.jl/dev/reference/Problem/#Dionysos.Problem.ProblemType) objects;
* the solver $\mathcal{O}$ implementents the abstract type [`AbstractOptimizer`](https://jump.dev/MathOptInterface.jl/stable/reference/models/#MathOptInterface.AbstractOptimizer) of [`MathOptInterface`](https://github.com/jump-dev/MathOptInterface.jl).

So a control problem $(\mathcal{S},\Sigma)$ can be solved by $\mathcal{O}$ via the [`JuMP`](https://github.com/jump-dev/JuMP.jl) interface, with Dionysos inheriting JuMP's powerful and practical optimization framework.

# Overview of the code structure
Description of the core of the Dionysos.jl package, the [`src`](https://github.com/dionysos-dev/Dionysos.jl/tree/master/src) folder:

| Subfolder        | Description |
| :--------------- | :---------- |
| [`utils`](https://github.com/dionysos-dev/Dionysos.jl/tree/master/src/utils) | Contains useful functions, data structures, classic search algorithms, file management, ... |
| [`domain`](https://github.com/dionysos-dev/Dionysos.jl/tree/master/src/domain) | Contains structures defining the domain of a system |
| [`system`](https://github.com/dionysos-dev/Dionysos.jl/tree/master/src/system) | Contains a description of specific systems |
| [`problem`](https://github.com/dionysos-dev/Dionysos.jl/tree/master/src/problem) | Contains control problems that can be solved by Dionysos solvers |
| [`symbolic`](https://github.com/dionysos-dev/Dionysos.jl/tree/master/src/symbolic) | Contains the data structures needed to encode the abstractions |
| [`optim`](https://github.com/dionysos-dev/Dionysos.jl/tree/master/src/optim) | Contains the solvers |


# Systems

The system types supported in Dionysos.jl are:
* [`MathematicalSystems`](https://juliareach.github.io/MathematicalSystems.jl/latest/lib/types/#MathematicalSystems.AbstractSystem), which proposes generic and flexible system definitions (e.g.     discrete-time/continuous-time, constrained, noisy systems), such that, for example, the system [`MathematicalSystems.NoisyConstrainedAffineControlDiscreteSystem`](https://juliareach.github.io/MathematicalSystems.jl/latest/lib/types/#MathematicalSystems.NoisyConstrainedAffineControlDiscreteSystem) of the form $$x(k+1) = A x(k) + B u(k) + c + D w(k), \ x(k)\in\mathcal{X}, \ u(k)\in\mathcal{U},\ w(k)\in\mathcal{W}\ \forall k$$
where $\mathcal{X}$ is the state constraint, $\mathcal{U}$ is the input constraint and $\mathcal{W}$ is the noise constraint.
* [`HybridSystems`](https://blegat.github.io/HybridSystems.jl/stable/lib/types/#HybridSystems.AbstractHybridSystem), which extends the class of systems of [`MathematicalSystems`](https://juliareach.github.io/MathematicalSystems.jl/latest/lib/types/#MathematicalSystems.AbstractSystem) to hybrid systems.


# Problems

The problem types supported in Dionysos.jl are:

| Type         | Description |
| :--------------- | :---------- |
| [`Reach-avoid optimal control problem`](https://dionysos-dev.github.io/Dionysos.jl/dev/reference/Problem/#Dionysos.Problem.OptimalControlProblem) | $(\mathcal{S},\mathcal{I},\mathcal{T},\mathcal{V},\mathcal{C}, T)$, where $\mathcal{S}$ is the system, $\mathcal{I}$ is the initial set, $\mathcal{T}$ is the target set, $\mathcal{V}:\mathcal{X}\rightarrow \mathbb{R}$ is the cost state function, $\mathcal{C}:\mathcal{X}\times \mathcal{U}\rightarrow \mathbb{R}$ is the transition cost function, $T$ is the time limit to satisfy the specification. |
| [`Safety control problem`](https://dionysos-dev.github.io/Dionysos.jl/dev/reference/Problem/#Dionysos.Problem.SafetyProblem) | $(\mathcal{S},\mathcal{I},\mathcal{S}, T)$, where $\mathcal{S}$ is the system, $\mathcal{I}$ is the initial set, $\mathcal{S}$ is the safe set,  $T$ is the time during which safety must be ensured. |

Extensions for linear temporal logic (LTL) specifications are currently being implemented.

# Solvers

The following tables summarize the different solvers. 


**Abstraction-based solver types implemented in Dionysos.jl:**

| Type          | Full vs partial discretization | Partition vs Cover | Shape | Local controller | Abstraction | System | Reference | 
| :--------------- | :---------- | :---------- | :---------- | :---------- | :---------- | :---------- | :---------- |
| [`Uniform grid abstraction`](https://dionysos-dev.github.io/Dionysos.jl/dev/reference/Optim/#Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer) | Full | Partition | Hyperrectangle | Piece-wise constant | Non-determinisitic | Continuous-time | [`SCOTS: A Tool for the Synthesis of Symbolic Controllers`](https://dl.acm.org/doi/abs/10.1145/2883817.2883834) |
| [`Lazy abstraction`](https://dionysos-dev.github.io/Dionysos.jl/dev/reference/Optim/#Dionysos.Optim.Abstraction.LazyAbstraction.Optimizer) | Partial | Partition | Hyperrectangle | Piece-wise constant | Non-determinisitic | Continuous-time | [`Alternating simulation on hierarchical abstractions`](https://ieeexplore.ieee.org/abstract/document/9683448/?casa_token=AXyECYU9FdwAAAAA:ERfvlbkORIbfGLbDd42d2K5K9SLVOjl-8kRs9pfp7lMa4QZEv0W4VgzlP8FshBlxXQF4ZQrDuzk) |
| [`Hierarchical abstraction`](https://dionysos-dev.github.io/Dionysos.jl/dev/reference/Optim/#Dionysos.Optim.Abstraction.HierarchicalAbstraction.Optimizer) | Partial | Partition | Hyperrectangle | Piece-wise constant | Non-determinisitic | Continuous-time | [`Alternating simulation on hierarchical abstractions`](https://ieeexplore.ieee.org/abstract/document/9683448/?casa_token=AXyECYU9FdwAAAAA:ERfvlbkORIbfGLbDd42d2K5K9SLVOjl-8kRs9pfp7lMa4QZEv0W4VgzlP8FshBlxXQF4ZQrDuzk) |
| [`Ellipsoid abstraction`](https://dionysos-dev.github.io/Dionysos.jl/dev/reference/Optim/#Dionysos.Optim.Abstraction.EllipsoidsAbstraction.Optimizer) | Full | Cover | Ellipsoid | Piece-wise affine | Determinisitic | Discrete-time affine | [`State-feedback Abstractions for Optimal Control of Piecewise-affine Systems`](https://arxiv.org/abs/2204.00315) |
| [`Lazy ellipsoid abstraction`](https://dionysos-dev.github.io/Dionysos.jl/dev/reference/Optim/#Dionysos.Optim.Abstraction.LazyEllipsoidsAbstraction.Optimizer) | Partial | Cover | Ellipsoid | Piece-wise affine | Determinisitic | Discrete-time non-linear | Not yet published |

**Non abstraction-based solver types implemented in Dionysos.jl:**

| Type          | Description | Reference |
| :--------------- | :---------- | :---------- |
| [`Bemporad Morari`](https://dionysos-dev.github.io/Dionysos.jl/dev/reference/Optim/#Dionysos.Optim.BemporadMorari.Optimizer) | Optimal control of hybrid systems via a predictive control scheme using mixed integer quadratic programming (MIQP) online optimization procedures. | [`Control of systems integrating logic, dynamics, and constraints`](https://www.sciencedirect.com/science/article/abs/pii/S0005109898001782)
| [`BranchAndBound`](https://dionysos-dev.github.io/Dionysos.jl/dev/reference/Optim/#Dionysos.Optim.BranchAndBound.Optimizer) | Optimal control of hybrid systems via a predictive control scheme combining a branch and bound algorithm that can refine Q-functions using Lagrangian duality. | [`Abstraction-based branch and bound approach to Q-learning for hybrid optimal control`](https://proceedings.mlr.press/v144/legat21a.html)


**Solver interface**

Each solver is defined by a module which must implement the abstract type [`AbstractOptimizer`](https://jump.dev/MathOptInterface.jl/stable/reference/models/#MathOptInterface.AbstractOptimizer) and the [`Optimize!`](https://jump.dev/MathOptInterface.jl/stable/reference/models/#MathOptInterface.optimize!) function.
For example, for the UniformGridAbstraction solver, this structure and function are defined as follows

```julia
using JuMP

mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    concrete_problem::Union{Nothing, PR.ProblemType}
    abstract_problem::Union{Nothing, PR.OptimalControlProblem, PR.SafetyProblem}
    abstract_system::Union{Nothing, SY.SymbolicModelList}
    abstract_controller::Union{Nothing, UT.SortedTupleSet{2, NTuple{2, Int}}}
    concrete_controller::Any
    state_grid::Union{Nothing, DO.Grid}
    input_grid::Union{Nothing, DO.Grid}
    function Optimizer{T}() where {T}
        return new{T}(nothing, nothing, nothing, nothing, nothing, nothing, nothing)
    end
end
```

and

```julia
function MOI.optimize!(optimizer::Optimizer)
    # Build the abstraction
    abstract_system = build_abstraction(
        optimizer.concrete_problem.system,
        optimizer.state_grid,
        optimizer.input_grid,
    )
    optimizer.abstract_system = abstract_system
    # Build the abstract problem
    abstract_problem = build_abstract_problem(optimizer.concrete_problem, abstract_system)
    optimizer.abstract_problem = abstract_problem
    # Solve the abstract problem
    abstract_controller = solve_abstract_problem(abstract_problem)
    optimizer.abstract_controller = abstract_controller
    # Solve the concrete problem
    optimizer.concrete_controller =
        solve_concrete_problem(abstract_system, abstract_controller)
    return
end
```

# Running an example
In this section, we outline how to define and solve a control problem with Dionsysos.
For an executable version of this example, see [`Solvers: Path planning problem`](https://dionysos-dev.github.io/Dionysos.jl/dev/generated/Path%20planning/#Example:-Path-planning-problem) in the documentation.

Define a control problem, i.e., the system and the specification of the desired closed loop behaviour.
To do this, you can define new ones yourself or directly load an existing benchmark, for example

```julia
concrete_problem = PathPlanning.problem(; simple = true, approx_mode = PathPlanning.GROWTH);
concrete_system = concrete_problem.system;
```

Choose the solver you wish to use
```julia
using JuMP
optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
```

Define the solver's meta-parameters
```julia
x0 = SVector(0.0, 0.0, 0.0);
hx = SVector(0.2, 0.2, 0.2);
state_grid = DO.GridFree(x0, hx);
u0 = SVector(0.0, 0.0);
hu = SVector(0.3, 0.3);
input_grid = DO.GridFree(u0, hu);
```

Set the solver's meta-parameters
```julia
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
```

Solve the control problem
```julia
MOI.optimize!(optimizer)
```

Get the results
```julia
abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
```

In Dionysos, all the structures that could be relevant to plot (such as trajectories, state-space discretization, specifications, obstacles, etc.) 
have an associated [`@recipe`](https://github.com/JuliaPlots/RecipesBase.jl) function, which makes it very easy to plot all the results using the single common [`plot`](https://docs.juliaplots.org/latest/generated/unitfulext_plots/) function of [`Plots.jl`](https://github.com/JuliaPlots/Plots.jl).
For example
```julia
using Plots

plot!(concrete_system.X; color = :yellow, opacity = 0.5);
plot!(abstract_system.Xdom; color = :blue, opacity = 0.5);
plot!(concrete_problem.initial_set; color = :green, opacity = 0.2);
plot!(concrete_problem.target_set; dims = [1, 2], color = :red, opacity = 0.2);
plot!(control_trajectory; ms = 0.5)
```



