# Abstraction-based control

Given a mathematical description of the system dynamics and the specifications describing the desired closed-loop behavior of the system, 
*abstraction-based control* techniques involve synthesizing a *correct-by-construction* controller through a *systematic* three-step procedure. 
First, both the original system and the specifications are transposed into an abstract domain, resulting in an abstract system and corresponding abstract specifications. 
This step generally involves a complete discretization of the state and input spaces, typically with uniform hyperrectangles.
Next, an abstract controller is synthesized to solve this abstract control problem. Finally, the third step, referred to as the *refinement procedure*, involves deducing a controller for the original control problem from the abstract controller. The value of this approach lies in the substitution of the original system (often an infinite system) with a finite system, which enables it to leverage powerful control tools in the domain of symbolic control. This three steps procedure is illustrated on the following figure.

![The three steps of abstraction-based control](./assets/Abstraction-procedure.png)

Although this approach offers a safety-critical framework, it suffers from the curse of dimensionality due to the exponential growth of the number of states with respect to the dimension.
In order to render these techniques practical, it is necessary to construct *smart abstractions*, i.e., they differ from classical techniques in that the partitioning is designed smartly, using optimization-based design techniques, and computed iteratively, unlike the classical approach which uses an a priori defined approach, sub-optimal and subject to the curse of dimensionality.
To this end, we introduce solvers called *lazy solvers* (i.e. postponing heavier numerical operations) that co-design the abstraction and the controller to reduce the computed part of the abstraction.


# Standard form problem

Dionysos aims to design a controller for a system $\mathcal{S}$ so that the closed-loop system satisfies the specification $\Sigma$ where:
* the system $\mathcal{S}$ is specified by 
[`MathematicalSystems`](https://juliareach.github.io/MathematicalSystems.jl/latest/lib/types/#MathematicalSystems.AbstractSystem) or 
[`HybridSystems`](https://blegat.github.io/HybridSystems.jl/stable/lib/types/#HybridSystems.AbstractHybridSystem) objects;
* the specification $\Sigma$ is specified by [`ProblemType`]() objects;
* the solver $\mathcal{O}$ implementents the abstract type [`AbstractOptimizer`](https://jump.dev/MathOptInterface.jl/stable/reference/models/#MathOptInterface.AbstractOptimizer) of [`MathOptInterface`](https://github.com/jump-dev/MathOptInterface.jl).

So a control problem $(\mathcal{S},\Sigma)$ can be addressed by a solver $\mathcal{O}$ via the [`JuMP`](https://github.com/jump-dev/JuMP.jl) interface, with Dionysos inheriting JuMP's powerful and practical optimization framework.

# Overview of the code structure
Description of the core of the Dionysos.jl package, the [`src`](https://github.com/dionysos-dev/Dionysos.jl/tree/master/src) folder:
| Subfolder        | Description |
| :--------------- | :---------- |
| [`utils`](https://github.com/dionysos-dev/Dionysos.jl/tree/master/src/utils) | Contains useful functions, data structures, classic search algorithms, file management, ... |
| [`domain`](https://github.com/dionysos-dev/Dionysos.jl/tree/master/src/domain) | Contains structures defining the domain of a system |
| [`system`](https://github.com/dionysos-dev/Dionysos.jl/tree/master/src/system) | Contains specific systems and controllers |
| [`problem`](https://github.com/dionysos-dev/Dionysos.jl/tree/master/src/problem) | Contains specifications |
| [`symbolic`](https://github.com/dionysos-dev/Dionysos.jl/tree/master/src/symbolic) | Contains the data structures needed to encode the abstractions |
| [`optim`](https://github.com/dionysos-dev/Dionysos.jl/tree/master/src/optim) | Contains the solvers |


# Systems

The system types supported in Dionysos.jl are:
* [`MathematicalSystems`](https://juliareach.github.io/MathematicalSystems.jl/latest/lib/types/#MathematicalSystems.AbstractSystem),
    which proposes generic and flexible system definitions (e.g. discrete-time/continuous-time, constrained, noisy systems), such that, for example, the system
    [`MathematicalSystems.NoisyConstrainedAffineControlDiscreteSystem`](https://juliareach.github.io/MathematicalSystems.jl/latest/lib/types/#MathematicalSystems.NoisyConstrainedAffineControlDiscreteSystem)
    of the form $$x(k+1) = A x(k) + B u(k) + c + D w(k), \ x(k)\in\mathcal{X}, \ u(k)\in\mathcal{U},\ w(k)\in\mathcal{W}\ \forall k$$
    where 
    * $\mathcal{X}$ is the state constraints;
    * $\mathcal{U}$ is the input constraints;
    * $\mathcal{W}$ is the noise constraints.
* [`HybridSystems`](https://blegat.github.io/HybridSystems.jl/stable/lib/types/#HybridSystems.AbstractHybridSystem), which extends the class of systems of [`MathematicalSystems`](https://juliareach.github.io/MathematicalSystems.jl/latest/lib/types/#MathematicalSystems.AbstractSystem) to hybrid systems.



# Specifications

The specification types implemented in Dionysos.jl are:
| Type         | Description |
| :--------------- | :---------- |
| [`Reach-avoid optimal control problem`](https://dionysos-dev.github.io/Dionysos.jl/dev/reference/Problem/#Dionysos.Problem.OptimalControlProblem) | $(\mathcal{S},\mathcal{I},\mathcal{T},\mathcal{V},\mathcal{C}, T)$, where $\mathcal{S}$ is the system, $\mathcal{I}$ is the initial set, $\mathcal{T}$ is the target set, $\mathcal{V}:\mathcal{X}\rightarrow \mathbb{R}$ is the cost state function, $\mathcal{C}:\mathcal{X}\times \mathcal{U}\rightarrow \mathbb{R}$ is the transition cost function, $T$ is the time limit to satisfy the specification. |
| [`Safety control problem`](https://dionysos-dev.github.io/Dionysos.jl/dev/reference/Problem/#Dionysos.Problem.SafetyProblem) | $(\mathcal{S},\mathcal{I},\mathcal{S}, T)$, where $\mathcal{S}$ is the system, $\mathcal{I}$ is the initial set, $\mathcal{S}$ is the safe set,  $T$ is the time during which safety must be ensured. |

Extensions for linear temporal logic (LTL) specifications are currently being implemented.

# Solvers

The following tables summarize the different solver types in abbreviated form. 


**Abstraction-based solver types implemented in Dionysos.jl:**
| Type          | Full vs partial discretization | Partition vs Cover | Shape | Local controller | Abstraction | System | Reference | 
| :--------------- | :---------- | :---------- | :---------- | :---------- | :---------- | :---------- | :---------- |
| [`SCOTS`](@ref) | Full | Partition | Hyperrectangle | Piece-wise constant | Non-determinisitic | Continuous-time | [`SCOTS: A Tool for the Synthesis of Symbolic Controllers`](https://dl.acm.org/doi/abs/10.1145/2883817.2883834) |
| [`Lazy abstraction`](@ref) | Partial | Partition | Hyperrectangle | Piece-wise constant | Non-determinisitic | Continuous-time | [`Alternating simulation on hierarchical abstractions`](https://ieeexplore.ieee.org/abstract/document/9683448/?casa_token=AXyECYU9FdwAAAAA:ERfvlbkORIbfGLbDd42d2K5K9SLVOjl-8kRs9pfp7lMa4QZEv0W4VgzlP8FshBlxXQF4ZQrDuzk) |
| [`Hierarchical`](@ref) | Partial | Partition | Hyperrectangle | Piece-wise constant | Non-determinisitic | Continuous-time | [`Alternating simulation on hierarchical abstractions`](https://ieeexplore.ieee.org/abstract/document/9683448/?casa_token=AXyECYU9FdwAAAAA:ERfvlbkORIbfGLbDd42d2K5K9SLVOjl-8kRs9pfp7lMa4QZEv0W4VgzlP8FshBlxXQF4ZQrDuzk) |
| [`Ellipsoid absraction`](@ref) | Full | Cover | Ellipsoid | Piece-wise affine | Determinisitic | Discrete-time affine | [`State-feedback Abstractions for Optimal Control of Piecewise-affine Systems`](https://arxiv.org/abs/2204.00315) |
| [`Lazy ellipsoid absraction`](@ref) | Partial | Cover | Ellipsoid | Piece-wise affine | Determinisitic | Discrete-time non-linear | Not yet published |

**Non abstraction-based solver types implemented in Dionysos.jl:**
| Type          | Description | Reference |
| :--------------- | :---------- | :---------- |
| [`Bemporad Morari`](@ref) |  | [`Control of systems integrating logic, dynamics, and constraints`](https://www.sciencedirect.com/science/article/abs/pii/S0005109898001782)
| [`BranchAndBound`](@ref) |  | [`Abstraction-based branch and bound approach to Q-learning for hybrid optimal control`](https://proceedings.mlr.press/v144/legat21a.html)


**Solver interface**

Each solver is defined by a module which must define the structure [`AbstractOptimizer`](https://jump.dev/MathOptInterface.jl/stable/reference/models/#MathOptInterface.AbstractOptimizer) and implement the [`Optimize!`](https://jump.dev/MathOptInterface.jl/stable/reference/models/#MathOptInterface.optimize!) function.
For example, for the SCOTS solver, this structure and function are defined as follows

```
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

```
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
For an executable version of this example, see [`Example: Path planning problem`](https://dionysos-dev.github.io/Dionysos.jl/dev/generated/Path%20planning/#Example:-Path-planning-problem) in the documentation.

First you need to define a control problem, i.e., the system and the specification of the desired closed loop behaviour.
To do this, you can define new ones yourself or directly load an existing benchmark, for example

```
concrete_problem = PathPlanning.problem(; simple = true, approx_mode = "growth");
concrete_system = concrete_problem.system;
```

Choose the solver you wish to use
```
using JuMP
optimizer = MOI.instantiate(AB.SCOTSAbstraction.Optimizer)
```

Define the solver's meta-parameters
```
x0 = SVector(0.0, 0.0, 0.0);
hx = SVector(0.2, 0.2, 0.2);
state_grid = DO.GridFree(x0, hx);
u0 = SVector(0.0, 0.0);
hu = SVector(0.3, 0.3);
input_grid = DO.GridFree(u0, hu);
```

Set the solver's meta-parameters
```
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
```

Solve the control problem
```
MOI.optimize!(optimizer)
```

Get the results
```
abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
```

In Dionysos, all the structures that could be relevant to plot (such as trajectories, state-space discretization, specifications, obstacles, etc.) 
have an associated [`@recipe`](https://github.com/JuliaPlots/RecipesBase.jl) function, which makes it very easy to plot all the results using the single common [`plot`](https://docs.juliaplots.org/latest/generated/unitfulext_plots/) function of [`Plots.jl`](https://github.com/JuliaPlots/Plots.jl).
For example
```
using Plots

plot!(concrete_system.X; color = :yellow, opacity = 0.5);
plot!(abstract_system.Xdom; color = :blue, opacity = 0.5);
plot!(concrete_problem.initial_set; color = :green, opacity = 0.2);
plot!(concrete_problem.target_set; dims = [1, 2], color = :red, opacity = 0.2);
plot!(UT.DrawTrajectory(x_traj); ms = 0.5)
```



