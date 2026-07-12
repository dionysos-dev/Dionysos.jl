# System

The `System` module represents and manipulates dynamical systems and their approximations. It extends
[`MathematicalSystems`](https://github.com/JuliaReach/MathematicalSystems.jl) (and
[`HybridSystems`](https://github.com/blegat/HybridSystems.jl) for hybrid automata) rather than
redefining systems from scratch.

## Concrete systems

The systems we control are ordinary `MathematicalSystems` objects, e.g. a continuous-time system:

```julia
concrete_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
    dynamics,  # f(x, u)
    n_X,       # state dimension
    n_U,       # input dimension
    _X_,       # state constraints
    _U_,       # input constraints
)
```

## System approximations

To reason about a system's behaviour during abstraction, `System` introduces **approximations** of
the evolution map, split by time domain into
[`DiscreteTimeSystemApproximation`](@ref Dionysos.System.DiscreteTimeSystemApproximation) and
[`ContinuousTimeSystemApproximation`](@ref Dionysos.System.ContinuousTimeSystemApproximation) (both
subtypes of [`SystemApproximation`](@ref Dionysos.System.SystemApproximation)), and by soundness into:

- **Underapproximations** — every returned trajectory is feasible under the true dynamics (used to
  *find* trajectories).
- **Overapproximations** — the true evolution is *contained* in the returned set (used for sound,
  robust abstraction). Growth bounds ([`ContinuousTimeGrowthBound`](@ref Dionysos.System.ContinuousTimeGrowthBound))
  and linearizations ([`ContinuousTimeLinearized`](@ref Dionysos.System.ContinuousTimeLinearized)) are
  the two overapproximation strategies behind the `GROWTH` / `LINEARIZED` abstraction modes.

A continuous-time approximation is turned into a discrete-time one with `discretize(approx, tstep)`.

## Affine approximation

Local affine approximations of nonlinear dynamics (linearization plus Lipschitz error bounds),
consumed by the lazy-ellipsoids abstraction and the ellipsoidal backward trajectory certifier. Built
via [`build_affine_approximation`](@ref Dionysos.System.build_affine_approximation) and its
providers.

## Controllers

Controllers follow a trait-based protocol (argument order `(controller, memory, measurement)`):
`controller_kind`, `output_control`, `is_defined`, and — for dynamic controllers — `initial_state` /
`update_state`, with static controllers inheriting memoryless defaults. Controllers are **plain data**
(tables, sets, struct callables), so they can be serialized to JLD2 and reloaded. See also the small
PID controllers in [`PIDControllers`](@ref Dionysos.System.PIDControllers).

## Transition synthesis

Synthesis of affine controllers certifying ellipsoid-to-ellipsoid transitions of affine systems (via
S-procedure LMIs), used by the ellipsoidal abstraction solvers and the backward trajectory certifier:
[`solve_transition`](@ref Dionysos.System.solve_transition).

## Trajectories

Discrete and continuous trajectory types plus the closed-loop simulation engine
([`get_closed_loop_trajectory`](@ref Dionysos.System.get_closed_loop_trajectory)).

## API reference

```@autodocs
Modules = [Dionysos.System]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```

## PID controllers

```@autodocs
Modules = [Dionysos.System.PIDControllers]
Filter  = _is_public
Order   = [:module, :type, :function, :constant]
```
