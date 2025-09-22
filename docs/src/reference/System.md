# System Representations and Approximations

This module defines tools to represent and manipulate dynamical systems and their approximations.

## Concrete Systems

The systems we aim to control are defined using types from external packages such as:

- [`MathematicalSystems.jl`](https://github.com/JuliaReach/MathematicalSystems.jl) for standard (e.g., continuous/discrete-time) control systems.
- [`HybridSystems.jl`](https://github.com/blegat/HybridSystems.jl) for hybrid automata.

For example, a continuous-time system might be defined as:

```julia
concrete_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
    dynamics(),  # system dynamics (function)
    n_X,         # state dimension
    n_U,         # input dimension
    _X_,         # state constraints
    _U_          # input constraints
)
```

## Symbolic/Abstract Systems

Symbolic models (used for abstraction-based control) are constructed separately using the symbolic abstraction module, typically resulting in a [`SymbolicModelList`](@ref Dionysos.Symbolic.SymbolicModelList).

## System Approximations

To reason about the system’s behavior during abstraction, we introduce **approximations** of the system's evolution. These are grouped into:

- [`DiscreteTimeSystemApproximation`](@ref Dionysos.System.DiscreteTimeSystemApproximation)
- [`ContinuousTimeSystemApproximation`](@ref Dionysos.System.ContinuousTimeSystemApproximation)

Both are subtypes of [`SystemApproximation`](@ref Dionysos.System.SystemApproximation), and can represent either **underapproximations** or **overapproximations** of the system dynamics.

```@docs
Dionysos.System.SystemApproximation
Dionysos.System.DiscreteTimeSystemApproximation
Dionysos.System.ContinuousTimeSystemApproximation
```

The following functions define the `SystemApproximation` interface:

- `get_system(approx::SystemApproximation)`: Returns the underlying concrete system.
- `is_continuous_time(approx)`: Returns `true` if the approximation is continuous-time, i.e., a [`ContinuousTimeSystemApproximation`](@ref Dionysos.System.ContinuousTimeSystemApproximation).
- `is_over_approximation(approx::SystemApproximation)`: Return `true` if `approx` is a [`DiscreteTimeSystemOverApproximation`](@ref Dionysos.System.DiscreteTimeSystemOverApproximation) or a [`ContinuousTimeSystemOverApproximation`](@ref Dionysos.System.ContinuousTimeSystemOverApproximation).
- `discretize(approx::ContinuousTimeSystemApproximation, tstep)::DiscreteTimeSystemApproximation`: Returns a  [`DiscreteTimeSystemApproximation`](@ref Dionysos.System.DiscreteTimeSystemApproximation) with given time step.
- `get_system_map(approx)`: Returns the map representing the system's evolution.

### Underapproximations

These approximations guarantee that all returned trajectories are feasible under the system dynamics.

```@docs
Dionysos.System.DiscreteTimeSystemUnderApproximation
Dionysos.System.ContinuousTimeSystemUnderApproximation
```

The following function define the `underapproximation` interface:

```@docs
Dionysos.System.get_under_approximation_map
```

### Overapproximations

These approximations guarantee that the true system evolution is **contained** in the returned set, making them useful for safety and robust control.

```@docs
Dionysos.System.DiscreteTimeSystemOverApproximation
Dionysos.System.ContinuousTimeSystemOverApproximation
```

The following function define the `overapproximation` interface:

```@docs
Dionysos.System.get_over_approximation_map
```

### Concrete implementations of abstract approximation types
```@docs
Dionysos.System.DiscreteTimeCenteredSimulation
Dionysos.System.ContinuousTimeCenteredSimulation
```

```@docs
Dionysos.System.DiscreteTimeRandomSimulation
Dionysos.System.ContinuousTimeRandomSimulation
```

```@docs
Dionysos.System.DiscreteTimeOverApproximationMap
Dionysos.System.ContinuousTimeSystemOverApproximationMap
```

```@docs
Dionysos.System.DiscreteTimeGrowthBound
Dionysos.System.ContinuousTimeGrowthBound
```

```@docs
Dionysos.System.DiscreteTimeLinearized
Dionysos.System.ContinuousTimeLinearized
```

## Trajectories 
```@docs
Dionysos.System.wrap_coord
Dionysos.System.DiscreteTrajectory
Dionysos.System.ContinuousTrajectory
Dionysos.System.HybridTrajectory
Dionysos.System.Trajectory
Dionysos.System.Control_trajectory
Dionysos.System.Cost_control_trajectory
```