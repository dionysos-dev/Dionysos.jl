# System

Represents and manipulates dynamical systems and their approximations, extending
[`MathematicalSystems`](https://github.com/JuliaReach/MathematicalSystems.jl) and
[`HybridSystems`](https://github.com/blegat/HybridSystems.jl). It covers system approximations
(over-/under-approximations of the dynamics, growth bounds and linearizations), local affine
approximations, the trait-based controller protocol, ellipsoidal transition synthesis
([`solve_transition`](@ref Dionysos.System.solve_transition)), and trajectories with the closed-loop
simulation engine ([`get_closed_loop_trajectory`](@ref Dionysos.System.get_closed_loop_trajectory)).
Controllers are plain data, so they can be serialized and reloaded.

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
