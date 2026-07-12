"""
    AbstractLift

A composable, optional generalization applied to a base [`AbstractSymbolicModel`](@ref):
it adds one factor (an extra axis) to the abstraction. Concrete lifts:

- [`ClockLift`](@ref) — adds a monotone time axis `t`;
- a mode lift (planned) — adds the hybrid mode axis `k`.

The four abstractions `x`, `(x,t)`, `(x,k)`, `(x,k,t)` are obtained by applying
zero, one, or both lifts to a base model. Each lift implements `lift(l, model)`
(and, in a later phase, `lift(l, spec)` so specifications extend the same way).
"""
abstract type AbstractLift end

"""
    lift(l::AbstractLift, base::AbstractSymbolicModel) -> AbstractSymbolicModel

Apply the lift `l` to `base`, returning a symbolic model with one additional factor.
"""
lift(l::AbstractLift, base::AbstractSymbolicModel) =
    error("implement `lift` for $(typeof(l)) on $(typeof(base))")
