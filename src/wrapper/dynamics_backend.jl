# ----------------------------------------------------------------------------------------
# Dynamics backends: how the parsed dynamics *expressions* become a callable `f(x, u)`.
#
# This is the only part of the front-end that needs a heavyweight optional dependency, which
# is why it is a swappable protocol rather than a hardwired path.
# ----------------------------------------------------------------------------------------

"""
    AbstractDynamicsBackend

Strategy for turning the dynamics expressions of a [`ModelIR`](@ref) into a callable
`f(x, u)`. Implement [`compile_dynamics`](@ref) for a new backend.

Available:
- [`SymbolicADBackend`](@ref) — traces the expressions with Symbolics and emits a fused
  function; the production path, provided by `DionysosMathOptSymbolicADExt`.
"""
abstract type AbstractDynamicsBackend end

"""
    SymbolicADBackend <: AbstractDynamicsBackend

Compile the dynamics by tracing them symbolically and building a fused Julia function with
common-subexpression elimination. Requires `Symbolics` and `MathOptSymbolicAD` to be loaded;
the implementation lives in `DionysosMathOptSymbolicADExt`.
"""
struct SymbolicADBackend <: AbstractDynamicsBackend end

"""
    compile_dynamics(backend, ir::ModelIR) -> f

Compile the dynamics of `ir` into a callable `f(x, u)` returning the state derivative
(continuous time) or the successor state (discrete time). `x` and `u` follow the orders given
by [`state_indices`](@ref) and [`input_indices`](@ref).
"""
function compile_dynamics(backend::AbstractDynamicsBackend, ::ModelIR)
    # Only a fallback: each backend's real method is more specific. `SymbolicADBackend` is
    # handled here rather than by its own stub method, because the extension defines exactly
    # that signature and Julia forbids overwriting it during precompilation.
    if backend isa SymbolicADBackend
        error(
            "The JuMP front-end needs Symbolics.jl and MathOptSymbolicAD.jl to compile the " *
            "dynamics; load them with `using Symbolics, MathOptSymbolicAD`.",
        )
    end
    return error("implement `compile_dynamics` for $(typeof(backend))")
end

"""
    default_dynamics_backend()

The backend used when the user sets none.
"""
default_dynamics_backend() = SymbolicADBackend()
