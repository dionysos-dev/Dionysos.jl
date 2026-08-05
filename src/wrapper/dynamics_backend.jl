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
    compile_dynamics(backend, ir::ModelIR, dynamics = ir.dynamics) -> f

Compile `dynamics` — one expression per variable, `nothing` for the non-states — into a
callable `f(x, u)` returning the state derivative (continuous time) or the successor state
(discrete time). `x` and `u` follow the orders given by [`state_indices`](@ref) and
[`input_indices`](@ref).

The dynamics are passed explicitly rather than read from `ir`, because a hybrid model compiles
one function per mode, and a transition's reset map goes through the same path.
"""
compile_dynamics(backend::AbstractDynamicsBackend, ir::ModelIR) =
    compile_dynamics(backend, ir, ir.dynamics)

function compile_dynamics(backend::AbstractDynamicsBackend, ::ModelIR, ::AbstractVector)
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
    NonlinearEvaluatorBackend <: AbstractDynamicsBackend

Evaluate the dynamics expressions directly with `MOI.Nonlinear`, needing no optional
dependency at all.

This is the **test and fallback** path, not the production one: evaluation is interpreted and
allocates on every call, and the abstraction calls it millions of times. Its purpose is to make
the front-end loadable, testable and documentable without Symbolics —
[`SymbolicADBackend`](@ref) is what you want for real problems.
"""
struct NonlinearEvaluatorBackend <: AbstractDynamicsBackend end

function compile_dynamics(
    ::NonlinearEvaluatorBackend,
    ir::ModelIR,
    dynamics::AbstractVector,
)
    x_idx = state_indices(ir)
    u_idx = input_indices(ir)

    nlp_model = MOI.Nonlinear.Model()
    for i in x_idx
        MOI.Nonlinear.add_constraint(nlp_model, dynamics[i], MOI.EqualTo(0.0))
    end
    vars = MOI.VariableIndex.(eachindex(ir.variables))
    evaluator = MOI.Nonlinear.Evaluator(nlp_model, MOI.Nonlinear.SparseReverseMode(), vars)
    MOI.initialize(evaluator, Symbol[])

    nvars = length(ir.variables)
    nstates = length(x_idx)
    return function (x, u)
        point = zeros(nvars)
        for (j, i) in enumerate(x_idx)
            point[i] = x[j]
        end
        for (j, i) in enumerate(u_idx)
            point[i] = u[j]
        end
        out = zeros(nstates)
        MOI.eval_constraint(evaluator, out, point)
        return out
    end
end

"""
    default_dynamics_backend()

The backend used when the user sets none.
"""
default_dynamics_backend() = SymbolicADBackend()
