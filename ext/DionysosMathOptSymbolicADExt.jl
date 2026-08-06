module DionysosMathOptSymbolicADExt

# The only Symbolics-dependent part of the JuMP front-end: compiling the parsed dynamics
# expressions into a callable `f(x, u)`. Everything else (parsing, roles, lowering, the
# optimizer itself) lives in `src/wrapper/` and needs no optional dependency.

import Dionysos
const Wrapper = Dionysos.Wrapper

import MathOptInterface as MOI
import Symbolics
import MathOptSymbolicAD

# One dynamics expression, re-expressed over the Symbolics variables `xu` with its stored
# parameters substituted back in.
function _symbolic(func::MathOptSymbolicAD._Function, xu)
    Symbolics.@variables(p[1:length(func.data)])
    e = MathOptSymbolicAD._expr_to_symbolics(
        func.model,
        func.expr,
        p,
        xu[func.ordered_variables],
    )
    return Symbolics.substitute(e, Dict(p[i] => func.data[i] for i in eachindex(p)))
end

function Wrapper.compile_dynamics(
    ::Wrapper.SymbolicADBackend,
    ir::Wrapper.ModelIR,
    dynamics::AbstractVector,
)
    x_idx = Wrapper.state_indices(ir)
    u_idx = Wrapper.input_indices(ir)

    nlp_model = MOI.Nonlinear.Model()
    for i in x_idx
        # The set does not matter; the constraint is only a carrier for the expression.
        MOI.Nonlinear.add_constraint(nlp_model, dynamics[i], MOI.EqualTo(0.0))
    end

    vars = MOI.VariableIndex.(eachindex(ir.variables))
    evaluator = MOI.Nonlinear.Evaluator(nlp_model, MathOptSymbolicAD.DefaultBackend(), vars)
    MOI.initialize(evaluator, Symbol[])

    # System eq x' = F_sys(x, u)
    Symbolics.@variables(x[1:length(x_idx)], u[1:length(u_idx)])
    xu = Vector{eltype(x)}(undef, length(ir.variables))
    xu[x_idx] = x
    xu[u_idx] = u
    expr = [_symbolic(func, xu) for func in evaluator.backend.constraints]

    # The second returned function is the in-place version, which we don't want.
    # `expression = Val{false}` yields a Julia function directly; `cse = true` detects common
    # sub-expressions (like `α` in the path planning example), cf.
    # https://discourse.julialang.org/t/detecting-function-composition-in-symbolics-jl/115885/6
    F_sys, _ = Symbolics.build_function(
        expr,
        collect(x),
        collect(u);
        expression = Val{false},
        cse = true,
    )
    return F_sys
end

end
