# ----------------------------------------------------------------------------------------
# The JuMP-level vocabulary: the nonlinear operators used to write dynamics and
# specifications, plus the obstacle syntax.
# ----------------------------------------------------------------------------------------

export ∂, Δ, final, start

function _diff end
"""
    ∂(x)

Continuous-time derivative of the state `x`, used to declare dynamics:
`@constraint(model, ∂(x[1]) == u[1])`. A variable carrying a `∂` (or [`Δ`](@ref))
constraint is a **state**; see [`infer_roles!`](@ref).
"""
const ∂ = JuMP.NonlinearOperator(_diff, :∂)

function _delta end
"""
    Δ(x)

Discrete-time successor of the state `x`: `@constraint(model, Δ(x[1]) == x[1] + u[1])`.
Mixing `Δ` and [`∂`](@ref) in one model is an error.
"""
const Δ = JuMP.NonlinearOperator(_delta, :Δ)

function _final end
"""
    final(x)

The value of `x` at the end of the trajectory, used to declare a reach target:
`@constraint(model, final(x[1]) in MOI.Interval(a, b))`.
"""
const final = JuMP.NonlinearOperator(_final, :final)

function _start end
"""
    start(x)

The value of `x` at the beginning of the trajectory, used to declare an initial set:
`@constraint(model, start(x[1]) in MOI.Interval(a, b))`. The `@variable(model, x, start = v)`
keyword is the singleton special case.
"""
const start = JuMP.NonlinearOperator(_start, :start)

"""
    OuterSet{S} <: MOI.AbstractVectorSet

Marks the complement of `inner`: `@constraint(model, x ∉ O)` wraps `O` in an `OuterSet`,
and the lowering carves it out of the state space `X`.
"""
struct OuterSet{S <: MOI.AbstractVectorSet} <: MOI.AbstractVectorSet
    inner::S
end

MOI.dimension(set::OuterSet) = MOI.dimension(set.inner)
Base.copy(set::OuterSet) = OuterSet(copy(set.inner))

# Deliberate type piracy: JuMP has no parsing rule for `∉`; this method makes the public obstacle
# syntax `@constraint(model, x ∉ obstacle)` work by wrapping the set in `OuterSet`. Purely additive
# (JuMP errors on `∉` otherwise). Remove once JuMP exposes an extension point for it.
function JuMP.parse_constraint_call(
    error_fn::Function,
    vectorized::Bool,
    ::Val{:∉},
    lhs,
    rhs,
)
    @assert !vectorized
    f, parse_code1 = JuMP._rewrite_expression(lhs)
    set, parse_code2 = JuMP._rewrite_expression(rhs)
    parse_code = quote
        $parse_code1
        $parse_code2
    end
    build_call = :(build_constraint($error_fn, $f, $(OuterSet)($set)))
    return parse_code, build_call
end
