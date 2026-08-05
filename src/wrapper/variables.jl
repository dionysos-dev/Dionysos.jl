# ----------------------------------------------------------------------------------------
# Variables and their roles.
#
# The wrapper never asks the user to declare what a variable *is*; the role is inferred from
# how the variable is used across the complete model. See `src/wrapper/README.md` §3.
# ----------------------------------------------------------------------------------------

"""
    VariableRole

What a declared JuMP variable turns out to be, inferred by [`infer_roles!`](@ref):

- `STATE` — carries `∂`/`Δ` dynamics;
- `INPUT` — appears on some dynamics right-hand side and is not a state;
- `CLOCK` — a state whose dynamics is exactly `∂(t) == 1` (running) or `∂(t) == 0` (frozen);
- `PARAMETER` — fixed and inlined rather than discretized;
- `DISTURBANCE` — an exogenous input, only ever set explicitly.
"""
@enum VariableRole STATE INPUT CLOCK PARAMETER DISTURBANCE

"""
    VariableInfo

Everything the front-end knows about one JuMP variable: its bounds, the initial/target
intervals declared for it, its inferred [`VariableRole`](@ref), its position within that role
group, and its JuMP name (used to make error messages nameable).
"""
mutable struct VariableInfo
    name::String
    lower::Float64
    upper::Float64
    start::MOI.Interval{Float64}
    target::MOI.Interval{Float64}
    role::Union{Nothing, VariableRole}
    index::Int
end

function VariableInfo()
    return VariableInfo(
        "",
        -Inf,
        Inf,
        MOI.Interval(-Inf, Inf),
        MOI.Interval(-Inf, Inf),
        nothing,
        0,
    )
end

# A variable the user constrained in some way other than through the dynamics.
_has_start(v::VariableInfo) = isfinite(v.start.lower) || isfinite(v.start.upper)
_has_target(v::VariableInfo) = isfinite(v.target.lower) || isfinite(v.target.upper)

"""
    describe(v::VariableInfo, i::Int) -> String

How to refer to variable `i` in a user-facing message: its JuMP name when JuMP gave it one,
otherwise its MOI index.
"""
describe(v::VariableInfo, i::Int) = isempty(v.name) ? "variable #$i" : "`$(v.name)`"

# ----------------------------------------------------------------------------------------
# Collecting the variables referenced by an MOI function
# ----------------------------------------------------------------------------------------

_collect_variables!(acc::Set{Int}, ::Any) = acc
_collect_variables!(acc::Set{Int}, v::MOI.VariableIndex) = (push!(acc, v.value); acc)

function _collect_variables!(acc::Set{Int}, f::MOI.ScalarNonlinearFunction)
    for arg in f.args
        _collect_variables!(acc, arg)
    end
    return acc
end

function _collect_variables!(acc::Set{Int}, f::MOI.ScalarAffineFunction)
    for t in f.terms
        push!(acc, t.variable.value)
    end
    return acc
end

function _collect_variables!(acc::Set{Int}, f::MOI.ScalarQuadraticFunction)
    for t in f.affine_terms
        push!(acc, t.variable.value)
    end
    for t in f.quadratic_terms
        push!(acc, t.variable_1.value)
        push!(acc, t.variable_2.value)
    end
    return acc
end
