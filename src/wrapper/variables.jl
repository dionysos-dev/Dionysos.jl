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
    # Set by `set_role!` when inference cannot see the truth — typically because the dynamics
    # were supplied as a Julia function rather than written as equations.
    declared_role::Union{Nothing, VariableRole}
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

"""
    set_role!(x, role)

Declare what `x` is, instead of letting the wrapper infer it. `x` may be a single variable or
an array of them, and `role` one of `Dionysos.STATE`, `Dionysos.INPUT`, `Dionysos.CLOCK`,
`Dionysos.DISTURBANCE`.

Needed in two situations where inference cannot see the truth. When the dynamics are supplied
as a Julia function there are no expressions to infer from, so the states have to be named:

```julia
set_role!(x, Dionysos.STATE)
set_attribute(model, "dynamics", (x, u) -> [x[2], -sin(x[1]) + u[1]])
```

And a **disturbance** is indistinguishable from an input by usage — both merely appear on a
dynamics right-hand side — so it must always be said. Declaring it changes who owns the signal:
the synthesized controller must then work for *every* value the disturbance takes in its bounds,
rather than being free to choose it:

```julia
@variable(model, -0.1 <= w <= 0.1)
set_role!(w, Dionysos.DISTURBANCE)
@constraint(model, ∂(x) == -x + u + w)
```

Variables left undeclared are inputs.
"""
function set_role!(x::JuMP.GenericVariableRef, role::VariableRole)
    role in (STATE, INPUT, CLOCK, DISTURBANCE) || error(
        "`set_role!` takes STATE, INPUT, CLOCK or DISTURBANCE; $(role) has no meaning for " *
        "the lowering yet.",
    )
    JuMP.set_attribute(JuMP.owner_model(x), "role[$(JuMP.index(x).value)]", role)
    return role
end

function set_role!(x::AbstractArray, role::VariableRole)
    for v in x
        set_role!(v, role)
    end
    return role
end

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
