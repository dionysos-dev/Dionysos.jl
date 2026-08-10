# ----------------------------------------------------------------------------------------
# Parsing constraints written on a mode or a transition.
#
# They arrive with their set wrapped in a `ScopedSet`/`ScopedVectorSet` (see `modes.jl`), so
# the routing is by set type and the scope says where the content belongs. Guard versus reset
# needs no heuristic: a `ScalarNonlinearFunction` headed by `∂`/`Δ` is a reset map, an affine
# function is a guard or a bound.
# ----------------------------------------------------------------------------------------

# Interval implied by a scalar set, before the affine function is inverted.
_set_bounds(set::MOI.LessThan) = (-Inf, set.upper)
_set_bounds(set::MOI.GreaterThan) = (set.lower, Inf)
_set_bounds(set::MOI.EqualTo) = (set.value, set.value)
_set_bounds(set::MOI.Interval) = (set.lower, set.upper)

# `coef * x + constant ∈ set` as bounds on `x`, for the single-variable case.
function _variable_bounds(func::MOI.ScalarAffineFunction, set)
    term = func.terms[]
    iszero(term.coefficient) && error("A bound or guard must have a non-zero coefficient.")
    lo, hi = _set_bounds(set)
    a = (lo - func.constant) / term.coefficient
    b = (hi - func.constant) / term.coefficient
    return term.variable.value, min(a, b), max(a, b)
end

# `Σ aᵢ xᵢ + constant ∈ set` with more than one variable. The constant moves to the bounds so
# what is kept is exactly the coefficients and the interval they lie in.
function _affine_guard(func::MOI.ScalarAffineFunction, set)
    coefficients = Dict{Int, Float64}()
    for term in func.terms
        coefficients[term.variable.value] =
            get(coefficients, term.variable.value, 0.0) + term.coefficient
    end
    lo, hi = _set_bounds(set)
    return AffineGuard(coefficients, lo - func.constant, hi - func.constant)
end

function _record_bounds!(lower::Dict{Int, Float64}, upper::Dict{Int, Float64}, i, lo, hi)
    isfinite(lo) && (lower[i] = haskey(lower, i) ? max(lower[i], lo) : lo)
    isfinite(hi) && (upper[i] = haskey(upper, i) ? min(upper[i], hi) : hi)
    return
end

# ---- affine constraints: per-mode bounds, or a transition guard ----

const _BoundSet = Union{MOI.LessThan, MOI.GreaterThan, MOI.EqualTo, MOI.Interval}

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{<:MOI.ScalarAffineFunction},
    ::Type{<:ScopedSet{<:_BoundSet}},
)
    return true
end

function MOI.add_constraint(
    model::Optimizer,
    func::MOI.ScalarAffineFunction,
    set::ScopedSet{<:_BoundSet},
)
    if length(func.terms) == 1
        i, lo, hi = _variable_bounds(func, set.inner)
        _apply_affine!(model.ir, set.scope, i, lo, hi)
    else
        _apply_multi_affine!(model.ir, set.scope, _affine_guard(func, set.inner))
    end
    return next_constraint_index!(model.ir, func, set)
end

function _apply_affine!(ir::ModelIR, scope::ModeScope, i, lo, hi)
    m = mode!(ir, scope.id)
    return _record_bounds!(m.lower, m.upper, i, lo, hi)
end

function _apply_affine!(ir::ModelIR, scope::TransitionScope, i, lo, hi)
    t = transition!(ir, scope)
    return _record_bounds!(t.guard_lower, t.guard_upper, i, lo, hi)
end

# A mode's state and input sets are boxes the abstraction discretizes coordinate by coordinate,
# so a multi-variable constraint has nowhere to go there. On a transition it becomes a
# half-space, which the guard — an arbitrary set — can carry.
function _apply_multi_affine!(::ModelIR, scope::ModeScope, ::AffineGuard)
    return error(
        "A mode's bounds are per-coordinate, so a constraint over several variables cannot " *
        "be written on mode $(scope.id); the state and input sets are boxes. Write one " *
        "constraint per variable. (On a *transition* a multi-variable constraint is allowed: " *
        "it becomes a half-space of the guard.)",
    )
end

function _apply_multi_affine!(ir::ModelIR, scope::TransitionScope, guard::AffineGuard)
    push!(transition!(ir, scope).affine_guards, guard)
    return
end

# ---- `[x, y] in Guard(S)`: a set-valued guard on a transition ----

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VectorOfVariables},
    ::Type{<:ScopedVectorSet{<:Guard}},
)
    return true
end

function MOI.add_constraint(
    model::Optimizer,
    func::MOI.VectorOfVariables,
    set::ScopedVectorSet{<:Guard},
)
    scope = set.scope
    scope isa TransitionScope || error(
        "`Guard` belongs to a transition, not to a mode (mode $(scope.id)). A mode's own " *
        "region is written with `Always`.",
    )
    push!(transition!(model.ir, scope).guard_sets, (copy(func.variables), set.inner.inner))
    return MOI.ConstraintIndex{typeof(func), typeof(set)}(
        next_constraint_index!(model.ir, func, set).value,
    )
end

# ---- nonlinear equalities: per-mode dynamics, or a transition reset map ----

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.ScalarNonlinearFunction},
    ::Type{<:ScopedSet{<:MOI.EqualTo}},
)
    return true
end

function MOI.add_constraint(
    model::Optimizer,
    func::MOI.ScalarNonlinearFunction,
    set::ScopedSet{<:MOI.EqualTo},
)
    ir = model.ir
    if func.head == :- && length(func.args) == 2
        lhs, rhs = func.args
        if lhs isa MOI.ScalarNonlinearFunction && length(lhs.args) == 1
            var = lhs.args[]
            if var isa MOI.VariableIndex
                if lhs.head == :∂ || lhs.head == :Δ
                    _apply_dynamics!(ir, set.scope, var.value, rhs, lhs.head)
                    return next_constraint_index!(ir, func, set)
                elseif lhs.head == :final
                    _apply_final!(ir, set.scope, var.value, MOI.Interval(rhs, rhs))
                    return next_constraint_index!(ir, func, set)
                end
            end
        end
    end
    dump(func)
    return error("Unsupported")
end

function _apply_dynamics!(ir::ModelIR, scope::ModeScope, i, expr, head::Symbol)
    _set_time_domain!(ir, head === :∂ ? CONTINUOUS : DISCRETE)
    mode!(ir, scope.id).dynamics[i] = expr
    return
end

# On a transition, `Δ(x) == …` is the reset map: where the state lands when the switch is
# taken. It is an instantaneous jump, so it says nothing about whether the *modes* evolve in
# continuous or discrete time — `∂` and `Δ` read the same here, and neither fixes the domain.
function _apply_dynamics!(ir::ModelIR, scope::TransitionScope, i, expr, ::Symbol)
    transition!(ir, scope).resets[i] = expr
    return
end

# ---- `final(x) in …` written on a mode ----

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.ScalarNonlinearFunction},
    ::Type{<:ScopedSet{<:MOI.Interval}},
)
    return true
end

function MOI.add_constraint(
    model::Optimizer,
    func::MOI.ScalarNonlinearFunction,
    set::ScopedSet{<:MOI.Interval},
)
    if length(func.args) == 1
        var = func.args[]
        if var isa MOI.VariableIndex && func.head in (:final, :start)
            _apply_final!(model.ir, set.scope, var.value, set.inner, func.head)
            return next_constraint_index!(model.ir, func, set)
        end
    end
    dump(func)
    return error("Unsupported")
end

function _apply_final!(ir::ModelIR, scope::ModeScope, i, interval, head::Symbol = :final)
    kind = head === :start ? START : FINAL
    m = mode!(ir, scope.id)
    push!(m.specs, SpecEntry(kind, interval, [MOI.VariableIndex(i)]))
    return
end

_apply_final!(::ModelIR, scope::TransitionScope, args...) = error(
    "A `final`/`start` constraint cannot be written on a transition (transition $(scope.id)); " *
    "put it on a mode.",
)

# ---- specification markers written on a mode ----

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VectorOfVariables},
    ::Type{<:ScopedVectorSet{<:SpecSet}},
)
    return true
end

function MOI.add_constraint(
    model::Optimizer,
    func::MOI.VectorOfVariables,
    set::ScopedVectorSet{<:SpecSet},
)
    scope = set.scope
    scope isa ModeScope || error(
        "A specification set belongs to a mode, not to a transition (transition $(scope.id)).",
    )
    entry = SpecEntry(
        spec_kind(set.inner),
        set.inner.inner,
        copy(func.variables),
        stay_on_first_entry(set.inner),
    )
    push!(mode!(model.ir, scope.id).specs, entry)
    return MOI.ConstraintIndex{typeof(func), typeof(set)}(
        next_constraint_index!(model.ir, func, set).value,
    )
end
