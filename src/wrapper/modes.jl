# ----------------------------------------------------------------------------------------
# Modes and transitions.
#
# A mode or a transition is a *scope* you write ordinary `@constraint`s on. Both subtype
# `JuMP.AbstractModel`, which is all JuMP needs for `@constraint(scope, …)` to work — no type
# piracy, since `_valid_model(::AbstractModel, ::Any) = nothing` already holds.
#
# A scoped constraint is forwarded to the parent model with its *set* wrapped in a
# [`ScopedSet`](@ref) carrying the scope. That matters: constraints written on a scope then
# travel through JuMP's caching and bridging layers like any other constraint, instead of
# being stashed in the optimizer behind their back. It also makes them self-describing — the
# whole hybrid topology is recoverable from the constraints alone.
# ----------------------------------------------------------------------------------------

export add_mode!, add_transition!, @mode, Guard

"""
    Guard(S)

Marks `S` as the guard of a transition — the states from which the switch may be taken:

```julia
add_transition!(model, a => b) do t
    @constraint(t, [x, y] in Guard(LazySets.Ball2([0.0, 0.0], 1.0)))
end
```

`S` may be any bounded `LazySet` and must span the state vector, in declaration order. The
marker is needed because JuMP has no parsing rule that would tell a bare set apart from a
bound; it plays the same role as [`Final`](@ref) and [`Always`](@ref) do on a mode.

A guard can equally be written as ordinary inequalities — `@constraint(t, x <= 1)` per
coordinate, or a multi-variable `@constraint(t, x + y <= 1)`, which becomes a half-space. All
of them intersect: a transition's guard is everything written on it, narrowed by the source
mode's own state set.
"""
struct Guard{S} <: MOI.AbstractVectorSet
    inner::S
    dim::Int
end

Guard(set) = Guard(set, LazySets.dim(set))
MOI.dimension(set::Guard) = set.dim
Base.copy(set::Guard) = set
JuMP.in_set_string(::MIME, set::Guard) = "in Guard($(nameof(typeof(set.inner))))"

"""
    ModeScope

Marks a constraint as belonging to mode `id`.
"""
struct ModeScope
    id::Int
end

"""
    TransitionScope

Marks a constraint as belonging to the transition `id` from mode `source` to mode `target`.
The switching type rides along so a transition is fully described by any one of its
constraints.
"""
struct TransitionScope{SW}
    id::Int
    source::Int
    target::Int
    switching::SW
end

const Scope = Union{ModeScope, TransitionScope}

"""
    ScopedSet(inner, scope)

A scalar set tagged with the mode or transition it was written on.
"""
struct ScopedSet{S <: MOI.AbstractScalarSet, C <: Scope} <: MOI.AbstractScalarSet
    inner::S
    scope::C
end
Base.copy(set::ScopedSet) = ScopedSet(copy(set.inner), set.scope)

"""
    ScopedVectorSet(inner, scope)

A vector set — a specification marker or an obstacle — tagged with its mode or transition.
"""
struct ScopedVectorSet{S <: MOI.AbstractVectorSet, C <: Scope} <: MOI.AbstractVectorSet
    inner::S
    scope::C
end
MOI.dimension(set::ScopedVectorSet) = MOI.dimension(set.inner)
Base.copy(set::ScopedVectorSet) = ScopedVectorSet(copy(set.inner), set.scope)

# ----------------------------------------------------------------------------------------
# The scopes themselves
# ----------------------------------------------------------------------------------------

"""
    Mode

One discrete mode of a hybrid model: a scope carrying its own dynamics, its own state and
input bounds, and its own specifications. Created with [`@mode`](@ref) or [`add_mode!`](@ref).
"""
struct Mode <: JuMP.AbstractModel
    model::JuMP.GenericModel
    id::Int
    name::Symbol
end

"""
    Transition

A switch from one mode to another: a scope carrying the **guard** (which states enable it) and
the **reset map** (where the state lands). Created with [`add_transition!`](@ref).
"""
mutable struct Transition{SW} <: JuMP.AbstractModel
    model::JuMP.GenericModel
    id::Int
    source::Int
    target::Int
    switching::SW
    nconstraints::Int
end

scope(m::Mode) = ModeScope(m.id)
scope(t::Transition) = TransitionScope(t.id, t.source, t.target, t.switching)

"""
    set_attribute(mode::Mode, name::String, value)

Set a solver option for this mode alone — each mode of a hybrid model is abstracted by its own
sub-solver, so `state_grid`, `input_grid` and `time_step` are per-mode:

```julia
set_attribute(off, "state_grid", MP.GridFree(SVector(0.0), SVector(0.1)))
```

Options set on the model itself apply to every mode unless a mode overrides them.
"""
function JuMP.set_attribute(m::Mode, name::String, value)
    return JuMP.set_attribute(m.model, "mode[$(m.id)].$(name)", value)
end

# `model_convert` broadcasts over the model for vector constraints, and JuMP only defines
# `broadcastable` for `GenericModel`.
Base.broadcastable(m::Mode) = Ref(m)
Base.broadcastable(t::Transition) = Ref(t)

# A scope subtypes `JuMP.AbstractModel` so that `@constraint` works on it, which also makes JuMP's
# generic `show` apply — and that one summarises a *model*, asking for `objective_sense`,
# `num_variables` and friends. A scope has none of those, so it prints itself.
#
# `print` needs its own method: JuMP defines `Base.print(::IO, ::AbstractModel)`, which is more
# specific than the generic `Base.print` that would otherwise fall back to `show`. Without it
# `println(mode)` and `"$mode"` still reach JuMP's model summary and throw.
Base.show(io::IO, m::Mode) = print(io, "Mode(:", m.name, ")")
Base.print(io::IO, m::Mode) = show(io, m)

# Constraint display. Without these, JuMP prints a scoped constraint as
# `… in Dionysos.Wrapper.ScopedSet{MathOptInterface.EqualTo{Float64}, …}(…)`, putting the
# wrapper's internal types in front of the reader; the scope is what they actually want to see.
_scope_suffix(s::ModeScope) = "  [mode $(s.id)]"
_scope_suffix(s::TransitionScope) = "  [transition $(s.source) → $(s.target)]"

JuMP.in_set_string(mode::MIME, set::ScopedSet) =
    JuMP.in_set_string(mode, set.inner) * _scope_suffix(set.scope)

JuMP.in_set_string(mode::MIME, set::ScopedVectorSet) =
    JuMP.in_set_string(mode, set.inner) * _scope_suffix(set.scope)

function Base.show(io::IO, t::Transition)
    plural = t.nconstraints == 1 ? "" : "s"
    return print(
        io,
        "Transition(mode ",
        t.source,
        " => mode ",
        t.target,
        ", ",
        t.nconstraints,
        " constraint",
        plural,
        ")",
    )
end

Base.print(io::IO, t::Transition) = show(io, t)

_record_constraint!(::Mode) = nothing
_record_constraint!(t::Transition) = (t.nconstraints += 1; nothing)

function JuMP.add_constraint(
    s::Union{Mode, Transition},
    con::JuMP.ScalarConstraint,
    name::String = "",
)
    _record_constraint!(s)
    scoped = JuMP.ScalarConstraint(con.func, ScopedSet(con.set, scope(s)))
    return JuMP.add_constraint(s.model, scoped, name)
end

function JuMP.add_constraint(
    s::Union{Mode, Transition},
    con::JuMP.VectorConstraint,
    name::String = "",
)
    _record_constraint!(s)
    scoped = JuMP.VectorConstraint(con.func, ScopedVectorSet(con.set, scope(s)), con.shape)
    return JuMP.add_constraint(s.model, scoped, name)
end

# ----------------------------------------------------------------------------------------
# Building them
# ----------------------------------------------------------------------------------------

# Scope ids are allocated in the JuMP model's `ext` dictionary rather than in the optimizer,
# because JuMP empties the optimizer before copying the model into it.
function _next_scope_id!(model::JuMP.GenericModel, key::Symbol)
    id = get(model.ext, key, 0) + 1
    model.ext[key] = id
    return id
end

"""
    add_mode!(model, name = :mode) -> Mode

Declare a discrete mode. Write its dynamics, bounds and specifications on the returned scope:

```julia
off = add_mode!(model, :off)
@constraint(off, ∂(T) == -α * (T - Ta))
@constraint(off, u == 0)
```

[`@mode`](@ref) is the same thing with the name taken from the variable.
"""
function add_mode!(model::JuMP.GenericModel, name::Symbol = :mode)
    return Mode(model, _next_scope_id!(model, :_dionysos_modes), name)
end

"""
    @mode(model, name)

Declare a mode. Like `@variable`, this **binds `name` in the calling scope** and registers the
mode in the model's object dictionary — no assignment needed:

```julia
@mode(model, off)
@constraint(off, ∂(T) == -α * (T - Ta))     # `off` is already bound
model[:off] === off                          # and registered
```

The macro also returns the mode, so `off = @mode(model, off)` works too; it is simply
redundant. [`add_mode!`](@ref) is the function form, for a name computed at runtime.
"""
macro mode(model, name)
    quote
        local mode = add_mode!($(esc(model)), $(QuoteNode(name)))
        $(esc(model))[$(QuoteNode(name))] = mode
        $(esc(name)) = mode
    end
end

"""
    add_transition!(f, model, source => target; switching = AutonomousSwitching())

Declare a switch from `source` to `target` and populate it with a do-block:

```julia
add_transition!(model, off => on) do t
    @constraint(t, T <= 19)        # guard: which states enable the switch
    @constraint(t, Δ(T) == T)      # reset map (identity if omitted)
end
```

A constraint containing `∂`/`Δ` is the **reset map**; anything else is a **guard**. Several
transitions may share the same `(source, target)` pair — each is its own object, so their
guards and resets never get mixed up.

A transition needs at least a guard; an always-enabled switch is written with a guard covering
the whole state set.

**A switch is an input the synthesis chooses.** The guard says *where* the switch is available,
not *that* it is taken: from a state inside the guard the abstraction offers both the switch and
the mode's own dynamics, and the controller picks. So a guard is a permission, and "stay put"
remains available everywhere.

Environment-forced switching — a fault, an impact, a threshold the physics obliges — is therefore
**not** expressible today: a controller synthesized here assumes it may decline any switch. Only
`HybridSystems.AutonomousSwitching()` is accepted, as the default placeholder value; that name is
`HybridSystems`' and does not describe the semantics implemented here.
"""
function add_transition!(
    f::Function,
    model::JuMP.GenericModel,
    pair::Pair{Mode, Mode};
    switching = HybridSystems.AutonomousSwitching(),
)
    # The value is carried into `HybridSystems` and never read again: the abstraction gives every
    # switch its own input symbol, so what it builds is a controller-chosen switch whatever this
    # says. Refusing the other value keeps the model from claiming a semantics nothing enforces.
    switching isa HybridSystems.AutonomousSwitching || error(
        "`switching = $(typeof(switching))` is not supported; pass " *
        "`HybridSystems.AutonomousSwitching()` (the default) or omit it. Note that the name is " *
        "`HybridSystems`': the abstraction gives each switch its own input symbol, so the " *
        "switch is one the synthesis *chooses* from inside the guard. Environment-forced " *
        "switching is not implemented.",
    )
    source, target = pair
    id = _next_scope_id!(model, :_dionysos_transitions)
    transition = Transition(model, id, source.id, target.id, switching, 0)
    f(transition)
    transition.nconstraints > 0 || error(
        "The transition $(source.name) → $(target.name) has no guard: nothing tells the " *
        "solver when the switch may be taken. Add one, e.g. `@constraint(t, x <= 1.0)`, or " *
        "a guard covering the whole state set for an always-enabled switch.",
    )
    return transition
end
