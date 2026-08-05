# ----------------------------------------------------------------------------------------
# Specifications.
#
# A specification set is written `@constraint(model, x in Final(S))`, wrapping the *set*
# rather than the variables. That is forced by JuMP: `final` is a `NonlinearOperator`, which
# only builds an expression when one of its arguments is a JuMP scalar, so `final(x)` on a
# whole vector calls the underlying function and raises a `MethodError`. Wrapping the set is
# also what lets `S` be any bounded `LazySet` instead of a box built coordinate by coordinate.
# ----------------------------------------------------------------------------------------

export Start, Final, Always, EventuallyAlways

"""
    SpecKind

Which temporal role a specification set plays: the initial set (`START`), a reach target
(`FINAL`), an invariant (`ALWAYS`), or a reach-and-stay target (`EVENTUALLY_ALWAYS`).
"""
@enum SpecKind START FINAL ALWAYS EVENTUALLY_ALWAYS

"""
    Start(S)

Marks `S` as the initial set: `@constraint(model, x in Start(S))`.

`S` may be any bounded `LazySet` and must span the whole state vector. The
`@variable(model, x, start = v)` keyword is the singleton special case.
"""
struct Start{S} <: MOI.AbstractVectorSet
    inner::S
    dim::Int
end

"""
    Final(S)

Marks `S` as a reach target — ◇, the trajectory must eventually enter `S`:
`@constraint(model, x in Final(S))`. Lowers to an
[`OptimalControlProblem`](@ref Dionysos.Problem.OptimalControlProblem).
"""
struct Final{S} <: MOI.AbstractVectorSet
    inner::S
    dim::Int
end

"""
    Always(S)

Marks `S` as an invariant — □, the trajectory must never leave `S`:
`@constraint(model, x in Always(S))`. On its own it lowers to a
[`SafetyProblem`](@ref Dionysos.Problem.SafetyProblem); together with a [`Final`](@ref) set it
becomes a reach-avoid problem, where it is folded into the state space.

This is *not* the same as `x ∉ O`: `∉` removes a region from the state space so it is never
abstracted, whereas `Always` keeps it representable so synthesis can actively avoid it.
"""
struct Always{S} <: MOI.AbstractVectorSet
    inner::S
    dim::Int
end

"""
    EventuallyAlways(S)

Marks `S` as a reach-and-stay target — ◇□, the trajectory must reach `S` and remain there:
`@constraint(model, x in EventuallyAlways(S))`. Lowers to a
[`ReachAndStayProblem`](@ref Dionysos.Problem.ReachAndStayProblem).
"""
struct EventuallyAlways{S} <: MOI.AbstractVectorSet
    inner::S
    dim::Int
end

"""
    SpecSet

Any of the specification markers [`Start`](@ref), [`Final`](@ref), [`Always`](@ref),
[`EventuallyAlways`](@ref).
"""
const SpecSet = Union{Start, Final, Always, EventuallyAlways}

Start(set) = Start(set, LazySets.dim(set))
Final(set) = Final(set, LazySets.dim(set))
Always(set) = Always(set, LazySets.dim(set))
EventuallyAlways(set) = EventuallyAlways(set, LazySets.dim(set))

spec_kind(::Start) = START
spec_kind(::Final) = FINAL
spec_kind(::Always) = ALWAYS
spec_kind(::EventuallyAlways) = EVENTUALLY_ALWAYS

MOI.dimension(set::SpecSet) = set.dim
Base.copy(set::SpecSet) = set

"""
    SpecEntry

One specification set as written by the user: its temporal role, the set itself, and the
variables it constrains.
"""
struct SpecEntry{S}
    kind::SpecKind
    set::S
    variables::Vector{MOI.VariableIndex}
end
