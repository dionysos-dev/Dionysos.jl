# ----------------------------------------------------------------------------------------
# Specifications.
#
# A specification set is written `@constraint(model, x in Final(S))`, wrapping the *set*
# rather than the variables. That is forced by JuMP: `final` is a `NonlinearOperator`, which
# only builds an expression when one of its arguments is a JuMP scalar, so `final(x)` on a
# whole vector calls the underlying function and raises a `MethodError`. Wrapping the set is
# also what lets `S` be any bounded `LazySet` instead of a box built coordinate by coordinate.
# ----------------------------------------------------------------------------------------

export Start, Final, Always, EventuallyAlways, Label, @specification

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

# Keep the marker, not its type signature, in front of the reader when a constraint is printed.
_marker_name(::Start) = "Start"
_marker_name(::Final) = "Final"
_marker_name(::Always) = "Always"
_marker_name(::EventuallyAlways) = "EventuallyAlways"

JuMP.in_set_string(
    ::MIME,
    set::SpecSet,
) = "in $(_marker_name(set))($(nameof(typeof(set.inner))))"

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

# ----------------------------------------------------------------------------------------
# The general layer: a temporal formula over named regions
# ----------------------------------------------------------------------------------------

"""
    Label(S; semantics = Mapping.INNER)

Name a region of the state space so a temporal formula can refer to it. The name is the
**constraint's own name**, so no separate registration call is needed:

```julia
@constraint(model, goal,   x in Label(target))
@constraint(model, hazard, x in Label(obstacle; semantics = MP.OUTER))
@specification(model, ltl"F(goal) & G(!hazard)")
```

`semantics` says how the region is discretized: `INNER` keeps only cells fully inside it — the
conservative reading for something you must reach — while `OUTER` keeps every cell touching it,
the conservative reading for something you must avoid.
"""
struct Label{S} <: MOI.AbstractVectorSet
    inner::S
    dim::Int
    semantics::UT.INCL_MODE
end

Label(set; semantics::UT.INCL_MODE = UT.INNER) = Label(set, LazySets.dim(set), semantics)

MOI.dimension(set::Label) = set.dim
Base.copy(set::Label) = set

JuMP.in_set_string(
    ::MIME,
    set::Label,
) = "in Label($(nameof(typeof(set.inner))), $(set.semantics))"

"""
    LabelEntry

One named region: the atomic proposition it defines, the set, how it is discretized, and the
variables it is written over. `name` is filled in when JuMP forwards the constraint's name.
"""
mutable struct LabelEntry
    name::String
    set::Any
    semantics::UT.INCL_MODE
    variables::Vector{MOI.VariableIndex}
end

"""
    @specification(model, formula)

Attach the temporal formula the trajectory must satisfy, over the regions named with
[`Label`](@ref):

```julia
@specification(model, ltl"F(goal) & G(!hazard)")
```

`formula` may be a `Spot.SpotFormula` (with `using Spot`) or any
`Optim.DiscreteSystems.AbstractSpecStepper` — a hand-written monitor, for instance. A model
carrying a formula lowers to a
[`CoSafeLTLProblem`](@ref Dionysos.Problem.CoSafeLTLProblem).
"""
macro specification(model, formula)
    return quote
        JuMP.set_attribute($(esc(model)), "specification", $(esc(formula)))
    end
end

"""
    to_stepper(specification)

Turn what the user attached with [`@specification`](@ref) into the automaton the co-safe solver
steps. A stepper passes through unchanged; `DionysosSpotExt` adds the method that compiles a
`Spot.SpotFormula`.
"""
to_stepper(spec::OPDS.AbstractSpecStepper) = spec

function to_stepper(spec)
    return error(
        "Cannot turn a $(typeof(spec)) into a specification monitor. Pass a Spot formula " *
        "(`using Spot`, then `ltl\"F(goal)\"`), or an " *
        "`Optim.DiscreteSystems.AbstractSpecStepper` such as a `FunctionMonitor`.",
    )
end
