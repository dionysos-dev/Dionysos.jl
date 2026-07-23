# ----------------------------
# Lifted specifications
# ----------------------------
#
# A specification is composed by the same optional lifts as the abstraction it is
# evaluated against (see `Dionysos.Symbolic.AbstractLift`): a base state
# constraint over `x`, optionally wrapped by a time window (the clock lift) and/or
# a per-mode dispatch (the mode lift). Concrete-point membership is `satisfies`;
# the abstract-state evaluation lives in `Symbolic` as `states_satisfying`.

"""
    AbstractSpecification

Root of the composable specification hierarchy. A spec constrains the (possibly
augmented) state `x` / `(x,t)` / `(x,t,k)` and is built by wrapping a base
[`StateSpec`](@ref) with the optional lifts [`TimedSpec`](@ref) (time) and
[`HybridSpec`](@ref) (mode), mirroring the abstraction it is evaluated against.
"""
abstract type AbstractSpecification end

"""
    StateSpec{S} <: AbstractSpecification

Base specification: the state `x` must lie in `set` under inclusion mode
`incl_mode`.
"""
struct StateSpec{S} <: AbstractSpecification
    set::S
    incl_mode::UT.INCL_MODE
end
StateSpec(set) = StateSpec(set, UT.INNER)

"""
    TimedSpec{B} <: AbstractSpecification

Clock lift of a spec: the base spec `base` must hold *and* the time value must lie
in `[tmin, tmax]`.
"""
struct TimedSpec{B <: AbstractSpecification} <: AbstractSpecification
    base::B
    tmin::Float64
    tmax::Float64
end

"""
    HybridSpec{S} <: AbstractSpecification

Mode lift of a spec: `per_mode[k]` is the spec that must hold in mode `k`; modes
absent from the mapping are not part of the specification.
"""
struct HybridSpec{S <: AbstractSpecification} <: AbstractSpecification
    per_mode::Dict{Int, S}
end

"""
    hybrid_reach_spec(state_sets, time_sets, mode_ids; incl_mode = UT.INNER) -> HybridSpec

Build a mode-indexed timed specification from parallel `(state_sets, time_sets,
mode_ids)`: mode `mode_ids[i]` requires `x ∈ state_sets[i]` within the time window
of `time_sets[i]`.
"""
function hybrid_reach_spec(
    state_sets,
    time_sets,
    mode_ids;
    incl_mode::UT.INCL_MODE = UT.INNER,
)
    @assert length(state_sets) == length(time_sets) == length(mode_ids) "state/time/mode lengths must match"
    @assert allunique(mode_ids) "mode_ids must be unique (one spec per mode); got $mode_ids"
    pairs = map(eachindex(mode_ids)) do i
        base = StateSpec(state_sets[i], incl_mode)
        timed = TimedSpec(
            base,
            LazySets.low(time_sets[i], 1),
            LazySets.high(time_sets[i], 1),
        )
        return mode_ids[i] => timed
    end
    return HybridSpec(Dict(pairs))
end

# ---- Concrete-point membership ----

"Whether the concrete state `x` satisfies the base spec (time-agnostic)."
satisfies(s::StateSpec, x) = x ∈ s.set
satisfies(s::StateSpec, x, _t) = x ∈ s.set

"Whether `(x, t)` satisfies the timed spec."
satisfies(s::TimedSpec, x, t) = (s.tmin <= t <= s.tmax) && satisfies(s.base, x)

"Whether the augmented `(x, t, k)` satisfies the mode-indexed spec (clock-lifted mode)."
function satisfies(s::HybridSpec, x, t, k)
    sub = get(s.per_mode, k, nothing)
    sub === nothing && return false
    return satisfies(sub, x, t)
end

"Whether the augmented `(x, k)` satisfies the mode-indexed spec (time-free mode)."
function satisfies(s::HybridSpec, x, k::Integer)
    sub = get(s.per_mode, k, nothing)
    sub === nothing && return false
    return satisfies(sub, x)
end
