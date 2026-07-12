"""
    ClockLift{C} <: AbstractLift

The time lift: augments a base abstraction with a monotone clock
[`ClockAbstraction`](@ref), turning a model of `x` into a model of `(x, t)`.
The same lift applies to a plain continuous model or (per mode) to a hybrid mode,
so it is the single reusable mechanism for time-dependent abstractions.
"""
struct ClockLift{C <: ClockAbstraction} <: AbstractLift
    clock::C
end

"""
    ClockLiftedSymbolicModel{B, C, A} <: AbstractSymbolicModel

Result of `lift(::ClockLift, base)`: the synchronized product of a base model with
a clock. Abstract states are the flattened `(base_state_id, time_id)` pairs; the
concrete coordinate is the augmented vector `[x; t]` (time is the trailing
dimension). For an active clock a base transition `q → q′` is replicated as
`(q, p) → (q′, p+1)` across time steps; a frozen clock keeps `p` fixed.

Inputs pass through to `base` unchanged (the clock adds no inputs).
"""
struct ClockLiftedSymbolicModel{B <: AbstractSymbolicModel, C <: ClockAbstraction, A} <:
       AbstractSymbolicModel
    base::B
    clock::C
    automaton::A
    flat::FlatIndex{Tuple{Int, Int}}
end

function lift(l::ClockLift, base::AbstractSymbolicModel)
    clock = l.clock
    ntime = length(clock.tsteps)
    base_states = collect(enum_states(base))

    # Flatten the (base_state, time_id) product, time-major so that a mode's states
    # are grouped by time slice.
    keys = Vector{Tuple{Int, Int}}(undef, length(base_states) * ntime)
    idx = 0
    for p in 1:ntime
        for q in base_states
            idx += 1
            keys[idx] = (q, p)
        end
    end
    flat = FlatIndex(keys)

    automaton = IndexedAutomatonList(n_flat(flat), get_n_input(base))
    for (target, source, input) in enum_transitions(base)
        if clock.is_active
            @inbounds for p in 1:(ntime - 1)
                s = flat_id(flat, (source, p))
                t = flat_id(flat, (target, p + 1))
                add_transition!(automaton, s, t, input)
            end
        else
            s = flat_id(flat, (source, 1))
            t = flat_id(flat, (target, 1))
            add_transition!(automaton, s, t, input)
        end
    end

    return ClockLiftedSymbolicModel(base, clock, automaton, flat)
end

# ---- AbstractSymbolicModel surface (states/automaton) ----

get_automaton(m::ClockLiftedSymbolicModel) = m.automaton
get_n_state(m::ClockLiftedSymbolicModel) = n_flat(m.flat)
enum_states(m::ClockLiftedSymbolicModel) = 1:get_n_state(m)
get_state_dim(m::ClockLiftedSymbolicModel) = get_state_dim(m.base) + 1

# The lift owns no transition metadata (the automaton is rebuilt from the base's),
# so a "light" copy for controller saving is the model itself.
without_metadata(m::ClockLiftedSymbolicModel) = m

# ---- Inputs pass through to the base ----

get_n_input(m::ClockLiftedSymbolicModel) = get_n_input(m.base)
enum_inputs(m::ClockLiftedSymbolicModel) = enum_inputs(m.base)
get_concrete_input(m::ClockLiftedSymbolicModel, u::Int) = get_concrete_input(m.base, u)
get_abstract_input(m::ClockLiftedSymbolicModel, u) = get_abstract_input(m.base, u)

# ---- Concrete ↔ abstract over the augmented coordinate [x; t] ----

"Concretize a flattened id to the augmented coordinate `[x; t]`."
function get_concrete_state(m::ClockLiftedSymbolicModel, id::Int)
    (q, p) = flat_key(m.flat, id)
    x = get_concrete_state(m.base, q)
    return vcat(x, SVector(int2time(m.clock, p)))
end

"""
    base_state_and_time(m::ClockLiftedSymbolicModel, id) -> (x, t)

Decompose a flattened id into the base concrete state `x` (unaugmented, keeping the
base's coordinate type) and the clock value `t`. Used by the mode composition to
build a hybrid coordinate without re-splitting the `[x; t]` vector.
"""
function base_state_and_time(m::ClockLiftedSymbolicModel, id::Int)
    (q, p) = flat_key(m.flat, id)
    return get_concrete_state(m.base, q), int2time(m.clock, p)
end

"""
    abstract_switch_target(mode_model, reset_coord) -> Int

Local abstract state of a guarded-switch target in `mode_model`, from the concrete
coordinate produced by the reset map (`0` if out of range). For a clock-lifted mode
the coordinate is `[x; t]` and the time is rounded *up* (the smallest clock step
`≥ t`); for a time-free mode it is just `x`.
"""
function abstract_switch_target(m::ClockLiftedSymbolicModel, xt)
    x = xt[1:(end - 1)]
    t = xt[end]
    q = get_abstract_state(m.base, x)
    (q === nothing || q <= 0) && return 0
    p = ceil_time2int(m.clock, t)
    (p <= 0 || p > length(m.clock.tsteps)) && return 0
    return flat_id(m.flat, (q, p))
end

function abstract_switch_target(m::SymbolicModel, x)
    q = get_abstract_state(m, x)
    return (q === nothing || q <= 0) ? 0 : q
end

"Abstract an augmented coordinate `[x; t]` to its flattened id (`0` if absent)."
function get_abstract_state(m::ClockLiftedSymbolicModel, xt)
    x = xt[1:(end - 1)]
    t = xt[end]
    q = get_abstract_state(m.base, x)
    (q === nothing || q <= 0) && return 0
    p = floor_time2int(m.clock, t)
    p <= 0 && return 0
    return flat_id(m.flat, (q, p))
end

"""
    get_states_from_set(m::ClockLiftedSymbolicModel, set, incl_mode)

Flattened states whose spatial part lies in the projection of `set` onto the base
coordinates and whose time index lies in the projection onto the trailing (time)
dimension. `set` is a box over `[x; t]`.
"""
function get_states_from_set(m::ClockLiftedSymbolicModel, set, incl_mode::MP.INCL_MODE)
    base_set = UT.box(LazySets.low(set)[1:(end - 1)], LazySets.high(set)[1:(end - 1)])
    d = LazySets.dim(set)
    tmin = LazySets.low(set, d)
    tmax = LazySets.high(set, d)

    base_states = get_states_from_set(m.base, base_set, incl_mode)
    pmin = ceil_time2int(m.clock, tmin)
    pmax = floor_time2int(m.clock, tmax)

    states = Int[]
    for q in base_states
        for p in pmin:pmax
            id = flat_id(m.flat, (q, p))
            id > 0 && push!(states, id)
        end
    end
    return states
end
