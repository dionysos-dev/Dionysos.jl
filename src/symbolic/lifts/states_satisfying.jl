# Evaluate a lifted specification (`Problem.AbstractSpecification`) against a
# (possibly lifted) symbolic model, returning the abstract states that satisfy it.
# The spec's lift stack is walked in lockstep with the model's: a base model meets
# a `StateSpec`, a clock-lifted model meets a `TimedSpec`, a hybrid model meets a
# `HybridSpec`.

"Abstract states of a base model satisfying a base state spec."
states_satisfying(m::SymbolicModel, s::PR.StateSpec) =
    get_states_from_set(m, s.set, s.incl_mode)

"Clock-lifted model + base spec: matching base states at every time index."
function states_satisfying(m::ClockLiftedSymbolicModel, s::PR.StateSpec)
    base_states = states_satisfying(m.base, s)
    return _at_time_indices(m, base_states, 1:length(m.clock.tsteps))
end

"Clock-lifted model + timed spec: matching base states within the time window."
function states_satisfying(m::ClockLiftedSymbolicModel, s::PR.TimedSpec)
    base_states = states_satisfying(m.base, s.base)
    pmin = ceil_time2int(m.clock, s.tmin)
    pmax = floor_time2int(m.clock, s.tmax)
    return _at_time_indices(m, base_states, pmin:pmax)
end

function _at_time_indices(m::ClockLiftedSymbolicModel, base_states, prange)
    result = Int[]
    for q in base_states
        for p in prange
            id = flat_id(m.flat, (q, p))
            id > 0 && push!(result, id)
        end
    end
    return result
end

"Hybrid model + mode spec: per-mode local states mapped to global ids."
function states_satisfying(m::HybridSymbolicModel, s::PR.HybridSpec)
    result = Int[]
    for (mode_id, sub) in s.per_mode
        (1 <= mode_id <= length(m.mode_models)) || continue
        for local_id in states_satisfying(m.mode_models[mode_id], sub)
            gid = flat_id(m.flat, (local_id, mode_id))
            gid > 0 && push!(result, gid)
        end
    end
    return result
end
