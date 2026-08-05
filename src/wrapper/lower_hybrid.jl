# ----------------------------------------------------------------------------------------
# Lowering a hybrid model: `ModelIR` with modes → `HybridSystems.HybridSystem` + a problem
# over augmented states.
#
# Every mode is a plain physical system here, so the augmented state is `(x, mode)`. Clocks —
# and with them `(x, t, mode)` — are a later phase.
# ----------------------------------------------------------------------------------------

# The bounds of `i` inside mode `m`: the model-level declaration, narrowed by any per-mode
# constraint.
function _mode_bound(ir::ModelIR, m::ModeIR, i::Int)
    v = ir.variables[i]
    return get(m.lower, i, v.lower), get(m.upper, i, v.upper)
end

function _mode_box(ir::ModelIR, m::ModeIR, idx::Vector{Int})
    lo = Float64[]
    hi = Float64[]
    for i in idx
        l, h = _mode_bound(ir, m, i)
        push!(lo, l)
        push!(hi, h)
    end
    return UT.box(SVector{length(idx)}(lo...), SVector{length(idx)}(hi...))
end

# The dynamics of one mode as the dense vector the backend expects: the mode's own equation
# for each state, falling back to the model-level one.
function _mode_dynamics(ir::ModelIR, m::ModeIR)
    dynamics =
        Vector{Union{Nothing, MOI.ScalarNonlinearFunction}}(nothing, length(ir.variables))
    for i in eachindex(ir.variables)
        dynamics[i] = get(m.dynamics, i, ir.dynamics[i])
    end
    return dynamics
end

function _build_mode_system(ir::ModelIR, m::ModeIR, backend, x_idx, u_idx)
    dynamics = _mode_dynamics(ir, m)
    missing_state = findfirst(i -> dynamics[i] === nothing, x_idx)
    missing_state === nothing || error(
        "Mode $(m.id) has no dynamics for $(describe(ir.variables[x_idx[missing_state]], x_idx[missing_state])). " *
        "Give every state an equation in every mode, or declare it at the model level.",
    )

    f = compile_dynamics(backend, ir, dynamics)
    X = _mode_box(ir, m, x_idx)
    U = _mode_box(ir, m, u_idx)

    if ir.time_domain == CONTINUOUS
        return MS.ConstrainedBlackBoxControlContinuousSystem(
            f,
            LazySets.dim(X),
            LazySets.dim(U),
            X,
            U,
        )
    else
        return MS.ConstrainedBlackBoxControlDiscreteSystem(
            (x, u) -> f(x, u),
            LazySets.dim(X),
            LazySets.dim(U),
            X,
            U,
        )
    end
end

# A guard is a box over the state coordinates: the source mode's state set, narrowed by the
# bounds written on the transition. An unconstrained coordinate stays free within the mode.
function _guard_box(ir::ModelIR, t::TransitionIR, x_idx::Vector{Int})
    source = ir.modes[t.source]
    lo = Float64[]
    hi = Float64[]
    for i in x_idx
        l, h = _mode_bound(ir, source, i)
        push!(lo, max(l, get(t.guard_lower, i, -Inf)))
        push!(hi, min(h, get(t.guard_upper, i, Inf)))
    end
    return UT.box(SVector{length(x_idx)}(lo...), SVector{length(x_idx)}(hi...))
end

# The reset applied when the switch is taken. With no `Δ(x) == …` on the transition this is
# the identity, which is what switching mode usually does to the state.
function _reset_function(ir::ModelIR, t::TransitionIR, backend, x_idx::Vector{Int})
    isempty(t.resets) && return identity

    dynamics =
        Vector{Union{Nothing, MOI.ScalarNonlinearFunction}}(nothing, length(ir.variables))
    for i in x_idx
        # A state with no reset keeps its value; `x + 0` expresses that as an expression the
        # dynamics backend can compile alongside the others.
        dynamics[i] = get(t.resets, i) do
            return MOI.ScalarNonlinearFunction(:+, Any[MOI.VariableIndex(i), 0.0])
        end
    end

    reset = compile_dynamics(backend, ir, dynamics)
    # The compiled function takes `(x, u)`; a reset map is applied to the state alone.
    u_idx = input_indices(ir)
    u0 = SVector{length(u_idx)}(zeros(length(u_idx))...)
    return state -> reset(state, u0)
end

"""
    build_hybrid_system(ir, backend) -> HybridSystems.HybridSystem

Assemble the hybrid automaton: one system per mode, and one
[`GuardedResetMap`](@ref Dionysos.System.GuardedResetMap) per transition.
"""
function build_hybrid_system(ir::ModelIR, backend)
    ids = mode_ids(ir)
    x_idx = state_indices(ir)
    u_idx = input_indices(ir)

    mode_systems = [_build_mode_system(ir, ir.modes[k], backend, x_idx, u_idx) for k in ids]

    transitions = [ir.transitions[k] for k in sort!(collect(keys(ir.transitions)))]
    automaton = HybridSystems.GraphAutomaton(length(ids))
    reset_maps = ST.GuardedResetMap[]
    switchings = Any[]
    for (i, t) in enumerate(transitions)
        HybridSystems.add_transition!(automaton, t.source, t.target, i)
        push!(
            reset_maps,
            ST.GuardedResetMap(
                _guard_box(ir, t, x_idx),
                _reset_function(ir, t, backend, x_idx),
            ),
        )
        push!(switchings, t.switching)
    end

    return HybridSystems.HybridSystem(automaton, mode_systems, reset_maps, switchings)
end

# ---- specifications over modes ----

# The specification set of one mode and one kind, as a box over the state coordinates: either a
# marker set given over the whole state, or per-coordinate `final(x[i])` intervals.
function _mode_spec_set(ir::ModelIR, m::ModeIR, kind::SpecKind, x_idx::Vector{Int})
    entries = filter(e -> e.kind === kind, m.specs)
    isempty(entries) && return nothing

    whole = filter(e -> length(e.variables) == length(x_idx), entries)
    if !isempty(whole)
        length(whole) == 1 || error("Mode $(m.id) has several `$kind` sets; keep one.")
        return whole[].set
    end

    # Per-coordinate intervals, defaulting to the mode's own bounds.
    lo = Float64[]
    hi = Float64[]
    for i in x_idx
        l, h = _mode_bound(ir, m, i)
        entry = findfirst(e -> e.variables == [MOI.VariableIndex(i)], entries)
        if entry !== nothing
            interval = entries[entry].set
            isfinite(interval.lower) && (l = interval.lower)
            isfinite(interval.upper) && (h = interval.upper)
        end
        push!(lo, l)
        push!(hi, h)
    end
    return UT.box(SVector{length(x_idx)}(lo...), SVector{length(x_idx)}(hi...))
end

function _hybrid_spec(ir::ModelIR, kind::SpecKind, x_idx::Vector{Int})
    pairs = Pair{Int, PR.StateSpec}[]
    for k in mode_ids(ir)
        set = _mode_spec_set(ir, ir.modes[k], kind, x_idx)
        set === nothing || push!(pairs, k => PR.StateSpec(set))
    end
    isempty(pairs) && return nothing
    return PR.HybridSpec(Dict(pairs))
end

# The augmented initial state `(x, mode)`. The continuous part is the centre of the declared
# initial box — a hybrid problem starts from a point, not a region — and the mode is the one
# carrying a `start` constraint, defaulting to the first.
function _hybrid_initial_state(ir::ModelIR, x_idx::Vector{Int})
    box = _coordinate_box(ir, x_idx, :start)
    x0 = SVector{length(x_idx)}(LazySets.center(box)...)

    starting = [k for k in mode_ids(ir) if any(e -> e.kind === START, ir.modes[k].specs)]
    length(starting) <= 1 || error(
        "Several modes carry a `start` constraint ($(starting)); the model can only begin in one.",
    )
    return (x0, isempty(starting) ? first(mode_ids(ir)) : starting[])
end

"""
    build_hybrid_problem(ir, backend; time_step = nothing) -> Problem.ProblemType

Assemble the hybrid control problem: a `HybridSystem`, an augmented initial state `(x, mode)`,
and a mode-indexed specification.
"""
function build_hybrid_problem(ir::ModelIR, backend; time_step = nothing)
    x_idx = state_indices(ir)
    system = build_hybrid_system(ir, backend)
    initial_state = _hybrid_initial_state(ir, x_idx)

    target = _hybrid_spec(ir, FINAL, x_idx)
    safe = _hybrid_spec(ir, ALWAYS, x_idx)

    if target !== nothing
        return PR.OptimalControlProblem(
            system,
            initial_state,
            target,
            nothing,
            nothing,
            _horizon(ir, false, time_step),
        )
    elseif safe !== nothing
        return PR.SafetyProblem(system, initial_state, safe, _horizon(ir, true, time_step))
    end
    return error(
        "A hybrid model needs a specification on at least one mode: add a `Final` or `Always` " *
        "set, or a `final(x[i]) in …` constraint, to the mode you care about.",
    )
end
