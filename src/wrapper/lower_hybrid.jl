# ----------------------------------------------------------------------------------------
# Lowering a hybrid model: `ModelIR` with modes → `HybridSystems.HybridSystem` + a problem
# over augmented states.
#
# The augmented state is `(x, mode)`, or `(x, t, mode)` when the model declares a clock (see
# `clock.jl`). The clock is what decides the arity of guards and reset maps, so it is consulted
# throughout rather than handled as a special case at the end.
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
    # The element type is spelled out because `idx` may be empty: a mode whose only control is
    # the switch itself has no continuous input, and `SVector{0}()` would infer `Union{}`.
    return LazySets.Hyperrectangle(;
        low = SVector{length(idx), Float64}(lo...),
        high = SVector{length(idx), Float64}(hi...),
    )
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
    supplied = m.user_dynamics === nothing ? ir.user_dynamics : m.user_dynamics

    if supplied === nothing
        missing_state = findfirst(i -> dynamics[i] === nothing, x_idx)
        missing_state === nothing || error(
            "Mode $(m.id) has no dynamics for $(describe(ir.variables[x_idx[missing_state]], x_idx[missing_state])). " *
            "Give every state an equation in every mode, or declare it at the model level.",
        )
    end

    f = dynamics_function(ir, backend, dynamics, supplied)
    X = _mode_box(ir, m, x_idx)
    U = _mode_box(ir, m, u_idx)

    physical = if ir.time_domain == CONTINUOUS
        MS.ConstrainedBlackBoxControlContinuousSystem(f, LazySets.dim(X), LazySets.dim(U), X, U)
    else
        MS.ConstrainedBlackBoxControlDiscreteSystem(
            (x, u) -> f(x, u),
            LazySets.dim(X),
            LazySets.dim(U),
            X,
            U,
        )
    end

    # A clock turns the mode into the physical system paired with its time axis, which is what
    # makes the augmented state `(x, t, mode)` instead of `(x, mode)`.
    clock = clock_index(ir)
    clock === nothing && return physical
    return ST.VectorContinuousSystem([physical, clock_system(ir, m, clock)])
end

# A guard is a box over the *augmented* state — `x`, plus the clock when there is one, since
# that is the space `HybridSystems` applies the transition in. It is the source mode's own set,
# narrowed by the bounds written on the transition; an unconstrained coordinate stays free.
function _guard_box(ir::ModelIR, t::TransitionIR, x_idx::Vector{Int})
    source = ir.modes[t.source]
    clock = clock_index(ir)
    coords = clock === nothing ? x_idx : vcat(x_idx, clock)
    lo = Float64[]
    hi = Float64[]
    for i in coords
        l, h = _mode_bound(ir, source, i)
        lo_i = max(l, get(t.guard_lower, i, -Inf))
        hi_i = min(h, get(t.guard_upper, i, Inf))
        # Crossed bounds would reach `Hyperrectangle` as a negative radius and assert from inside
        # LazySets, naming neither the transition nor the variable.
        lo_i <= hi_i || error(
            "Transition $(t.id) (mode $(t.source) → mode $(t.target)) has an empty guard: on " *
            "$(describe(ir.variables[i], i)) it asks for $(lo_i) ≤ x ≤ $(hi_i), which no " *
            "state satisfies, so the switch could never be taken. Note the guard is " *
            "intersected with the source mode's own bounds ([$(l), $(h)]), so a guard outside " *
            "that range is empty even if its own bounds are consistent.",
        )
        push!(lo, lo_i)
        push!(hi, hi_i)
    end
    return LazySets.Hyperrectangle(;
        low = SVector{length(coords)}(lo...),
        high = SVector{length(coords)}(hi...),
    )
end

# A multi-variable affine guard `lo ≤ a'x ≤ hi` as the half-spaces bounding it. `a` is spread
# over the state coordinates in declaration order.
function _guard_half_spaces(
    ir::ModelIR,
    t::TransitionIR,
    g::AffineGuard,
    x_idx::Vector{Int},
)
    a = zeros(length(x_idx))
    for (i, c) in g.coefficients
        j = findfirst(==(i), x_idx)
        j === nothing && error(
            "A guard constrains the state at the moment of the switch, but this one mentions " *
            "$(describe(ir.variables[i], i)), which is not a state.",
        )
        a[j] = c
    end

    spaces = LazySets.HalfSpace{Float64, Vector{Float64}}[]

    # Coefficients that cancel (`x - x <= 1`) leave a constraint mentioning no state at all.
    # It is then either vacuous or unsatisfiable, and neither is a half-space — LazySets would
    # assert on the zero normal vector without naming the transition.
    if all(iszero, a)
        g.lower <= 0.0 <= g.upper && return spaces      # vacuous: constrains nothing
        return error(
            "Transition $(t.id) (mode $(t.source) → mode $(t.target)) has a guard constraint " *
            "whose coefficients cancel out, leaving $(g.lower) ≤ 0 ≤ $(g.upper) — which " *
            "nothing satisfies, so the switch could never be taken.",
        )
    end

    isfinite(g.upper) && push!(spaces, LazySets.HalfSpace(copy(a), g.upper))
    isfinite(g.lower) && push!(spaces, LazySets.HalfSpace(-a, -g.lower))
    return spaces
end

function _check_guard_set(vars, S, x_idx::Vector{Int}, t::TransitionIR)
    idx = [v.value for v in vars]
    idx == x_idx && LazySets.dim(S) == length(x_idx) || error(
        "The `Guard` set on transition $(t.id) must be given over the whole state vector, in " *
        "declaration order: got a $(nameof(typeof(S))) of dimension $(LazySets.dim(S)) over " *
        "$(length(idx)) variable(s), against $(length(x_idx)) state(s).",
    )
    return
end

# The full enabling condition: the box of per-coordinate bounds, narrowed by every half-space
# and every `Guard` set written on the transition.
#
# The intersections are built *binary and lazily*. `IntersectionArray` would be the obvious
# container, but the grid discretisation resolves it by computing a concrete `intersection`,
# which has no method for most set pairs; nested binary `Intersection` is evaluated through
# support functions and membership instead, and enumerates correctly for every combination.
function _guard_set(ir::ModelIR, t::TransitionIR, x_idx::Vector{Int})
    box = _guard_box(ir, t, x_idx)
    has_set_guard(t) || return box

    # With a clock the guard lives in `[x; t]`, so an `x`-only set would have to be extruded
    # across the time axis as a Cartesian product — which the `INNER` discretisation cannot
    # enumerate. Per-coordinate bounds already span the clock correctly, so they still work.
    clock_index(ir) === nothing || error(
        "Transition $(t.id) carries a guard that is not a box, but the model has a clock. " *
        "A non-box guard would have to be extruded across the time axis, which the " *
        "discretization does not support. Use per-coordinate bounds on a clocked model, or " *
        "drop the clock.",
    )

    guard = box
    for g in t.affine_guards
        for half_space in _guard_half_spaces(ir, t, g, x_idx)
            guard = LazySets.Intersection(guard, half_space)
        end
    end
    for (vars, S) in t.guard_sets
        _check_guard_set(vars, S, x_idx, t)
        guard = LazySets.Intersection(guard, S)
    end
    return guard
end

# The reset of the physical coordinates, or `nothing` when the transition leaves them alone.
function _physical_reset(ir::ModelIR, t::TransitionIR, backend, x_idx::Vector{Int})
    any(i -> haskey(t.resets, i), x_idx) || return nothing

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
    u0 = SVector{length(u_idx), Float64}(zeros(length(u_idx))...)
    return x -> reset(x, u0)
end

# The reset applied when the switch is taken. With nothing written on the transition this is
# the identity, which is what switching mode usually does to the state.
function _reset_function(ir::ModelIR, t::TransitionIR, backend, x_idx::Vector{Int})
    isempty(t.resets) && return identity
    physical = _physical_reset(ir, t, backend, x_idx)
    clock = clock_index(ir)

    if clock === nothing
        return physical === nothing ? identity : physical
    end

    # With a clock the map runs on `[x; t]`. A clock reset must be a constant — "restart the
    # timer at 0" — because the time axis is discretized independently of the dynamics.
    reset_expr = get(t.resets, clock, nothing)
    clock_value = reset_expr === nothing ? nothing : _clock_rate(reset_expr)
    reset_expr === nothing ||
        clock_value !== nothing ||
        error("The clock reset on transition $(t.id) must be a constant, e.g. `Δ(t) == 0`.")

    physical === nothing && clock_value === nothing && return identity
    n = length(x_idx)
    return function (state)
        x = @view state[1:n]
        next_x = physical === nothing ? collect(x) : collect(physical(x))
        next_t = clock_value === nothing ? state[end] : clock_value
        return vcat(next_x, next_t)
    end
end

"""
    build_hybrid_system(ir, backend) -> HybridSystems.HybridSystem

Assemble the hybrid automaton: one system per mode, and one
[`GuardedResetMap`](@ref Dionysos.System.GuardedResetMap) per transition.
"""
function build_hybrid_system(ir::ModelIR, backend)
    check_time_domain(ir)
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
                _guard_set(ir, t, x_idx),
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
    clock = clock_index(ir)
    entries = filter(m.specs) do e
        e.kind === kind || return false
        # The clock carries the *time window* of the spec, not one of its state coordinates.
        return clock === nothing || !(MOI.VariableIndex(clock) in e.variables)
    end
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
    return LazySets.Hyperrectangle(;
        low = SVector{length(x_idx)}(lo...),
        high = SVector{length(x_idx)}(hi...),
    )
end

# The time window of a mode's spec: the interval written on the clock with `final(t) in …`,
# defaulting to the mode's whole time domain.
function _mode_time_window(ir::ModelIR, m::ModeIR, kind::SpecKind, clock::Int)
    lo, hi = _mode_bound(ir, m, clock)
    for e in m.specs
        e.kind === kind || continue
        e.variables == [MOI.VariableIndex(clock)] || continue
        interval = e.set
        isfinite(interval.lower) && (lo = interval.lower)
        isfinite(interval.upper) && (hi = interval.upper)
    end
    return lo, hi
end

# One mode's entry in a `HybridSpec`: the set, wrapped in the mode's time window when the model
# has a clock.
function _mode_spec(ir::ModelIR, m::ModeIR, kind::SpecKind, set, clock)
    spec = PR.StateSpec(set)
    clock === nothing && return spec
    lo, hi = _mode_time_window(ir, m, kind, clock)
    return PR.TimedSpec(spec, lo, hi)
end

# What a mode that carries no set of this kind means. For `FINAL` it is simply not a goal, so it
# is left out of the spec. For `ALWAYS` leaving it out would say "this mode is forbidden
# entirely" — writing one `Always` would silently outlaw every mode it was not written on — so
# the neutral reading is the mode's own state set.
_absent_mode_is_whole(kind::SpecKind) = kind === ALWAYS

function _hybrid_spec(ir::ModelIR, kind::SpecKind, x_idx::Vector{Int})
    clock = clock_index(ir)
    explicit = Dict{Int, Any}()
    for k in mode_ids(ir)
        set = _mode_spec_set(ir, ir.modes[k], kind, x_idx)
        set === nothing || (explicit[k] = set)
    end
    isempty(explicit) && return nothing

    fill_absent = _absent_mode_is_whole(kind)
    pairs = Pair{Int, PR.AbstractSpecification}[]
    for k in mode_ids(ir)
        m = ir.modes[k]
        set = get(explicit, k, nothing)
        if set === nothing
            fill_absent || continue
            set = _mode_box(ir, m, x_idx)
        end
        push!(pairs, k => _mode_spec(ir, m, kind, set, clock))
    end
    return PR.HybridSpec(Dict(pairs))
end

# The whole state space, mode by mode: the safe set of a reach-and-stay model that constrains
# only where it must end up, not where it may pass through.
function _whole_state_spec(ir::ModelIR, x_idx::Vector{Int})
    clock = clock_index(ir)
    pairs = map(mode_ids(ir)) do k
        m = ir.modes[k]
        return k => _mode_spec(ir, m, ALWAYS, _mode_box(ir, m, x_idx), clock)
    end
    return PR.HybridSpec(Dict(pairs))
end

# `_hybrid_spec` hands back the sets, so the `EventuallyAlways` option is read separately. The
# flag is written per mode; requiring agreement is better than letting one mode's reading win.
function _hybrid_stay_on_first_entry(ir::ModelIR)
    flags = Bool[]
    for k in mode_ids(ir), e in ir.modes[k].specs
        e.kind === EVENTUALLY_ALWAYS && push!(flags, e.stay_on_first_entry)
    end
    isempty(flags) && return false
    allequal(flags) || error(
        "The modes disagree on `stay_on_first_entry`: it is part of what `◇□` means, so every " *
        "`EventuallyAlways` in one model must use the same reading.",
    )
    return first(flags)
end

# Where the clock starts: its declared `start`, otherwise the beginning of its domain — a timer
# normally starts at zero rather than in the middle of its range.
function _clock_start(ir::ModelIR, clock::Int)
    v = ir.variables[clock]
    isfinite(v.start.lower) && return v.start.lower
    return v.lower
end

# Where the model starts, as a mode-indexed set. It used to be the *centre* of the declared box,
# so `Start(S)` meant "start at the middle of `S`" and a model reported `OPTIMAL` for one point
# while the states around it — inside the set the user declared — went unsynthesized.
#
# `OUTER` because a controller has to handle every state the initial set touches, including the
# cells it only clips. That is the reading the flat path uses too (`get_states_from_set(…,
# MP.OUTER)`), and it is what keeps a `start = v` point, whose box is degenerate, from
# discretizing to nothing.
#
# More than one mode may carry a `start`: with a set-valued initial condition, beginning in
# either of two modes is a meaningful thing to ask for.
function _hybrid_initial_spec(ir::ModelIR, x_idx::Vector{Int})
    clock = clock_index(ir)
    starting = [k for k in mode_ids(ir) if any(e -> e.kind === START, ir.modes[k].specs)]
    isempty(starting) && (starting = [first(mode_ids(ir))])

    pairs = map(starting) do k
        m = ir.modes[k]
        set = _mode_spec_set(ir, m, START, x_idx)
        set === nothing && (set = _coordinate_box(ir, x_idx, :start))
        spec = PR.StateSpec(set, UT.OUTER)
        if clock !== nothing
            t0 = _clock_start(ir, clock)
            spec = PR.TimedSpec(spec, t0, t0)
        end
        return k => spec
    end
    return PR.HybridSpec(Dict(pairs))
end

# Every specification kind the hybrid lowering consumes. A kind absent from this list would be
# parsed onto a mode and then silently dropped on the way to the problem — which is exactly how
# `EventuallyAlways` came to be accepted and ignored. Adding a `SpecKind` should fail here until
# the lowering handles it, rather than producing a controller for a different specification.
const _HYBRID_SPEC_KINDS = (START, FINAL, ALWAYS, EVENTUALLY_ALWAYS)

function _reject_unhandled_specs(ir::ModelIR)
    for k in mode_ids(ir), e in ir.modes[k].specs
        e.kind in _HYBRID_SPEC_KINDS && continue
        error(
            "Mode $k carries a `$(e.kind)` specification, which the hybrid lowering does not " *
            "handle. It would otherwise be dropped and the controller synthesized for the " *
            "rest of the model.",
        )
    end
    return
end

"""
    build_hybrid_problem(ir, backend; time_step = nothing) -> Problem.ProblemType

Assemble the hybrid control problem: a `HybridSystem`, an augmented initial state `(x, mode)`,
and a mode-indexed specification.
"""
function build_hybrid_problem(ir::ModelIR, backend; time_step = nothing)
    _reject_unhandled_specs(ir)
    x_idx = state_indices(ir)
    system = build_hybrid_system(ir, backend)
    initial_state = _hybrid_initial_spec(ir, x_idx)

    target = _hybrid_spec(ir, FINAL, x_idx)
    stay = _hybrid_spec(ir, EVENTUALLY_ALWAYS, x_idx)
    safe = _hybrid_spec(ir, ALWAYS, x_idx)

    target === nothing ||
        stay === nothing ||
        error(
            "The model asks both to reach a set (`Final`) and to reach-and-stay in one " *
            "(`EventuallyAlways`). Keep one: `◇□ T` already reaches `T`.",
        )

    if stay !== nothing
        # `ReachAndStayProblem` always carries a safe set, so a model that says only where to
        # settle is safe everywhere on the way.
        return PR.ReachAndStayProblem(
            system,
            initial_state,
            stay,
            safe === nothing ? _whole_state_spec(ir, x_idx) : safe,
            _horizon(ir, true, time_step);
            stay_on_first_entry = _hybrid_stay_on_first_entry(ir),
        )
    elseif target !== nothing
        # With both markers this is reach-avoid: before `safe_set` existed the `Always` set was
        # dropped here, since the hybrid path has no state space to fold it into.
        return PR.OptimalControlProblem(
            system,
            initial_state,
            target,
            nothing,
            nothing,
            _horizon(ir, false, time_step);
            safe_set = safe,
        )
    elseif safe !== nothing
        return PR.SafetyProblem(system, initial_state, safe, _horizon(ir, true, time_step))
    end
    return error(
        "A hybrid model needs a specification on at least one mode: add a `Final`, `Always` or " *
        "`EventuallyAlways` set, or a `final(x[i]) in …` constraint, to the mode you care about.",
    )
end
