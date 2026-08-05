# ----------------------------------------------------------------------------------------
# Lowering, part 2: `ModelIR` → a `Problem.ProblemType`.
#
# The problem *type* is inferred from which specification markers the user wrote, in the same
# spirit as inferring a variable's role from how it is used. See `src/wrapper/README.md` §5.
# ----------------------------------------------------------------------------------------

# The per-coordinate `start`/`final` intervals as a box over the state coordinates. A
# coordinate the user left unconstrained falls back to the *variable's own bounds*; before
# Phase 3 it contributed ±Inf, from which `UT.box` built a NaN radius and threw deep inside
# LazySets (plan.md, L7).
function _coordinate_box(ir::ModelIR, x_idx::Vector{Int}, field::Symbol)
    n = length(ir.variables)
    lb = Vector{Float64}(undef, n)
    ub = Vector{Float64}(undef, n)
    for (i, v) in enumerate(ir.variables)
        interval = getproperty(v, field)
        lb[i] = isfinite(interval.lower) ? interval.lower : v.lower
        ub[i] = isfinite(interval.upper) ? interval.upper : v.upper
    end
    return UT.box(_svec(lb, x_idx), _svec(ub, x_idx))
end

# The single specification set of kind `kind`, or `nothing`. A specification set must span
# exactly the state vector: a set over a subset of coordinates would need extruding, which the
# per-coordinate `final(x[i]) in …` form already expresses.
function _unique_spec(ir::ModelIR, kind::SpecKind, x_idx::Vector{Int})
    entries = filter(e -> e.kind === kind, ir.specs)
    isempty(entries) && return nothing
    length(entries) == 1 || error(
        "Several `$kind` specification sets were given; keep one (intersect them yourself " *
        "if that is what you meant).",
    )
    entry = entries[]
    vars = [v.value for v in entry.variables]
    vars == x_idx || error(
        "A `$kind` set must be given over the whole state vector, in declaration order. " *
        "Got $(length(vars)) variable(s) against $(length(x_idx)) state(s). Constrain single " *
        "coordinates with `final(x[i]) in MOI.Interval(a, b)` instead.",
    )
    return entry.set
end

# `horizon` is in seconds for a continuous-time model and in steps for a discrete-time one.
# The conversion happens here because no abstraction solver calls `Problem.discretize_problem`
# — they consume `problem.time` verbatim (plan.md, D7).
function _horizon(ir::ModelIR, round_up::Bool, time_step)
    ir.horizon === nothing && return PR.Infinity()
    ir.time_domain === DISCRETE && return round(Int, ir.horizon)
    time_step === nothing && error(
        "Setting `horizon` on a continuous-time model also needs `time_step`, which converts " *
        "the horizon (in seconds) into a number of abstraction steps.",
    )
    return PR.discretize_time(ir.horizon, time_step; round_up = round_up)
end

"""
    build_problem(ir::ModelIR, backend; time_step = nothing) -> Problem.ProblemType

Lower the model, compiling its dynamics through `backend`. A model with modes goes to
[`build_hybrid_problem`](@ref); everything else is monolithic.
"""
function build_problem(ir::ModelIR, backend::AbstractDynamicsBackend; time_step = nothing)
    is_hybrid(ir) && return build_hybrid_problem(ir, backend; time_step = time_step)
    f = compile_dynamics(backend, ir, ir.dynamics)
    return build_problem(ir, f; time_step = time_step)
end

"""
    build_problem(ir::ModelIR, f; time_step = nothing) -> Problem.ProblemType

Assemble the control problem. Which `ProblemType` is built follows from the specification
markers present:

| `Final` | `Always` | `EventuallyAlways` | problem |
| :-: | :-: | :-: | :--- |
| ✓ | | | `OptimalControlProblem` (reach) |
| ✓ | ✓ | | `OptimalControlProblem`, avoid-set folded into `X` |
| | ✓ | | `SafetyProblem` |
| | opt | ✓ | `ReachAndStayProblem` |
| | | | `AlternatingSimulationProblem` — build the abstraction only |

A per-coordinate `final(x[i]) in …` counts as a `Final` set.
"""
function build_problem(ir::ModelIR, f; time_step = nothing)
    x_idx = state_indices(ir)

    start_set = _unique_spec(ir, START, x_idx)
    final_set = _unique_spec(ir, FINAL, x_idx)
    always_set = _unique_spec(ir, ALWAYS, x_idx)
    stay_set = _unique_spec(ir, EVENTUALLY_ALWAYS, x_idx)

    initial_set = start_set === nothing ? _coordinate_box(ir, x_idx, :start) : start_set

    has_coordinate_target = any(_has_target, ir.variables)
    target_set = if final_set !== nothing
        final_set
    elseif has_coordinate_target
        _coordinate_box(ir, x_idx, :target)
    else
        nothing
    end

    if stay_set !== nothing
        safe_set = always_set === nothing ? state_box(ir, x_idx) : always_set
        return PR.ReachAndStayProblem(
            build_system(ir, f),
            initial_set,
            stay_set,
            safe_set,
            _horizon(ir, true, time_step),
        )
    elseif target_set !== nothing
        return PR.OptimalControlProblem(
            build_system(ir, f; restrict_to = always_set),
            initial_set,
            target_set,
            nothing,
            nothing,
            _horizon(ir, false, time_step),
        )
    elseif always_set !== nothing
        return PR.SafetyProblem(
            build_system(ir, f),
            initial_set,
            always_set,
            _horizon(ir, true, time_step),
        )
    else
        # No specification at all: the task is the abstraction itself.
        return PR.AlternatingSimulationProblem(build_system(ir, f), nothing)
    end
end
