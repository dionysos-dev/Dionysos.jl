# ----------------------------------------------------------------------------------------
# Lowering, part 2: `ModelIR` → a `Problem.ProblemType`.
# ----------------------------------------------------------------------------------------

# The per-coordinate `start`/`final` intervals, assembled into a box over the state
# coordinates. Coordinates the user left unconstrained currently contribute ±Inf; Phase 3
# replaces that with the variable's own bounds (plan.md, L7).
function _interval_box(intervals::Vector{MOI.Interval{Float64}}, x_idx::Vector{Int})
    lb = _svec([s.lower for s in intervals], x_idx)
    ub = _svec([s.upper for s in intervals], x_idx)
    return UT.box(lb, ub)
end

"""
    build_problem(ir::ModelIR, f) -> Problem.ProblemType

Assemble the control problem: the system from [`build_system`](@ref) plus the initial and
target sets declared with `start`/`final`.
"""
function build_problem(ir::ModelIR, f)
    x_idx = state_indices(ir)
    system = build_system(ir, f)

    initial_set = _interval_box([v.start for v in ir.variables], x_idx)
    target_set = _interval_box([v.target for v in ir.variables], x_idx)

    return PR.OptimalControlProblem(
        system,
        initial_set,
        target_set,
        nothing,
        nothing,
        PR.Infinity(),
    )
end
