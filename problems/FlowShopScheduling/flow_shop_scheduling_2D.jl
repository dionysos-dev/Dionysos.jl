module FlowShopScheduling2D

using StaticArrays

import HybridSystems as HS
import MathematicalSystems as MS
import LazySets

using LinearAlgebra
using Plots

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem

# Reset on a guarded switch: reinitialize the 2-D state to `x_init` and push the clock
# forward to at least `t_min` (the state passed in is the augmented `[x1, x2, t]`).
# On task completion the state restarts at `x_init`, and the clock is held at `t_min` when the
# task finished early. Guard and reset live in the augmented `[x; t]` space.
FlowShopResetMap(domain, x_init, t_min) =
    ST.GuardedResetMap(domain, state -> vcat(x_init, max(t_min, state[end])))

function get_matrices()
    A1 = SMatrix{2, 2}(0.95, 0.01, 0.0, 0.98)
    A2 = SMatrix{2, 2}(0.90, 0.05, -0.01, 0.92)
    A3 = SMatrix{2, 2}(0.85, 0.10, 0.0, 0.88)
    return A1, A2, A3
end

# 3-task flowshop with 2-D linear dynamics per task plus a per-task clock. Each task
# has its own state/input constraints and time window; a guarded switch resets the
# state and advances the clock, and the final target is a region in (x1, x2) within a
# time window. This is the 2-D analogue of `flow_shop_scheduling_3D.jl`, small enough
# to abstract and solve on a laptop.
function system()
    A1, A2, A3 = get_matrices()

    B1 = SMatrix{2, 2}(1.0, 0.1, 0.0, 0.8)
    B2 = SMatrix{2, 2}(0.9, 0.2, -0.001, 0.7)

    task1_dynamics(x, u) = A1 * x + B1 * u
    task2_dynamics(x, u) = A2 * x + B2 * u
    task3_dynamics(x, u) = A3 * x .+ u

    # State and input spaces (2-D)
    X1 = LazySets.Hyperrectangle(; low = [-1.0, -1.0], high = [5.0, 5.0])
    U1 = LazySets.Hyperrectangle(; low = [-1.0, -1.0], high = [4.0, 4.0])
    X2 = LazySets.Hyperrectangle(; low = [-2.5, -2.5], high = [2.5, 2.5])
    U2 = LazySets.Hyperrectangle(; low = [-1.2, -1.2], high = [9.2, 6.2])
    X3 = LazySets.Hyperrectangle(; low = [-1.0, -1.0], high = [5.0, 3.0])
    U3 = LazySets.Hyperrectangle(; low = [-1.5, -1.5], high = [7.5, 5.5])

    task1_system =
        MS.ConstrainedBlackBoxControlContinuousSystem(task1_dynamics, 2, 2, X1, U1)
    task2_system =
        MS.ConstrainedBlackBoxControlContinuousSystem(task2_dynamics, 2, 2, X2, U2)
    task3_system =
        MS.ConstrainedBlackBoxControlContinuousSystem(task3_dynamics, 2, 2, X3, U3)

    # Time systems (clocks) for each task
    task_1_time_system = MS.ConstrainedLinearContinuousSystem(
        [1.0;;],
        LazySets.Hyperrectangle(; low = [0.0], high = [2.0]),
    )
    task_2_time_system = MS.ConstrainedLinearContinuousSystem(
        [1.0;;],
        LazySets.Hyperrectangle(; low = [1.5], high = [4.0]),
    )
    task_3_time_system = MS.ConstrainedLinearContinuousSystem(
        [1.0;;],
        LazySets.Hyperrectangle(; low = [5.0], high = [7.0]),
    )

    modes_systems = [
        ST.VectorContinuousSystem([task1_system, task_1_time_system]),
        ST.VectorContinuousSystem([task2_system, task_2_time_system]),
        ST.VectorContinuousSystem([task3_system, task_3_time_system]),
    ]

    # Guards (acceptance regions) over [x1, x2, t]
    task1_target = LazySets.Hyperrectangle(; low = [0.5, 1.5, 0.0], high = [5.0, 5.0, 2.0])
    task2_target = LazySets.Hyperrectangle(; low = [1.0, 0.0, 1.5], high = [2.5, 2.5, 4.0])

    # Reset maps for each transition. The reset's `t_min` advances the clock to the
    # start of the next task's time window (task 2: 1.5, task 3: 5.0), so the switched
    # state lands inside the target task's clock domain.
    t1_t2_reset_map = FlowShopResetMap(task1_target, [0.0, 0.0], 2.0)
    t2_t3_reset_map = FlowShopResetMap(task2_target, [1.0, 1.0], 5.0)
    reset_maps = [t1_t2_reset_map, t2_t3_reset_map]

    automaton = HS.GraphAutomaton(3)
    HS.add_transition!(automaton, 1, 2, 1)
    HS.add_transition!(automaton, 2, 3, 2)
    switchings = [HS.AutonomousSwitching(), HS.AutonomousSwitching()]

    return HS.HybridSystem(automaton, modes_systems, reset_maps, switchings)
end

function create_jacobian(A_matrix)
    L = MMatrix{2, 2, Float64}(A_matrix)
    n = size(L, 1)
    for i in 1:n
        for j in 1:n
            if i != j
                L[i, j] = abs(L[i, j])
            end
        end
    end
    return L
end

function jacobian_bounds()
    A1, A2, A3 = get_matrices()
    return [u -> create_jacobian(A1), u -> create_jacobian(A2), u -> create_jacobian(A3)]
end

function problem()
    hybrid_system = system()

    initial_state = ([0.0, 0.0], 0.0, 1)   # (x1, x2, t, mode) -> ([x1,x2], t, mode)

    target_set = PR.hybrid_reach_spec(
        [LazySets.Hyperrectangle(; low = [0.5, 0.5], high = [5.0, 3.0])],
        [LazySets.Hyperrectangle(; low = [5.0], high = [7.0])],
        [3],
    )

    transition_cost_function = (aug_state, u) -> 1.0

    return PR.OptimalControlProblem(
        hybrid_system,
        initial_state,
        target_set,
        nothing,
        transition_cost_function,
    )
end

# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------

const TASK_COLORS = [:steelblue, :darkorange, :seagreen, :crimson, :purple]
_task_color(k) = TASK_COLORS[mod1(k, length(TASK_COLORS))]

# Draw a filled/dashed rectangle [tlo, thi] × [xlo, xhi].
function _rectangle!(fig, tlo, thi, xlo, xhi; kwargs...)
    plot!(fig, [tlo, thi, thi, tlo, tlo], [xlo, xlo, xhi, xhi, xlo]; kwargs...)
    return fig
end

"""
    system_plot!(; problem, xlims = (-2.5, 5.2), ylims = (-2.5, 5.2))

Return a per-frame drawer `(fig, x, u, k) -> fig` for
[`Dionysos.animate_trajectory_dashboard`](@ref). Given the continuous state `x = [x1, x2]`,
input `u`, and active task `k`, it draws — in the `(x1, x2)` plane — every task's guard
acceptance region and the final target, then the current point coloured by the active task
and annotated with the task and input.

The hybrid closed-loop result is a channelled `Trajectory` (via
`HybridSystemAbstraction.channelled_trajectory`), so the dashboard reads the continuous
state, feeds the task `k` here, and adds the task-vs-time panel.
"""
function system_plot!(; problem, xlims = (-2.5, 5.2), ylims = (-2.5, 5.2))
    hs = problem.system
    transitions = collect(HS.transitions(hs.automaton))
    task_modes = collect(HS.states(hs.automaton))
    target = problem.target_set.per_mode[task_modes[end]]
    xtarget = target.base.set

    return function (fig, x, u, k)
        # Guard acceptance regions and the final target, projected to the (x1, x2) plane.
        for tr in transitions
            guard = HS.guard(hs, tr)
            guard === nothing && continue
            src = HS.source(hs.automaton, tr)
            _rectangle!(
                fig,
                LazySets.low(guard, 1),
                LazySets.high(guard, 1),   # x1 range
                LazySets.low(guard, 2),
                LazySets.high(guard, 2);   # x2 range
                linecolor = _task_color(src),
                linestyle = :dash,
                linewidth = 2,
                fill = (0, 0.10),
                fillcolor = _task_color(src),
                label = "",
            )
        end
        _rectangle!(
            fig,
            LazySets.low(xtarget, 1),
            LazySets.high(xtarget, 1),
            LazySets.low(xtarget, 2),
            LazySets.high(xtarget, 2);
            linecolor = :magenta,
            linewidth = 3,
            fill = (0, 0.18),
            fillcolor = :magenta,
            label = "",
        )

        # Current physical state, coloured by the active task.
        scatter!(
            fig,
            [x[1]],
            [x[2]];
            color = _task_color(k),
            markersize = 9,
            marker = :circle,
            label = "",
        )

        u_str =
            u isa AbstractString ? u :
            "[" * join(round.(Float64.(u); digits = 2), ", ") * "]"
        annotate!(
            fig,
            xlims[1] + 0.04 * (xlims[2] - xlims[1]),
            ylims[2] - 0.06 * (ylims[2] - ylims[1]),
            text("task $k\nu = $u_str", 9, :left),
        )
        xlims!(fig, xlims...)
        ylims!(fig, ylims...)
        return fig
    end
end

end # module
