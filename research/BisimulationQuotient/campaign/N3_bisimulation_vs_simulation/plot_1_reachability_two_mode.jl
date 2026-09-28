# The figure: the grid's certified region growing toward the bisimulation's, never past it.
#
# One panel per grid resolution, then the bisimulation. Green is what the controller wins and red
# what it does not, as everywhere else in this folder, so the panels can be read side by side
# against the other experiments' figures. The target region is an outline in a colour disjoint from
# both, since a filled region would tint the cells underneath and blur the very thing the figure is
# about.
#
# Both abstractions discretize the SAME set: the grid's state set is the certificate's outermost
# sublevel set, not the working set X. Without that the grid covers X while the quotient covers a
# rotated box reaching well beyond it, the panels have different extents, and the quotient's region
# *looks* smaller than the grid's while actually being larger.
#
# What the eye should see: green spreading across the panels as h shrinks, and the last panel --
# the exact answer -- containing all of it. Measured pointwise in `grid_vs_quotient.jl`:
#
#   h=0.40     4 of 53 probes certified,   209 cells      h=0.10   31 of 53,   3 839 cells
#   h=0.20    27 of 53,                    919 cells      h=0.05   43 of 53,  15 701 cells
#   bisimulation                                          53 of 53,     379 cells
#
# and the grid certifies nothing the bisimulation does not, at any resolution.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))

using LinearAlgebra
import MathematicalSystems

const SY = DI.Symbolic
const OPDS = OP.DiscreteSystems

const A = [[0.70 0.10; 0.00 0.65], [0.60 -0.15; 0.10 0.55]]

mode_of(u) = clamp(round(Int, u[1]), 1, length(A))
step_map(x, u) = A[mode_of(u)] * x

function exact_image_box(rect, u)
    c = LazySets.center(rect)
    r = LazySets.radius_hyperrectangle(rect)
    M = A[mode_of(u)]
    return LazySets.Hyperrectangle(M * c, abs.(M) * r)
end

"""
Controllable and uncontrollable grid cells for reaching `target`, as plottable boxes.
"""
function grid_cells(X, target, h)
    system = MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem(
        step_map,
        2,
        1,
        X,
        LazySets.Hyperrectangle(; low = SVector(1.0), high = SVector(2.0)),
    )
    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("concrete_problem"),
        PR.AlternatingSimulationProblem(system, X),
    )
    grid = MP.GridFree(SVector(0.0, 0.0), SVector(h, h))
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), grid)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("input_grid"),
        MP.GridFree(SVector(1.0), SVector(1.0)),
    )
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("approx_mode"),
        AB.UniformGridAbstraction.USER_DEFINED,
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("overapproximation_map"), exact_image_box)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.optimize!(optimizer)

    abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    automaton = SY.get_automaton(abstract_system)
    mapping = SY.get_state_mapping(abstract_system)
    targets = collect(MP.get_states_from_set(mapping, target, MP.INNER))
    _, controllable, _, _ = OPDS.compute_worst_case_cost_controller(automaton, targets)

    won = Set(controllable)
    boxes_win, boxes_lose = LazySets.Hyperrectangle[], LazySets.Hyperrectangle[]
    r = SVector(h / 2, h / 2)
    for q in MP.enum_states(mapping)
        c = MP.get_coord_by_pos(grid, MP.get_pos_by_state(mapping, q))
        push!(q in won ? boxes_win : boxes_lose, LazySets.Hyperrectangle(c, r))
    end
    return boxes_win, boxes_lose, SY.get_n_state(abstract_system)
end

# Every panel is clipped to the working set. Without this the quotient panel spans the PCLF's
# sublevel sets -- a rotated box reaching well beyond X -- and its certified region *looks* smaller
# than the grid's while actually being larger, which is the opposite of what the figure shows.
const PANEL_LIMS = (-5.0, 5.0)

function grid_panel(boxes_win, boxes_lose, problem, title)
    fig = plot(; aspect_ratio = :equal, legend = false, title = title)
    for b in boxes_lose
        plot!(fig, b; color = LOSING_COLOR, fillalpha = 1.0, linewidth = 0)
    end
    for b in boxes_win
        plot!(fig, b; color = WINNING_COLOR, fillalpha = 1.0, linewidth = 0)
    end
    plot!(
        fig,
        problem;
        plot_region = false,
        observation_region_alpha = 0.0,
        observation_colors = OBSERVATION_COLORS,
        observation_linewidth = 2.0,
    )
    xlims!(fig, PANEL_LIMS...)
    ylims!(fig, PANEL_LIMS...)
    return fig
end

function quotient_winning(quotient, target_obs::Int)
    win = Set{Int}(id for (id, state) in quotient.states if state.obs == target_obs)
    changed = true
    while changed
        changed = false
        for (id, state) in quotient.states
            id in win && continue
            if any(t -> t[2] in win, state.next)
                push!(win, id)
                changed = true
            end
        end
    end
    return win
end

gr()

(; f, problem, X, R1, R2) = two_mode_problem()

pclf = PCLF.compute_symmetric_2n_faces_polyhedral_pieces_pclf(
    f,
    PCLF.generate_DeBruijn_edges(2, 0),
    JuMP.optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => 1000);
    Gmats = rotation_templates([1]; θ = π / 6, mode = :rotation),
    MLF = true,
    verbose = false,
)
(; optimizer, quotient) = build_quotient(
    problem,
    pclf;
    atol = 1e-3,
    max_levels = 100,
    max_slices = 8,
    print_level = 0,
)
win = quotient_winning(quotient, 1)
lose = setdiff(Set(keys(quotient.states)), win)

# The grid discretizes the certificate's outermost sublevel set, which is exactly what the quotient
# covers, so the panels show the same region and the cell counts in the titles are comparable.
domain = PCLF.get_sublevel_set(pclf.pieces[first(keys(pclf.pieces))], maximum(optimizer.Γ))

panels = Any[]
for h in [0.4, 0.2, 0.1, 0.05]
    boxes_win, boxes_lose, n = grid_cells(domain, R1, h)
    push!(panels, grid_panel(boxes_win, boxes_lose, problem, "grid h=$h  ($n cells)"))
    println("h = ", h, ": ", n, " cells, ", length(boxes_win), " controllable")
end

exact = plot(;
    aspect_ratio = :equal,
    legend = false,
    title = "bisimulation ($(length(quotient.states)) cells)",
)
_plot_winning_losing!(exact, quotient, collect(win), collect(lose), nothing; linewidth = 0)
plot!(
    exact,
    problem;
    plot_region = false,
    observation_region_alpha = 0.0,
    observation_colors = OBSERVATION_COLORS,
    observation_linewidth = 2.0,
)
xlims!(exact, PANEL_LIMS...)
ylims!(exact, PANEL_LIMS...)
push!(panels, exact)

fig = plot(panels...; layout = (1, length(panels)), size = (420 * length(panels), 420))
savefig(fig, joinpath(@__DIR__, "certified_sets.png"))
println("\nwrote certified_sets.png")
println("quotient: ", length(quotient.states), " cells, ", length(win), " controllable")
