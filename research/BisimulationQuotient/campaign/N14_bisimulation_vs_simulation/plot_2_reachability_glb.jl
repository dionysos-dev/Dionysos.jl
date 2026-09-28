# Gol-Lazar-Belta with PLAIN REACHABILITY: the companion figure to `plot_glb_cosafe.jl`.
#
# Same system, same domain, same grid, same 4.61-way branching -- only the specification differs.
# Here the grid works: the certified region grows monotonically with refinement and converges toward
# the exact answer from below, which is what an over-approximation is supposed to do.
#
#   reach R1     h=0.5   23.6 % of the domain      h=0.125   43.7 %
#                h=0.25  36.7 %
#   their φ      h=0.25   0.71 %                   h=0.125    0.31 %   (and falling)
#
# Put side by side with the co-safe figure, the pair isolates the cause: 4.61-way spurious branching
# is survivable for reachability and fatal for a formula that constrains the ORDER of events. The
# adversary cannot stop the system reaching a set -- the contraction guarantees that -- but every
# spurious branch is another chance to enter R2 early, or to hit R1 immediately after R3.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))

using LinearAlgebra
import MathematicalSystems

const SY = DI.Symbolic
const OPDS = OP.DiscreteSystems

(; f, problem, X, R1, R2, R3) = gol_lazar_belta_problem()
const A = UT.mode_matrices(f)

mode_of(u) = clamp(round(Int, u[1]), 1, length(A))
step_map(x, u) = A[mode_of(u)] * x

function exact_image_box(rect, u)
    c = LazySets.center(rect)
    r = LazySets.radius_hyperrectangle(rect)
    M = A[mode_of(u)]
    return LazySets.Hyperrectangle(M * c, abs.(M) * r)
end

const CACHE = joinpath(dirname(dirname(@__DIR__)), "examples", "gol_lazar_belta_pclf.jld2")
isfile(CACHE) || error("run examples/gol_lazar_belta_pclf.jl first to write the cache")

opt0 = import_optimizer_jld2(CACHE)
quotient = MOI.get(opt0, MOI.RawOptimizerAttribute("bisimulation_quotient"))
sets = [PCLF.get_sublevel_set(p, maximum(opt0.Γ)) for (_, p) in opt0.pclf.pieces]
domain = sets[argmax([LazySets.volume(s) for s in sets])]

"""
Certified and uncertified grid cells for reaching `target`, as plottable boxes.
"""
function grid_cells(target, h)
    system = MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem(
        step_map,
        2,
        1,
        domain,
        LazySets.Hyperrectangle(; low = SVector(1.0), high = SVector(2.0)),
    )
    opt = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(
        opt,
        MOI.RawOptimizerAttribute("concrete_problem"),
        PR.AlternatingSimulationProblem(system, domain),
    )
    grid = MP.GridFree(SVector(0.0, 0.0), SVector(h, h))
    MOI.set(opt, MOI.RawOptimizerAttribute("state_grid"), grid)
    MOI.set(
        opt,
        MOI.RawOptimizerAttribute("input_grid"),
        MP.GridFree(SVector(1.0), SVector(1.0)),
    )
    MOI.set(
        opt,
        MOI.RawOptimizerAttribute("approx_mode"),
        AB.UniformGridAbstraction.USER_DEFINED,
    )
    MOI.set(opt, MOI.RawOptimizerAttribute("overapproximation_map"), exact_image_box)
    MOI.set(opt, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.optimize!(opt)

    abstract_system = MOI.get(opt, MOI.RawOptimizerAttribute("abstract_system"))
    mapping = SY.get_state_mapping(abstract_system)
    targets = collect(MP.get_states_from_set(mapping, target, MP.INNER))
    _, controllable, _, _ =
        OPDS.compute_worst_case_cost_controller(SY.get_automaton(abstract_system), targets)

    won = Set(controllable)
    r = SVector(h / 2, h / 2)
    win, lose = LazySets.Hyperrectangle[], LazySets.Hyperrectangle[]
    for q in MP.enum_states(mapping)
        c = MP.get_coord_by_pos(grid, MP.get_pos_by_state(mapping, q))
        push!(q in won ? win : lose, LazySets.Hyperrectangle(c, r))
    end
    return win, lose, SY.get_n_state(abstract_system)
end

"""
Backward reachability on the quotient: a cell wins when some mode leads into the winning set.
"""
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

const LIMS = (-12.0, 12.0)

function panel(win, lose, title)
    fig = plot(; aspect_ratio = :equal, legend = false, title = title)
    for b in lose
        plot!(fig, b; color = LOSING_COLOR, fillalpha = 1.0, linewidth = 0)
    end
    for b in win
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
    xlims!(fig, LIMS...)
    ylims!(fig, LIMS...)
    return fig
end

gr()

# `two_mode_problem` labels R1 as observation 1; so does `gol_lazar_belta_problem`.
qwin = collect(quotient_winning(quotient, 1))
qlose = collect(setdiff(Set(keys(quotient.states)), Set(qwin)))
println("quotient: ", length(qwin), " of ", length(quotient.states), " cells reach R1")

panels = Any[]
for h in [0.5, 0.25, 0.125]
    win, lose, n = grid_cells(R1, h)
    area = length(win) * h^2
    println(
        "h = ",
        h,
        ": ",
        n,
        " cells, ",
        length(win),
        " certified (area ",
        round(area; digits = 2),
        ", ",
        round(100 * area / 424.12; digits = 1),
        " %)",
    )
    push!(
        panels,
        panel(win, lose, "grid h=$h  ($n cells, $(round(100*area/424.12; digits=1))%)"),
    )
end

exact = plot(;
    aspect_ratio = :equal,
    legend = false,
    title = "bisimulation ($(length(quotient.states)) cells)",
)
_plot_winning_losing!(exact, quotient, qwin, qlose, nothing; linewidth = 0)
plot!(
    exact,
    problem;
    plot_region = false,
    observation_region_alpha = 0.0,
    observation_colors = OBSERVATION_COLORS,
    observation_linewidth = 2.0,
)
xlims!(exact, LIMS...)
ylims!(exact, LIMS...)
push!(panels, exact)

fig = plot(panels...; layout = (1, length(panels)), size = (420 * length(panels), 420))
savefig(fig, joinpath(@__DIR__, "glb_reach_sets.png"))
println("\nwrote glb_reach_sets.png")
