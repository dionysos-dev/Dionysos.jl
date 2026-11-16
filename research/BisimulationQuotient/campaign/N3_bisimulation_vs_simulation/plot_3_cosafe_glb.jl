# The Gol-Lazar-Belta figure: grid against bisimulation, under their co-safe LTL formula.
#
#     φ = (!R2 U D) & F(R1) & ((R3 -> X(!R1)) U D)
#
# Same convention as everywhere in this folder: green is what the controller wins, red what it does
# not, observation regions as outlines in a disjoint colour. Both abstractions discretize the same
# set, the certificate's outermost sublevel set.
#
# READ THIS BEFORE USING THE FIGURE. The grid's certified region **shrinks** as the grid refines --
# 13.0, 10.0, 3.0 in area at h = 1.0, 0.5, 0.25 -- which a sound over-approximation must not do. The
# anomaly is under investigation (see plan.md); the panels are drawn precisely so it can be seen
# rather than argued about. The quotient panel is trustworthy: it reproduces E8's recorded 8 794 of
# 10 611 cells.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))

using Spot
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
D = MOI.get(opt0, MOI.RawOptimizerAttribute("D"))

sets = [PCLF.get_sublevel_set(p, maximum(opt0.Γ)) for (_, p) in opt0.pclf.pieces]
domain = sets[argmax([LazySets.volume(s) for s in sets])]

φ = ltl"((!R2 U D) & F(R1) & ((R3 -> X(!R1)) U D))"
regions = Dict(:D => D, :R1 => R1, :R2 => R2, :R3 => R3)
x0 = SVector(-4.0, -7.0)

"""
Certified and uncertified grid cells at step `h`, as plottable boxes.
"""
function grid_cells(h)
    system = MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem(
        step_map,
        2,
        1,
        domain,
        LazySets.Hyperrectangle(; low = SVector(1.0), high = SVector(2.0)),
    )
    concrete = PR.CoSafeLTLProblem(
        system,
        LazySets.Hyperrectangle(; low = [x0[1], x0[2]], high = [x0[1], x0[2]]),
        Dionysos.spot_stepper(φ),
        regions,
        Dict{Symbol, Any}(
            :D => MP.INNER,
            :R1 => MP.INNER,
            :R2 => MP.OUTER,
            :R3 => MP.OUTER,
        ),
    )
    opt = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(opt, MOI.RawOptimizerAttribute("concrete_problem"), concrete)
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
    won = Set(opt.control_solver.abstract_optimizer.controllable_set)

    r = SVector(h / 2, h / 2)
    win, lose = LazySets.Hyperrectangle[], LazySets.Hyperrectangle[]
    for q in MP.enum_states(mapping)
        c = MP.get_coord_by_pos(grid, MP.get_pos_by_state(mapping, q))
        push!(q in won ? win : lose, LazySets.Hyperrectangle(c, r))
    end
    return win, lose, SY.get_n_state(abstract_system)
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

println("solving their formula on the quotient ...")
qres = synthesize_cosafe_ltl(
    f,
    quotient,
    Dionysos.spot_stepper(φ),
    regions,
    Dict(:D => -1, :R1 => 1, :R2 => 2, :R3 => 3),
    x0;
    print_level = 0,
)
qwin = collect(qres.controllable_set)
qlose = collect(setdiff(Set(keys(quotient.states)), Set(qwin)))
println("quotient: ", length(qwin), " of ", length(quotient.states), " cells certified")

panels = Any[]
for h in [1.0, 0.5, 0.25]
    win, lose, n = grid_cells(h)
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
        ")",
    )
    push!(panels, panel(win, lose, "grid h=$h  ($n cells, area $(round(area; digits=1)))"))
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
savefig(fig, joinpath(@__DIR__, "glb_cosafe_sets.png"))
println("\nwrote glb_cosafe_sets.png")
