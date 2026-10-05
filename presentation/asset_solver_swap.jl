# The solver-interface slide: the same JuMP model solved twice, with one
# attribute changed. Everything else — the dynamics, the specification, the
# grids, the closed-loop call — is byte-identical between the two runs.
#
# Output: presentation/gallery/solver_swap.png (overwrites)
#
# Run: julia --project=test presentation/asset_solver_swap.jl

ENV["GKSwstype"] = "100"

using StaticArrays, JuMP, Plots
import LazySets
using Symbolics, MathOptSymbolicAD
using Dionysos
const DI = Dionysos
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const AB = DI.Optim.Abstraction

const WALLS =
    [([1.0, 0.0], [1.2, 9.0]), ([2.2, 0.0], [2.4, 5.0]), ([2.2, 6.0], [2.4, 10.0])]

function jacobian_bound(u)
    β = abs(u[1] / cos(atan(tan(u[2]) / 2)))
    return SMatrix{3, 3}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, β, β, 0.0)
end

function solve_with(mode)
    model = Model(Dionysos.Optimizer)
    x_low, x_upp = [0.0, 0.0, -pi - 0.4], [4.0, 10.0, pi + 0.4]
    @variable(model, x_low[i] <= x[i = 1:3] <= x_upp[i], start = [0.4, 0.4, 0.0][i])
    @variable(model, -1 <= u[1:2] <= 1)
    @expression(model, α, atan(tan(u[2]) / 2))
    @constraint(model, ∂(x[1]) == u[1] * cos(α + x[3]) * sec(α))
    @constraint(model, ∂(x[2]) == u[1] * sin(α + x[3]) * sec(α))
    @constraint(model, ∂(x[3]) == u[1] * tan(u[2]))
    @constraint(model, final(x[1]) in MOI.Interval(3.0, 3.6))
    @constraint(model, final(x[2]) in MOI.Interval(0.3, 0.8))
    for (lo, hi) in WALLS
        @constraint(model, x[1:2] ∉ MOI.HyperRectangle(lo, hi))
    end

    set_attribute(model, "jacobian_bound", jacobian_bound)
    set_attribute(model, "approx_mode", mode)          ## the one line that differs
    set_attribute(model, "time_step", 0.3)
    set_attribute(
        model,
        "state_grid",
        MP.GridFree(SVector(0.0, 0.0, 0.0), SVector(0.2, 0.2, 0.2)),
    )
    set_attribute(model, "input_grid", MP.GridFree(SVector(0.0, 0.0), SVector(0.3, 0.3)))
    set_attribute(model, "print_level", 0)

    t = @elapsed optimize!(model)
    abs_sys = get_attribute(model, "abstract_system")
    traj = Dionysos.simulate(model, SVector(0.4, 0.4, 0.0); nsteps = 100)
    return (;
        status = termination_status(model),
        time = t,
        nstates = SY.get_n_state(abs_sys),
        ntrans = SY.get_n_transitions(abs_sys),
        traj = traj,
        problem = get_attribute(model, "concrete_problem"),
    )
end

results = Dict{Symbol, Any}()
for (name, mode) in (
    (:growth, AB.UniformGridAbstraction.GROWTH),
    (:linearized, AB.UniformGridAbstraction.LINEARIZED),
)
    println("— $name —")
    r = solve_with(mode)
    results[name] = r
    println(
        "  ",
        r.status,
        ", ",
        r.nstates,
        " states, ",
        r.ntrans,
        " transitions, ",
        round(r.time; digits = 1),
        " s, ",
        length(ST.states(r.traj)),
        " steps",
    )
end

function panel(r, title)
    fig = plot(; aspect_ratio = :equal, title = title, legend = false)
    plot!(fig, r.problem.system.X; color = :grey, opacity = 1.0, label = "")
    for (lo, hi) in WALLS
        plot!(
            fig,
            LazySets.Hyperrectangle(; low = lo, high = hi);
            color = :black,
            opacity = 1.0,
            label = "",
        )
    end
    plot!(
        fig,
        LazySets.Hyperrectangle(; low = [3.0, 0.3], high = [3.6, 0.8]);
        color = :green,
        opacity = 0.7,
        label = "",
    )
    xs = collect(ST.states(r.traj))
    plot!(
        fig,
        [x[1] for x in xs],
        [x[2] for x in xs];
        color = :orange,
        linewidth = 3,
        label = "",
    )
    return fig
end

g, l = results[:growth], results[:linearized]
fig = plot(
    panel(g, "GROWTH"),
    panel(l, "LINEARIZED");
    layout = (1, 2),
    size = (900, 620),
    dpi = 160,
)
out = joinpath(@__DIR__, "gallery", "solver_swap.png")
savefig(fig, out)
println("saved ", out)
println("GROWTH     : $(g.ntrans) transitions, $(round(g.time; digits=1)) s")
println("LINEARIZED : $(l.ntrans) transitions, $(round(l.time; digits=1)) s")
