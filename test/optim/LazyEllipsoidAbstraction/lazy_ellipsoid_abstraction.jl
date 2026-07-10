module TestMain

using Test
using StaticArrays
using Random
using MathematicalSystems, HybridSystems
using JuMP
using Clarabel
import MathOptInterface as MOI
using Plots

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const OP = DI.Optim
const AB = OP.Abstraction

# Don’t plot on CI (keeps tests fast + headless-safe)
const _NO_PLOT = get(ENV, "CI", "false") == "true"

include("../../../problems/non_linear.jl")

@testset "LazyEllipsoidsAbstraction (NonLinear) end-to-end" begin
    # Deterministic randomness
    Random.seed!(0)

    concrete_problem = NonLinear.problem()
    concrete_system = concrete_problem.system

    sdp_opt = optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => true)

    # Parameters (kept identical to your script)
    maxδx = 100
    maxδu = 10 * 2
    λ = 0.01
    k1 = 1
    k2 = 1
    RRTstar = false
    continues = false
    maxIter = 100

    optimizer = MOI.instantiate(AB.LazyEllipsoidsAbstraction.Optimizer)
    AB.LazyEllipsoidsAbstraction.set_optimizer!(
        optimizer,
        concrete_problem,
        sdp_opt,
        maxδx,
        maxδu,
        λ,
        k1,
        k2,
        RRTstar,
        continues,
        maxIter,
    )

    # Build abstraction + solve
    MOI.optimize!(optimizer)

    abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
    concrete_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
    abstract_lyap_fun = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_lyap_fun"))
    concrete_lyap_fun = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_lyap_fun"))

    @test abstract_system !== nothing
    @test abstract_problem !== nothing
    @test concrete_controller !== nothing
    @test abstract_lyap_fun !== nothing
    @test concrete_lyap_fun !== nothing

    # ----------------------------
    # Simulation
    # ----------------------------
    cost_eval(x, u) = concrete_problem.transition_cost(x, u)
    reached(x) = x ∈ concrete_problem.target_set
    nstep = typeof(concrete_problem.time) == PR.Infinity ? 100 : concrete_problem.time

    # initial state: the center of the (ellipsoidal) initial set
    x0 = UT.get_center(concrete_problem.initial_set)

    x_traj, u_traj = ST.get_closed_loop_trajectory(
        concrete_system,
        concrete_controller,
        x0,
        nstep;
        stopping = reached,
        f_map_override = (x, u) -> concrete_system.f_eval(x, u, [0, 0]),
    )

    xs = collect(ST.enum_elems(x_traj))
    @test !isempty(xs)
    @test any(reached, xs)

    _, cost_true = ST.get_cost_trajectory(x_traj, u_traj, cost_eval)
    cost_bound = concrete_lyap_fun(x0)

    @test isfinite(cost_true)
    @test isfinite(cost_bound)
    @test cost_true <= cost_bound + 1e-9

    # ----------------------------
    # Plots (optional)
    # ----------------------------
    if !_NO_PLOT
        fig1 = plot(;
            aspect_ratio = :equal,
            xtickfontsize = 10,
            ytickfontsize = 10,
            guidefontsize = 16,
            titlefontsize = 14,
        )
        xlabel!("\$x_1\$")
        ylabel!("\$x_2\$")
        title!("Specifications and domains")

        plot!(fig1, concrete_system.X; color = :yellow, opacity = 0.5)
        if hasproperty(concrete_system, :obstacles)
            for obs in concrete_system.obstacles
                plot!(fig1, obs; color = :black)
            end
        end
        plot!(fig1, abstract_system; with_arrows = false, cost = false)
        plot!(fig1, concrete_problem.initial_set; color = :green)
        plot!(fig1, concrete_problem.target_set; color = :red)

        fig2 = plot(;
            aspect_ratio = :equal,
            xtickfontsize = 10,
            ytickfontsize = 10,
            guidefontsize = 16,
            titlefontsize = 14,
        )
        title!("Abstractions")
        plot!(fig2, abstract_system; with_arrows = true, cost = false)

        fig3 = plot(;
            aspect_ratio = :equal,
            xtickfontsize = 10,
            ytickfontsize = 10,
            guidefontsize = 16,
            titlefontsize = 14,
        )
        xlabel!("\$x_1\$")
        ylabel!("\$x_2\$")
        title!("Trajectory and Lyapunov-like Fun.")

        if hasproperty(concrete_system, :obstacles)
            for obs in concrete_system.obstacles
                plot!(fig3, obs; color = :black)
            end
        end
        plot!(fig3, abstract_system; with_arrows = false, cost = true)
        plot!(fig3, x_traj; color = :black)

        @test isa(fig1, Plots.Plot{Plots.GRBackend})
        @test isa(fig2, Plots.Plot{Plots.GRBackend})
        @test isa(fig3, Plots.Plot{Plots.GRBackend})
    end
end

end # module TestMain
