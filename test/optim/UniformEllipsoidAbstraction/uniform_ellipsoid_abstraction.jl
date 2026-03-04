module TestMain

using Test
using StaticArrays
using LinearAlgebra
using JuMP
using Clarabel
import MathOptInterface as MOI
import CDDLib
import Random
using Plots

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

# ------------------------------------------------------------------
# Optional plotting gate
# ------------------------------------------------------------------
const _NO_PLOT = get(ENV, "CI", "false") == "true"

# ------------------------------------------------------------------
# Load PWA system
# ------------------------------------------------------------------
include("../../../problems/pwa_sys.jl")

@testset "UniformEllipsoidAbstraction on PWAsys (end-to-end)" begin
    Random.seed!(0)

    # ----------------------------
    # Problem + system
    # ----------------------------
    lib = CDDLib.Library()

    Usz = 70
    Wsz = 3
    dt = 0.01

    concrete_problem =
        PWAsys.problem(; lib = lib, dt = dt, Usz = Usz, Wsz = Wsz, simple = true)
    concrete_system = concrete_problem.system

    # ----------------------------
    # Build EmptyProblem (abstraction-only phase)
    # ----------------------------
    @test haskey(concrete_system.ext, :X)
    empty_problem = PR.EmptyProblem(concrete_system, concrete_system.ext[:X])

    # ----------------------------
    # Ellipsoidal grid definition
    # ----------------------------
    n_step = 3
    origin = SVector(0.0, 0.0)
    h = SVector(1.0 / n_step, 1.0 / n_step)

    nx = size(concrete_system.resetmaps[1].A, 1)
    @test nx == 2  # this test is written for the 2D PWAsys instance

    P = (1 / nx) * diagm((h ./ 2) .^ (-2))
    Pm = P
    R = h ./ 2

    state_grid = MP.GridEllipsoidalRectangular(origin, h, P)

    # ----------------------------
    # SDP solver
    # ----------------------------
    opt_sdp = optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => true)

    # ----------------------------
    # Instantiate abstraction optimizer
    # ----------------------------
    optimizer = MOI.instantiate(AB.UniformEllipsoidAbstraction.Optimizer)

    MOI.set(optimizer, MOI.RawOptimizerAttribute("empty_problem"), empty_problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("incl_mode"), MP.INNER)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("P"), P)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("Pm"), Pm)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("R"), R)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("sdp_solver"), opt_sdp)

    # Q_aug
    nx, nu = 2, 2
    naug = nx + nu + 1
    Q_aug = Matrix{Float64}(I, naug, naug) * (dt^2)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("Q_aug"), Q_aug)

    # ----------------------------
    # Phase 1: Build abstraction only
    # ----------------------------
    MOI.optimize!(optimizer)

    abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    @test abstract_system !== nothing

    nX = SY.get_n_state(abstract_system)
    @test nX > 0

    @test length(optimizer.abstraction_solver.transitionCost) ≥ 0

    # ----------------------------
    # Phase 2: Solve control problem on the abstraction
    # ----------------------------
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    MOI.optimize!(optimizer)

    abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
    abstract_system2 = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    concrete_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
    concrete_value_function =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_value_function"))
    abstract_value_function =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_value_function"))

    @test abstract_problem !== nothing
    @test abstract_system2 !== nothing
    @test concrete_controller !== nothing
    @test concrete_value_function !== nothing
    @test abstract_value_function !== nothing

    # ----------------------------
    # Closed-loop simulation
    # ----------------------------

    # Return PWA mode for a given x
    get_mode(x) = findfirst(m -> (x ∈ m.X), concrete_system.resetmaps)

    function f_eval1(x, u)
        # choose the state among intersected cells with minimal abstract value
        states = SY.get_abstract_states(abstract_system2, x)
        @test !isempty(states)  # should map to at least one cell (if within X)

        smin = argmin(s -> abstract_value_function(s), states)
        c = SY.get_concrete_state(abstract_system2, smin)
        m = get_mode(c)

        W = concrete_system.ext[:W]
        w = (2 * (Random.rand(2) .^ (1 / 4)) .- 1) .* W[:, 1]

        return concrete_system.resetmaps[m].A * x +
               concrete_system.resetmaps[m].B * u +
               concrete_system.resetmaps[m].c +
               w
    end

    cost_eval(x, u) = UT.function_value(concrete_problem.transition_cost[1][1], x, u)

    nstep = typeof(concrete_problem.time) == PR.Infinity ? 100 : concrete_problem.time

    function reached(x)
        states = SY.get_abstract_states(abstract_system2, x)
        return !isempty(states ∩ abstract_problem.target_set)
    end

    x0 = concrete_problem.initial_set

    x_traj, u_traj = ST.get_closed_loop_trajectory(
        concrete_system,
        concrete_controller,
        x0,
        nstep;
        stopping = reached,
        f_map_override = f_eval1,
    )

    xs = collect(ST.enum_elems(x_traj))
    us = collect(ST.enum_elems(u_traj))

    @test !isempty(xs)
    @test length(xs) ≤ nstep + 1
    @test length(us) == length(xs) - 1 || length(us) == length(xs)

    c_traj, cost_true = ST.get_cost_trajectory(x_traj, u_traj, cost_eval)
    cost_bound = concrete_value_function(x0)

    # ----------------------------
    # Assertions (robust ones)
    # ----------------------------
    @test isfinite(cost_true)
    @test isfinite(cost_bound)
    @test cost_true ≤ cost_bound + 1e-8
    @test any(reached, xs)

    # ----------------------------
    # Plotting
    # ----------------------------
    if !_NO_PLOT
        fig = plot(; aspect_ratio = :equal)
        X = concrete_system.ext[:X]
        plot!(X; color = :grey, opacity = 1.0, label = "")
        plot!(abstract_system2; value_function = abstract_value_function)
        Xmap = SY.get_state_mapping(abstract_system2)
        plot!(
            (
                SY.get_state_set_from_states(
                    abstract_system2,
                    abstract_problem.initial_set,
                ),
                Xmap,
            );
            color = :green,
            efficient = false,
            opacity = 0.6,
        )
        plot!(
            (
                SY.get_state_set_from_states(abstract_system2, abstract_problem.target_set),
                Xmap,
            );
            color = :red,
            efficient = false,
            opacity = 0.6,
        )
        plot!(UT.DrawPoint(concrete_problem.initial_set); color = :green, opacity = 1.0)
        plot!(UT.DrawPoint(concrete_problem.target_set); color = :red, opacity = 1.0)
        plot!(x_traj; ms = 2.0, arrows = false, color = :blue)
        @test isa(fig, Plots.Plot{Plots.GRBackend})
    end
end

end # module TestMain
