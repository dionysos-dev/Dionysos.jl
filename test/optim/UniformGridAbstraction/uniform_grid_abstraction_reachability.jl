module TestMain

using Test
using StaticArrays, MathematicalSystems, Plots
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

using JuMP
import MathOptInterface as MOI

include("../../../problems/PathPlanning/path_planning.jl")

@testset "UniformGridAbstraction PathPlanning (simple)" begin
    concrete_problem = PathPlanning.problem(; simple = true)
    concrete_system = concrete_problem.system

    # grids
    x0_grid = SVector(0.0, 0.0, 0.0)
    hx = SVector(0.2, 0.2, 0.2)
    state_grid = MP.GridFree(x0_grid, hx)

    u0_grid = SVector(0.0, 0.0)
    hu = SVector(0.3, 0.3)
    input_grid = MP.GridFree(u0_grid, hu)

    # solve abstraction
    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.3)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("jacobian_bound"),
        PathPlanning.jacobian_bound(),
    )

    MOI.optimize!(optimizer)

    abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
    abstract_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
    concrete_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
    discrete_time_system =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system"))

    # --- core sanity checks ---
    @test abstract_system !== nothing
    @test abstract_problem !== nothing
    @test abstract_controller !== nothing
    @test concrete_controller !== nothing
    @test discrete_time_system !== nothing

    # Mapping sanity
    Xmap = SY.get_state_mapping(abstract_system)
    @test Xmap !== nothing
    @test MP.get_n_state(Xmap) > 0

    # closed-loop simulation
    nstep = 100
    reached(x) = (x ∈ concrete_problem.target_set)

    x0 = SVector(0.4, 0.4, 0.0)

    # controller should return something at x0
    u0 = ST.output_control(concrete_controller, ST.initial_state(concrete_controller), x0)
    @test u0 !== nothing

    traj = ST.get_closed_loop_trajectory(
        discrete_time_system,
        concrete_controller,
        x0,
        nstep;
        stopping = reached,
    )

    x_traj = traj.x
    u_traj = traj.u

    # --- trajectory checks ---
    xs = collect(ST.enum_elems(x_traj))
    us = collect(ST.enum_elems(u_traj))

    @test !isempty(xs)
    @test length(xs) <= nstep + 1
    @test length(us) == length(xs) - 1 || length(us) == length(xs)

    # Must reach target at some time (stop condition)
    @test any(reached, xs)

    # Stronger: last state is in target if stopping works by termination
    @test reached(xs[end])

    no_plot = false
    @static if get(ENV, "CI", "false") == "false" &&
               (isdefined(@__MODULE__, :no_plot) && no_plot == false)
        fig = plot(; aspect_ratio = :equal)
        plot!(concrete_system.X; color = :yellow, opacity = 0.5)
        plot!(abstract_system; color = :blue, opacity = 0.5)
        plot!(concrete_problem.initial_set; color = :green, opacity = 0.2)
        plot!(concrete_problem.target_set; dims = [1, 2], color = :red, opacity = 0.2)

        plot!(
            (
                SY.get_state_set_from_states(abstract_system, abstract_problem.initial_set),
                Xmap,
            );
            color = :green,
        )
        plot!(
            (
                SY.get_state_set_from_states(abstract_system, abstract_problem.target_set),
                Xmap,
            );
            color = :red,
        )
        plot!(x_traj; ms = 0.5)

        @test isa(fig, Plots.Plot{Plots.GRBackend})
    end
end

end # module TestMain
