module TestMain

using Test
using StaticArrays, MathematicalSystems, Plots
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

using JuMP
import MathOptInterface as MOI

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "dc_dc.jl"))

@testset "UniformGridAbstraction DCDC (growth mode)" begin
    concrete_problem = DCDC.problem()
    concrete_system = concrete_problem.system

    # grids
    x0_grid = SVector(0.0, 0.0)
    hx = SVector(2.0 / 4.0e3, 2.0 / 4.0e3)
    state_grid = MP.GridFree(x0_grid, hx)

    u0_grid = SVector(1)
    hu = SVector(1)
    input_grid = MP.GridFree(u0_grid, hu)

    # solve abstraction
    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("approx_mode"),
        Dionysos.Optim.Abstraction.UniformGridAbstraction.GROWTH,
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), DCDC.jacobian_bound())
    MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.5)

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

    # Optional: mapping sanity
    Xmap = SY.get_state_mapping(abstract_system)
    @test Xmap !== nothing
    @test MP.get_n_state(Xmap) > 0

    # --- closed-loop trajectory sanity ---
    nstep = 300
    x0 = SVector(1.2, 5.6)

    # controller should return an input at x0
    u0 = ST.output_control(concrete_controller, ST.initial_state(concrete_controller), x0)
    @test u0 !== nothing

    x_traj, u_traj =
        ST.get_closed_loop_trajectory(discrete_time_system, concrete_controller, x0, nstep)

    xs = collect(ST.enum_elems(x_traj))
    us = collect(ST.enum_elems(u_traj))

    @test !isempty(xs)
    @test length(xs) <= nstep + 1
    @test length(us) == length(xs) - 1 || length(us) == length(xs)

    no_plot = false
    @static if get(ENV, "CI", "false") == "false" &&
               (isdefined(@__MODULE__, :no_plot) && no_plot == false)
        fig = plot(; aspect_ratio = :equal)
        plot!(concrete_system.X)
        plot!(x_traj)
        @test isa(fig, Plots.Plot{Plots.GRBackend})
    end
end

end # module TestMain
