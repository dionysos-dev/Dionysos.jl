# Regression test for the USER_DEFINED approximation mode: the optimizer used to
# reference non-existent type names (`ST.ContinuousTimeOverapproximationMap`), so any
# run with a custom overapproximation map failed with an UndefVarError.
module TestMain

using Test
using StaticArrays
using MathematicalSystems
import LazySets
import MathOptInterface as MOI
using Dionysos

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

sleep(0.1)
println("Started USER_DEFINED approximation-mode tests")

_X_ = UT.box(SVector(-2.0, -2.0), SVector(2.0, 2.0))
_U_ = UT.box(SVector(-1.0), SVector(1.0))

state_grid = MP.GridFree(SVector(0.0, 0.0), SVector(0.5, 0.5))
input_grid = MP.GridFree(SVector(0.0), SVector(0.5))

function build_abstraction(concrete_system, overapproximation_map; time_step = nothing)
    problem = DI.Problem.AlternatingSimulationProblem(concrete_system, concrete_system.X)

    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("approx_mode"),
        AB.UniformGridAbstraction.USER_DEFINED,
    )
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("overapproximation_map"),
        overapproximation_map,
    )
    if time_step !== nothing
        MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), time_step)
    end
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.optimize!(optimizer)
    return optimizer
end

@testset "USER_DEFINED overapproximation (continuous-time system)" begin
    F_sys(x, u) = SVector(u[1], -x[1])
    concrete_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        F_sys,
        2,
        1,
        _X_,
        _U_,
    )

    tstep = 0.3
    # Exact flow of the center, inflated by a fixed margin: a valid custom map.
    over_map = (rect, u, h) -> begin
        x = LazySets.center(rect)
        r = LazySets.radius_hyperrectangle(rect)
        Fx = ST.runge_kutta4(F_sys, x, u, h, 5)
        return LazySets.Hyperrectangle(Fx, r .+ 0.1)
    end

    optimizer = build_abstraction(concrete_system, over_map; time_step = tstep)

    abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    @test abstract_system !== nothing
    @test SY.ntransitions(SY.get_automaton(abstract_system)) > 0

    approx =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system_approximation"))
    @test approx isa ST.DiscreteTimeOverApproximationMap
end

@testset "USER_DEFINED overapproximation (discrete-time system)" begin
    fd(x, u) = SVector(0.9 * x[1] + 0.1 * u[1], 0.9 * x[2])
    concrete_system =
        MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem(fd, 2, 1, _X_, _U_)

    over_map = (rect, u) -> begin
        x = LazySets.center(rect)
        r = LazySets.radius_hyperrectangle(rect)
        return LazySets.Hyperrectangle(fd(x, u), r .+ 0.05)
    end

    optimizer = build_abstraction(concrete_system, over_map)

    abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    @test abstract_system !== nothing
    @test SY.ntransitions(SY.get_automaton(abstract_system)) > 0

    approx =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system_approximation"))
    @test approx isa ST.DiscreteTimeOverApproximationMap
end

sleep(0.1)
println("End test")

end # module TestMain
