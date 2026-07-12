module TestClockLiftedContinuous

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

using MathOptInterface
const MOI = MathOptInterface

@testset "Continuous quantified-time via the reusable ClockLift" begin
    sys = single_integrator(; n = 1, xbound = 2.0, ubound = 1.0)   # ẋ = u

    # Specification over (x, t): reach x ∈ [0.5, 2] within t ∈ [1, 2], from x ≈ 0 at t = 0.
    initial_set = UT.box([-0.1, 0.0], [0.1, 0.0])
    target_set = UT.box([0.5, 1.0], [2.0, 2.0])
    concrete_problem = PR.OptimalControlProblem(
        sys,
        initial_set,
        target_set,
        nothing,
        nothing,
        PR.Infinity(),
    )

    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("state_grid"),
        MP.GridFree(SVector(0.0), SVector(0.25)),
    )
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("input_grid"),
        MP.GridFree(SVector(0.0), SVector(0.5)),
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.5)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("approx_mode"),
        AB.UniformGridAbstraction.GROWTH,
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), u -> SMatrix{1, 1}(0.0))
    # The reusable time lift: the SAME ClockLift the hybrid solver applies per mode,
    # here wrapping a plain continuous abstraction into an (x, t) model.
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("clock"),
        SY.ClockAbstraction(UT.box([0.0], [2.0]), 0.5),   # tsteps [0,.5,1,1.5,2]
    )

    MOI.optimize!(optimizer)

    abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    concrete_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

    # The continuous abstraction is now an (x, t) product built by the clock lift.
    @test abstract_system isa SY.ClockLiftedSymbolicModel
    @test SY.get_state_dim(abstract_system) == 2
    q1 = first(SY.enum_states(abstract_system))
    @test length(SY.get_concrete_state(abstract_system, q1)) == 2
    @test SY.get_n_state(abstract_system) == SY.get_n_state(abstract_system.base) * 5

    # Synthesis over (x, t) succeeded, and the concretized controller is queryable
    # at an (x, t) point.
    @test MOI.get(optimizer, MOI.RawOptimizerAttribute("success"))
    @test !isnothing(concrete_controller)
    @test ST.output_control(concrete_controller, nothing, SVector(0.0, 0.0)) !== nothing
end

end # module
