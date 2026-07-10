module TestOptimizerCommon

using Test

import MathOptInterface as MOI
import Dionysos
const PR = Dionysos.Problem
const OP = Dionysos.Optim
const UGA = OP.Abstraction.UniformGridAbstraction

@testset "AbstractDionysosOptimizer (leaf contract)" begin
    leaf = MOI.instantiate(UGA.OptimizerSafetyProblem)

    @test MOI.supports(leaf, MOI.RawOptimizerAttribute("anything"))

    MOI.set(leaf, MOI.RawOptimizerAttribute("print_level"), 2)
    @test MOI.get(leaf, MOI.RawOptimizerAttribute("print_level")) == 2

    @test_throws MOI.UnsupportedAttribute MOI.set(
        leaf,
        MOI.RawOptimizerAttribute("nonexistent_attr"),
        1,
    )
    @test_throws MOI.UnsupportedAttribute MOI.get(
        leaf,
        MOI.RawOptimizerAttribute("nonexistent_attr"),
    )

    # MOI.Silent maps onto the print_level field.
    @test MOI.supports(leaf, MOI.Silent())
    MOI.set(leaf, MOI.Silent(), true)
    @test MOI.get(leaf, MOI.RawOptimizerAttribute("print_level")) == 0
    @test MOI.get(leaf, MOI.Silent())
    MOI.set(leaf, MOI.Silent(), false)
    @test MOI.get(leaf, MOI.RawOptimizerAttribute("print_level")) == 1
    @test !MOI.get(leaf, MOI.Silent())

    # SolveTimeSec: explicit override backed by abstract_problem_time_sec.
    @test MOI.get(leaf, MOI.SolveTimeSec()) == 0.0
end

@testset "CompositeDionysosOptimizer (forwarding + problem dispatch)" begin
    opt = MOI.instantiate(UGA.Optimizer)
    @test MOI.is_empty(opt)

    # Setting a sub-solver attribute lazily creates the abstraction solver and forwards.
    MOI.set(opt, MOI.RawOptimizerAttribute("time_step"), 0.5)
    @test !MOI.is_empty(opt)
    @test opt.abstraction_solver.time_step == 0.5
    @test MOI.get(opt, MOI.RawOptimizerAttribute("time_step")) == 0.5

    # Own fields win over sub-solver fields with the same name.
    MOI.set(opt, MOI.RawOptimizerAttribute("print_level"), 0)
    @test opt.print_level == 0
    @test opt.abstraction_solver.print_level == 1

    @test_throws MOI.UnsupportedAttribute MOI.set(
        opt,
        MOI.RawOptimizerAttribute("nonexistent_attr"),
        1,
    )
    @test_throws MOI.UnsupportedAttribute MOI.get(
        opt,
        MOI.RawOptimizerAttribute("nonexistent_attr"),
    )

    # "concrete_problem" dispatches on the problem type to pick the control solver.
    problem = PR.SafetyProblem(nothing, nothing, nothing)
    MOI.set(opt, MOI.RawOptimizerAttribute("concrete_problem"), problem)
    @test opt.control_solver isa UGA.OptimizerSafetyProblem
    @test MOI.get(opt, MOI.RawOptimizerAttribute("concrete_problem")) === problem
    @test opt.abstraction_solver.alternating_simulation_problem isa
          PR.AlternatingSimulationProblem

    # The abstraction, once attached, is kept when the control task changes.
    asp = opt.abstraction_solver.alternating_simulation_problem
    reach = PR.OptimalControlProblem(nothing, nothing, nothing, nothing, nothing, 10)
    MOI.set(opt, MOI.RawOptimizerAttribute("concrete_problem"), reach)
    @test opt.control_solver isa UGA.OptimizerOptimalControlProblem
    @test opt.abstraction_solver.alternating_simulation_problem === asp

    # An unsupported problem type errors instead of silently misconfiguring.
    @test_throws ErrorException MOI.set(
        opt,
        MOI.RawOptimizerAttribute("concrete_problem"),
        42,
    )

    # An AlternatingSimulationProblem selects the abstraction-only path.
    MOI.set(
        opt,
        MOI.RawOptimizerAttribute("concrete_problem"),
        PR.AlternatingSimulationProblem(nothing, nothing),
    )
    @test opt.control_solver === nothing

    # The composite reads SolveTimeSec and Silent from the shared base.
    @test MOI.get(opt, MOI.SolveTimeSec()) == 0.0
    MOI.set(opt, MOI.Silent(), true)
    @test opt.print_level == 0
end

end # module
