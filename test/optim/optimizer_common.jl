module TestOptimizerCommon

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import MathOptInterface as MOI
import LazySets
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

@testset "solution status is reported by the solver itself" begin
    # The status attributes live on the solver, not on the JuMP front-end, so a direct-MOI user
    # gets the same answers. `ẋ = u` on [-1, 1] with a reachable target.
    function reach_optimizer(;
        U = LazySets.Hyperrectangle(; low = SVector(-1.0), high = SVector(1.0)),
        target = LazySets.Hyperrectangle(; low = SVector(-0.5), high = SVector(0.5)),
    )
        f(x, u) = u
        X = LazySets.Hyperrectangle(; low = SVector(-1.0), high = SVector(1.0))
        system =
            MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(f, 1, 1, X, U)
        problem = PR.OptimalControlProblem(
            system,
            LazySets.Hyperrectangle(; low = SVector(-0.9), high = SVector(-0.8)),
            target,
            nothing,
            nothing,
            PR.Infinity(),
        )

        opt = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
        MOI.set(opt, MOI.RawOptimizerAttribute("concrete_problem"), problem)
        MOI.set(opt, MOI.RawOptimizerAttribute("time_step"), 1.0)
        MOI.set(
            opt,
            MOI.RawOptimizerAttribute("approx_mode"),
            AB.UniformGridAbstraction.CENTER_SIMULATION,
        )
        MOI.set(
            opt,
            MOI.RawOptimizerAttribute("state_grid"),
            MP.GridFree(SVector(0.0), SVector(0.25)),
        )
        MOI.set(
            opt,
            MOI.RawOptimizerAttribute("input_grid"),
            MP.GridFree(SVector(0.0), SVector(0.5)),
        )
        MOI.set(opt, MOI.RawOptimizerAttribute("print_level"), 0)
        return opt
    end

    solvable = reach_optimizer()
    @test MOI.get(solvable, MOI.TerminationStatus()) == MOI.OPTIMIZE_NOT_CALLED
    @test MOI.get(solvable, MOI.PrimalStatus()) == MOI.NO_SOLUTION
    @test MOI.get(solvable, MOI.ResultCount()) == 0

    MOI.optimize!(solvable)
    @test MOI.get(solvable, MOI.TerminationStatus()) == MOI.OPTIMAL
    @test MOI.get(solvable, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
    @test MOI.get(solvable, MOI.ResultCount()) == 1
    @test MOI.get(solvable, MOI.DualStatus()) == MOI.NO_SOLUTION
    @test occursin("controller", MOI.get(solvable, MOI.RawStatusString()))

    # Every admissible input drives the state *down*, so a target above the initial set cannot
    # be reached. Failure is LOCALLY_INFEASIBLE — never INFEASIBLE, which would claim more than
    # a sound-but-incomplete abstraction can prove.
    unreachable = reach_optimizer(;
        U = LazySets.Hyperrectangle(; low = SVector(-1.0), high = SVector(-0.5)),
        target = LazySets.Hyperrectangle(; low = SVector(0.7), high = SVector(0.9)),
    )
    MOI.optimize!(unreachable)
    @test MOI.get(unreachable, MOI.TerminationStatus()) == MOI.LOCALLY_INFEASIBLE
    @test MOI.get(unreachable, MOI.PrimalStatus()) == MOI.NO_SOLUTION
    @test MOI.get(unreachable, MOI.ResultCount()) == 0
    @test occursin("finer", MOI.get(unreachable, MOI.RawStatusString()))

    # An abstraction-only problem has nothing to fail: building it *is* the task.
    abstraction = reach_optimizer()
    MOI.set(
        abstraction,
        MOI.RawOptimizerAttribute("concrete_problem"),
        PR.AlternatingSimulationProblem(single_integrator(), nothing),
    )
    @test MOI.get(abstraction, MOI.TerminationStatus()) == MOI.OPTIMIZE_NOT_CALLED
    MOI.optimize!(abstraction)
    @test MOI.get(abstraction, MOI.TerminationStatus()) == MOI.OPTIMAL
    @test MOI.get(abstraction, MOI.PrimalStatus()) == MOI.NO_SOLUTION
    @test occursin("no control objective", MOI.get(abstraction, MOI.RawStatusString()))
end

end # module
