module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import MathOptInterface as MOI

# Fast end-to-end check of the UniformGridAbstraction reach-and-stay optimizer
# (`OptimizerReachAndStayProblem`), which the lifted-control template drives.
@testset "ControllerReachAndStay" begin
    # Tiny 1D single integrator: ẋ = u. Drive into the target and stay (u = 0).
    # Use center-simulation (deterministic transitions) so the coarse abstraction
    # is exact and the test stays fast and robust.
    F_sys(x, u) = SVector(u[1])

    _X_ = UT.box(SVector(-2.0), SVector(2.0))
    _U_ = UT.box(SVector(-1.0), SVector(1.0))

    concrete_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        F_sys,
        1,
        1,
        _X_,
        _U_,
    )

    _I_ = UT.box(SVector(-1.9), SVector(-1.7))   # start near the left edge
    _T_ = UT.box(SVector(-0.4), SVector(0.4))    # target around the origin
    _S_ = _X_                                    # safe set = whole domain

    concrete_problem = PR.ReachAndStayProblem(concrete_system, _I_, _T_, _S_)

    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("state_grid"),
        MP.GridFree(SVector(0.0), SVector(0.2)),
    )
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("input_grid"),
        MP.GridFree(SVector(0.0), SVector(0.5)),
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 1.0)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("approx_mode"),
        AB.UniformGridAbstraction.CENTER_SIMULATION,
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 0)

    MOI.optimize!(optimizer)

    winning_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("winning_set"))
    concrete_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
    success = MOI.get(optimizer, MOI.RawOptimizerAttribute("success"))
    solve_time = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem_time_sec"))

    @test winning_set !== nothing
    @test concrete_controller !== nothing
    @test success == true                         # initial set is inside the winning set
    @test solve_time >= 0.0

    # The controller must produce an admissible input from a state in the target.
    x_in_target = SVector(0.0)
    @test ST.is_defined(concrete_controller, nothing, x_in_target)
    u = ST.output_control(concrete_controller, nothing, x_in_target)
    @test u !== nothing
end

# `stay_on_first_entry` swaps the nested μ/ν persistence fixed point for safety-then-
# reachability: the invariant core of the target is computed once, and the run is steered into
# it. The property gained is that the run never leaves the target after arriving; the price is
# a winning set that can only shrink.
@testset "ControllerReachAndStay: stay_on_first_entry" begin
    F_sys(x, u) = SVector(u[1])
    _X_ = UT.box(SVector(-2.0), SVector(2.0))
    _U_ = UT.box(SVector(-1.0), SVector(1.0))
    concrete_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        F_sys,
        1,
        1,
        _X_,
        _U_,
    )
    _I_ = UT.box(SVector(-1.9), SVector(-1.7))
    _T_ = UT.box(SVector(-0.4), SVector(0.4))

    function solve(stay)
        problem = PR.ReachAndStayProblem(
            concrete_system,
            _I_,
            _T_,
            _X_;
            stay_on_first_entry = stay,
        )
        opt = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
        MOI.set(opt, MOI.RawOptimizerAttribute("concrete_problem"), problem)
        MOI.set(
            opt,
            MOI.RawOptimizerAttribute("state_grid"),
            MP.GridFree(SVector(0.0), SVector(0.2)),
        )
        MOI.set(
            opt,
            MOI.RawOptimizerAttribute("input_grid"),
            MP.GridFree(SVector(0.0), SVector(0.5)),
        )
        MOI.set(opt, MOI.RawOptimizerAttribute("time_step"), 1.0)
        MOI.set(
            opt,
            MOI.RawOptimizerAttribute("approx_mode"),
            AB.UniformGridAbstraction.CENTER_SIMULATION,
        )
        MOI.set(opt, MOI.RawOptimizerAttribute("print_level"), 0)
        MOI.optimize!(opt)
        return opt
    end

    loose, strict = solve(false), solve(true)

    @test MOI.get(strict, MOI.RawOptimizerAttribute("success"))

    # The stronger requirement cannot certify more states than the weaker one. Both runs build
    # the same abstraction from the same grid, so the abstract state ids are comparable.
    function win(o)
        ws = MOI.get(o, MOI.RawOptimizerAttribute("winning_set"))
        abs_sys = MOI.get(o, MOI.RawOptimizerAttribute("abstract_system"))
        return Set(MP.enum_states(ws, SY.get_state_mapping(abs_sys)))
    end
    @test issubset(win(strict), win(loose))

    # The property itself: from the initial state the run enters the target and never leaves.
    controller = MOI.get(strict, MOI.RawOptimizerAttribute("concrete_controller"))
    system = MOI.get(strict, MOI.RawOptimizerAttribute("discrete_time_system"))
    traj = ST.get_closed_loop_trajectory(
        system,
        controller,
        SVector(-1.8),
        40;
        stopping = _ -> false,
    )
    xs = ST.states(traj)
    k = findfirst(x -> x ∈ _T_, xs)
    @test k !== nothing
    @test all(x -> x ∈ _T_, xs[k:end])
end

end # module TestMain
