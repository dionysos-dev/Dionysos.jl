module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import LazySets
import LinearAlgebra as LA
import MathOptInterface as MOI
using Symbolics
import Clarabel

include("../../problems/ToyProblem/toy_problem.jl")

const EB = AB.EllipsoidalBackwardTrajectoryCertifier

# Contract test for the trajectory certifiers. Builds a tiny reach trajectory with the
# optimizer-based generator, then certifies it with both the uniform-grid certifier and
# the ellipsoidal backward certifier. Adapted from the certification case-study script.
@testset "trajectory certifiers" begin
    Δt = 0.3
    _X_ = UT.box(SVector(-2.0, -2.0), SVector(4.0, 4.0))
    _U_ = UT.box(SVector(-1.0, -1.0), SVector(1.0, 1.0))
    concrete_system = ToyProblem.system(; _X_ = _X_, _U_ = _U_)
    jacobian_bound = ToyProblem.jacobian_bound()

    build_optimizer = function ()
        opt = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
        MOI.set(
            opt,
            MOI.RawOptimizerAttribute("concrete_problem"),
            PR.AlternatingSimulationProblem(concrete_system, concrete_system.X),
        )
        MOI.set(
            opt,
            MOI.RawOptimizerAttribute("state_grid"),
            MP.GridFree(SVector(-2.0, -2.0), SVector(0.2, 0.2)),
        )
        MOI.set(
            opt,
            MOI.RawOptimizerAttribute("input_grid"),
            MP.GridFree(SVector(-1.0, -1.0), SVector(0.5, 0.5)),
        )
        MOI.set(opt, MOI.RawOptimizerAttribute("time_step"), Δt)
        MOI.set(
            opt,
            MOI.RawOptimizerAttribute("approx_mode"),
            AB.UniformGridAbstraction.GROWTH,
        )
        MOI.set(opt, MOI.RawOptimizerAttribute("jacobian_bound"), jacobian_bound)
        MOI.set(opt, MOI.RawOptimizerAttribute("n_samples"), 1)
        MOI.set(opt, MOI.Silent(), true)
        MOI.optimize!(opt)
        return opt
    end

    _I_ = UT.box(SVector(-1.7, -1.7), SVector(-1.6, -1.6))
    target_set = UT.set_union([
        UT.box(SVector(-1.0, 3.0), SVector(-0.3, 3.7)),
        UT.box(SVector(1.0, 2.0), SVector(3.0, 3.7)),
    ])
    concrete_problem =
        PR.OptimalControlProblem(concrete_system, _I_, target_set, nothing, nothing, 20)

    # A reach trajectory to certify.
    optimizer = build_optimizer()
    gen = AB.OptimizerTrajectoryGenerator.TrajectoryGenerator(
        optimizer;
        initial_state = SVector(-1.65, -1.65),
        concrete = false,
        nstep = 20,
    )
    AB.set_problem!(gen, concrete_problem)
    AB.generate!(gen)
    traj = AB.get_trajectory(gen)
    @test traj !== nothing

    # 1) Uniform-grid trajectory certifier: builds a tube around the trajectory and
    #    re-synthesizes on the restricted region.
    ug_optimizer = build_optimizer()
    ug_cert = AB.UniformGridTrajectoryCertifier.TrajectoryCertifier(;
        optimizer = ug_optimizer,
        radius = 0.3,
        n_between = 0,
    )
    AB.set_problem!(ug_cert, concrete_problem)
    AB.set_trajectory!(ug_cert, traj)
    AB.certify!(ug_cert)
    @test AB.get_success(ug_cert) isa Bool
    @test AB.get_solve_time(ug_cert) >= 0.0

    # 2) Ellipsoidal backward certifier: SDP-based tube certification.
    Symbolics.@variables x[1:2]
    Symbolics.@variables u[1:2]
    Symbolics.@variables w[1:2]
    fsymbolic = [x[1] + Δt * (u[1] + w[1]), x[2] + Δt * (u[2] + w[2])]
    Wformat = UT.box(SVector(0.0, 0.0), SVector(0.0, 0.0))
    provider = ST.SymbolicAffineApproximationProvider(
        fsymbolic,
        collect(x),
        collect(u),
        collect(w),
        [0.0, 0.0],
        ST.format_input_set(_U_),
        ST.format_noise_set(Wformat),
    )
    adaptive_opts = EB.AdaptiveLinearizationBoxOptions(
        false,
        [0.2, 0.2],
        [0.01, 0.01],
        [1.0, 1.0],
        [0.5, 0.5],
        [0.01, 0.01],
        [1.0, 1.0],
        2.0,
        1.05,
        1,
        1e-8,
        false,
        [1.0],
        :first_consistent,
        true,
    )
    ellip_opts = EB.EllipsoidalBackwardOptions(;
        maxδx = 30,
        maxδu = 1.0,
        λ = 0.05,
        terminal_shape = Matrix{Float64}(LA.I, 2, 2) * 0.5^2,
        linearization_δx = [0.2, 0.2],
        linearization_δu = [0.5, 0.5],
        adaptive_boxes = adaptive_opts,
        use_log_det = false,
    )
    eb_cert = EB.TrajectoryCertifier(provider, Clarabel.Optimizer, ellip_opts)
    AB.set_problem!(eb_cert, concrete_problem)
    AB.set_trajectory!(eb_cert, traj)
    AB.certify!(eb_cert)
    @test AB.get_success(eb_cert) isa Bool
    @test EB.get_result(eb_cert) !== nothing
end

end # module TestMain
