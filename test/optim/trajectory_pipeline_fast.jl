# Fast pipeline coverage (kept OUT of :slow): the driver loop, the adaptive
# linearization-box search, the soundness gates as units, and the small library
# pieces the demos lean on — on an exactly-linear integrator whose SDPs solve in
# milliseconds.

module TestTrajectoryPipelineFast

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

using LinearAlgebra
using Random
import Serialization
using Clarabel, JuMP
import LazySets
import MathOptInterface as MOI
using Symbolics

const EB = AB.EllipsoidalTrajectoryCertifier
const MPPI = AB.MPPITrajectoryGenerator

sdp = optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => true)

# 2-D single integrator ẋ = u — exact under discretization, zero Hessian.
Δt = 0.25
_X_ = LazySets.Hyperrectangle(; low = SVector(-2.0, -2.0), high = SVector(2.0, 2.0))
_U_ = LazySets.Hyperrectangle(; low = SVector(-1.0, -1.0), high = SVector(1.0, 1.0))
_I_ = LazySets.Hyperrectangle(; low = SVector(-1.05, -1.05), high = SVector(-0.95, -0.95))
_T_ = LazySets.Hyperrectangle(; low = SVector(-0.3, -0.3), high = SVector(0.3, 0.3))
sys = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
    (x, u) -> SVector(u[1], u[2]),
    2,
    2,
    _X_,
    _U_,
)
problem =
    PR.discretize_problem(PR.OptimalControlProblem(sys, _I_, _T_, nothing, nothing, 10), Δt)
f = MathematicalSystems.mapping(problem.system)

x0 = SVector(-1.0, -1.0)
us = [SVector(0.5, 0.5) for _ in 1:8]
xs = [x0]
for u in us
    push!(xs, f(xs[end], u))
end
traj = ST.Trajectory(xs; inputs = us)

Symbolics.@variables xv[1:2] uv[1:2] wv[1:2]
fsymbolic = [xv[1] + Δt * (uv[1] + wv[1]), xv[2] + Δt * (uv[2] + wv[2])]
Wset = LazySets.Hyperrectangle(; low = SVector(0.0, 0.0), high = SVector(0.0, 0.0))
provider = ST.SymbolicAffineApproximationProvider(
    fsymbolic,
    collect(xv),
    collect(uv),
    collect(wv),
    [0.0, 0.0],
    ST.format_input_set(_U_),
    ST.format_noise_set(Wset),
)

adaptive(objective) = EB.AdaptiveLinearizationBoxOptions(;
    ΔX_initial = [0.2, 0.2],
    ΔX_min = [0.01, 0.01],
    ΔX_max = [2.0, 2.0],
    ΔU_initial = [0.5, 0.5],
    ΔU_min = [0.01, 0.01],
    ΔU_max = [2.0, 2.0],
    search_scales = [0.5, 1.0, 1.5],
    objective = objective,
)

@testset "adaptive box search certifies (both search objectives)" begin
    for box_objective in (:first_consistent, :max_volume)
        # domain_cap keeps the size-maximizing funnels inside X by construction
        # (without it the maximin chain grows past the domain near the start and
        # the state-domain gate rejects a posteriori — the documented failure
        # mode the cap exists for).
        opts = EB.ChainOptions(;
            maxδx = 5.0,
            maxδu = 2.0,
            terminal_shrink = 0.9,
            linearization_δu = [0.5, 0.5],
            adaptive_boxes = adaptive(box_objective),
            objective = :maximin,
            domain_cap = true,
        )
        cert = EB.BackwardCertifier(provider, sdp, opts)
        AB.set_problem!(cert, problem)
        AB.set_trajectory!(cert, traj)
        AB.certify!(cert)
        res = EB.get_result(cert)
        @test res.success
        @test res.failed_k === nothing
        @test res.terminal_contained == true
        @test length(res.lmi_data.ellipsoids) == length(us) + 1
        @test res.lmi_data.reason === nothing
        @test all(rec.status == :ok for rec in res.steps)
        @test AB.get_controller(cert) isa ST.FunnelController
    end
end

@testset "two-step rescue stays inert on a healthy chain" begin
    opts = EB.ChainOptions(;
        maxδx = 5.0,
        maxδu = 2.0,
        linearization_δu = [0.5, 0.5],
        adaptive_boxes = adaptive(:first_consistent),
        objective = :maximin,
        domain_cap = true,
        two_step_rescue = true,
    )
    cert = EB.BackwardCertifier(provider, sdp, opts)
    AB.set_problem!(cert, problem)
    AB.set_trajectory!(cert, traj)
    AB.certify!(cert)
    @test AB.get_success(cert)
end

@testset "soundness gates as units" begin
    E_in = LazySets.Ellipsoid([0.0, 0.0], Matrix(0.01 * I, 2, 2))
    E_out = LazySets.Ellipsoid([1.95, 0.0], Matrix(0.25 * I, 2, 2))

    rec_ok = EB.StepRecord(1, :ok, 0.0, E_in, nothing, (;))
    rec_out = EB.StepRecord(1, :ok, 0.0, E_out, nothing, (;))

    @test EB.collapse_gate(rec_ok, 0.0) === nothing
    @test EB.collapse_gate(rec_ok, 0.5) isa String            # 0.1-radius < 0.5 floor
    @test EB.state_domain_gate(rec_ok, _X_) === nothing
    @test EB.state_domain_gate(rec_out, _X_) isa String       # sticks out of X

    # Endpoint gate: the terminal/entry ellipsoid never enters a StepRecord.
    @test EB.endpoint_gate(E_in, _X_, 0.0, true, "terminal") === nothing
    @test EB.endpoint_gate(E_out, _X_, 0.0, true, "terminal") isa String
    @test EB.endpoint_gate(E_in, _X_, 0.5, false, "entry") isa String

    # Obstacle domain: provable disjointness required.
    hole =
        LazySets.Hyperrectangle(; low = SVector(-0.05, -0.05), high = SVector(0.05, 0.05))
    Xobst = UT.set_minus(_X_, hole)
    @test EB.state_domain_gate(rec_ok, Xobst) isa String      # overlaps the hole
    @test EB.endpoint_gate(E_out, _X_, 0.0, false, "entry") === nothing  # domain check off
end

@testset "driver: round log, best-so-far, MOI statuses" begin
    gen = MPPI.TrajectoryGenerator(;
        rng = Random.MersenneTwister(1),
        seed_trajectory = traj,
        nstep = 8,
        nsamples = 30,
        niter = 3,
        noise = MPPI.GaussianMPPINoise(SVector(0.1, 0.1)),
        project_input = u -> SVector(clamp(u[1], -1.0, 1.0), clamp(u[2], -1.0, 1.0)),
        cost = AB.CompositeCost(
            AB.ReachObjectiveCost(problem.target_set),
            AB.InputEffortCost(0.01),
        ),
    )
    opts = EB.ChainOptions(;
        maxδx = 5.0,
        maxδu = 2.0,
        linearization_δu = [0.5, 0.5],
        adaptive_boxes = adaptive(:first_consistent),
        objective = :maximin,
        domain_cap = true,
    )
    cert = EB.BackwardCertifier(provider, sdp, opts)

    driver = AB.TrajectoryCertificationOptimizer.Optimizer(gen, cert)
    @test MOI.get(driver, MOI.TerminationStatus()) == MOI.OPTIMIZE_NOT_CALLED
    MOI.set(driver, MOI.RawOptimizerAttribute("concrete_problem"), problem)
    MOI.set(driver, MOI.RawOptimizerAttribute("max_rounds"), 2)
    MOI.optimize!(driver)

    rounds = MOI.get(driver, MOI.RawOptimizerAttribute("rounds"))
    round_log = MOI.get(driver, MOI.RawOptimizerAttribute("round_log"))
    @test length(round_log) == rounds
    @test all(haskey(r, :gen_success) && haskey(r, :cert_success) for r in round_log)
    if MOI.get(driver, MOI.RawOptimizerAttribute("success"))
        @test MOI.get(driver, MOI.TerminationStatus()) == MOI.LOCALLY_SOLVED
        @test MOI.get(driver, MOI.RawOptimizerAttribute("controller")) isa
              ST.FunnelController
    else
        @test MOI.get(driver, MOI.TerminationStatus()) == MOI.LOCALLY_INFEASIBLE
    end
    @test MOI.get(driver, MOI.RawOptimizerAttribute("trajectory")) isa ST.Trajectory
end

@testset "generator interface verbs" begin
    gen = MPPI.TrajectoryGenerator(;
        rng = Random.MersenneTwister(1),
        seed_trajectory = traj,
        nstep = 8,
        nsamples = 10,
        niter = 1,
        noise = MPPI.GaussianMPPINoise(SVector(0.1, 0.1)),
        project_input = u -> SVector(clamp(u[1], -1.0, 1.0), clamp(u[2], -1.0, 1.0)),
        cost = AB.ReachObjectiveCost(problem.target_set),
    )
    @test AB.set_stop_on_success!(gen, false) == true    # returns the previous value
    @test gen.stop_on_success == false
    AB.set_horizon!(gen, 5)
    @test gen.nstep == 5
    @test AB.get_seed(gen) === traj
end

@testset "TerminalPullCost (periodic endpoint pull)" begin
    pull = AB.TerminalPullCost([π, 0.0], [0.5, 1.0]; w = 2.0, periods = [2π, nothing])
    # Endpoint at π + 2π: the periodic difference is 0 — no pull.
    @test AB.cost_final(pull, 0.0, [3π, 0.0]) ≈ 0.0 atol = 1e-12
    # Endpoint half a radius off in θ: w · (0.25/0.5)² = 0.5.
    @test AB.cost_final(pull, 0.0, [π + 0.25, 0.0]) ≈ 0.5
end

@testset "FunnelController serialization round-trip" begin
    E1 = LazySets.Ellipsoid([0.0, 0.0], Matrix(0.1 * I, 2, 2))
    E2 = LazySets.Ellipsoid([0.5, 0.0], Matrix(0.05 * I, 2, 2))
    κ = MathematicalSystems.AffineMap([1.0 0.0], [0.1])
    ctrl = ST.FunnelController([κ], [E1, E2])

    io = IOBuffer()
    Serialization.serialize(io, ctrl)
    seekstart(io)
    ctrl2 = Serialization.deserialize(io)
    @test ctrl2 isa ST.FunnelController
    @test length(ctrl2.kappas) == 1
    @test Matrix(ctrl2.kappas[1].A) == Matrix(κ.A)
    @test LazySets.center(ctrl2.ellipsoids[2]) == LazySets.center(E2)
end

end # module
