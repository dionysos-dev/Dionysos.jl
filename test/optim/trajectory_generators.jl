module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import LazySets
import LinearAlgebra as LA
import MathOptInterface as MOI
import MathematicalSystems as MS
using Random

include("../../problems/Integrator/integrator.jl")

# Contract test for the three trajectory generators (optimizer-based, MPPI, composite).
# Tiny Integrator instance; each generator must run and expose the expected result API.
@testset "trajectory generators (optimizer / MPPI / composite)" begin
    _X_ = LazySets.Hyperrectangle(; low = SVector(-2.0, -2.0), high = SVector(4.0, 4.0))
    _U_ = LazySets.Hyperrectangle(; low = SVector(-1.0, -1.0), high = SVector(1.0, 1.0))
    concrete_system = Integrator.system(; _X_ = _X_, _U_ = _U_)
    jacobian_bound = Integrator.jacobian_bound()

    # Abstraction (alternating simulation) on a coarse grid.
    asp = PR.AlternatingSimulationProblem(concrete_system, concrete_system.X)
    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), asp)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("state_grid"),
        MP.GridFree(SVector(-2.0, -2.0), SVector(0.2, 0.2)),
    )
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("input_grid"),
        MP.GridFree(SVector(-1.0, -1.0), SVector(0.5, 0.5)),
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.3)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("approx_mode"),
        AB.UniformGridAbstraction.GROWTH,
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), jacobian_bound)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("n_samples"), 1)
    MOI.set(optimizer, MOI.Silent(), true)
    MOI.optimize!(optimizer)

    _I_ = LazySets.Hyperrectangle(; low = SVector(-1.7, -1.7), high = SVector(-1.6, -1.6))
    target_set = UT.set_union([
        LazySets.Hyperrectangle(; low = SVector(-1.0, 3.0), high = SVector(-0.3, 3.7)),
        LazySets.Hyperrectangle(; low = SVector(1.0, 2.0), high = SVector(3.0, 3.7)),
    ])
    concrete_problem =
        PR.OptimalControlProblem(concrete_system, _I_, target_set, nothing, nothing, 20)

    # 1) Optimizer-based generator: rolls out the abstract controller.
    opt_gen = AB.OptimizerTrajectoryGenerator.TrajectoryGenerator(
        optimizer;
        initial_state = SVector(-1.65, -1.65),
        concrete = false,
        nstep = 20,
    )
    AB.set_problem!(opt_gen, concrete_problem)
    AB.generate!(opt_gen)
    seed = AB.get_trajectory(opt_gen)
    @test seed !== nothing
    @test AB.get_success(opt_gen) == true      # reach must succeed to seed MPPI
    @test !isempty(ST.states(seed))
    @test !isempty(ST.inputs(seed))

    # 1b) Same generator, concrete rollout with a default initial state: with
    #     `initial_state = nothing` the start is picked as the center of the problem
    #     initial set (exercises `select_initial_state` on a `LazySet`), and
    #     `concrete = true` simulates the concrete closed loop.
    conc_gen = AB.OptimizerTrajectoryGenerator.TrajectoryGenerator(
        optimizer;
        initial_state = nothing,
        concrete = true,
        nstep = 20,
    )
    AB.set_problem!(conc_gen, concrete_problem)
    AB.generate!(conc_gen)
    conc_traj = AB.get_trajectory(conc_gen)
    @test conc_traj !== nothing
    @test !isempty(ST.states(conc_traj))
    @test AB.get_success(conc_gen) isa Bool

    # 2) MPPI generator, seeded by the optimizer trajectory (small sample budget).
    Δt = 0.3
    discrete_problem = PR.discretize_problem(concrete_problem, Δt)
    noise_sampler = (rng, u, k) -> SVector(0.3 * randn(rng), 0.3 * randn(rng))
    project_input = u -> SVector(clamp(u[1], -1.0, 1.0), clamp(u[2], -1.0, 1.0))
    trajectory_cost = (problem, traj) -> sum(LA.norm(u)^2 for u in ST.inputs(traj))
    mppi_gen = AB.MPPITrajectoryGenerator.TrajectoryGenerator(;
        rng = Random.MersenneTwister(0),
        seed_trajectory = seed,
        nstep = 20,
        nsamples = 40,
        niter = 3,
        λ = 1.0,
        noise_sampler = noise_sampler,
        project_input = project_input,
        trajectory_cost = trajectory_cost,
        hard_constraint = false,
    )
    AB.set_problem!(mppi_gen, discrete_problem)
    AB.generate!(mppi_gen)
    @test AB.get_trajectory(mppi_gen) !== nothing
    @test AB.get_success(mppi_gen) isa Bool

    # 3) Composite generator: optimizer seed refined by MPPI.
    combo_gen = AB.CompositeTrajectoryGenerator.TrajectoryGenerator(
        opt_gen,
        mppi_gen;
        Δt = 0.3,
        num_substeps = 5,
    )
    AB.set_problem!(combo_gen, concrete_problem)
    AB.generate!(combo_gen)
    @test AB.get_trajectory(combo_gen) !== nothing
end

# The P2 MPPI mechanics (plan.md §3): penalty rollouts, ESS-adaptive temperature,
# structured noise, determinism, stage costs vs closures, CEM mode.
@testset "MPPI mechanics" begin
    _X_ = LazySets.Hyperrectangle(; low = SVector(-2.0, -2.0), high = SVector(4.0, 4.0))
    _U_ = LazySets.Hyperrectangle(; low = SVector(-1.0, -1.0), high = SVector(1.0, 1.0))
    concrete_system = Integrator.system(; _X_ = _X_, _U_ = _U_)
    _I_ = LazySets.Hyperrectangle(; low = SVector(-1.7, -1.7), high = SVector(-1.6, -1.6))
    _T_ = LazySets.Hyperrectangle(; low = SVector(1.0, 2.0), high = SVector(3.0, 3.7))
    problem = PR.discretize_problem(
        PR.OptimalControlProblem(concrete_system, _I_, _T_, nothing, nothing, 25),
        0.3,
    )
    f = MS.mapping(problem.system)

    x0 = SVector(-1.65, -1.65)
    zero_seed = begin
        xs = [x0]
        for _ in 1:20
            push!(xs, f(xs[end], SVector(0.0, 0.0)))
        end
        ST.Trajectory(xs; inputs = fill(SVector(0.0, 0.0), 20))
    end

    make_gen(; kwargs...) = AB.MPPITrajectoryGenerator.TrajectoryGenerator(;
        rng = Random.MersenneTwister(7),
        seed_trajectory = zero_seed,
        nstep = 20,
        nsamples = 150,
        niter = 15,
        project_input = u -> SVector(clamp(u[1], -1.0, 1.0), clamp(u[2], -1.0, 1.0)),
        kwargs...,
    )

    # (a) Penalty rollouts: a control sequence that leaves X must cost *more* than one
    # that stays, under hard_constraint, even for a pure input-effort cost (the old
    # truncation semantics made it strictly cheaper).
    let
        us_leave = fill(SVector(-1.0, -1.0), 20)   # drives past the lower-left corner
        us_stay = fill(SVector(0.05, 0.05), 20)
        cost = AB.CompositeCost(AB.InputEffortCost(1.0))
        wrapf = x -> x
        X = MS.stateset(problem.system)
        c_leave, frozen_leave = AB.rollout_cost(f, x0, us_leave, wrapf, cost, X, true, 1e6)
        c_stay, frozen_stay = AB.rollout_cost(f, x0, us_stay, wrapf, cost, X, true, 1e6)
        @test frozen_leave > 0
        @test frozen_stay == 0
        @test c_leave > c_stay
    end

    # (b) Stage-cost terms match an equivalent trajectory closure on the same rollout.
    let
        us = [SVector(0.3 * sin(k), 0.2 * cos(k)) for k in 1:20]
        terms = AB.CompositeCost(
            AB.InputEffortCost(0.5),
            AB.DomainPenaltyCost(problem.system.X; w = 123.0),
        )
        X = MS.stateset(problem.system)
        c_terms, _ = AB.rollout_cost(f, x0, us, x -> x, terms, X, false, 0.0)
        traj, _ = AB.rollout_trajectory(f, x0, us, x -> x, X, false)
        closure =
            0.5 * sum(LA.norm(u)^2 for u in ST.inputs(traj)) +
            123.0 * count(x -> !(x ∈ problem.system.X), ST.states(traj))
        @test c_terms ≈ closure rtol = 1e-12
    end

    # (c) ESS-adaptive temperature hits its target where a fixed λ collapses.
    let
        costs = collect(range(0.0, 1000.0; length = 200))  # wide spread
        _, ess_fixed, _ = AB.MPPITrajectoryGenerator._softmin_weights(costs, 1.0)
        @test ess_fixed < 2.0                     # weight collapse at λ = 1
        λ = AB.MPPITrajectoryGenerator._ess_lambda(costs, 20.0)
        _, ess_adaptive, _ = AB.MPPITrajectoryGenerator._softmin_weights(costs, λ)
        @test isapprox(ess_adaptive, 20.0; rtol = 0.15)
    end

    # (d) Structured noise + ESS + antithetic end-to-end, deterministic given the rng.
    reach_cost = AB.CompositeCost(
        AB.ReachObjectiveCost(problem.target_set),
        AB.InputEffortCost(0.01),
        AB.DomainPenaltyCost(problem.system.X),
    )
    results = map(1:2) do _
        gen = make_gen(;
            noise = AB.MPPITrajectoryGenerator.GaussianMPPINoise(SVector(0.3, 0.3)),
            cost = reach_cost,
            ess_target = 0.1,
            antithetic = true,
        )
        AB.set_problem!(gen, problem)
        AB.generate!(gen)
        (
            AB.get_success(gen),
            AB.MPPITrajectoryGenerator.get_diagnostics(gen).best_cost,
            collect(ST.states(AB.get_trajectory(gen))[end]),
        )
    end
    @test results[1] == results[2]                # seeded determinism
    @test results[1][1] == true                   # reaches the target
    diag_gen = make_gen(;
        noise = AB.MPPITrajectoryGenerator.GaussianMPPINoise(SVector(0.3, 0.3)),
        cost = reach_cost,
        ess_target = 0.1,
    )
    AB.set_problem!(diag_gen, problem)
    AB.generate!(diag_gen)
    d = AB.MPPITrajectoryGenerator.get_diagnostics(diag_gen)
    @test d.ess > 1.0 && isfinite(d.weight_entropy) && d.λ_used > 0.0

    # (e) CEM update mode solves the same problem.
    let
        gen = make_gen(;
            noise = AB.MPPITrajectoryGenerator.GaussianMPPINoise(SVector(0.3, 0.3)),
            cost = reach_cost,
            update_rule = :cem,
            elite_frac = 0.2,
        )
        AB.set_problem!(gen, problem)
        AB.generate!(gen)
        @test AB.get_success(gen) == true
    end
end

@testset "truncate_at_target picks the DEEPEST hit, not the first" begin
    T = LazySets.Hyperrectangle(; low = SVector(-0.3, -0.3), high = SVector(0.3, 0.3))
    prob = PR.OptimalControlProblem(
        MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem(
            (x, u) -> x,
            2,
            2,
            nothing,
            nothing,
        ),
        nothing,
        T,
        nothing,
        nothing,
    )
    xs = [
        SVector(-1.0, -1.0),
        SVector(-0.5, -0.5),
        SVector(0.25, 0.25),   # first hit (shallow: 83% of the half-width out)
        SVector(0.0, 0.0),     # deepest hit
        SVector(0.5, 0.5),     # back outside
    ]
    tr = ST.Trajectory(xs; inputs = fill(SVector(0.0, 0.0), 4))
    truncated = AB.MPPITrajectoryGenerator.truncate_at_target(prob, tr)
    states_t = collect(ST.states(truncated))
    @test states_t[end] == SVector(0.0, 0.0)
    @test length(states_t) == 4
end

end # module TestMain
