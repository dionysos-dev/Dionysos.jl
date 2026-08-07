module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import LazySets
import LinearAlgebra as LA
import MathOptInterface as MOI
using Random

include("../../problems/Integrator/integrator.jl")

# Contract test for the three trajectory generators (optimizer-based, MPPI, composite).
# Tiny Integrator instance; each generator must run and expose the expected result API.
@testset "trajectory generators (optimizer / MPPI / composite)" begin
    _X_ = UT.box(SVector(-2.0, -2.0), SVector(4.0, 4.0))
    _U_ = UT.box(SVector(-1.0, -1.0), SVector(1.0, 1.0))
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

    _I_ = UT.box(SVector(-1.7, -1.7), SVector(-1.6, -1.6))
    target_set = UT.set_union([
        UT.box(SVector(-1.0, 3.0), SVector(-0.3, 3.7)),
        UT.box(SVector(1.0, 2.0), SVector(3.0, 3.7)),
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

end # module TestMain
