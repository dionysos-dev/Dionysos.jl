# C1 — generator ablation (plan.md §6b): does the P2 machinery earn its place?
# Integrator reach task from a zero-input seed. Compares the legacy closure-cost
# MPPI against the typed stage-cost path, ESS-adaptive temperature, antithetic
# sampling, and the CEM update, over seeded runs.
#
#     julia --project=bench research/TrajectoryCertificationOptimizer/campaigns/c1_generator_ablation.jl

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const AB = DI.Optim.Abstraction
const MPPI = AB.MPPITrajectoryGenerator

import LazySets
import LinearAlgebra as LA
import MathematicalSystems as MS
using StaticArrays
using Random

include(
    joinpath(dirname(dirname(pathof(Dionysos))), "problems", "Integrator", "integrator.jl"),
)
include(joinpath(@__DIR__, "runner.jl"))
import .CampaignRunner

Δt = 0.3
_X_ = LazySets.Hyperrectangle(; low = SVector(-2.0, -2.0), high = SVector(4.0, 4.0))
_U_ = LazySets.Hyperrectangle(; low = SVector(-1.0, -1.0), high = SVector(1.0, 1.0))
concrete_system = Integrator.system(; _X_ = _X_, _U_ = _U_)
_I_ = LazySets.Hyperrectangle(; low = SVector(-1.7, -1.7), high = SVector(-1.6, -1.6))
_T_ = LazySets.Hyperrectangle(; low = SVector(1.0, 2.0), high = SVector(3.0, 3.7))
problem = PR.discretize_problem(
    PR.OptimalControlProblem(concrete_system, _I_, _T_, nothing, nothing, 25),
    Δt,
)
f = MS.mapping(problem.system)
x0 = SVector(-1.65, -1.65)
nstep = 20

zero_seed = begin
    xs = [x0]
    for _ in 1:nstep
        push!(xs, f(xs[end], SVector(0.0, 0.0)))
    end
    ST.Trajectory(xs; inputs = fill(SVector(0.0, 0.0), nstep))
end

stage_cost = AB.CompositeCost(
    AB.ReachObjectiveCost(problem.target_set),
    AB.InputEffortCost(0.01),
    AB.DomainPenaltyCost(problem.system.X),
)

closure_cost = function (prob, tr)
    xs = ST.states(tr)
    best = minimum(LA.norm(x - LazySets.center(prob.target_set)) for x in xs)
    hit = findfirst(x -> x ∈ prob.target_set, xs)
    bonus = hit === nothing ? 0.0 : -1000.0 / hit
    violation = sum(x -> x ∈ prob.system.X ? 0.0 : 1000.0, xs)
    return 100.0 * best +
           0.01 * sum(LA.norm(u)^2 for u in ST.inputs(tr)) +
           violation +
           bonus
end

project = u -> SVector(clamp(u[1], -1.0, 1.0), clamp(u[2], -1.0, 1.0))

function run_one(config, rng)
    gen = MPPI.TrajectoryGenerator(;
        rng = rng,
        seed_trajectory = zero_seed,
        nstep = nstep,
        nsamples = 150,
        niter = 15,
        project_input = project,
        config.gen_kwargs...,
    )
    AB.set_problem!(gen, problem)
    stats = @timed AB.generate!(gen)
    d = MPPI.get_diagnostics(gen)
    return (;
        success = AB.get_success(gen),
        cost = d.best_cost,
        ess = d.ess,
        time_s = stats.time,
    )
end

gauss = MPPI.GaussianMPPINoise(SVector(0.3, 0.3))
configs = [
    (;
        label = "closure_lam1",
        gen_kwargs = (;
            noise_sampler = (rng, u, k) -> SVector(0.3 * randn(rng), 0.3 * randn(rng)),
            trajectory_cost = closure_cost,
            λ = 1.0,
        ),
    ),
    (; label = "stage_lam1", gen_kwargs = (; noise = gauss, cost = stage_cost, λ = 1.0)),
    (;
        label = "stage_ess",
        gen_kwargs = (; noise = gauss, cost = stage_cost, ess_target = 0.1),
    ),
    (;
        label = "stage_ess_anti",
        gen_kwargs = (;
            noise = gauss,
            cost = stage_cost,
            ess_target = 0.1,
            antithetic = true,
        ),
    ),
    (;
        label = "cem",
        gen_kwargs = (; noise = gauss, cost = stage_cost, update_rule = :cem),
    ),
]

CampaignRunner.run_campaign(;
    name = "c1_generator_ablation",
    configs = configs,
    run_one = run_one,
    nseeds = 10,
)
