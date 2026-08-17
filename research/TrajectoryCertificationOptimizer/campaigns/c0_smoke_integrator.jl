# C0 — smoke campaign (plan.md §6b, P0): proves the multi-seed campaign runner on the
# cheapest possible case. MPPI-only reach on the 2-D integrator, seeded with a
# zero-input rollout; two sample budgets and two temperatures. Runs in seconds.
#
#     julia --project=bench research/TrajectoryCertificationOptimizer/campaigns/c0_smoke_integrator.jl

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const AB = DI.Optim.Abstraction

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
concrete_problem = PR.OptimalControlProblem(concrete_system, _I_, _T_, nothing, nothing, 25)
discrete_problem = PR.discretize_problem(concrete_problem, Δt)

x0 = SVector(-1.65, -1.65)

# Zero-input seed rollout: MPPI only needs an input sequence to perturb.
function zero_seed(nstep)
    f = MS.mapping(discrete_problem.system)
    u0 = SVector(0.0, 0.0)
    xs = [x0]
    for _ in 1:nstep
        push!(xs, f(xs[end], u0))
    end
    return ST.Trajectory(xs; inputs = fill(u0, nstep))
end

trajectory_cost = function (problem, tr)
    xs = ST.states(tr)
    us = ST.inputs(tr)
    best = minimum(LA.norm(x - LazySets.center(problem.target_set)) for x in xs)
    hit_idx = findfirst(x -> x ∈ problem.target_set, xs)
    bonus = hit_idx === nothing ? 0.0 : -1000.0 / hit_idx
    violation = sum(x -> x ∈ problem.system.X ? 0.0 : 1000.0, xs)
    return 100.0 * best + 0.01 * sum(LA.norm(u)^2 for u in us) + violation + bonus
end

function run_one(config, rng)
    nstep = 25
    gen = AB.MPPITrajectoryGenerator.TrajectoryGenerator(;
        rng = rng,
        seed_trajectory = zero_seed(nstep),
        nstep = nstep,
        nsamples = config.nsamples,
        niter = 10,
        λ = config.λ,
        noise_sampler = (rng, u, k) -> SVector(0.3 * randn(rng), 0.3 * randn(rng)),
        project_input = u -> SVector(clamp(u[1], -1.0, 1.0), clamp(u[2], -1.0, 1.0)),
        trajectory_cost = trajectory_cost,
        hard_constraint = false,
    )
    AB.set_problem!(gen, discrete_problem)
    stats = @timed AB.generate!(gen)
    traj = AB.get_trajectory(gen)
    return (;
        success = AB.get_success(gen),
        cost = traj === nothing ? NaN : trajectory_cost(discrete_problem, traj),
        time_s = stats.time,
    )
end

configs = [
    (; label = "ns50_lam1", nsamples = 50, λ = 1.0),
    (; label = "ns200_lam1", nsamples = 200, λ = 1.0),
    (; label = "ns50_lam50", nsamples = 50, λ = 50.0),
    (; label = "ns200_lam50", nsamples = 200, λ = 50.0),
]

CampaignRunner.run_campaign(;
    name = "c0_smoke_integrator",
    configs = configs,
    run_one = run_one,
    nseeds = 5,
)
