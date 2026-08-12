# C2 — certifier ablation (plan.md §6b): one factor at a time around the mainline
# backward configuration (fixed global normalization, :vertices remainder,
# :maximin objective) on the certified pendulum swing-up. Reproduces Florentin's
# normalization sweep on mainline (none / fixed / trajectory-std) and ranks the
# remainder models and the SDP size objectives. Generation is cached per seed —
# every config certifies the same trajectory, so comparisons are paired.
#
#     julia --project=bench research/TrajectoryCertificationOptimizer/campaigns/c2_certifier_ablation.jl

include(joinpath(@__DIR__, "pendulum_pipeline.jl"))
include(joinpath(@__DIR__, "runner.jl"))
import .CampaignRunner

const FAILED_ROW = (;
    success = false,
    gen_success = false,
    V0 = NaN,
    Vmin = NaN,
    Vmed = NaN,
    coverage = NaN,
    chain_frac = 0.0,
    cert_time_s = NaN,
)

function run_one(config, rng)
    gen = cached_pendulum(rng)
    gen === nothing && return FAILED_ROW

    t =
        config.norm === :std ? trajectory_std(gen.traj) :
        config.norm === :none ? [1.0, 1.0] : T_FIXED
    opts = pend_backward_options(
        t;
        objective = config.objective,
        remainder_model = config.remainder,
    )
    cert = EB.BackwardCertifier(zprovider(t), sdp, opts)
    AB.set_problem!(cert, zproblem(t))
    AB.set_trajectory!(cert, ztraj(gen.traj, t))
    AB.certify!(cert)
    res = EB.get_result(cert)

    K = length(ST.inputs(gen.traj))
    areas = res.success ? funnel_areas(res, t) : Float64[]
    sorted = sort(areas)
    return (;
        success = res.success,
        gen_success = gen.gen_success,
        V0 = isempty(areas) ? NaN : areas[1],
        Vmin = isempty(sorted) ? NaN : sorted[1],
        Vmed = isempty(sorted) ? NaN : sorted[max(1, div(length(sorted), 2))],
        coverage = res.initial_coverage === nothing ? NaN : Float64(res.initial_coverage),
        chain_frac = res.failed_k === nothing ? 1.0 : (K - res.failed_k) / K,
        cert_time_s = AB.get_solve_time(cert),
    )
end

configs = [
    (;
        label = "fixed_vert_maximin",
        norm = :fixed,
        remainder = :vertices,
        objective = :maximin,
    ),
    (;
        label = "none_vert_maximin",
        norm = :none,
        remainder = :vertices,
        objective = :maximin,
    ),
    (;
        label = "std_vert_maximin",
        norm = :std,
        remainder = :vertices,
        objective = :maximin,
    ),
    (;
        label = "fixed_ball_maximin",
        norm = :fixed,
        remainder = :ball,
        objective = :maximin,
    ),
    (;
        label = "fixed_john_maximin",
        norm = :fixed,
        remainder = :john_ball,
        objective = :maximin,
    ),
    (;
        label = "fixed_vert_logdet",
        norm = :fixed,
        remainder = :vertices,
        objective = :logdet,
    ),
    (;
        label = "fixed_vert_trace",
        norm = :fixed,
        remainder = :vertices,
        objective = :trace,
    ),
]

CampaignRunner.run_campaign(;
    name = "c2_certifier_ablation",
    configs = configs,
    run_one = run_one,
    nseeds = parse(Int, get(ENV, "CAMPAIGN_NSEEDS", "8")),
)
