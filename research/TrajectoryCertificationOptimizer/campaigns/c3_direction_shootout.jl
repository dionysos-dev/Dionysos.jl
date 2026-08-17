# C3 — direction shoot-out (plan.md §6b): backward vs forward :fixed vs forward
# :free vs the bidirectional handoff, per system — the fully actuated linear
# integrator (B = Δt·I, both directions should close) against the underactuated
# pendulum swing-up (the regime where backward is expected to dominate). All
# directions share the mainline step configuration (:vertices, :maximin) and
# certify the SAME cached per-seed trajectory, so comparisons are paired.
#
# Reported per run: certified entry area (physical frame), fraction of the chain
# certified, handoff position (bidirectional), certification time.
#
#     julia --project=bench research/TrajectoryCertificationOptimizer/campaigns/c3_direction_shootout.jl

include(joinpath(@__DIR__, "pendulum_pipeline.jl"))
include(joinpath(@__DIR__, "runner.jl"))
import .CampaignRunner

include(joinpath(PROBLEMS, "Integrator", "integrator.jl"))

# ------------------------------------------------------------
# Integrator leg: C1's reach task and its CEM winner, plus an exact symbolic
# provider (RK4 of ẋ = u is x + Δt·u — zero Hessian, so remainders vanish and
# fixed generous linearization boxes are free).
# ------------------------------------------------------------

const int_Δt = 0.3
const int_X = LazySets.Hyperrectangle(; low = SVector(-2.0, -2.0), high = SVector(4.0, 4.0))
const int_U = LazySets.Hyperrectangle(; low = SVector(-1.0, -1.0), high = SVector(1.0, 1.0))
const int_system = Integrator.system(; _X_ = int_X, _U_ = int_U)
const int_I =
    LazySets.Hyperrectangle(; low = SVector(-1.7, -1.7), high = SVector(-1.6, -1.6))
const int_T = LazySets.Hyperrectangle(; low = SVector(1.0, 2.0), high = SVector(3.0, 3.7))
const int_problem = PR.discretize_problem(
    PR.OptimalControlProblem(int_system, int_I, int_T, nothing, nothing, 25),
    int_Δt,
)
const int_x0 = SVector(-1.65, -1.65)
const int_nstep = 20

const int_seed = begin
    f = MS.mapping(int_problem.system)
    xs = [int_x0]
    for _ in 1:int_nstep
        push!(xs, f(xs[end], SVector(0.0, 0.0)))
    end
    ST.Trajectory(xs; inputs = fill(SVector(0.0, 0.0), int_nstep))
end

const int_cost = AB.CompositeCost(
    AB.ReachObjectiveCost(int_problem.target_set),
    AB.InputEffortCost(0.01),
    AB.DomainPenaltyCost(int_problem.system.X),
)

function generate_integrator(rng)
    gen = MPPI.TrajectoryGenerator(;
        rng = rng,
        seed_trajectory = int_seed,
        nstep = int_nstep,
        nsamples = 150,
        niter = 15,
        noise = MPPI.GaussianMPPINoise(SVector(0.3, 0.3)),
        project_input = u -> SVector(clamp(u[1], -1.0, 1.0), clamp(u[2], -1.0, 1.0)),
        cost = int_cost,
        update_rule = :cem,
    )
    AB.set_problem!(gen, int_problem)
    AB.generate!(gen)
    traj = AB.get_trajectory(gen)
    traj === nothing && return nothing
    return (; traj = traj, gen_success = AB.get_success(gen), gen_time = 0.0)
end

const _INT_CACHE = Dict{UInt64, Any}()
cached_integrator(rng) = get!(() -> generate_integrator(rng), _INT_CACHE, rand(rng, UInt64))

const int_provider = begin
    Symbolics.@variables x1 x2 u1 u2 v1 v2 T
    f_disc = ST.runge_kutta4((x, u) -> [u[1], u[2]], [x1, x2], [u1, u2], T, 1)
    fsymbolic =
        Symbolics.substitute([f_disc[1] + v1, f_disc[2] + v2], Dict(T => int_Δt))
    Wset = LazySets.Hyperrectangle(; low = SVector(0.0, 0.0), high = SVector(0.0, 0.0))
    ST.SymbolicAffineApproximationProvider(
        fsymbolic,
        [x1, x2],
        [u1, u2],
        [v1, v2],
        [0.0, 0.0],
        ST.format_input_set(int_U),
        ST.format_noise_set(Wset),
    )
end

function int_backward_options()
    return EB.ChainOptions(;
        maxδx = 12.0,
        maxδu = 2.0,
        λ = 0.001,
        terminal_shrink = 0.95,
        linearization_δx = [4.0, 4.0],
        linearization_δu = [2.0, 2.0],
        adaptive_boxes = nothing,
        objective = :maximin,
        remainder_model = :vertices,
        domain_cap = true,
        check_state_domain = true,
    )
end

# The two forward modes need OPPOSITE λ regimes (measured, k=1 ladder):
# :fixed's α is dimensionless, and a size-dominated objective (λ ≪ 1) polishes
# the min-α solution onto the strict-PSD razor edge where the solver tolerance
# eats the ε margin and validation (correctly) rejects a "solved" step — λ = 0.5
# regularizes it (full 20/20 chain). :free's trace term is in absolute tube
# units, so λ = 0.5 lets the cost term buy input effort by drifting the tube off
# the nominal (×2.5 growth per step, input death at k=4) — the shipped λ = 0.01
# is right there, PROVIDED q_min floors the tube near the entry scale (4e-3 vs
# entry 5e-3 chains 20/20; 1e-3 dies of source conditioning at k=5; the 1e-9
# default needle-collapses at k=1).
function int_forward_options(; target_mode = :fixed)
    return EB.ForwardOptions(;
        target_mode = target_mode,
        maxδu = 2.0,
        λ = target_mode === :fixed ? 0.5 : 0.01,
        q_min = 4e-3,
        linearization_δu = [2.0, 2.0],
        linearization_δx_margin = 1.1,
        check_state_domain = true,
        remainder_model = :vertices,
    )
end

# ------------------------------------------------------------
# The shoot-out
# ------------------------------------------------------------

_entry_area(res, t) = isempty(res.lmi_data.ellipsoids) ? NaN : funnel_areas(res, t)[1]

# Certified transitions out of K: backward dies AT failed_k (certified above it),
# forward dies at failed_k (certified below it); failed_k == K + 1 is a complete
# chain whose terminal gate failed.
_backward_frac(res, K) = res.failed_k === nothing ? 1.0 : (K - res.failed_k) / K
_forward_frac(res, K) = res.failed_k === nothing ? 1.0 : (res.failed_k - 1) / K

function _certify_direction(direction, make_backward, make_forward, prob, traj, t)
    K = length(ST.inputs(traj))
    if direction === :backward
        cert = make_backward()
        AB.set_problem!(cert, prob)
        AB.set_trajectory!(cert, traj)
        AB.certify!(cert)
        res = EB.get_result(cert)
        return (;
            success = res.success,
            entry_area = res.success ? _entry_area(res, t) : NaN,
            certified_frac = _backward_frac(res, K),
            handoff_frac = NaN,
            cert_time_s = res.solve_time_sec,
        )
    elseif direction === :bidirectional
        fwd, bwd = make_forward(:fixed), make_backward()
        out = EB.bidirectional_certify!(fwd, bwd, prob, traj)
        fres, bres = out.forward_result, out.backward_result
        return (;
            success = out.success,
            entry_area = out.success ? _entry_area(fres, t) : NaN,
            certified_frac = out.success ? 1.0 :
                             max(_forward_frac(fres, K), _backward_frac(bres, K)),
            handoff_frac = out.k_handoff === nothing ? NaN : out.k_handoff / (K + 1),
            cert_time_s = fres.solve_time_sec + bres.solve_time_sec,
        )
    else
        mode = direction === :forward_fixed ? :fixed : :free
        cert = make_forward(mode)
        AB.set_problem!(cert, prob)
        AB.set_trajectory!(cert, traj)
        AB.certify!(cert)
        res = EB.get_result(cert)
        return (;
            success = res.success,
            entry_area = _entry_area(res, t),
            certified_frac = _forward_frac(res, K),
            handoff_frac = NaN,
            cert_time_s = res.solve_time_sec,
        )
    end
end

const FAILED_ROW = (;
    success = false,
    gen_success = false,
    entry_area = NaN,
    certified_frac = 0.0,
    handoff_frac = NaN,
    cert_time_s = NaN,
)

function run_one(config, rng)
    if config.system === :integrator
        gen = cached_integrator(rng)
        gen === nothing && return FAILED_ROW
        t = [1.0, 1.0]
        out = _certify_direction(
            config.direction,
            () -> EB.BackwardCertifier(int_provider, sdp, int_backward_options()),
            mode -> EB.ForwardCertifier(
                int_provider,
                sdp,
                int_forward_options(; target_mode = mode),
            ),
            int_problem,
            gen.traj,
            t,
        )
        return (; out..., gen_success = gen.gen_success)
    else
        gen = cached_pendulum(rng)
        gen === nothing && return FAILED_ROW
        t = T_FIXED
        escale = get(config, :entry_scale, nothing)
        eshape = escale === nothing ? nothing : pend_entry_shape(t, escale)
        out = _certify_direction(
            config.direction,
            () -> EB.BackwardCertifier(zprovider(t), sdp, pend_backward_options(t)),
            mode -> EB.ForwardCertifier(
                zprovider(t),
                sdp,
                pend_forward_options(; target_mode = mode, entry_shape = eshape),
            ),
            zproblem(t),
            ztraj(gen.traj, t),
            t,
        )
        return (; out..., gen_success = gen.gen_success)
    end
end

configs = [
    (; label = "int_backward", system = :integrator, direction = :backward),
    (; label = "int_fwd_fixed", system = :integrator, direction = :forward_fixed),
    (; label = "int_fwd_free", system = :integrator, direction = :forward_free),
    (; label = "int_bidir", system = :integrator, direction = :bidirectional),
    (; label = "pend_backward", system = :pendulum, direction = :backward),
    (; label = "pend_fwd_fixed", system = :pendulum, direction = :forward_fixed),
    (; label = "pend_fwd_free", system = :pendulum, direction = :forward_free),
    (; label = "pend_bidir", system = :pendulum, direction = :bidirectional),
    # Entry ladder: can a forward tube defend a QUARTER of the initial set?
    # Separates "the full-set tube is indefensible" (physics) from tool failure,
    # and gives the handoff a nestable forward tube.
    (;
        label = "pend_fwd_fix_s25",
        system = :pendulum,
        direction = :forward_fixed,
        entry_scale = 0.25,
    ),
    (;
        label = "pend_bidir_s25",
        system = :pendulum,
        direction = :bidirectional,
        entry_scale = 0.25,
    ),
]

CampaignRunner.run_campaign(;
    name = "c3_direction_shootout",
    configs = configs,
    run_one = run_one,
    nseeds = parse(Int, get(ENV, "CAMPAIGN_NSEEDS", "8")),
)
