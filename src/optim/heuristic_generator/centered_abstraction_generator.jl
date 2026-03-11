export CenteredAbstractionConfig, CenteredAbstractionGenerator

import Dionysos
import MathOptInterface as MOI

const DI = Dionysos
const ST = DI.System
const OP = DI.Optim
const AB = OP.Abstraction
# j'aime pas trop cette logique (à discuter avec julien (c'est long pour pas grand chose))
struct CenteredAbstractionConfig{TH, TU, TJ, TP, FX0} 
    Δt::Float64
    hx::TH
    Udom::TU
    jacobian_bound::TJ
    periodicity::TP
    nstep::Int
    num_substeps::Int
    x0_provider::FX0
end

function CenteredAbstractionConfig(
    Δt::Real,
    hx,
    Udom,
    jacobian_bound,
    periodicity,
    nstep::Integer,
    x0_provider,
)
    return CenteredAbstractionConfig{ # je modifierai les typages une bonne fois pour toutes
        typeof(hx),
        typeof(Udom),
        typeof(jacobian_bound),
        typeof(periodicity),
        typeof(x0_provider),
    }(
        Float64(Δt),
        hx,
        Udom,
        jacobian_bound,
        periodicity,
        Int(nstep),
        5,
        x0_provider,
    )
end

function CenteredAbstractionConfig(
    Δt::Real,
    hx,
    Udom,
    jacobian_bound,
    periodicity,
    nstep::Integer,
    num_substeps::Integer,
    x0_provider,
)
    return CenteredAbstractionConfig{
        typeof(hx),
        typeof(Udom),
        typeof(jacobian_bound),
        typeof(periodicity),
        typeof(x0_provider),
    }(
        Float64(Δt),
        hx,
        Udom,
        jacobian_bound,
        periodicity,
        Int(nstep),
        Int(num_substeps),
        x0_provider,
    )
end

mutable struct CenteredAbstractionGenerator{P, C, O, CT} <: AbstractHeuristicGenerator
    problem::Union{Nothing, P}
    config::C
    optimizer::Union{Nothing, O}
    candidate::Union{Nothing, CT}
    success::Bool
    solve_time_sec::Float64
end

function _periodicity(cfg::CenteredAbstractionConfig) # je sais pas trop comment gérer ça 
    p = cfg.periodicity
    p === nothing && return nothing
    p isa NamedTuple || error("cfg.periodicity must be nothing or a NamedTuple")

    ks = Tuple(keys(p))
    required = (:with_period, :periodic_dims, :periodic_periods)
    allowed = (:with_period, :periodic_dims, :periodic_periods, :periodic_start)

    for k in required
        k in ks || error("cfg.periodicity is missing required key: $(k)")
    end
    for k in ks
        k in allowed || error("cfg.periodicity has unsupported key: $(k)")
    end
    p.with_period isa Bool || error("cfg.periodicity.with_period must be Bool")

    return p
end

function _configure_optimizer!(optimizer, problem, cfg::CenteredAbstractionConfig, p) # litéralement ce que l'on fait dans dio, y'a pas une constrution plus smart ? 
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("h"), cfg.hx)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("UMapping"), cfg.Udom)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), cfg.jacobian_bound)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), cfg.Δt)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("approx_mode"),
        AB.UniformGridAbstraction.CENTER_SIMULATION,
    )

    if p !== nothing
        MOI.set(optimizer, MOI.RawOptimizerAttribute("use_periodic_mapping"), p.with_period)
        if p.with_period
            MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_dims"), p.periodic_dims)
            MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_periods"), p.periodic_periods)
            if hasproperty(p, :periodic_start)
                MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_start"), p.periodic_start)
            end
        end
    end

    return optimizer
end

function set_problem!(gen::CenteredAbstractionGenerator, prob) # j'aime bien
    gen.problem = prob
    gen.optimizer = nothing
    gen.candidate = nothing
    gen.success = false
    gen.solve_time_sec = 0.0
    return gen
end

function generate!(gen::CenteredAbstractionGenerator) #j'aime bien
    cfg = gen.config
    @assert gen.problem !== nothing
    @assert cfg.Δt > 0.0
    @assert cfg.nstep >= 1
    @assert cfg.num_substeps >= 1
    p = _periodicity(cfg)
    t0 = time()

    optimizer =
        gen.optimizer === nothing ? MOI.instantiate(AB.UniformGridAbstraction.Optimizer) :
        gen.optimizer

    _configure_optimizer!(optimizer, gen.problem, cfg, p)
    MOI.optimize!(optimizer)

    concrete_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

    concrete_system = gen.problem.system
    disc_system =
        ST.discretize_continuous_system(
            concrete_system,
            cfg.Δt;
            num_substeps = cfg.num_substeps,
        )

    x0 = cfg.x0_provider(gen.problem)

    wrap = identity
    if p !== nothing && p.with_period
        if hasproperty(p, :periodic_start)
            wrap = ST.get_periodic_wrapper(
                p.periodic_dims,
                p.periodic_periods;
                start = p.periodic_start,
            )
        else
            wrap = ST.get_periodic_wrapper(p.periodic_dims, p.periodic_periods)
        end
    end

    traj = ST.get_closed_loop_trajectory(
        disc_system,
        concrete_controller,
        x0,
        cfg.nstep;
        stopping = x -> false,
        wrap = wrap,
    )

    x_traj = traj.x
    u_traj = traj.u

    candidate = CandidateTrajectory(
        x_traj,
        u_traj;
        Ts = cfg.Δt,
        source = :centered,
        metadata = (; hx = cfg.hx, nstep = cfg.nstep),
    )

    success = length(x_traj) > 0
    if success && hasproperty(gen.problem, :target_set)
        last_state = last(ST.enum_elems(x_traj))
        success = last_state ∈ gen.problem.target_set
    end

    gen.optimizer = optimizer
    gen.candidate = candidate
    gen.success = success
    gen.solve_time_sec = time() - t0
    return gen
end

get_trajectory(gen::CenteredAbstractionGenerator) = gen.candidate
get_success(gen::CenteredAbstractionGenerator) = gen.success
get_solve_time(gen::CenteredAbstractionGenerator) = gen.solve_time_sec
