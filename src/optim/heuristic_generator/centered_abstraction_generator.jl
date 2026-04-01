export CenteredAbstractionConfig, CenteredAbstractionGenerator

import Dionysos
import MathOptInterface as MOI

const DI = Dionysos
const ST = DI.System
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

struct CenteredAbstractionConfig{TH, TU, TJ, TP, FX0}
    Δt::Float64
    hx::TH
    Udom::TU
    jacobian_bound::TJ
    periodicity::TP
    nstep::Int
    num_substeps::Int
    x0_provider::FX0
    # Choix de la trajectoire candidate :
    # - :closed_loop   -> simulation avec le controleur concret
    # - :abstract_traj -> reconstruction depuis les etats abstraits
    trajectory_mode::Symbol
end

function CenteredAbstractionConfig(
    Δt::Real,
    hx,
    Udom,
    jacobian_bound,
    periodicity,
    nstep::Integer,
    x0_provider;
    num_substeps::Integer = 5,
    trajectory_mode::Symbol = :closed_loop,
)
    return CenteredAbstractionConfig(
        Float64(Δt),
        hx,
        Udom,
        jacobian_bound,
        periodicity,
        Int(nstep),
        Int(num_substeps),
        x0_provider,
        trajectory_mode,
    )
end

# Compatibilite avec l'ancienne arite positionnelle qui expose num_substeps.
function CenteredAbstractionConfig(
    Δt::Real,
    hx,
    Udom,
    jacobian_bound,
    periodicity,
    nstep::Integer,
    num_substeps::Integer,
    x0_provider;
    trajectory_mode::Symbol = :closed_loop,
)
    return CenteredAbstractionConfig(
        Δt,
        hx,
        Udom,
        jacobian_bound,
        periodicity,
        nstep,
        x0_provider;
        num_substeps = num_substeps,
        trajectory_mode = trajectory_mode,
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

function CenteredAbstractionGenerator(problem, config)
    P = problem === nothing ? DI.Problem.ProblemType : typeof(problem)
    return CenteredAbstractionGenerator{
        P,
        typeof(config),
        MOI.AbstractOptimizer,
        CandidateTrajectory,
    }(
        problem,
        config,
        nothing,
        nothing,
        false,
        0.0,
    )
end

function _periodicity(cfg::CenteredAbstractionConfig)
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

function _configure_optimizer!(optimizer, problem, cfg::CenteredAbstractionConfig, p)
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

function _build_wrap_function(p)
    if p === nothing || !p.with_period
        return identity
    end

    if hasproperty(p, :periodic_start)
        return ST.get_periodic_wrapper(
            p.periodic_dims,
            p.periodic_periods;
            start = p.periodic_start,
        )
    end

    return ST.get_periodic_wrapper(p.periodic_dims, p.periodic_periods)
end

function _build_closed_loop_candidate(problem, optimizer, cfg, p)
    # on simule le systeme discret avec le controleur concret.
    concrete_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

    concrete_system = problem.system
    disc_system = ST.discretize_continuous_system(
        concrete_system,
        cfg.Δt;
        num_substeps = cfg.num_substeps,
    )

    x0 = cfg.x0_provider(problem)
    wrap = _build_wrap_function(p)

    target_set = hasproperty(problem, :target_set) ? problem.target_set : nothing
    stopfun = target_set === nothing ? (_ -> false) : (x -> (wrap(x) ∈ target_set))

    traj = ST.get_closed_loop_trajectory(
        disc_system,
        concrete_controller,
        x0,
        cfg.nstep;
        stopping = stopfun,
        wrap = wrap,
    )

    length(traj.x) >= 2 || return nothing

    return CandidateTrajectory(
        traj.x,
        traj.u;
        Ts = cfg.Δt,
        source = :centered_closed_loop,
        metadata = (
            ;
            hx = cfg.hx,
            nstep = cfg.nstep,
            num_substeps = cfg.num_substeps,
            trajectory_mode = :closed_loop,
        ),
    )
end

function _select_best_abstract_step(abs_sys, k_abs, q::Int, value_fun_tab)
    # On choisit le couple (u, q_next) qui minimise la valeur abstraite du successeur. (je sais pas trop si ça a du sens, à terme je devrais généraliser potentiellement renvoyer plusieurs abstract traj)
    u_candidates = k_abs(q)
    u_candidates === nothing && return nothing, nothing

    inputs = u_candidates isa AbstractVector ? u_candidates : (u_candidates,)
    best_u = nothing
    best_q = nothing
    best_cost = Inf

    for u_sym in inputs
        for q_next in SY.post(abs_sys, q, u_sym)
            cost = value_fun_tab[q_next]
            if isfinite(cost) && cost < best_cost
                best_cost = cost
                best_u = u_sym
                best_q = q_next
            end
        end
    end

    return best_u, best_q
end

function _build_abstract_candidate(problem, optimizer, cfg, p)
    # Mode abstrait : on reconstruit une trajectoire concrete a partir des centres -> le centre du premier point de la traj c'est le centre de la cellule qui lui est assigné 
    # des cellules abstraites et des entrees abstraites choisies par le controleur.
    abs_sys = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    abs_ctrl = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
    value_fun_tab = MOI.get(optimizer, MOI.RawOptimizerAttribute("value_fun_tab"))

    k_abs = abs_ctrl.h
    wrap = _build_wrap_function(p)
    x0 = wrap(cfg.x0_provider(problem))

    q = SY.get_abstract_state(abs_sys, x0)
    q === nothing && return nothing
    SY.is_allowed_state(abs_sys, q) || return nothing
    isfinite(value_fun_tab[q]) || return nothing

    qs = Int[q]
    u_syms = Int[]
    xs = [wrap(SY.get_concrete_state(abs_sys, q))]

    for _ in 1:cfg.nstep
        # Une valeur abstraite nulle signifie que la cellule cible est atteinte.
        iszero(value_fun_tab[q]) && break

        u_sym, q_next = _select_best_abstract_step(abs_sys, k_abs, q, value_fun_tab)
        (u_sym === nothing || q_next === nothing) && break

        push!(u_syms, u_sym)
        push!(qs, q_next)
        push!(xs, wrap(SY.get_concrete_state(abs_sys, q_next)))

        q = q_next
    end

    length(xs) >= 2 || return nothing
    us = [SY.get_concrete_input(abs_sys, u_sym) for u_sym in u_syms]

    return CandidateTrajectory(
        ST.Trajectory(xs),
        ST.Trajectory(us);
        Ts = cfg.Δt,
        source = :centered_abstract,
        metadata = (
            ;
            hx = cfg.hx,
            nstep = cfg.nstep,
            num_substeps = cfg.num_substeps,
            trajectory_mode = :abstract_traj,
            q_traj = qs,
            u_sym_traj = u_syms,
        ),
    )
end

function _build_candidate(problem, optimizer, cfg, p)
    if cfg.trajectory_mode == :closed_loop
        return _build_closed_loop_candidate(problem, optimizer, cfg, p)
    elseif cfg.trajectory_mode == :abstract_traj
        return _build_abstract_candidate(problem, optimizer, cfg, p)
    end

    error("Unsupported trajectory_mode: $(cfg.trajectory_mode)")
end

function _candidate_reaches_target(problem, candidate, optimizer, cfg)
    candidate === nothing && return false

    # En mode abstrait, le dernier etat abstrait doit etre une cellule cible.
    if cfg.trajectory_mode == :abstract_traj &&
       candidate.metadata isa NamedTuple &&
       hasproperty(candidate.metadata, :q_traj)
        value_fun_tab = MOI.get(optimizer, MOI.RawOptimizerAttribute("value_fun_tab"))
        last_q = last(candidate.metadata.q_traj)
        return isfinite(value_fun_tab[last_q]) && iszero(value_fun_tab[last_q])
    end

    hasproperty(problem, :target_set) || return true
    last_state = last(ST.enum_elems(candidate.x_traj))
    return last_state ∈ problem.target_set
end

function set_problem!(gen::CenteredAbstractionGenerator, prob)
    gen.problem = prob
    gen.optimizer = nothing
    gen.candidate = nothing
    gen.success = false
    gen.solve_time_sec = 0.0
    return gen
end

function generate!(gen::CenteredAbstractionGenerator)
    cfg = gen.config
    @assert gen.problem !== nothing "Call set_problem!(gen, problem) first."
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

    candidate = _build_candidate(gen.problem, optimizer, cfg, p)
    success = _candidate_reaches_target(gen.problem, candidate, optimizer, cfg)

    gen.optimizer = optimizer
    gen.candidate = candidate
    gen.success = success
    gen.solve_time_sec = time() - t0
    return gen
end

get_trajectory(gen::CenteredAbstractionGenerator) = gen.candidate
get_success(gen::CenteredAbstractionGenerator) = gen.success
get_solve_time(gen::CenteredAbstractionGenerator) = gen.solve_time_sec
