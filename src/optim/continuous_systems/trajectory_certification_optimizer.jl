module TrajectoryCertificationOptimizer

import MathOptInterface as MOI
import Dionysos

const DI = Dionysos
const ST = DI.System
const PR = DI.Problem
const OP = DI.Optim
const AB = OP.Abstraction

export Optimizer

mutable struct Optimizer{TG, TC, T} <: OP.AbstractDionysosOptimizer
    concrete_problem::Union{Nothing, PR.ProblemType}
    trajectory_generator::TG
    trajectory_certifier::TC

    # Generate ⇄ certify loop (plan.md §6): on a failed certification, `replan!`
    # (a `(generator, certifier) -> Nothing` hook, e.g. re-seeding or retargeting
    # the prefix at the certified suffix's entry ellipsoid) reconfigures the
    # generator and another round runs, up to `max_rounds` in total. Even with no
    # hook, rounds differ: the generator's rng advances between rounds.
    max_rounds::Int
    replan!::Any

    # Optional `traj -> traj` transform applied between generation and
    # certification — e.g. `ST.unwrap_trajectory` for periodic states, which the
    # certifiers (linearizing in ℝⁿ) require (plan.md §4.2-H).
    prepare_trajectory::Any

    # Optional problem override for the certifier: generation and certification may
    # consume different framings of the same task (discrete-time vs continuous,
    # globally normalized coordinates, seam-split target sets). `nothing` means the
    # certifier gets `concrete_problem` like the generator.
    certifier_problem::Any

    trajectory::Any
    controller::Any

    success::Bool
    rounds::Int
    # One entry per round: (; round, gen_success, gen_time, cert_run,
    # cert_success, cert_time) — the generate-vs-certify split per round.
    round_log::Vector{NamedTuple}
    solve_time_sec::T
    print_level::Int
end

function Optimizer(trajectory_generator, trajectory_certifier)
    return Optimizer{typeof(trajectory_generator), typeof(trajectory_certifier), Float64}(
        nothing,
        trajectory_generator,
        trajectory_certifier,
        1,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        false,
        0,
        NamedTuple[],
        0.0,
        1,
    )
end

MOI.is_empty(opt::Optimizer) = opt.concrete_problem === nothing

function MOI.set(opt::Optimizer, param::MOI.RawOptimizerAttribute, value)
    if Symbol(param.name) === :concrete_problem
        # Storage only — `optimize!` propagates the problem to the generator and
        # certifier (propagating here too would run e.g. `discretize_problem` twice).
        opt.concrete_problem = value
        return
    end
    return OP.set_field_attribute!(opt, param, value)
end

function MOI.optimize!(opt::Optimizer)
    @assert opt.concrete_problem !== nothing "Set concrete_problem first."

    t0 = time()

    opt.trajectory = nothing
    opt.controller = nothing
    opt.success = false
    opt.rounds = 0
    empty!(opt.round_log)
    opt.solve_time_sec = 0.0

    AB.set_problem!(opt.trajectory_generator, opt.concrete_problem)
    AB.set_problem!(
        opt.trajectory_certifier,
        opt.certifier_problem === nothing ? opt.concrete_problem : opt.certifier_problem,
    )

    for round in 1:max(opt.max_rounds, 1)
        opt.rounds = round

        AB.generate!(opt.trajectory_generator)

        gen_success = AB.get_success(opt.trajectory_generator)
        # Keep the best-so-far candidate: a later round whose generator produced
        # nothing must not erase an earlier trajectory.
        traj = AB.get_trajectory(opt.trajectory_generator)
        traj === nothing || (opt.trajectory = traj)

        cert_run = false
        cert_success = false
        if gen_success && traj !== nothing
            cert_traj =
                opt.prepare_trajectory === nothing ? traj : opt.prepare_trajectory(traj)
            AB.set_trajectory!(opt.trajectory_certifier, cert_traj)
            AB.certify!(opt.trajectory_certifier)
            cert_run = true
            cert_success = AB.get_success(opt.trajectory_certifier)
        end

        push!(
            opt.round_log,
            (;
                round,
                gen_success,
                gen_time = AB.get_solve_time(opt.trajectory_generator),
                cert_run,
                cert_success,
                cert_time = cert_run ? AB.get_solve_time(opt.trajectory_certifier) : NaN,
            ),
        )

        if cert_success
            opt.controller = AB.get_controller(opt.trajectory_certifier)
            opt.success = true
            break
        end

        round >= opt.max_rounds && break
        opt.replan! === nothing ||
            opt.replan!(opt.trajectory_generator, opt.trajectory_certifier)
    end

    opt.solve_time_sec = time() - t0

    return
end

function MOI.get(opt::Optimizer, ::MOI.TerminationStatus)
    opt.rounds == 0 && return MOI.OPTIMIZE_NOT_CALLED
    return opt.success ? MOI.LOCALLY_SOLVED : MOI.LOCALLY_INFEASIBLE
end

function MOI.get(opt::Optimizer, ::MOI.RawStatusString)
    opt.rounds == 0 && return "optimize! not called"
    opt.success && return "certified in $(opt.rounds) round(s)"
    return "no certified trajectory after $(opt.rounds) round(s)"
end

end # module
