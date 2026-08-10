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
    # generator and another round runs, up to `max_rounds` in total.
    max_rounds::Int
    replan!::Any

    trajectory::Any
    controller::Any

    success::Bool
    rounds::Int
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
        false,
        0,
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
    opt.solve_time_sec = 0.0

    AB.set_problem!(opt.trajectory_generator, opt.concrete_problem)
    AB.set_problem!(opt.trajectory_certifier, opt.concrete_problem)

    for round in 1:max(opt.max_rounds, 1)
        opt.rounds = round

        AB.generate!(opt.trajectory_generator)

        gen_success = AB.get_success(opt.trajectory_generator)
        opt.trajectory = AB.get_trajectory(opt.trajectory_generator)

        if gen_success
            AB.set_trajectory!(opt.trajectory_certifier, opt.trajectory)
            AB.certify!(opt.trajectory_certifier)

            if AB.get_success(opt.trajectory_certifier)
                opt.controller = AB.get_controller(opt.trajectory_certifier)
                opt.success = true
                break
            end
        end

        round == opt.max_rounds && break
        opt.replan! === nothing ||
            opt.replan!(opt.trajectory_generator, opt.trajectory_certifier)
    end

    opt.solve_time_sec = time() - t0

    return
end

end # module
