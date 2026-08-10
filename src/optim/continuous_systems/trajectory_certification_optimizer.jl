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

    trajectory::Any
    controller::Any

    success::Bool
    solve_time_sec::T
    print_level::Int
end

function Optimizer(trajectory_generator, trajectory_certifier)
    return Optimizer{typeof(trajectory_generator), typeof(trajectory_certifier), Float64}(
        nothing,
        trajectory_generator,
        trajectory_certifier,
        nothing,
        nothing,
        false,
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
    opt.solve_time_sec = 0.0

    AB.set_problem!(opt.trajectory_generator, opt.concrete_problem)
    AB.set_problem!(opt.trajectory_certifier, opt.concrete_problem)

    AB.generate!(opt.trajectory_generator)

    gen_success = AB.get_success(opt.trajectory_generator)
    opt.trajectory = AB.get_trajectory(opt.trajectory_generator)

    if !gen_success
        opt.solve_time_sec = time() - t0
        return
    end

    AB.set_trajectory!(opt.trajectory_certifier, opt.trajectory)
    AB.certify!(opt.trajectory_certifier)

    cert_success = AB.get_success(opt.trajectory_certifier)

    opt.controller = AB.get_controller(opt.trajectory_certifier)
    opt.success = gen_success && cert_success
    opt.solve_time_sec = time() - t0

    return
end

end # module
