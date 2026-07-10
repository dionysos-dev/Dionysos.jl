module OptimizerTrajectoryGenerator

import MathOptInterface as MOI

import ..AbstractTrajectoryGenerator
import ..set_problem!
import ..generate!
import ..get_trajectory
import ..get_success
import ..get_solve_time

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const SY = DI.Symbolic
const PR = DI.Problem

mutable struct TrajectoryGenerator <: AbstractTrajectoryGenerator
    # Inputs
    problem::Union{Nothing, PR.ProblemType}
    optimizer::MOI.AbstractOptimizer

    # Parameters
    initial_state::Any
    # Initial state used for closed-loop simulation.
    # If `nothing`, one representative state is selected from the problem initial set.

    concrete::Bool
    # If true, simulate the concrete closed-loop system.
    # If false, simulate the abstract closed-loop system and concretize both states and inputs.

    nstep::Int
    # Maximum number of closed-loop simulation steps.

    # Outputs
    trajectory::Union{Nothing, ST.ClosedLoopTrajectory}
    success::Bool
    solve_time_sec::Float64
end

function TrajectoryGenerator(
    optimizer;
    initial_state = nothing,
    concrete::Bool = true,
    nstep::Int,
)
    return TrajectoryGenerator(
        nothing,
        optimizer,
        initial_state,
        concrete,
        nstep,
        nothing,
        false,
        NaN,
    )
end

function set_problem!(gen::TrajectoryGenerator, problem::PR.ProblemType)
    gen.problem = problem
    return gen
end

function generate!(gen::TrajectoryGenerator)
    problem = gen.problem
    problem === nothing && error("No problem attached. Call set_problem! first.")

    t0 = time()

    MOI.set(gen.optimizer, MOI.RawOptimizerAttribute("concrete_problem"), problem)
    MOI.optimize!(gen.optimizer)

    gen.solve_time_sec = time() - t0

    gen.trajectory =
        gen.concrete ? _generate_concrete_trajectory(gen) :
        _generate_abstract_trajectory(gen)

    gen.success = PR.trajectory_success(problem, gen.trajectory.x)

    return gen
end

get_trajectory(gen::TrajectoryGenerator) = gen.trajectory
get_success(gen::TrajectoryGenerator) = gen.success
get_solve_time(gen::TrajectoryGenerator) = gen.solve_time_sec

select_initial_state(initial_set::UT.AbstractSetNode) = UT.get_center(initial_set)
select_initial_state(initial_set::Vector{Int}) = first(initial_set)

function _initial_state(gen::TrajectoryGenerator, problem)
    gen.initial_state !== nothing && return gen.initial_state
    return select_initial_state(problem.initial_set)
end

function _generate_concrete_trajectory(gen::TrajectoryGenerator)
    optimizer = gen.optimizer
    discrete_time_system =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system"))
    concrete_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

    x0 = _initial_state(gen, gen.problem)

    traj = ST.get_closed_loop_trajectory(
        discrete_time_system,
        concrete_controller,
        x0,
        gen.nstep;
        trajectory_success = xtraj -> PR.trajectory_success(gen.problem, xtraj),
    )

    return ST.ClosedLoopTrajectory(traj.x, traj.u)
end

function _generate_abstract_trajectory(gen::TrajectoryGenerator)
    optimizer = gen.optimizer
    abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
    autom = abstract_problem.system
    abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    abstract_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))

    q0 = SY.get_abstract_state(abstract_system, _initial_state(gen, gen.problem))

    traj_abs = SY.get_closed_loop_trajectory(
        autom,
        abstract_controller,
        q0,
        gen.nstep;
        trajectory_success = qtraj -> PR.trajectory_success(abstract_problem, qtraj),
    )

    xs = [SY.get_concrete_state(abstract_system, q) for q in traj_abs.x.seq]

    us = [SY.get_concrete_input(abstract_system, u) for u in traj_abs.u.seq]

    return ST.ClosedLoopTrajectory(ST.Trajectory(xs), ST.Trajectory(us))
end

end # end module
