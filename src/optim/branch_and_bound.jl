export BranchAndBound

module BranchAndBound

import Dionysos
const DI = Dionysos
const CO = DI.Control
const PR = DI.Problem
const OP = DI.Optim

using DataStructures # for BinaryMinHeap
using JuMP, Polyhedra
using HybridSystems

mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    continuous_solver::Any
    mixed_integer_solver::Any
    Q_function_init::Any
    Q_function::Any
    lower_bound::Any
    horizon::Int
    indicator::Bool
    log_level::Int
    log_iter::AbstractVector{Int}
    max_iter::Int
    max_time::Float64
    feasible_solution_callback::Union{Nothing, Function}
    problem::Union{Nothing, PR.OptimalControlProblem}
    status::MOI.TerminationStatusCode
    upper_bound::T
    solve_time::Float64
    num_total::Int
    num_iter::Int
    num_done::Int
    num_pruned_bound::Int
    num_pruned_inf::Int
    best_traj::Union{Nothing, CO.ContinuousTrajectory}
    function Optimizer{T}() where {T}
        return new{T}(
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            0,
            false,
            1,
            0:100:typemax(Int),
            1000,
            60.0,
            nothing,
            nothing,
            MOI.OPTIMIZE_NOT_CALLED,
            NaN,
            NaN,
            0,
            0,
            0,
            0,
            0,
            nothing,
        )
    end
end

MOI.is_empty(optimizer::Optimizer) = optimizer.problem === nothing

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    if param.name == "problem"
        if !(value isa PR.OptimalControlProblem)
            throw(MOI.UnsupportedAttribute(param, "$(typeof(value)) not supported"))
        end
    end
    return setproperty!(model, Symbol(param.name), value)
end
function MOI.get(model::Optimizer, param::MOI.RawOptimizerAttribute)
    return getproperty(model, Symbol(param.name))
end

struct Candidate{T, TT}
    lower_bound::T
    upper_bound::T
    traj::CO.DiscreteTrajectory{TT}
end

function Base.isless(a::Candidate, b::Candidate)
    if length(a.traj) == length(b.traj)
        return isless(a.lower_bound, b.lower_bound)
    else
        return length(a.traj) > length(b.traj)
    end
end

function number_of_nodes(prob)
    num = JuMP.Containers.@container([0:(prob.time), modes(prob.system)], 0)
    num[0, prob.target_set] = 1
    for i in 1:(prob.time)
        for mode in modes(prob.system)
            for t in out_transitions(prob.system, mode)
                num[i, mode] += num[i - 1, target(prob.system, t)]
            end
            if num[i, mode] > 0
                # This is feasible so we count this mode as a node at time `i`.
                num[i, mode] += 1
            end
        end
    end
    return num
end

function candidate(prob, algo::Optimizer{T}, Q_function, traj) where {T}
    sub_algo = optimizer_with_attributes(
        OP.BemporadMorari.Optimizer{T},
        "continuous_solver" => algo.continuous_solver,
        "mixed_integer_solver" => algo.mixed_integer_solver,
        "indicator" => algo.indicator,
        "log_level" => algo.log_level,
    )
    modes = [[target(prob.system, t)] for t in traj.transitions]
    state_cost = collect(prob.state_cost[1:length(traj)])
    left = prob.time - length(traj)
    terminal_cost = OP.value_function(Q_function, left, CO.last_mode(prob.system, traj))
    @assert isempty(state_cost) == iszero(length(traj))
    if !isempty(state_cost)
        # We use `Ref` as `Base.broadcastable` is not defined for functions that
        # are in `terminal_cost`.
        state_cost = [state_cost[1:(end - 1)]; [state_cost[end] .+ Ref(terminal_cost)]]
    end
    cont_traj_prob = PR.OptimalControlProblem(
        prob.system,
        prob.initial_set,
        CO.last_mode(prob.system, traj),
        state_cost,
        prob.transition_cost[1:length(traj)],
        length(traj),
    )
    start_optimizer = MOI.instantiate(sub_algo)
    MOI.set(start_optimizer, MOI.RawOptimizerAttribute("modes"), modes)
    MOI.set(start_optimizer, MOI.RawOptimizerAttribute("problem"), cont_traj_prob)
    MOI.optimize!(start_optimizer)
    if MOI.get(start_optimizer, MOI.TerminationStatus()) != MOI.OPTIMAL
        return
    end
    left = prob.time - length(traj)
    start_cost = MOI.get(start_optimizer, MOI.ObjectiveValue())
    start_sol = MOI.get(start_optimizer, CO.ContinuousTrajectoryAttribute())
    if isempty(state_cost)
        lb = start_cost + CO.function_value(terminal_cost, prob.initial_set[2])
    else
        lb = start_cost
        start_cost -= CO.function_value(terminal_cost, start_sol.x[end])
    end
    if left <= algo.horizon || false # prob.allow_less_iterations
        horizon_prob = PR.OptimalControlProblem(
            prob.system,
            (
                CO.last_mode(prob.system, traj),
                length(traj) == 0 ? prob.initial_set[2] : start_sol.x[end],
            ),
            prob.target_set,
            prob.state_cost[(length(traj) + 1):end],
            prob.transition_cost[(length(traj) + 1):end],
            min(left, algo.horizon),
        )
        horizon_optimizer = MOI.instantiate(sub_algo)
        MOI.set(horizon_optimizer, MOI.RawOptimizerAttribute("problem"), horizon_prob)
        MOI.optimize!(horizon_optimizer)
        status_horizon = MOI.get(horizon_optimizer, MOI.PrimalStatus())
    else
        status_horizon = MOI.NO_SOLUTION
    end
    if status_horizon == MOI.NO_SOLUTION
        ub = typemax(T) # ∞
        sol_traj = nothing
    else
        ub = start_cost + MOI.get(horizon_optimizer, MOI.ObjectiveValue())
        sol_horizon = MOI.get(horizon_optimizer, CO.ContinuousTrajectoryAttribute())
        sol_traj = CO.ContinuousTrajectory(
            [start_sol.x; sol_horizon.x],
            [start_sol.u; sol_horizon.u],
        )
    end
    return Candidate(lb, ub, traj), sol_traj
end

using Printf

const NUM_NODES = ["#nodes", "#queued", "#done", "#pruned", "#infeas", "#left"]

# Inspired from Pavito's `printgap`
function print_info(optimizer::Optimizer, last_iter::Bool, start_time, num_queued)
    num_nodes = [
        optimizer.num_total,
        num_queued,
        optimizer.num_done,
        optimizer.num_pruned_bound,
        optimizer.num_pruned_inf,
    ]
    push!(num_nodes, 2optimizer.num_total - sum(num_nodes))
    @assert num_nodes[end] >= 0
    if optimizer.log_level >= 1
        len = max(maximum(length, NUM_NODES), length(string(num_nodes[1])))
        if optimizer.num_iter == 0 || optimizer.log_level >= 2
            @printf "\n%-5s | %-14s | %-11s" "Iter." "Best feasible" "Time (s)"
            for s in NUM_NODES
                print(" | ", lpad(s, len))
            end
            println()
        end
        if last_iter || (optimizer.num_iter in optimizer.log_iter)
            @printf "%5d | %+14.6e | %11.3e" optimizer.num_iter optimizer.upper_bound (
                time() - start_time
            )
            for num in num_nodes
                print(" | ", lpad(string(num), len))
            end
            println()
        end
        flush(stdout)
        flush(stderr)
    end
    return
end

function MOI.optimize!(optimizer::Optimizer{T}) where {T}
    start_time = time()
    prob = optimizer.problem
    if iszero(prob.time)
        if optimizer.problem.target_set == optimizer.problem.initial_set[1]
            optimizer.status = MOI.OPTIMAL
            optimizer.best_traj = CO.ContinuousTrajectory(Vector{T}[], Vector{T}[])
            optimizer.upper_bound = zero(T)
        else
            optimizer.status = MOI.INFEASIBLE
        end
        return
    end
    num_nodes = number_of_nodes(prob)
    optimizer.num_total = num_nodes[prob.time, prob.initial_set[1]]
    if iszero(optimizer.num_total)
        optimizer.status = MOI.INFEASIBLE
        return
    end
    if optimizer.lower_bound === nothing
        optimizer.lower_bound = OP.DiscreteLowerBoundAlgo{T}(optimizer.continuous_solver)
    end
    optimizer.Q_function = OP.instantiate(prob, optimizer.lower_bound)
    if optimizer.Q_function_init !== nothing
        optimizer.Q_function = OP.q_merge(optimizer.Q_function, optimizer.Q_function_init)
    end
    candidate_0, optimizer.best_traj = candidate(
        prob,
        optimizer,
        optimizer.Q_function,
        CO.DiscreteTrajectory{transitiontype(prob.system)}(prob.initial_set[1]),
    )
    candidate_0 === nothing && return
    optimizer.upper_bound = candidate_0.upper_bound
    optimizer.num_iter = 0
    candidates = BinaryMinHeap([candidate_0])
    optimizer.num_done = 0
    optimizer.num_pruned_bound = 0
    optimizer.num_pruned_inf = 0
    while true
        if isempty(candidates)
            optimizer.status =
                optimizer.best_traj === nothing ? MOI.INFEASIBLE : MOI.OPTIMAL
        elseif optimizer.num_iter >= optimizer.max_iter
            optimizer.status = MOI.ITERATION_LIMIT
        elseif time() - start_time >= optimizer.max_time
            optimizer.status = MOI.TIME_LIMIT
        end
        stop = optimizer.status != MOI.OPTIMIZE_NOT_CALLED
        print_info(optimizer, stop, start_time, length(candidates))
        stop && break
        optimizer.num_iter += 1
        cur_candidate = pop!(candidates)
        if cur_candidate.lower_bound < optimizer.upper_bound
            for t in
                out_transitions(prob.system, CO.last_mode(prob.system, cur_candidate.traj))
                new_traj = CO.append(cur_candidate.traj, t)
                # Shortcuts avoiding to solve the QP in case, it's infeasible.
                nnodes = num_nodes[
                    prob.time - length(new_traj),
                    CO.last_mode(prob.system, new_traj),
                ]
                iszero(nnodes) && continue
                new_candidate_traj =
                    candidate(prob, optimizer, optimizer.Q_function, new_traj)
                if new_candidate_traj === nothing
                    optimizer.num_pruned_inf += nnodes - 1
                    optimizer.num_done += 1
                else
                    new_candidate, sol_traj = new_candidate_traj
                    if sol_traj !== nothing
                        # TODO Instead of `new_traj`, we should complete it
                        # with the discrete solution of the last part.
                        # If `horizon == 0`, there is no difference but otherwise,
                        # it does matter
                        OP.learn(
                            optimizer.Q_function,
                            prob,
                            new_traj,
                            sol_traj,
                            optimizer.lower_bound,
                        )
                    end
                    if new_candidate.upper_bound < optimizer.upper_bound
                        if optimizer.feasible_solution_callback !== nothing
                            optimizer.feasible_solution_callback(new_candidate, sol_traj)
                        end
                        optimizer.best_traj = sol_traj
                        optimizer.upper_bound = new_candidate.upper_bound
                    end
                    if length(new_candidate.traj) < prob.time
                        if new_candidate.lower_bound < optimizer.upper_bound
                            push!(candidates, new_candidate)
                        else
                            optimizer.num_pruned_bound += nnodes - 1
                            optimizer.num_done += 1
                        end
                    else
                        optimizer.num_done += 1
                    end
                end
            end
        else
            left = prob.time - length(cur_candidate.traj)
            nnodes = num_nodes[left, CO.last_mode(prob.system, cur_candidate.traj)]
            optimizer.num_pruned_bound += nnodes - 1
        end
        optimizer.num_done += 1
    end
end

function MOI.get(optimizer::Optimizer, ::CO.ContinuousTrajectoryAttribute)
    return optimizer.best_traj
end
function MOI.get(optimizer::Optimizer, ::MOI.ObjectiveValue)
    return optimizer.upper_bound
end
function MOI.get(optimizer::Optimizer, ::MOI.TerminationStatus)
    return optimizer.status
end

end
