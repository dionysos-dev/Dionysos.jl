export BranchAndBound

using DataStructures # for BinaryMinHeap

struct BranchAndBound{C, M, V<:AbstractVector{Int}}
    continuous_solver::C
    mixed_integer_solver::M
    horizon::Int
    indicator::Bool
    log_level::Int
    log_iter::V
    max_iter::Int
    max_time::Float64
    rel_gap::Float64
    feasible_solution_callback::Union{Nothing, Function}
end
function BranchAndBound(continuous_solver, mixed_integer_solver; horizon=0, max_iter=1000, max_time=60.0, rel_gap=1e-3, indicator::Bool = false, log_level = 1, log_iter = 0:100:max_iter,
                        feasible_solution_callback = nothing)
    return BranchAndBound(continuous_solver, mixed_integer_solver, horizon, indicator, log_level, log_iter, max_iter, max_time, rel_gap, feasible_solution_callback)
end

struct Candidate{T, TT}
    lower_bound::T
    upper_bound::T
    traj::DiscreteTrajectory{TT}
end

function Base.isless(a::Candidate, b::Candidate)
    if length(a.traj) == length(b.traj)
        return isless(a.lower_bound, b.lower_bound)
    else
        return length(a.traj) > length(b.traj)
    end
end

# TODO It assumes that the cost is `Fill{...}`, otherwise, we should compute a different
#      cost for each time
function minimum_transition_cost(prob, transition, solver, log_level = 0)
    model = Model(solver)
    from = source(prob.system, transition)
    to = target(prob.system, transition)
    @variable(model, x0[1:statedim(prob.system, from)] in hrep(stateset(prob.system, from)))
    @variable(model, x1[1:statedim(prob.system, to)] in hrep(stateset(prob.system, to)))
    @variable(model, u[1:inputdim(resetmap(prob.system, transition))])
    # We set `indicator` to `true`, there is no integer variables anyway so it does not matter.
    algo = BemporadMorari(solver, solver, true, 0)
    # We use `1` as we asssume the cost is the same along time
    t = 1
    δ_mode = IndicatorVariables([to], t)
    state_cost, δ_mode = hybrid_cost(model, fillify(prob.state_cost[t][[to]]), x1, u, δ_mode)
    symbols = [symbol(prob.system, transition)]
    δ_trans = IndicatorVariables(symbols, t)
    δ_trans = hybrid_constraints(model, fillify(prob.system.resetmaps[symbols]), x0, x1, u, algo, δ_trans)
    trans_cost, δ_trans = hybrid_cost(model, fillify(prob.transition_cost[t][symbols]), x0, u, δ_trans)
    @objective(model, Min, state_cost + trans_cost)
    optimize!(model)
    if termination_status(model) == MOI.OPTIMAL
        return objective_value(model)
    elseif termination_status(model) == MOI.INFEASIBLE
        return Inf
    else
        @error("Termination status: $(termination_status(model)), raw status: $(raw_status(model))")
    end
end

function distances(prob, solver)
    syst = prob.system
    dists = JuMP.Containers.@container([0:prob.number_of_time_steps, modes(syst)], Inf)
    dists[0, prob.q_T] = 0.0
    transition_cost = HybridSystems.transition_property(syst, Float64)
    for t in transitions(syst)
        transition_cost[t] = minimum_transition_cost(prob, t, solver)
    end
    for i in 1:prob.number_of_time_steps
        for mode in modes(syst)
            for t in out_transitions(syst, mode)
                dists[i, mode] = min(dists[i, mode], transition_cost[t] + dists[i - 1, target(syst, t)])
            end
        end
    end
    return dists
end

function candidate(prob, algo, discrete_lb, traj)
    sub_algo = BemporadMorari(algo.continuous_solver, algo.mixed_integer_solver, algo.indicator, algo.log_level)
    modes = [[target(prob.system, t)] for t in traj.transitions]
    cont_traj_prob = OptimalControlProblem(
        prob.system, prob.q_0, prob.x_0,
        prob.state_cost[1:length(traj)], prob.transition_cost[1:length(traj)],
        last_mode(prob.system, traj), length(traj)
    )
    sol = optimal_control(cont_traj_prob, sub_algo, modes)
    sol === nothing && return
    left = prob.number_of_time_steps - length(traj)
    lb = sol[2] + discrete_lb[left, last_mode(prob.system, traj)]
    horizon_prob = OptimalControlProblem(
        prob.system, last_mode(prob.system, traj), length(traj) == 0 ? prob.x_0 : sol[1].x[end],
        prob.state_cost[(length(traj) + 1):end], prob.transition_cost[(length(traj) + 1):end],
        prob.q_T, min(left, algo.horizon)
    )
    sol_horizon = optimal_control(horizon_prob, sub_algo)
    if sol_horizon === nothing
        ub = Inf
        sol_traj = nothing
    else
        ub = sol[2] + sol_horizon[2]
        sol_traj = ContinuousTrajectory(
            [sol[1].x; sol_horizon[1].x],
            [sol[1].u; sol_horizon[1].u]
        )
    end
    return Candidate(lb, ub, traj), sol_traj
end

using Printf

# Inspired from Pavito's `printgap`
function print_info(log_level, log_iter, num_iter, ub, start_time)
    if log_level >= 1
        if num_iter == 1 || log_level >= 2
            @printf "\n%-5s | %-14s | %-11s\n" "Iter." "Best feasible" "Time (s)"
        end
        if num_iter in log_iter
            @printf "%5d | %+14.6e | %11.3e\n" num_iter ub (time() - start_time)
        end
        flush(stdout)
        flush(stderr)
    end
    return
end

function optimal_control(
    prob::OptimalControlProblem,
    algo::BranchAndBound
)
    start_time = time()
    iszero(prob.number_of_time_steps) && return _zero_steps(prob)
    discrete_lb = distances(prob, algo.continuous_solver)
    candidate_0, best_traj = candidate(prob, algo, discrete_lb, DiscreteTrajectory{transitiontype(prob.system)}(prob.q_0))
    candidate_0 === nothing && return
    ub = candidate_0.upper_bound
    num_iter = 0
    candidates = BinaryMinHeap([candidate_0])
    while !isempty(candidates) && num_iter < algo.max_iter && time() - start_time < algo.max_time
        num_iter += 1
        cur_candidate = pop!(candidates)
        if cur_candidate.lower_bound < ub
            for t in out_transitions(prob.system, last_mode(prob.system, cur_candidate.traj))
                new_candidate_traj = candidate(prob, algo, discrete_lb, append(cur_candidate.traj, t))
                if new_candidate_traj !== nothing
                    new_candidate, sol_traj = new_candidate_traj
                    if new_candidate.upper_bound < ub
                        if algo.feasible_solution_callback !== nothing
                            algo.feasible_solution_callback(new_candidate, sol_traj)
                        end
                        best_traj = sol_traj
                        ub = new_candidate.upper_bound
                    end
                    if length(new_candidate.traj) < prob.number_of_time_steps
                        push!(candidates, new_candidate)
                    end
                end
            end
        end
        print_info(algo.log_level, algo.log_iter, num_iter, ub, start_time)
    end
    best_traj === nothing && return
    return best_traj, ub
end
