module BranchAndBound
using DataStructures, JuMP
using Printf

export Abstract_BB_Problem

mutable struct Node
    elem::Any
    depth::Any
    lower_bound::Any
    upper_bound::Any
    sol::Any
    parent::Union{Nothing, Node}
    ext::Any
end

function Node(
    elem,
    depth;
    lower_bound = -Inf,
    upper_bound = Inf,
    sol = nothing,
    parent = nothing,
    ext = nothing,
)
    return Node(elem, depth, lower_bound, upper_bound, sol, parent, ext)
end

# function used in the priority queue of the BB
# highest depth and best lower bound (optimistic DFS)
function Base.isless(a::Node, b::Node)
    #return isless(a.lower_bound, b.lower_bound)
    if a.depth == b.depth
        return isless(a.lower_bound, b.lower_bound)
    else
        return a.depth > b.depth
    end
end

"""
    Abstract_BB_Problem #abstract branch and bound problem
Sould implement the methods below
"""
abstract type Abstract_BB_Problem end

# function to check if the problem is feasible at the beginning
function check_trivial_infeasibility(prob::Abstract_BB_Problem)
    return false
end
# return the initial instance of the problem
function get_first_instance(prob::Abstract_BB_Problem) end
# compute an upper bound, update the field upper_bound and sol of node
function compute_upper_bound!(prob::Abstract_BB_Problem, node::Node) end
# compute an lower bound, update the field lower_bound of node
function compute_lower_bound!(prob::Abstract_BB_Problem, node::Node) end
# return the children nodes of node
function expand(prob::Abstract_BB_Problem, node::Node) end

mutable struct Optimizer <: MOI.AbstractOptimizer
    problem::Abstract_BB_Problem
    lower_bound::Float64
    upper_bound::Float64
    best_sol::Any
    solve_time::Float64
    num_iter::Int64
    max_iter::Int
    max_time::Float64
    status::MOI.TerminationStatusCode
    num_total::Int         # number of node created
    num_pruned_inf::Int    # number of node pruned (by infeasibility)
    num_pruned_bound::Int  # number of node pruned (by optimality)
    log_level::Int
end
function Optimizer(problem, max_iter, max_time; log_level = 0)
    return Optimizer(
        problem,
        0.0,
        Inf,
        nothing,
        0.0,
        0,
        max_iter,
        max_time,
        MOI.OPTIMIZE_NOT_CALLED,
        0,
        0,
        0,
        log_level,
    )
end

has_solution(optimizer::Optimizer) = optimizer.best_sol !== nothing

function MOI.optimize!(optimizer::Optimizer)
    start_time = time()
    prob = optimizer.problem

    if check_trivial_infeasibility(prob)
        optimizer.status = MOI.INFEASIBLE
        return
    end

    node_0 = Node(get_first_instance(prob), 0)
    optimizer.num_total += 1
    compute_lower_bound!(prob, node_0)

    candidates = BinaryMinHeap([node_0])
    while true
        if isempty(candidates)
            optimizer.status = optimizer.best_sol === nothing ? MOI.INFEASIBLE : MOI.OPTIMAL
        elseif optimizer.num_iter >= optimizer.max_iter
            optimizer.status = MOI.ITERATION_LIMIT
        elseif time() - start_time >= optimizer.max_time
            optimizer.status = MOI.TIME_LIMIT
        end
        stop = optimizer.status != MOI.OPTIMIZE_NOT_CALLED
        print_info(optimizer, stop, start_time, length(candidates))
        stop && break
        optimizer.num_iter += 1
        cur_node = pop!(candidates)
        println(cur_node.elem)
        if cur_node.lower_bound < optimizer.upper_bound
            compute_upper_bound!(prob, cur_node)
            if cur_node.upper_bound == -Inf
                optimizer.num_pruned_inf += 1
            else
                if cur_node.upper_bound < optimizer.upper_bound
                    optimizer.upper_bound = cur_node.upper_bound
                    optimizer.best_sol = cur_node.sol
                end
                children = expand(prob, cur_node)
                optimizer.num_total += length(children)
                for child in children # necessary for the PQ
                    compute_lower_bound!(prob, child)
                    if child.lower_bound == Inf
                        optimizer.num_pruned_inf += 1
                    end
                end

                for child in children # (not really necessary)
                    if child.lower_bound < optimizer.upper_bound
                        push!(candidates, child)
                    else
                        optimizer.num_pruned_bound += 1
                    end
                end
            end
        else
            optimizer.num_pruned_bound += 1
        end
    end
end

const NUM_NODES = ["#nodes", "#pruned", "#infeas", "#queued"]

function print_info(optimizer::Optimizer, last_iter::Bool, start_time, num_queued)
    num_nodes = [
        optimizer.num_total,
        optimizer.num_pruned_bound,
        optimizer.num_pruned_inf,
        num_queued,
    ]
    if optimizer.log_level >= 1
        len = max(maximum(length, NUM_NODES), length(string(num_nodes[1])))
        if optimizer.num_iter == 0 || optimizer.log_level >= 2
            @printf "\n%-5s | %-14s | %-11s" "Iter." "Best feasible" "Time (s)"
            for s in NUM_NODES
                print(" | ", lpad(s, len))
            end
            println()
        end
        if optimizer.log_level >= 2
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
end # end module
