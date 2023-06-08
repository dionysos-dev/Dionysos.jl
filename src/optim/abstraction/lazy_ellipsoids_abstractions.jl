export LazyEllipsoidsAbstractions

module LazyEllipsoidsAbstractions

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const CO = DI.Control
const SY = DI.Symbolic
const PR = DI.Problem

using JuMP

mutable struct Optimizer <: MOI.AbstractOptimizer
    problem #::Union{Nothing, PR.OptimalControlProblem}
    distance
    rand_state
    new_conf
    keep
    stop_crit
    RRTstar::Bool
    compute_transition
    maxIter::Int    
    maxδx
    maxδu
    λ
    sdp_opt
    k1
    k2
    tree::Union{Nothing, UT.Tree} #later we could create a symmodel with overelapping cells

    # sdp_solver::Union{Nothing, MOI.OptimizerWithAttributes}
end

function build_OptimizerLazyEllipsoids(problem, distance, rand_state, new_conf, keep, stop_crit, RRTstar, compute_transition, maxIter, maxδx, maxδu, λ, sdp_opt, k1, k2)
    return Optimizer(problem, distance, rand_state, new_conf, keep, stop_crit, RRTstar, compute_transition, maxIter, maxδx, maxδu, λ, sdp_opt, k1, k2, nothing)
end

Optimizer() = Optimizer{Float64}()

MOI.is_empty(optimizer::Optimizer) = optimizer.problem === nothing

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    setproperty!(model, Symbol(param.name), value)
end
function MOI.get(model::Optimizer, param::MOI.RawOptimizerAttribute)
    getproperty(model, Symbol(param.name))
end
 
function MOI.optimize!(optimizer::Optimizer)
    problem = optimizer.problem
    Einit = problem.initial_set
    Etarget = problem.target_set

    distance = optimizer.distance
    rand_state = optimizer.rand_state
    new_conf = optimizer.new_conf
    keep = optimizer.keep
    stop_crit = optimizer.stop_crit
    maxIter = optimizer.maxIter
    RRTstar = optimizer.RRTstar
    compute_transition = optimizer.compute_transition
    k1 = optimizer.k1
    k2 = optimizer.k2
    tree = UT.RRT(Etarget, Einit, distance, rand_state, new_conf, keep, stop_crit, optimizer; maxIter=maxIter, RRTstar=RRTstar, compute_transition, k1=k1, k2=k2)
    
    optimizer.tree = tree 
    return 
end

end