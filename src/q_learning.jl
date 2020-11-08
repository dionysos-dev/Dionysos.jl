export DiscreteLowerBoundAlgo, HybridDualDynamicProgrammingAlgo

struct DiscreteLowerBoundAlgo{S}
    solver::S
end
struct DiscreteLowerBound{D}
    discrete_lb::D
end
function instantiate(prob::OptimalControlProblem, algo::DiscreteLowerBoundAlgo)
    syst = prob.system
    dists = JuMP.Containers.@container([0:prob.number_of_time_steps, modes(syst)], Inf)
    dists[0, prob.q_T] = 0.0
    transition_cost = HybridSystems.transition_property(syst, Float64)
    for t in transitions(syst)
        transition_cost[t] = minimum_transition_cost(prob, t, algo.solver)
    end
    for i in 1:prob.number_of_time_steps
        for mode in modes(syst)
            for t in out_transitions(syst, mode)
                dists[i, mode] = min(dists[i, mode], transition_cost[t] + dists[i - 1, target(syst, t)])
            end
        end
    end
    return DiscreteLowerBound(dists)
end
function value_function(Q::DiscreteLowerBound, left::Int, mode)
    return ConstantFunction(Q.discrete_lb[left, mode])
end

struct HybridDualDynamicProgrammingAlgo{S}
    solver::S
end
struct HybridDualDynamicProgramming{C, D}
    cuts::C
    discrete::DiscreteLowerBound{D}
end
function instantiate(prob::OptimalControlProblem, algo::HybridDualDynamicProgrammingAlgo)
    cuts = JuMP.Containers.@container(
        [0:prob.number_of_time_steps, mode = modes(prob.system)],
        AffineFunction{Float64}[]
    )
    discrete = instantiate(prob, DiscreteLowerBoundAlgo(algo.solver))
    return HybridDualDynamicProgramming(cuts, discrete)
end
function value_function(Q::HybridDualDynamicProgramming, left::Int, mode)
    return PolyhedralFunction(Q.discrete.discrete_lb[left, mode], Q.cuts[left, mode])
end
