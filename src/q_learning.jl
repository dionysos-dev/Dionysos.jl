using LinearAlgebra
using Polyhedra, CDDLib, ParameterJuMP

export DiscreteLowerBoundAlgo, HybridDualDynamicProgrammingAlgo

# TODO It assumes that the cost is `Fill{...}`, otherwise, we should compute a different
#      cost for each time
function minimum_transition_cost(prob, transition, solver, log_level = 0)
    model = Model(solver)
    from = source(prob.system, transition)
    to = target(prob.system, transition)
    @variable(model, x0[1:statedim(prob.system, from)] in hrep(stateset(prob.system, from)))
    @variable(model, x1[1:statedim(prob.system, to)] in hrep(stateset(prob.system, to)))
    @variable(model, u[1:inputdim(resetmap(prob.system, transition))])
    algo = optimizer_with_attributes(
        BemporadMorari.Optimizer,
        "continuous_solver" => solver,
        "log_level" => 0)
    # We use `1` as we asssume the cost is the same along time
    t = 1
    δ_mode = BemporadMorari.IndicatorVariables([to], t)
    state_cost, δ_mode = BemporadMorari.hybrid_cost(model, BemporadMorari.fillify(prob.state_cost[t][[to]]), x1, u, δ_mode)
    symbols = [symbol(prob.system, transition)]
    δ_trans = BemporadMorari.IndicatorVariables(symbols, t)
    δ_trans = BemporadMorari.hybrid_constraints(model, BemporadMorari.fillify(prob.system.resetmaps[symbols]), x0, x1, u, algo, δ_trans)
    trans_cost, δ_trans = BemporadMorari.hybrid_cost(model, BemporadMorari.fillify(prob.transition_cost[t][symbols]), x0, u, δ_trans)
    @objective(model, Min, state_cost + trans_cost)
    optimize!(model)
    if termination_status(model) == MOI.OPTIMAL
        return objective_value(model)
    elseif termination_status(model) in (MOI.INFEASIBLE,)#, MOI.INFEASIBLE_OR_UNBOUNDED)
        return Inf
    else
        @error("Candidate: Termination status: $(termination_status(model)), raw status: $(raw_status(model))")
    end
end



function _merge_with(combine, a, b)
    n = (axes(a, 1).stop, axes(b, 1).stop)
    number_of_time_steps = maximum(n)
    modes = axes(a, 2)
    return JuMP.Containers.@container(
        [i in 0:number_of_time_steps, mode in modes],
        if i > n[1]
            b[i, mode]
        elseif i > n[2]
            a[i, mode]
        else
            combine(a[i, mode], b[i, mode])
        end
    )
end

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
function q_merge(a::DiscreteLowerBound, b::DiscreteLowerBound)
    return DiscreteLowerBound(_merge_with(max, a.discrete_lb, b.discrete_lb))
end
function value_function(Q::DiscreteLowerBound, left::Int, mode)
    return ConstantFunction(Q.discrete_lb[left, mode])
end
function learn(::DiscreteLowerBound, prob, ::DiscreteTrajectory, ::ContinuousTrajectory, ::DiscreteLowerBoundAlgo) end

struct HybridDualDynamicProgrammingAlgo{S, T}
    solver::S
    new_cut_tol::T
end
struct HybridDualDynamicProgramming{C, D}
    cuts::C
    discrete::DiscreteLowerBound{D}
end
function _no_cuts(number_of_time_steps, modes)
    return JuMP.Containers.@container(
        [0:number_of_time_steps, modes],
        AffineFunction{Float64}[]
    )
end
function instantiate(prob::OptimalControlProblem, algo::HybridDualDynamicProgrammingAlgo)
    cuts = _no_cuts(prob.number_of_time_steps, modes(prob.system))
    discrete = instantiate(prob, DiscreteLowerBoundAlgo(algo.solver))
    return HybridDualDynamicProgramming(cuts, discrete)
end
function q_merge(a::DiscreteLowerBound, b::HybridDualDynamicProgramming)
    d = a.discrete_lb
    cuts = _no_cuts(axes(d, 1).stop, axes(d, 2))
    return q_merge(HybridDualDynamicProgramming(cuts, a), b)
end
function q_merge(a::HybridDualDynamicProgramming, b::HybridDualDynamicProgramming)
    return HybridDualDynamicProgramming(_merge_with(vcat, a.cuts, b.cuts), q_merge(a.discrete, b.discrete))
end
function value_function(Q::HybridDualDynamicProgramming, left::Int, mode)
    return PolyhedralFunction(Q.discrete.discrete_lb[left, mode], Q.cuts[left, mode])
end
function vertices(fs::Vector{PolyhedralFunction{T}}, X) where T
    h = hrep(X) * intersect(HalfSpace([-one(T)], -fs[1].lower_bound))
    cuts = HalfSpace{T, Vector{T}}[]
    for f in fs
        append!(cuts, [HalfSpace([p.a; -one(T)], -p.β) for p in f.pieces])
    end
    h_cut = h ∩ hrep(cuts, d = fulldim(h))
    p = polyhedron(h_cut, CDDLib.Library())
    removehredundancy!(p)
    return collect(points(vrep(p)))
end
function vertices(Q::HybridDualDynamicProgramming, prob, left, mode)
    return vertices([value_function(Q, i, mode) for i in left:prob.number_of_time_steps], stateset(prob.system, mode))
end
function learn(::HybridDualDynamicProgramming, prob, ::DiscreteTrajectory, ::ContinuousTrajectory,
               ::DiscreteLowerBoundAlgo)
end
function learn(Q::HybridDualDynamicProgramming, prob, dtraj::DiscreteTrajectory, ctraj::ContinuousTrajectory,
               algo::HybridDualDynamicProgrammingAlgo)
    for i in length(dtraj):-1:1
        x = i == 1 ? prob.x_0 : ctraj.x[i - 1]
        mode = source(prob.system, dtraj.transitions[i])
        left = prob.number_of_time_steps - i
        trans = filter(collect(out_transitions(prob.system, mode))) do t
            Q.discrete.discrete_lb[left, target(prob.system, t)] != Inf
        end
        verts = JuMP.Containers.@container(
            [t in trans],
            vertices(Q, prob, left, target(prob.system, t))
        )
        model = ModelWithParams(algo.solver)
        params = add_parameters(model, x)
        # FIXME assume they are all the same
        r = resetmap(prob.system, first(trans))
        @variable(model, u[1:inputdim(r)] in inputset(r))
        @variable(model, λ[t in trans, eachindex(verts[t])] ≥ 0)
        @constraint(model, sum(λ) == 1)
        @variable(model, θ)
        epi(i) = sum(λ[t, v] * verts[t][v][i] for t in trans for v in eachindex(verts[t]))
        x_next = r.A * params + r.B * u
        for j in eachindex(x_next)
            @constraint(model, x_next[j] == epi(j))
        end
        @constraint(model, θ >= epi(length(x_next) + 1))
        tos = [target(prob.system, t) for t in trans]
        δ_mode = JuMP.Containers.@container([mode in tos],
            sum(λ[t, v] for t in trans for v in eachindex(verts[t]) if target(prob.system, t) == mode)
        )
        state_cost, δ_mode = BemporadMorari.hybrid_cost(model, BemporadMorari.fillify(prob.state_cost[i][tos]), x_next, u, δ_mode)
        symbols = symbol.(prob.system, trans)
        δ_trans = JuMP.Containers.@container([s in symbols],
            sum(λ[t, v] for t in trans for v in eachindex(verts[t]) if symbol(prob.system, t) == s)
        )
        trans_cost, δ_trans = BemporadMorari.hybrid_cost(model, BemporadMorari.fillify(prob.transition_cost[i][symbols]), x, u, δ_trans)
        @objective(model, Min, trans_cost + state_cost + θ)
        optimize!(model)
        if termination_status(model) == MOI.OPTIMAL && dual_status(model) == MOI.FEASIBLE_POINT
            a = dual.(params)
            # OSQP does not support accessing `DualObjectiveValue`
            #β = dual_objective_value(model) - a ⋅ x
            β = objective_value(model) - a ⋅ x
            # rounding to avoid numerical inconsistencies with CDD
            #cut = AffineFunction(round.(a, digits=8), round(β, digits=8))
            cut = AffineFunction(a, β)
            V = value_function(Q, left + 1, mode)
            before = function_value(V, x)
            after = function_value(cut, x)
            # without the tolerance `1e-5`, CDDLib often throws `Numerically inconsistent`.
            if after < before + algo.new_cut_tol
                @info("Cut ignored: $cut, $after ≤ $before.")
            else
                @info("Cut added: $cut, $after > $before.")
                push!(Q.cuts[left + 1, mode], cut)
            end
        else
            @warn("No new cut generated as termination status is $(termination_status(model)), primal status is $(primal_status(model)), dual status is $(dual_status(model)): $(raw_status(model))")
        end
    end
end
