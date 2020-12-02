using LinearAlgebra

import MutableArithmetics
const MA = MutableArithmetics

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
    tight_tol::T
    log_level::Int
end
struct HybridDualDynamicProgramming{C,H<:JuMP.Containers.DenseAxisArray{Polyhedra.Intersection{Float64,Vector{Float64},Int}}, D}
    cuts::C
    domains::H
    discrete::DiscreteLowerBound{D}
end
function _no_cuts(number_of_time_steps, modes)
    return JuMP.Containers.@container(
        [0:number_of_time_steps, modes],
        AffineFunction{Float64}[]
    )
end
function _full_domains(number_of_time_steps, modes, d)
    return JuMP.Containers.@container(
        [0:number_of_time_steps, modes],
        hrep(HalfSpace{Float64, Vector{Float64}}[], d = d)
    )
end
function instantiate(prob::OptimalControlProblem, algo::HybridDualDynamicProgrammingAlgo)
    cuts = _no_cuts(prob.number_of_time_steps, modes(prob.system))
    domains = _full_domains(prob.number_of_time_steps, modes(prob.system), length(prob.x_0))
    discrete = instantiate(prob, DiscreteLowerBoundAlgo(algo.solver))
    return HybridDualDynamicProgramming(cuts, domains, discrete)
end
function q_merge(a::DiscreteLowerBound, b::HybridDualDynamicProgramming)
    d = a.discrete_lb
    cuts = _no_cuts(axes(d, 1).stop, axes(d, 2))
    domains = _full_domains(axes(d, 1).stop, axes(d, 2), fulldim(first(b.domains)))
    return q_merge(HybridDualDynamicProgramming(cuts, domains, a), b)
end
function q_merge(a::HybridDualDynamicProgramming, b::HybridDualDynamicProgramming)
    return HybridDualDynamicProgramming(
        _merge_with(vcat, a.cuts, b.cuts),
        _merge_with(intersect, a.domains, b.domains),
        q_merge(a.discrete, b.discrete)
    )
end
function value_function(Q::HybridDualDynamicProgramming, left::Int, mode)
    d = fulldim(Q.domains[left, mode])
    number_of_time_steps = axes(Q.cuts, 1).stop
    return PolyhedralFunction(
        Q.discrete.discrete_lb[left, mode],
        reduce(append!, (Q.cuts[i, mode] for i in left:number_of_time_steps),
               init = AffineFunction{Float64}[]),
        # TODO use intersect! once it is implemented in Polyhedra
        reduce(intersect, (Q.domains[i, mode] for i in left:number_of_time_steps),
               init = hrep(HalfSpace{Float64, Vector{Float64}}[], d = d))
    )
end
function vertices(f::PolyhedralFunction{T}, X) where T
    h = (hrep(X) ∩ f.domain) * intersect(HalfSpace([-one(T)], -f.lower_bound))
    cuts = [HalfSpace([p.a; -one(T)], -p.β) for p in f.pieces]
    h_cut = h ∩ hrep(cuts, d = fulldim(h))
    p = polyhedron(h_cut, CDDLib.Library())
    removehredundancy!(p)
    return collect(points(vrep(p)))
end
function vertices(Q::HybridDualDynamicProgramming, prob, left, mode)
    return vertices(value_function(Q, left, mode), stateset(prob.system, mode))
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
        U = inputset(r)
        hashyperplanes(U) && error("TODO: Q-learning with input set with hyperplanes")
        @variable(model, u[1:fulldim(U)])
        U_h = collect(halfspaces(U))
        @constraint(model, u_con[i in eachindex(U_h)], u ⋅ U_h[i].a <= U_h[i].β)
        @variable(model, λ[t in trans, eachindex(verts[t])] ≥ 0)
        @constraint(model, sum(λ) == 1)
        @variable(model, θ)
        epi(i) = sum(λ[t, v] * verts[t][v][i] for t in trans for v in eachindex(verts[t]))
        x_next = r.A * params + r.B * u
        @constraint(model, epi_con[j in eachindex(x_next)], x_next[j] == epi(j))
        @constraint(model, θ >= epi(length(x_next) + 1))
        tos = [target(prob.system, t) for t in trans]
        _sum(g) = reduce(+, g, init = zero(JuMP.VariableRef))
        δ_mode = JuMP.Containers.@container([mode in tos],
            _sum(λ[t, v] for t in trans for v in eachindex(verts[t]) if target(prob.system, t) == mode)
        )
        state_cost, δ_mode = BemporadMorari.hybrid_cost(model, BemporadMorari.fillify(prob.state_cost[i][tos]), x_next, u, δ_mode)
        symbols = symbol.(prob.system, trans)
        δ_trans = JuMP.Containers.@container([s in symbols],
            _sum(λ[t, v] for t in trans for v in eachindex(verts[t]) if symbol(prob.system, t) == s)
        )
        trans_cost, δ_trans = BemporadMorari.hybrid_cost(model, BemporadMorari.fillify(prob.transition_cost[i][symbols]), x, u, δ_trans)
        obj = trans_cost + state_cost + θ
        @objective(model, Min, obj)
        optimize!(model)
        if termination_status(model) != MOI.OPTIMAL || primal_status(model) != MOI.FEASIBLE_POINT
            if algo.log_level >= 1
                @warn("No new cut generated as termination status is $(termination_status(model)), primal status is $(primal_status(model)), dual status is $(dual_status(model)): $(raw_status(model))")
            end
            return
        end
        V = value_function(Q, left + 1, mode)
        before = function_value(V, x)
        # FIXME OSQP does not support accessing `DualObjectiveValue`
        #after = dual_objective_value(model)
        after = objective_value(model)
        if after < before + algo.new_cut_tol
            if algo.log_level >= 3
                @info("Cuts ignored: $after ≤ $before.")
            end
            return
        end

        dual_model = Model()
        @variable(dual_model, y_sum)
        @variable(dual_model, y[1:length(epi_con)])
        λ_cons = JuMP.Containers.@container(
            [t in trans, v in eachindex(verts[t])],
            y_sum - sum(verts[t][v][i] * y[i] for i in eachindex(y)) - verts[t][v][length(x_next) + 1]
        )
        u_value = value.(u)
        u_cons = r.B' * y
        for (xy, coef) in obj.terms
            if xy.a == xy.b
                u_idx = findfirst(isequal(xy.a), u)
                @assert u_idx !== nothing
                MA.mutable_operate!(MA.add_mul, u_cons[u_idx], -2coef * u_value[u_idx])
            else
                ua_idx = findfirst(isequal(xy.a), u)
                ub_idx = findfirst(isequal(xy.b), u)
                if ua_idx !== nothing && ub_idx !== nothing
                    MA.mutable_operate!(MA.add_mul, u_cons[ua_idx], -coef * u_value[ub_idx])
                    MA.mutable_operate!(MA.add_mul, u_cons[ub_idx], -coef * u_value[ua_idx])
                else
                    error("TODO")
                end
            end
        end
        for (x, coef) in obj.aff.terms
            x == θ && continue
            u_idx = findfirst(isequal(x), u)
            if u_idx !== nothing
                MA.mutable_operate!(MA.add_mul, u_cons[u_idx], -coef)
            else
                found = false
                for t in trans, v in eachindex(verts[t])
                    if x == λ[t, v]
                        @assert !found
                        found = true
                        MA.mutable_operate!(MA.add_mul, λ_cons[t, v], -coef)
                    end
                end
                @assert found
            end
        end

        for t in trans, v in eachindex(verts[t])
            if value(λ[t, v]) <= algo.tight_tol
                @constraint(dual_model, λ_cons[t, v] <= 0)
            else
                @constraint(dual_model, λ_cons[t, v] == 0)
            end
        end
        U_v = vrep(
            [zeros(length(u))],
            Line{Float64, Vector{Float64}}[],
            Ray{Float64, Vector{Float64}}[]
        )
        for i in eachindex(U_h)
            α = u_value ⋅ U_h[i].a
            if α >= U_h[i].β - (max(abs(α), abs(U_h[i].β)) + 1) * algo.tight_tol
                # TODO Use `convexhull!` once it is implemented in Polyhedra
                convexhull!(U_v, Ray(U_h.a))
            end
        end
        U_p = polyhedron(U_v, CDDLib.Library())
        removehredundancy!(U_p)
        U_h = hrep(U_p)
        @constraint(dual_model, u_cons in U_h)

        h = _LPHRep(backend(dual_model))
        names = dimension_names(h)
        p = polyhedron(h, CDDLib.Library())
        removevredundancy!(p)
        P = zeros(size(r.A, 2), size(r.A, 1) + 1)
        for i in eachindex(y)
            j = findfirst(isequal("y[$i]"), names)
            P[:, j] = r.A[i, :]
        end
        cuts = -P * vrep(p)

        for a in points(cuts)
            β = after - a ⋅ x
            cut = AffineFunction(a, β)
            # without some tolerance, CDDLib often throws `Numerically inconsistent`.
            if algo.log_level >= 2
                @info("Cut added: $cut, $after > $before.")
            end
            push!(Q.cuts[left + 1, mode], cut)
        end
        for r in rays(cuts)
            a = normalize(coord(r))
            cut = HalfSpace(a, a ⋅ x)
            if algo.log_level >= 2
                @info("Cut added: $cut, $after > $before.")
            end
            intersect!(Q.domains[left + 1, mode], cut)
        end
    end
end
