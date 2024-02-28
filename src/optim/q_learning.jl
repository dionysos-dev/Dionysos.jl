using LinearAlgebra

import MutableArithmetics
const MA = MutableArithmetics

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem

using HybridSystems
using Polyhedra
using JuMP

export DiscreteLowerBoundAlgo, HybridDualDynamicProgrammingAlgo

# TODO It assumes that the cost is `Fill{...}`, otherwise, we should compute a different
#      cost for each time
function minimum_transition_cost(
    prob,
    transition,
    solver,
    ::Type{T},
    log_level = 0,
) where {T}
    model = MOI.instantiate(solver; with_bridge_type = T)
    from = source(prob.system, transition)
    to = target(prob.system, transition)
    MOI.Bridges.add_bridge(model, Polyhedra.PolyhedraToLPBridge{T})
    x0, c0 = MOI.add_constrained_variables(
        model,
        Polyhedra.PolyhedraOptSet(hrep(stateset(prob.system, from))),
    )
    x1, c1 = MOI.add_constrained_variables(
        model,
        Polyhedra.PolyhedraOptSet(hrep(stateset(prob.system, to))),
    )
    u = MOI.add_variables(model, inputdim(resetmap(prob.system, transition)))
    algo = MOI.instantiate(
        optimizer_with_attributes(
            BemporadMorari.Optimizer{T},
            "continuous_solver" => solver,
            "log_level" => 0,
        ),
    )

    # We use `1` as we asssume the cost is the same along time
    t = 1
    δ_mode = BemporadMorari.IndicatorVariables([to], t)
    state_cost, δ_mode = BemporadMorari.hybrid_cost(
        model,
        BemporadMorari.fillify(prob.state_cost[t][[to]]),
        x1,
        u,
        δ_mode,
        T,
    )
    symbols = [symbol(prob.system, transition)]
    δ_trans = BemporadMorari.IndicatorVariables(symbols, t)
    δ_trans = BemporadMorari.hybrid_constraints(
        model,
        BemporadMorari.fillify(prob.system.resetmaps[symbols]),
        x0,
        x1,
        u,
        algo,
        δ_trans,
    )
    trans_cost, δ_trans = BemporadMorari.hybrid_cost(
        model,
        BemporadMorari.fillify(prob.transition_cost[t][symbols]),
        x0,
        u,
        δ_trans,
        T,
    )
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    obj = state_cost + trans_cost
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj)}(), obj)
    MOI.optimize!(model)
    status = MOI.get(model, MOI.TerminationStatus())
    if status == MOI.OPTIMAL
        return MOI.get(model, MOI.ObjectiveValue())
    elseif status in (MOI.INFEASIBLE,)#, MOI.INFEASIBLE_OR_UNBOUNDED)
        return typemax(T) # ∞
    else
        @error(
            "Candidate: Termination status: $status, raw status: $(MOI.get(model, MOI.RawStatusString()))"
        )
    end
end

function _merge_with(combine, a, b)
    n = (axes(a, 1).stop, axes(b, 1).stop)
    time = maximum(n)
    modes = axes(a, 2)
    return JuMP.Containers.@container(
        [i in 0:time, mode in modes],
        if i > n[1]
            b[i, mode]
        elseif i > n[2]
            a[i, mode]
        else
            combine(a[i, mode], b[i, mode])
        end
    )
end

struct DiscreteLowerBoundAlgo{T, S}
    solver::S
end
DiscreteLowerBoundAlgo{T}(solver) where {T} =
    DiscreteLowerBoundAlgo{T, typeof(solver)}(solver)
struct DiscreteLowerBound{D}
    discrete_lb::D
end
function instantiate(
    prob::PR.OptimalControlProblem,
    algo::DiscreteLowerBoundAlgo{T},
) where {T}
    syst = prob.system
    dists = JuMP.Containers.@container([0:(prob.time), modes(syst)], typemax(T))
    dists[0, prob.target_set] = 0.0
    transition_cost = HybridSystems.transition_property(syst, T)
    for t in transitions(syst)
        transition_cost[t] = minimum_transition_cost(prob, t, algo.solver, T)
    end
    for i in 1:(prob.time)
        for mode in modes(syst)
            for t in out_transitions(syst, mode)
                dists[i, mode] =
                    min(dists[i, mode], transition_cost[t] + dists[i - 1, target(syst, t)])
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
function learn(
    ::DiscreteLowerBound,
    prob,
    ::ST.DiscreteTrajectory,
    ::ST.ContinuousTrajectory,
    ::DiscreteLowerBoundAlgo,
) end

struct HybridDualDynamicProgrammingAlgo{T, S, P}
    solver::S
    polyhedra_library::P
    new_cut_tol::T
    tight_tol::T
    log_level::Int
end
struct HybridDualDynamicProgramming{
    T,
    C,
    H <: JuMP.Containers.DenseAxisArray{Polyhedra.Intersection{T, Vector{T}, Int}},
    D,
}
    cuts::C
    domains::H
    discrete::DiscreteLowerBound{D}
end
function _no_cuts(time, modes, ::Type{T}) where {T}
    return JuMP.Containers.@container([0:time, modes], UT.AffineFunction{T}[])
end
function _full_domains(time, modes, d, ::Type{T}) where {T}
    return JuMP.Containers.@container(
        [0:time, modes],
        hrep(HalfSpace{T, Vector{T}}[], d = d)
    )
end
function instantiate(
    prob::PR.OptimalControlProblem,
    algo::HybridDualDynamicProgrammingAlgo{T},
) where {T}
    cuts = _no_cuts(prob.time, modes(prob.system), T)
    domains = _full_domains(prob.time, modes(prob.system), length(prob.initial_set[2]), T)
    discrete = instantiate(prob, DiscreteLowerBoundAlgo{T}(algo.solver))
    return HybridDualDynamicProgramming(cuts, domains, discrete)
end
function q_merge(a::DiscreteLowerBound, b::HybridDualDynamicProgramming{T}) where {T}
    d = a.discrete_lb
    cuts = _no_cuts(axes(d, 1).stop, axes(d, 2), T)
    domains = _full_domains(axes(d, 1).stop, axes(d, 2), fulldim(first(b.domains)), T)
    return q_merge(HybridDualDynamicProgramming(cuts, domains, a), b)
end
function q_merge(a::HybridDualDynamicProgramming, b::HybridDualDynamicProgramming)
    return HybridDualDynamicProgramming(
        _merge_with(vcat, a.cuts, b.cuts),
        _merge_with(intersect, a.domains, b.domains),
        q_merge(a.discrete, b.discrete),
    )
end
function value_function(Q::HybridDualDynamicProgramming{T}, left::Int, mode) where {T}
    d = fulldim(Q.domains[left, mode])
    time = axes(Q.cuts, 1).stop
    return UT.PolyhedralFunction(
        Q.discrete.discrete_lb[left, mode],
        reduce(
            append!,
            (Q.cuts[i, mode] for i in left:time);
            init = UT.AffineFunction{T}[],
        ),
        # TODO use intersect! once it is implemented in Polyhedra
        reduce(
            intersect,
            (Q.domains[i, mode] for i in left:time);
            init = hrep(HalfSpace{T, Vector{T}}[]; d = d),
        ),
    )
end
function vertices(f::UT.PolyhedralFunction{T}, X, lib) where {T}
    h = (hrep(X) ∩ f.domain) * intersect(HalfSpace([-one(T)], -f.lower_bound))
    cuts = [HalfSpace([p.a; -one(T)], -p.β) for p in f.pieces]
    h_cut = h ∩ hrep(cuts; d = fulldim(h))
    p = polyhedron(h_cut, lib)
    removehredundancy!(p)
    return collect(points(vrep(p)))
end
function vertices(Q::HybridDualDynamicProgramming, prob, left, mode, lib)
    return vertices(value_function(Q, left, mode), stateset(prob.system, mode), lib)
end
function learn(
    ::HybridDualDynamicProgramming,
    prob,
    ::ST.DiscreteTrajectory,
    ::ST.ContinuousTrajectory,
    ::DiscreteLowerBoundAlgo,
) end
function learn(
    Q::HybridDualDynamicProgramming,
    prob,
    dtraj::ST.DiscreteTrajectory,
    ctraj::ST.ContinuousTrajectory,
    algo::HybridDualDynamicProgrammingAlgo{T},
) where {T}
    for i in length(dtraj):-1:1
        x = i == 1 ? prob.initial_set[2] : ctraj.x[i - 1]
        mode = source(prob.system, dtraj.transitions[i])
        left = prob.time - i
        trans = filter(collect(out_transitions(prob.system, mode))) do t
            return Q.discrete.discrete_lb[left, target(prob.system, t)] != typemax(T) # ∞
        end
        verts = JuMP.Containers.@container(
            [t in trans],
            vertices(Q, prob, left, target(prob.system, t), algo.polyhedra_library)
        )
        model = MOI.instantiate(algo.solver; with_bridge_type = T)
        params = x
        # FIXME assume they are all the same
        r = resetmap(prob.system, first(trans))
        U = inputset(r)
        hashyperplanes(U) && error("TODO: Q-learning with input set with hyperplanes")
        u = MOI.add_variables(model, fulldim(U))
        U_h = collect(halfspaces(U))
        u_con = [MOI.add_constraint(model, u ⋅ hs.a, MOI.LessThan(hs.β)) for hs in U_h]
        λ = Dict(
            (t, v) => MOI.add_constrained_variable(model, MOI.GreaterThan(zero(T)))[1]
            for t in trans for v in eachindex(verts[t])
        )
        #@variable(model, λ[t in trans, eachindex(verts[t])] ≥ 0)
        λs = collect(values(λ))
        MOI.add_constraint(model, BemporadMorari._sum(λs, T), MOI.EqualTo(one(T)))
        epi(i) = sum(
            λ[(t, v)] * convert(T, verts[t][v][i]) for t in trans for
            v in eachindex(verts[t])
        )
        x_next = r.A * params + r.B * u
        epi_con = [
            MOI.Utilities.normalize_and_add_constraint(
                model,
                x_next[j] - epi(j),
                MOI.EqualTo(zero(T)),
            ) for j in eachindex(x_next)
        ]
        θ = MOI.add_variable(model)
        MOI.add_constraint(model, θ - epi(length(x_next) + 1), MOI.GreaterThan(zero(T)))
        tos = [target(prob.system, t) for t in trans]
        δ_mode = JuMP.Containers.@container(
            [mode in tos],
            BemporadMorari._sum(
                (
                    λ[t, v] for t in trans for
                    v in eachindex(verts[t]) if target(prob.system, t) == mode
                ),
                T,
            )
        )
        state_cost, δ_mode = BemporadMorari.hybrid_cost(
            model,
            BemporadMorari.fillify(prob.state_cost[i][tos]),
            x_next,
            u,
            δ_mode,
            T,
        )
        symbols = symbol.(prob.system, trans)
        δ_trans = JuMP.Containers.@container(
            [s in symbols],
            BemporadMorari._sum(
                (
                    λ[t, v] for t in trans for
                    v in eachindex(verts[t]) if symbol(prob.system, t) == s
                ),
                T,
            )
        )
        trans_cost, δ_trans = BemporadMorari.hybrid_cost(
            model,
            BemporadMorari.fillify(prob.transition_cost[i][symbols]),
            x,
            u,
            δ_trans,
            T,
        )
        MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        obj = trans_cost + state_cost + θ
        MOI.set(model, MOI.ObjectiveFunction{typeof(obj)}(), obj)
        MOI.optimize!(model)
        if MOI.get(model, MOI.TerminationStatus()) != MOI.OPTIMAL ||
           MOI.get(model, MOI.PrimalStatus()) != MOI.FEASIBLE_POINT
            #if algo.log_level >= 1
            @warn(
                "No new cut generated as termination status is $(MOI.get(model, MOI.TerminationStatus())), primal status is $(MOI.get(model, MOI.PrimalStatus())), dual status is $(MOI.get(model, MOI.DualStatus())): $(MOI.get(model, MOI.RawStatusString()))"
            )
            #end
            return
        end
        V = value_function(Q, left + 1, mode)
        before = UT.function_value(V, x)
        # FIXME OSQP does not support accessing `DualObjectiveValue`
        #after = dual_objective_value(model)
        after = MOI.get(model, MOI.ObjectiveValue())
        if after < before + algo.new_cut_tol
            if algo.log_level >= 3
                @info("Cuts ignored: $after ≤ $before (difference is $(before - after)).")
            end
            return
        end

        dual_model = Model()
        @variable(dual_model, y_sum)
        @variable(dual_model, y[1:length(epi_con)])
        λ_cons = JuMP.Containers.@container(
            [t in trans, v in eachindex(verts[t])],
            y_sum - sum(verts[t][v][i] * y[i] for i in eachindex(y)) -
            verts[t][v][length(x_next) + 1]
        )
        u_value = MOI.get.(model, MOI.VariablePrimal(), u)
        u_cons = r.B' * y
        for term in obj.quadratic_terms
            if term.variable_1 == term.variable_2
                u_idx = findfirst(isequal(term.variable_1), u)
                @assert u_idx !== nothing
                # d(coef * u^2)/du = 2coef * u. In MOI, the quadratic term is stored as
                # `MOI.ScalarQuadraticTerm(2coef, u, u)` so we don't have to multiply by 2.
                MA.operate!(MA.add_mul, u_cons[u_idx], -term.coefficient * u_value[u_idx])
            else
                ua_idx = findfirst(isequal(term.variable_1), u)
                ub_idx = findfirst(isequal(term.variable_2), u)
                if ua_idx !== nothing && ub_idx !== nothing
                    MA.operate!(
                        MA.add_mul,
                        u_cons[ua_idx],
                        -term.coefficient * u_value[ub_idx],
                    )
                    MA.operate!(
                        MA.add_mul,
                        u_cons[ub_idx],
                        -term.coefficient * u_value[ua_idx],
                    )
                else
                    error("TODO")
                end
            end
        end
        for term in obj.affine_terms
            term.variable == θ && continue
            u_idx = findfirst(isequal(term.variable), u)
            if u_idx !== nothing
                MA.operate!(MA.add_mul, u_cons[u_idx], -term.coefficient)
            else
                found = false
                for t in trans, v in eachindex(verts[t])
                    if term.variable == λ[t, v]
                        @assert !found
                        found = true
                        MA.operate!(MA.add_mul, λ_cons[t, v], -term.coefficient)
                    end
                end
                @assert found
            end
        end

        for t in trans, v in eachindex(verts[t])
            if MOI.get(model, MOI.VariablePrimal(), λ[t, v]) <= algo.tight_tol
                @constraint(dual_model, λ_cons[t, v] <= 0)
            else
                @constraint(dual_model, λ_cons[t, v] == 0)
            end
        end

        pv = [zeros(length(u))]
        lv = Line{T, Vector{T}}[]
        rv = Ray{T, Vector{T}}[]
        U_v = Polyhedra.Hull(Polyhedra.FullDim_rec(pv, lv, rv), pv, lv, rv)

        for i in eachindex(U_h)
            α = u_value ⋅ U_h[i].a
            if α >= U_h[i].β - (max(abs(α), abs(U_h[i].β)) + 1) * algo.tight_tol
                # TODO Use `convexhull!` once it is implemented in Polyhedra
                convexhull!(U_v, Ray(U_h.a))
            end
        end
        U_p = polyhedron(U_v, algo.polyhedra_library)
        removehredundancy!(U_p)
        U_h = hrep(U_p)
        @constraint(dual_model, u_cons in U_h)

        h = hrep(dual_model)
        names = dimension_names(h)
        p = polyhedron(h, algo.polyhedra_library)
        removevredundancy!(p)
        P = zeros(size(r.A, 2), size(r.A, 1) + 1)
        for i in eachindex(y)
            j = findfirst(isequal("y[$i]"), names)
            P[:, j] = r.A[i, :]
        end
        #cuts = -P * removevredundancy(vrep(p), Polyhedra.default_solver(p))
        cuts = -P * vrep(p)

        for a in points(cuts)
            β = after - a ⋅ x
            cut = UT.AffineFunction{T}(a, β)
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
