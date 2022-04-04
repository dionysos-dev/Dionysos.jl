export BemporadMorari

module BemporadMorari

using Dionysos

import MutableArithmetics
const MA = MutableArithmetics

using FillArrays, MathematicalSystems, HybridSystems, JuMP, SemialgebraicSets, Polyhedra

@enum DiscretePresolveStatus OPTIMIZE_NOT_CALLED TRIVIAL FEASIBLE NO_MODE NO_TRANSITION

mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    continuous_solver
    mixed_integer_solver
    indicator::Bool
    log_level::Int
    modes::Union{Nothing, Vector{Vector{Int}}}
    problem::Union{Nothing, OptimalControlProblem}
    discrete_presolve_status::DiscretePresolveStatus
    inner
    x
    u
    δ_modes
    δ_transs
    function Optimizer{T}() where {T}
        return new{T}(
            nothing, nothing, true, 1, nothing, nothing, OPTIMIZE_NOT_CALLED,
            nothing, nothing, nothing, nothing, nothing
        )
    end
end

MOI.is_empty(optimizer::Optimizer) = optimizer.problem === nothing

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    setproperty!(model, Symbol(param.name), value)
end
function MOI.get(model::Optimizer, param::MOI.RawOptimizerAttribute)
    getproperty(model, Symbol(param.name))
end

function default_modes(system, q_T, N)
    return [
        t == N ? [q_T] : HybridSystems.modes(system)
        for t in 1:N
    ]
end

all_same(f, it) = all(el -> f(el) == f(first(it)), it)

struct IndicatorVariables{T, VT<:AbstractVector{T}}
    choices::VT
    time::Int
end
_name(::AbstractVector{<:Integer}) = "mode"
_name(::AbstractVector{<:HybridSystems.AbstractTransition}) = "trans"
_base_name(iv::IndicatorVariables) = "δ_$(_name(iv.choices))_$(iv.time)"

_Scalar = Union{JuMP.AbstractJuMPScalar, MOI.AbstractScalarFunction}

indicator_variables(model, δ::AbstractVector{<:_Scalar}, ::Type) = δ
_sum(g, T) = reduce(+, g, init = zero(MOI.ScalarAffineFunction{T}))
function indicator_variables(model, iv::IndicatorVariables, T::Type)
    δ = JuMP.Containers.@container([iv.choices], add_variable(model))
    add_constraint.(model, δ, MOI.ZeroOne())
    add_constraint(model, _sum(δ, T), MOI.EqualTo(one(T)))
    #@constraint(model, δ ∈ SOS1())
    return δ
end

function possible_transitions(system, froms, tos)
    set = Set{transitiontype(system)}()
    for from in froms
        for to in tos
            trs = transitions(system, from, to)
            union!(set, trs)
        end
    end
    return collect(set)
end

fillify(vector::Fill) = vector
function fillify(vector::AbstractVector)
    if length(vector) == 1 || all_same(identity, vector)
        return Fill(first(vector), length(vector))
    else
        return vector
    end
end
inner_vector(f, vector) = fillify(map(f, vector))

_value(model::MOI.ModelLike, v::MOI.VariableIndex) = MOI.get(model, MOI.VariablePrimal(), v)
_value(::JuMP.Model, v::JuMP.VariableRef) = value(v)
add_variable(model::MOI.ModelLike) = MOI.add_variable(model)
add_variable(model::JuMP.Model) = @variable(model)
add_constraint(model::MOI.ModelLike, func::AbstractVector, set) = MOI.add_constraint(model, MOI.Utilities.vectorize(func), set)
add_constraint(model::MOI.ModelLike, func, set) = MOI.add_constraint(model, func, set)
function add_constraint(model::MOI.ModelLike, func::MOI.AbstractScalarFunction, set::MOI.AbstractScalarSet)
    MOI.Utilities.normalize_and_add_constraint(model, func, set; allow_modify_function=true)
end
add_constraint(model::JuMP.Model, func, set) = JuMP.add_constraint(model, JuMP.build_constraint(error, func, set))
indicator_constraint(model::JuMP.Model, δ, func, set) = @constraint(model, δ => {func in set})

function hybrid_constraints(model, sets::Fill{<:Polyhedra.Rep}, x, algo::Optimizer, δ)
    set = first(sets)
    add_constraint(model, x, Polyhedra.PolyhedraOptSet(hrep(set)))
    return δ
end

using LinearAlgebra
function hybrid_constraints(model, sets::AbstractVector{<:Polyhedra.Rep}, x, algo::Optimizer{T}, δs) where {T}
    δs = indicator_variables(model, δs, T)
    for (δ, set) in zip(δs, sets)
        # TODO implement indicator with Polyhedra in Polyhedra.jl
        for h in halfspaces(hrep(set))
            if algo.indicator
                indicator_constraint(model, δ, x ⋅ h.a, MOI.LessThan(h.β))
            else
                M = maximum([support_function(h.a, set) for set in sets])
                add_constraint(model, x ⋅ h.a + (M - h.β) * δ, MOI.LessThan(M))
            end
        end
        for h in hyperplanes(hrep(set))
            indicator_constraint(model, δ, x ⋅ h.a, MOI.EqualTo(h.β))
        end
    end
    return δs
end

hybrid_constraints(model, ::AbstractVector{FullSpace}, x, algo::Optimizer, δ) = δ

function hybrid_constraints(model, systems::Fill{<:LinearControlMap}, x_prev, x, u, algo::Optimizer{T}, δ) where {T}
    system = first(systems)
    add_constraint.(model, x - system.A * x_prev - system.B * u, MOI.EqualTo(zero(T)))
    return δ
end

function hybrid_constraints(model, systems::AbstractVector{<:ConstrainedContinuousIdentitySystem}, x_prev, x, u, algo::Optimizer, δ)
    return hybrid_constraints(model, inner_vector(system -> system.X, systems), x, algo, δ)
end

function hybrid_constraints(model, systems::AbstractVector{<:ConstrainedLinearControlMap}, x_prev, x, u, algo::Optimizer, δ)
    δ = hybrid_constraints(model, inner_vector(system -> LinearControlMap(system.A, system.B), systems), x_prev, x, u, algo, δ)
    δ = hybrid_constraints(model, inner_vector(system -> system.X, systems), x, algo, δ)
    return hybrid_constraints(model, inner_vector(system -> system.U, systems), u, algo, δ)
end

# TODO move to MA
#Base.zero(::Type{MA.Zero}) = MA.Zero()

function hybrid_cost(model, costs::AbstractArray{ZeroFunction}, x, u, δ, ::Type{T}) where {T}
    return zero(T), δ
end
function hybrid_cost(model, costs::AbstractArray{<:ConstantFunction}, x, u, δ, ::Type{T}) where {T}
    δ = indicator_variables(model, δ, T)
    cost = [cost.value for cost in costs]
    return cost ⋅ δ, δ
end
function hybrid_cost(model, costs::Fill{<:ConstantFunction}, x, u, δ, ::Type{T}) where {T}
    return first(costs).value, δ
end
function hybrid_cost(model, costs::Fill{<:QuadraticControlFunction}, x, u, δ, ::Type{T}) where {T}
    cost = first(costs)
    return u' * cost.Q * u, δ
end
function hybrid_cost(model, costs::Fill{<:PolyhedralFunction}, x, u, δ, ::Type{T}) where {T}
    cost = first(costs)
    θ = add_variable(model)
    add_constraint(model, θ, MOI.GreaterThan(cost.lower_bound))
    for piece in cost.pieces
        add_constraint(model, θ - function_value(piece, x), MOI.GreaterThan(zero(T)))
    end
    add_constraint(model, x, Polyhedra.PolyhedraOptSet(cost.domain))
    return θ, δ
end

function transitions_constraints(model, system, modes_from, δ_from::IndicatorVariables, modes_to, δ_to::IndicatorVariables, trans, δ_trans::IndicatorVariables, ::Type)
    # Nothing to do, the impossible modes should have already been pruned
end
function transitions_constraints(model, system, modes_from, δ_from::IndicatorVariables, modes_to, δ_to::AbstractVector{<:_Scalar}, trans, δ_trans::IndicatorVariables, ::Type)
    # Nothing to do, the impossible modes should have already been pruned
end
function transitions_constraints(model, system, modes_from, δ_from::AbstractVector{<:_Scalar}, modes_to, δ_to::IndicatorVariables, trans, δ_trans::IndicatorVariables, ::Type)
    # Nothing to do, the impossible modes should have already been pruned
end
function transitions_constraints(model, system, modes_from, δ_from::AbstractVector{<:_Scalar}, modes_to, δ_to::AbstractVector{<:_Scalar}, trans, δ_trans::IndicatorVariables, T::Type)
    for (mode_from, from) in zip(modes_from, δ_from)
        for (mode_to, to) in zip(modes_to, δ_to)
            if !has_transition(system, mode_from, mode_to)
                # Logical constraint `!(δ_from && δ_to)`
                add_constraint(model, one(T) * from + to, MOI.LessThan(one(T)))
            end
        end
    end
end

function _sparse_set(δs, δ, t)
    if !(δ isa IndicatorVariables)
        if δs === nothing
            δs = JuMP.Containers.SparseAxisArray(Dict((t,) => δ))
        else
            δs[t] = δ
        end
    end
    return δs
end

function _zero_steps(optimizer)
    if optimizer.problem.q_T == optimizer.problem.q_0
        optimizer.discrete_presolve_status = TRIVIAL
    else
        optimizer.discrete_presolve_status = NO_MODE
    end
    return
end

function MOI.optimize!(optimizer::Optimizer{T}) where {T}
    prob = optimizer.problem
    if optimizer.modes === nothing
        optimizer.modes = default_modes(prob.system, prob.q_T, prob.number_of_time_steps)
    end
    iszero(prob.number_of_time_steps) && return _zero_steps(optimizer)
    modes = Vector{Vector{Int}}(undef, prob.number_of_time_steps)
    for t in 1:prob.number_of_time_steps
        modes_prev = t == 1 ? [prob.q_0] : modes[t - 1]
        modes[t] = filter(optimizer.modes[t]) do mode
            any(modes_prev) do mode_prev
                has_transition(prob.system, mode_prev, mode)
            end
        end
    end
    for t in (prob.number_of_time_steps - 1):-1:1
        modes[t] = filter(modes[t]) do mode
            any(modes[t + 1]) do mode_next
                has_transition(prob.system, mode, mode_next)
            end
        end
    end
    # TODO remove modes that are impossible
    if any(isempty, modes)
        optimizer.log_level >= 1 && @warn("`modes` is empty for some time step.")
        optimizer.discrete_presolve_status = NO_MODE
        return
    end

    model = MOIU.CachingOptimizer(MOIU.UniversalFallback(MOIU.Model{T}()), MOIU.MANUAL)
    # Use `model` instead of `optimizer.inner`
    # as it is type stable.
    optimizer.inner = model
    JuMP.Containers.@container(x[1:prob.number_of_time_steps, 1:statedim(prob.system, first(first(modes)))], add_variable(model))
    optimizer.x = x
    JuMP.Containers.@container(u[1:prob.number_of_time_steps, 1:inputdim(resetmap(prob.system, first(transitions(prob.system))))], add_variable(model))
    optimizer.u = u

    transs = [possible_transitions(prob.system, t == 1 ? [prob.q_0] : modes[t-1], modes[t])
              for t in 1:prob.number_of_time_steps]
    if any(isempty, transs)
        optimizer.log_level >= 1 && @warn("`transs` is empty for some time step.")
        optimizer.discrete_presolve_status = NO_TRANSITION
        return
    end

    optimizer.discrete_presolve_status = FEASIBLE

    total_cost = MA.Zero()

    δ_mode_prev = IndicatorVariables([prob.q_0], 0)
    δ_modes = nothing
    δ_transs = nothing

    for t in 1:prob.number_of_time_steps
        x_prev = t == 1 ? prob.x_0 : x[t - 1, :]
        xi = x[t, :]
        ui = u[t, :]
        δ_mode = IndicatorVariables(modes[t], t)
        δ_mode = hybrid_constraints(model, fillify(prob.system.modes[modes[t]]), x_prev, xi, ui, optimizer, δ_mode)
        symbols = symbol.(prob.system, transs[t])
        δ_trans = IndicatorVariables(transs[t], t)
        δ_trans = hybrid_constraints(model, fillify(prob.system.resetmaps[symbols]), x_prev, xi, ui, optimizer, δ_trans)
        state_cost, δ_mode = hybrid_cost(model, fillify(prob.state_cost[t][modes[t]]), xi, ui, δ_mode, T)
        total_cost = MA.operate!!(+, total_cost, state_cost)
        trans_cost, δ_trans = hybrid_cost(model, fillify(prob.transition_cost[t][symbols]), x_prev, ui, δ_trans, T)
        total_cost = MA.operate!!(+, total_cost, trans_cost)
        modes_prev = t == 1 ? [prob.q_0] : modes[t - 1]
        transitions_constraints(model, prob.system, modes_prev, δ_mode_prev, modes[t], δ_mode, transs[t], δ_trans, T)
        δ_mode_prev = δ_mode
        δ_modes = _sparse_set(δ_modes, δ_mode, t)
        δ_transs = _sparse_set(δ_transs, δ_trans, t)
    end
    optimizer.δ_modes = δ_modes
    optimizer.δ_transs = δ_transs

    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{typeof(total_cost)}(), total_cost)

    solver = if δ_modes === nothing && δ_transs === nothing
        optimizer.continuous_solver
    else
        optimizer.mixed_integer_solver
    end
    _optimizer = MOI.instantiate(solver, with_bridge_type=T)
    MOI.Bridges.add_bridge(_optimizer, Polyhedra.PolyhedraToLPBridge{T})
    MOIU.reset_optimizer(model, _optimizer)
    MOIU.attach_optimizer(model)

    if optimizer.log_level >= 2
        print(model)
    end
    MOI.optimize!(model)

    if !(MOI.get(model, MOI.TerminationStatus()) in [MOI.OPTIMAL, MOI.INFEASIBLE])
        if optimizer.log_level >= 1
            @warn("BemporadMorari: Termination status: $(MOI.get(model, MOI.TerminationStatus())), raw status: $(MOI.get(model, MOI.RawStatusString()))")
        end
    end
end

_rows(A::Matrix) = [A[i, :] for i in 1:size(A, 1)]
function MOI.get(optimizer::Optimizer, ::ContinuousTrajectoryAttribute)
    #if optimizer.log_level >= 1
    #    if δ_modes !== nothing
    #        @show (x -> value.(x)).(δ_modes)
    #    end
    #    if δ_transs !== nothing
    #        @show (x -> value.(x)).(δ_transs)
    #    end
    #end
    if optimizer.discrete_presolve_status == TRIVIAL
        return ContinuousTrajectory(Vector{Float64}[], Vector{Float64}[])
    else
        return ContinuousTrajectory(_rows(_value.(optimizer.inner, optimizer.x)), _rows(_value.(optimizer.inner, optimizer.u)))
    end
end

function MOI.get(optimizer::Optimizer, attr::Union{MOI.SolveTimeSec, MOI.ObjectiveValue})
    if optimizer.discrete_presolve_status == TRIVIAL
        return 0.0
    else
        return MOI.get(optimizer.inner, attr)
    end
end

function MOI.get(optimizer::Optimizer, attr::Union{MOI.TerminationStatus, MOI.PrimalStatus})
    if optimizer.discrete_presolve_status == OPTIMIZE_NOT_CALLED
        return attr isa MOI.TerminationStatus ? MOI.OPTIMIZE_NOT_CALLED : MOI.NO_SOLUTION
    elseif optimizer.discrete_presolve_status == FEASIBLE
        return MOI.get(optimizer.inner, attr)
    elseif optimizer.discrete_presolve_status == TRIVIAL
        return MOI.OPTIMAL
    else
        return attr isa MOI.TerminationStatus ? MOI.INFEASIBLE : MOI.NO_SOLUTION
    end
end

end
