module Opt

import MathOptInterface as MOI
import StaticArrays as SA
import MathematicalSystems
import Dionysos

@enum(VariableType, INPUT, STATE, MODE)

mutable struct Optimizer <: MOI.AbstractOptimizer
    inner
    nlp_model
    evaluator::Union{Nothing,MOI.Nonlinear.Evaluator}
    variable_values::Vector{Float64}
    derivative_values::Vector{Float64}
    variable_types::Vector{VariableType}
    variable_index::Vector{Int} # MOI.VariableIndex -> state/input/mode index
    lower::Vector{Float64}
    upper::Vector{Float64}
    start::Vector{Float64}
    target::Vector{MOI.Interval{Float64}}
    dynamic::Vector{Union{Nothing,MOI.ScalarNonlinearFunction}}
    obstacles::Vector{Tuple{Vector{MOI.VariableIndex},MOI.HyperRectangle}}
    objective_sense::MOI.OptimizationSense
    objective_function::MOI.AbstractScalarFunction
    nonlinear_index::Int
    function Optimizer()
        return new(
            MOI.instantiate(Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer),
            nothing,
            nothing,
            Float64[],
            Float64[],
            VariableType[],
            Int[],
            Float64[],
            Float64[],
            Float64[],
            MOI.Interval{Float64}[],
            Union{Nothing,MOI.ScalarNonlinearFunction}[],
            Tuple{Vector{MOI.VariableIndex},MOI.HyperRectangle}[],
            MOI.FEASIBILITY_SENSE,
            zero(MOI.ScalarAffineFunction{Float64}),
            0,
        )
    end
end

MOI.is_empty(model::Optimizer) = isempty(model.variable_types)

function MOI.empty!(model::Optimizer)
    empty!(model.variable_values)
    empty!(model.derivative_values)
    empty!(model.variable_types)
    empty!(model.variable_index)
    empty!(model.lower)
    empty!(model.upper)
    empty!(model.start)
    empty!(model.target)
    empty!(model.dynamic)
    empty!(model.obstacles)
    model.objective_sense = MOI.FEASIBILITY_SENSE
    model.objective_function = zero(MOI.ScalarAffineFunction{Float64})
    model.nonlinear_index = 0
    return
end

function MOI.add_variable(model::Optimizer)
    push!(model.variable_values, NaN)
    push!(model.variable_types, INPUT)
    push!(model.variable_index, 0)
    push!(model.lower, -Inf)
    push!(model.upper, Inf)
    push!(model.start, NaN)
    push!(model.target, MOI.Interval(-Inf, Inf))
    push!(model.dynamic, nothing)
    return MOI.VariableIndex(length(model.upper))
end

struct OuterSet{S<:MOI.AbstractVectorSet} <: MOI.AbstractVectorSet
    inner::S
end

MOI.dimension(set::OuterSet) = MOI.dimension(set.inner)
Base.copy(set::OuterSet) = OuterSet(copy(set.inner))

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VectorOfVariables},
    ::Type{OuterSet{MOI.HyperRectangle{Float64}}}
)
    return true
end

function MOI.add_constraint(
    model::Optimizer,
    func::MOI.VectorOfVariables,
    set::OuterSet{MOI.HyperRectangle{Float64}},
)
    push!(model.obstacles, (func.variables, set.inner))
    return MOI.ConstraintIndex{typeof(func),typeof(set)}(length(model.obstacles))
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VariableIndex},
    ::Type{<:Union{MOI.LessThan,MOI.GreaterThan}}
)
    return true
end

function MOI.add_constraint(
    model::Optimizer,
    func::MOI.VariableIndex,
    set::MOI.GreaterThan
)
    model.lower[func.value] = MOI.constant(set)
    return MOI.ConstraintIndex{typeof(func),typeof(set)}(func.value)
end

function MOI.add_constraint(
    model::Optimizer,
    func::MOI.VariableIndex,
    set::MOI.LessThan,
)
    model.upper[func.value] = MOI.constant(set)
    return MOI.ConstraintIndex{typeof(func),typeof(set)}(func.value)
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.ScalarNonlinearFunction},
    ::Type{<:MOI.Interval},
)
    return true
end

function MOI.add_constraint(
    model::Optimizer,
    func::MOI.ScalarNonlinearFunction,
    set::MOI.Interval,
)
    if length(func.args) == 1
        var = func.args[]
        if var isa MOI.VariableIndex && func.head == :final
            model.target[var.value] = set
            model.nonlinear_index += 1
            return MOI.ConstraintIndex{typeof(func),typeof(set)}(model.nonlinear_index)
        end
    end
    dump(func)
    error("Unsupported")
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.ScalarNonlinearFunction},
    ::Type{<:MOI.EqualTo},
)
    return true
end

function MOI.add_constraint(
    model::Optimizer,
    func::MOI.ScalarNonlinearFunction,
    set::MOI.EqualTo,
)
    if func.head == :- && length(func.args) == 2
        lhs, rhs = func.args
        if lhs isa MOI.ScalarNonlinearFunction && length(lhs.args) == 1
            var = lhs.args[]
            if var isa MOI.VariableIndex
                if lhs.head == :diff
                    model.dynamic[var.value] = rhs
                    model.nonlinear_index += 1
                    return MOI.ConstraintIndex{typeof(func),typeof(set)}(model.nonlinear_index)
                elseif lhs.head == :final
                    model.target[var.value] = MOI.Interval(rhs, rhs)
                    model.nonlinear_index += 1
                    return MOI.ConstraintIndex{typeof(func),typeof(set)}(model.nonlinear_index)
                end
            end
        end
    end
    dump(func)
    error("Unsupported")
end

function MOI.supports(
    ::Optimizer,
    ::Union{MOI.ObjectiveSense,MOI.ObjectiveFunction}
)
    return true
end

function MOI.supports(
    ::Optimizer,
    ::MOI.VariablePrimalStart,
    ::Type{MOI.VariableIndex},
)
    return true
end

function MOI.set(
    model::Optimizer,
    ::MOI.VariablePrimalStart,
    vi::MOI.VariableIndex,
    value,
)
    model.start[vi.value] = value
end

function MOI.set(
    model::Optimizer,
    ::MOI.ObjectiveSense,
    sense::MOI.OptimizationSense,
)
    model.objective_sense = sense
    return
end

function MOI.set(
    model::Optimizer,
    ::MOI.ObjectiveFunction,
    func::MOI.AbstractScalarFunction,
)
    model.objective_function = func
    return
end

MOI.supports_incremental_interface(::Optimizer) = true

function MOI.copy_to(model::Optimizer, src::MOI.ModelLike)
    return MOI.Utilities.default_copy_to(model, src)
end

_svec(vec, idx) = SA.SVector([vec[i] for i in idx]...)

function _full(model, def, vars, vals)
    v = fill(def, length(model.variable_index))
    for (var, val) in zip(vars, vals)
        v[var.value] = val
    end
    return v
end

function obstacles(model, x_idx)
    return [
        Dionysos.Utils.HyperRectangle(
            _svec(_full(model, -Inf, o[1], o[2].lower), x_idx),
            _svec(_full(model, Inf, o[1], o[2].upper), x_idx),
        )
        for o in model.obstacles
    ]
end

function dynamic!(model)
    MOI.eval_constraint(model.evaluator, model.derivative_values, model.variable_values)
    return
end

function dynamicofsystem(model)
    # System eq x' = F_sys(x, u)
    function F_sys(x, u)
        for i in eachindex(model.variable_values)
            if variable_type(model, MOI.VariableIndex(i)) == STATE
                model.variable_values[i] = x[model.variable_index[i]]
            else
                model.variable_values[i] = u[model.variable_index[i]]
            end
        end
        dynamic!(model)
        return model.derivative_values
    end

    # We define the growth bound function of $f$:
    ngrowthbound = 5
    function L_growthbound(u)
        β = abs(u[1] / cos(atan(tan(u[2]) / 2)))
        return SA.SMatrix{3, 3}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, β, β, 0.0)
    end

    return F_sys, L_growthbound, ngrowthbound
end

function system(
    model,
    _X_,
    _U_;
    sysnoise = SA.SVector(0.0, 0.0, 0.0),
    measnoise = SA.SVector(0.0, 0.0, 0.0),
    tstep = 0.3,
    nsys = 5,
)
    F_sys, L_growthbound, ngrowthbound = dynamicofsystem(model)
    contsys = Dionysos.System.NewControlSystemGrowthRK4(
        tstep,
        F_sys,
        L_growthbound,
        sysnoise,
        measnoise,
        nsys,
        ngrowthbound,
    )
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        contsys,
        Dionysos.Utils.get_dims(_X_),
        Dionysos.Utils.get_dims(_U_),
        _X_,
        _U_,
    )
end

function problem(model::Optimizer)
    x_idx = state_indices(model)
    u_idx = input_indices(model)
    _X_ = Dionysos.Utils.HyperRectangle(
        _svec(model.lower, x_idx),
        _svec(model.upper, x_idx),
    )
    _U_ = Dionysos.Utils.HyperRectangle(
        _svec(model.lower, u_idx),
        _svec(model.upper, u_idx),
    )
    _I_ = Dionysos.Utils.HyperRectangle(
        _svec(model.start, x_idx),
        _svec(model.start, x_idx),
    )
    _T_ = Dionysos.Utils.HyperRectangle(
        _svec([s.lower for s in model.target], x_idx),
        _svec([s.upper for s in model.target], x_idx),
    )

    _X_ = Dionysos.Utils.LazySetMinus(
        _X_,
        Dionysos.Utils.LazyUnionSetArray(obstacles(model, x_idx)),
    )

    sys = system(model, _X_, _U_)
    problem = Dionysos.Problem.OptimalControlProblem(
        sys, _I_, _T_, nothing, nothing, Dionysos.Problem.Infinity(),
    )
    return problem
end

function variable_type(model::Optimizer, vi::MOI.VariableIndex)
    if isnothing(model.dynamic[vi.value])
        return INPUT
    else
        return STATE
    end
end

function state_indices(model::Optimizer)
    return findall(eachindex(model.variable_index)) do i
        return variable_type(model, MOI.VariableIndex(i)) == STATE
    end
end

function input_indices(model::Optimizer)
    return findall(eachindex(model.variable_index)) do i
        return variable_type(model, MOI.VariableIndex(i)) == INPUT
    end
end

function setup!(model::Optimizer)
    input_index = 0
    state_index = 0
    model.nlp_model = MOI.Nonlinear.Model()
    for i in eachindex(model.variable_types)
        if variable_type(model, MOI.VariableIndex(i)) == STATE
            # The set does not matter
            MOI.Nonlinear.add_constraint(model.nlp_model, model.dynamic[i], MOI.EqualTo(0.0))
            state_index += 1
            model.variable_index[i] = state_index
        else
            input_index += 1
            model.variable_index[i] = input_index
        end
    end
    model.derivative_values = fill(NaN, state_index)
    backend = MOI.Nonlinear.SparseReverseMode()
    vars = MOI.VariableIndex.(eachindex(model.lower))
    model.evaluator = MOI.Nonlinear.Evaluator(model.nlp_model, backend, vars)
    MOI.initialize(model.evaluator, Symbol[])
    return
end

function MOI.optimize!(model::Optimizer)
    setup!(model)
    MOI.set(model.inner, MOI.RawOptimizerAttribute("concrete_problem"), problem(model))
    MOI.optimize!(model.inner)
    return
end

MOI.supports(::Optimizer, ::MOI.RawOptimizerAttribute) = true

function MOI.get(model::Optimizer, attr::MOI.RawOptimizerAttribute)
    return MOI.get(model.inner, attr)
end

function MOI.set(model::Optimizer, attr::MOI.RawOptimizerAttribute, value)
    return MOI.set(model.inner, attr, value)
end

end
