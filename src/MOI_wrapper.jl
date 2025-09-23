import MathOptInterface as MOI
import StaticArrays: SVector
import MathematicalSystems
import JuMP
import MathOptSymbolicAD
import Symbolics
import OSQP, HiGHS, Ipopt, Pavito, CDDLib
using Polyhedra
using FillArrays

export ∂,
    Δ,
    final,
    start,
    rem,
    TransitionAsModel,
    add_transition!,
    guard,
    resetmaps,
    themode,
    set_algorithm!

@enum(VariableType, INPUT, STATE, MODE)

@enum(TimeType, UNKNOWN, CONTINUOUS, DISCRETE, HYBRID)

DEFAULT_MODE = 1

# Define the class Variable to handle a variable value, type, index, lower, upper, start, and target.
mutable struct Variable
    value::Float64
    type::VariableType
    index::Int
    lower::Float64
    upper::Float64
    start::MOI.Interval{Float64}
    target::MOI.Interval{Float64}

    function Variable()
        return new(
            0.0,
            INPUT,
            0,
            -Inf,
            Inf,
            MOI.Interval(-Inf, Inf),
            MOI.Interval(-Inf, Inf),
        )
    end
end

# Define the class Objective to handle the sense and function.
mutable struct Objective
    sense::MOI.OptimizationSense
    func::MOI.AbstractScalarFunction

    function Objective()
        return new(MOI.FEASIBILITY_SENSE, zero(MOI.ScalarAffineFunction{Float64}))
    end
end

# Define the class Mode to handle nlp_model, evaluator, and dynamic.
# The dynamic is a dictionary with mode_name -> [dynamic_func for each state]).
# The evaluator is a dictionary with mode_name -> evaluator.
# The nlp_model is a dictionary with mode_name -> nlp_model.
# if no mode is specified, the default mode is used.
mutable struct Mode
    name::Int
    nlp_model::Any
    evaluator::Union{Nothing, MOI.Nonlinear.Evaluator}
    dynamic::Vector{Union{Nothing, MOI.ScalarNonlinearFunction}}

    function Mode(name, nvars)
        return new(name, MOI.Nonlinear.Model(), nothing, fill(nothing, nvars))
    end
end

# Define the class Transition to handle source, destination, guard, and resetmap.
struct Transition
    source::Int
    destination::Int
    guard::Vector{Union{Nothing, MOI.ScalarNonlinearFunction}}
    resetmap::Vector{Union{Nothing, MOI.ScalarNonlinearFunction}}

    function Transition(src, dst, nvars)
        return new(src, dst, fill(nothing, nvars), fill(nothing, nvars))
    end
end

# Define the class Optimizer to handle the inner, modes, time_type, variables, derivative_values, obstacles, objective, nonlinear_index, mode_variable, transitions, integer_variables, and indicator_constraints.
mutable struct Optimizer <: MOI.AbstractOptimizer
    inner::Any
    variables::Vector{Variable}
    modes::Dict{Int, Mode}
    time_type::TimeType
    derivative_values::Vector{Float64}
    obstacles::Vector{Tuple{Vector{MOI.VariableIndex}, MOI.HyperRectangle}}
    objective::Objective
    nonlinear_index::Int
    mode_variable::Union{Nothing, MOI.VariableIndex}
    transitions::Dict{String, Transition}
    integer_variables::Dict{MOI.VariableIndex, Union{MOI.Integer, MOI.ZeroOne}}
    indicator_constraints::Vector{
        Tuple{MOI.VariableIndex, Bool, MOI.ScalarNonlinearFunction},
    }

    function Optimizer()
        return new(
            MOI.instantiate(Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer),
            Vector{Variable}(),
            Dict{Int, Mode}(),
            UNKNOWN,
            Vector{Float64}(),
            Vector{Tuple{Vector{MOI.VariableIndex}, MOI.HyperRectangle}}(),
            Objective(),
            0,
            nothing,
            Dict{String, Transition}(),
            Dict{MOI.VariableIndex, Union{MOI.Integer, MOI.ZeroOne}}(),
            Vector{Tuple{MOI.VariableIndex, Bool, MOI.ScalarNonlinearFunction}}(),
        )
    end
end

MOI.is_empty(model::Optimizer) = isempty(model.variables)

function MOI.empty!(model::Optimizer)
    model.time_type = UNKNOWN
    empty!(model.variables)
    empty!(model.modes)
    empty!(model.derivative_values)
    empty!(model.obstacles)
    model.objective = Objective()
    model.nonlinear_index = 0
    model.mode_variable = nothing
    empty!(model.transitions)
    empty!(model.integer_variables)
    empty!(model.indicator_constraints)
    return
end

function MOI.add_variable(model::Optimizer)
    println(">>Adding variable")
    push!(model.variables, Variable())

    return MOI.VariableIndex(length(model.variables))
end

struct OuterSet{S <: MOI.AbstractVectorSet} <: MOI.AbstractVectorSet
    inner::S
end

MOI.dimension(set::OuterSet) = MOI.dimension(set.inner)
Base.copy(set::OuterSet) = OuterSet(copy(set.inner))

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VectorOfVariables},
    ::Type{OuterSet{MOI.HyperRectangle{Float64}}},
)
    return true
end

function MOI.add_constraint(
    model::Optimizer,
    func::MOI.VectorOfVariables,
    set::OuterSet{MOI.HyperRectangle{Float64}},
)
    push!(model.obstacles, (func.variables, set.inner))
    return MOI.ConstraintIndex{typeof(func), typeof(set)}(length(model.obstacles))
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VariableIndex},
    ::Type{<:Union{MOI.LessThan, MOI.GreaterThan}},
)
    return true
end

function MOI.add_constraint(model::Optimizer, func::MOI.VariableIndex, set::MOI.GreaterThan)
    model.variables[func.value].lower = MOI.constant(set)
    return MOI.ConstraintIndex{typeof(func), typeof(set)}(func.value)
end

function MOI.add_constraint(model::Optimizer, func::MOI.VariableIndex, set::MOI.LessThan)
    model.variables[func.value].upper = MOI.constant(set)
    return MOI.ConstraintIndex{typeof(func), typeof(set)}(func.value)
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.ScalarNonlinearFunction},
    ::Type{<:MOI.Interval},
)
    return true
end

function _new_nonlinear_index!(model, func, set)
    model.nonlinear_index += 1
    return MOI.ConstraintIndex{typeof(func), typeof(set)}(model.nonlinear_index)
end

function MOI.add_constraint(
    model::Optimizer,
    func::MOI.ScalarNonlinearFunction,
    set::MOI.Interval,
)
    if _is_guard_setup(func)
        println(">>Adding guard setup")
        _handle_hybrid_systems_table(model, func.args[1], func.args[2], set)
        return _new_nonlinear_index!(model, func, set)
    end

    if length(func.args) == 1
        var = func.args[]
        if var isa MOI.VariableIndex && func.head == :final
            model.variables[var.value].target = set
            return _new_nonlinear_index!(model, func, set)
        end

        if var isa MOI.VariableIndex && func.head == :start
            model.variables[var.value].start = set
            return _new_nonlinear_index!(model, func, set)
        end
    end
    dump(func)
    return error("Unsupported! Happened when adding belonging (in) constraint.")
end

# Adding support for nonlinear constraints with ==, >=, <=
function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.ScalarNonlinearFunction},
    ::Type{<:MOI.EqualTo},
)
    return true
end

function _handle_dynamic_error(head, type)
    if head == :∂ && type == DISCRETE
        error(
            "Cannot add constraint with `∂` since you already added a constraint with `Δ`.",
        )
    end

    if head == :Δ && type == CONTINUOUS
        error(
            "Cannot add constraint with `Δ` since you already added a constraint with `∂`.",
        )
    end
end

function _is_binary_expr(expr, head)
    return expr.head == head && length(expr.args) == 2
end

function _is_expr_var_equal_int(expr)
    return _is_binary_expr(expr, :(==)) &&
           expr.args[1] isa MOI.VariableIndex &&
           expr.args[2] isa Int
end

function _is_mode_setup(expr)
    return _is_binary_expr(expr, :(==)) && _is_expr_var_equal_int(expr.args[1])
end

function _is_dynamic_expr(expr)
    return (expr.head == :∂ || expr.head == :Δ) &&
           length(expr.args) == 1 &&
           expr.args[1] isa MOI.VariableIndex
end

function _is_continuous_dynamic(expr)
    return _is_binary_expr(expr, :-) &&
           expr.args[1].head == :∂ &&
           _is_dynamic_expr(expr.args[1])
end

function _is_discrete_dynamic(expr)
    return _is_binary_expr(expr, :-) &&
           expr.args[1].head == :Δ &&
           _is_dynamic_expr(expr.args[1])
end

function _is_resetmaps_setup(expr)
    return _is_binary_expr(expr, :(==)) &&
           _is_binary_expr(expr.args[1], :resetmaps) &&
           _is_expr_var_equal_int(expr.args[1].args[1]) &&
           _is_expr_var_equal_int(expr.args[1].args[2]) &&
           _is_discrete_dynamic(expr.args[2])
end

function _is_guard_setup(expr)
    return _is_binary_expr(expr, :(==)) &&
           _is_binary_expr(expr.args[1], :guard) &&
           _is_expr_var_equal_int(expr.args[1].args[1]) &&
           _is_expr_var_equal_int(expr.args[1].args[2])
end

function _update_dynamic_dict!(model, type, func, val, mode = DEFAULT_MODE)
    if !haskey(model.modes, mode)
        model.modes[mode] = Mode(mode, length(model.variables))
    end

    model.time_type = type

    return model.modes[mode].dynamic[val] = func
end

function _create_transition_if_not_exists!(model, src, dst)
    if isempty(model.transitions)
        model.transitions = Dict{String, Transition}()
    end

    transition = "$src->$dst"
    if !haskey(model.transitions, transition)
        model.transitions[transition] = Transition(src, dst, length(model.variables))
    end

    return transition
end

function _update_guard_table!(model, src, dst, func, val)
    transition = _create_transition_if_not_exists!(model, src, dst)

    return model.transitions[transition].guard[val] = func
end

function _update_resetmaps_table!(model, src, dst, func, val)
    transition = _create_transition_if_not_exists!(model, src, dst)

    return model.transitions[transition].resetmap[val] = func
end

function MOI.add_constraint(
    model::Optimizer,
    func::MOI.ScalarNonlinearFunction,
    set::MOI.EqualTo,
)
    if _is_mode_setup(func)
        println(">>Adding mode setup")
        _handle_hybrid_systems_table(model, func.args[1], func.args[2])
        return _new_nonlinear_index!(model, func, set)
    elseif _is_resetmaps_setup(func)
        println(">>Adding resetmaps setup")
        _handle_hybrid_systems_table(model, func.args[1], func.args[2])
        return _new_nonlinear_index!(model, func, set)
    elseif _is_guard_setup(func)
        println(">>Adding guard setup")
        _handle_hybrid_systems_table(model, func.args[1], func.args[2])
        return _new_nonlinear_index!(model, func, set)
    end

    if func.head == :- && length(func.args) == 2
        lhs, rhs = func.args
        if lhs isa MOI.ScalarNonlinearFunction && length(lhs.args) == 1
            var = lhs.args[]
            if var isa MOI.VariableIndex
                _handle_dynamic_error(lhs.head, model.time_type)

                if lhs.head == :∂
                    _update_dynamic_dict!(model, CONTINUOUS, rhs, var.value)
                    return _new_nonlinear_index!(model, func, set)
                elseif lhs.head == :Δ
                    _update_dynamic_dict!(model, DISCRETE, rhs, var.value)
                    return _new_nonlinear_index!(model, func, set)
                elseif lhs.head == :final
                    model.target[var.value] = MOI.Interval(rhs, rhs)
                    return _new_nonlinear_index!(model, func, set)
                end
            end
        end
    end

    dump(func)
    return error("Unsupported! Happened when adding equality constraint.")
end

function _handle_hybrid_systems_table(model::Optimizer, lhs, rhs, set = MOI.EqualTo)
    # Set hybrid systems tables
    if !(length(lhs.args) == 2)
        error(
            "Unsupported! If you are defining a mode the left-hand side must be something like mode == val. If you are defining a guard or resetmap, use add_transition.",
        )
    end

    # handling mode
    if lhs.head == :(==)
        lhs_var, lhs_val = lhs.args[1], lhs.args[2]
        if !(lhs_var isa MOI.VariableIndex && lhs_val isa Int)
            error(
                "Unsupported! If you are defining a mode the left-hand side must be something like mode == val.",
            )
        end

        if !(_is_continuous_dynamic(rhs))
            error(
                "Unsupported! If you are defining a mode the right-hand side must be a continuous dynamic.",
            )
        end
        rhs_rhs = rhs.args[2]
        rhs_var = rhs.args[1].args[1]

        if !isnothing(model.mode_variable) && model.mode_variable != lhs_var
            error(
                "The mode variable must be the same for all modes. We cannot have multiple mode variables.",
            )
        end

        model.mode_variable = lhs_var
        model.variables[lhs_var.value].type = MODE

        _update_dynamic_dict!(model, HYBRID, rhs_rhs, rhs_var.value, lhs_val)
    end

    # handling guard
    if lhs.head == :guard
        lhs_var, lhs_src, lhs_dst =
            lhs.args[1].args[1], lhs.args[1].args[2], lhs.args[2].args[2]

        # build a MOI.ScalarNonlinearFunction from rhs.args[1] a MOI.ScalarAffineFunction and set a MOI.Interval
        if rhs.args[1] isa MOI.ScalarAffineFunction && set isa MOI.Interval
            #TODO: tweak this to handle more cases
            var = rhs.args[1].terms[1].variable
            if set.lower == -Inf
                rhs = MOI.ScalarNonlinearFunction(:(>=), [var, set.upper])
            elseif set.upper == Inf
                rhs = MOI.ScalarNonlinearFunction(:(<=), [var, set.lower])
            elseif set.lower == set.upper
                rhs = MOI.ScalarNonlinearFunction(:(==), [var, set.lower])
            else
                error("Unsupported! The set must be an interval.")
            end
        end

        _update_guard_table!(model, lhs_src, lhs_dst, rhs, lhs_var.value)
    end

    # handling resetmap
    if lhs.head == :resetmaps
        lhs_var, lhs_src, lhs_dst =
            lhs.args[1].args[1], lhs.args[1].args[2], lhs.args[2].args[2]

        if !(_is_discrete_dynamic(rhs))
            error(
                "Unsupported! If you are defining a resetmap the right-hand side must be a discrete dynamic.",
            )
        end

        rhs_rhs = rhs.args[2]
        rhs_var = rhs.args[1].args[1]

        _update_resetmaps_table!(model, lhs_src, lhs_dst, rhs_rhs, rhs_var.value)
    end
end

# Adding support for objective sense
function MOI.supports(::Optimizer, ::Union{MOI.ObjectiveSense, MOI.ObjectiveFunction})
    return true
end

function MOI.set(model::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    model.objective.sense = sense
    return
end

function MOI.set(
    model::Optimizer,
    ::MOI.ObjectiveFunction,
    func::MOI.AbstractScalarFunction,
)
    model.objective.func = func
    return
end

# Adding support for primal start
function MOI.supports(::Optimizer, ::MOI.VariablePrimalStart, ::Type{MOI.VariableIndex})
    return true
end

function MOI.set(model::Optimizer, ::MOI.VariablePrimalStart, vi::MOI.VariableIndex, value)
    # create a MOI.Interval from the value if value is a scalar
    if !isa(value, MOI.Interval)
        value = MOI.Interval(value, value)
    end

    return model.variables[vi.value].start = value
end

# Adding supports for integer variables
function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VariableIndex},
    ::Type{<:Union{MOI.Integer, MOI.ZeroOne}},
)
    return true
end

function MOI.add_constraint(
    model::Optimizer,
    func::MOI.VariableIndex,
    set::Union{MOI.Integer, MOI.ZeroOne},
)
    model.integer_variables[func] = set
    return MOI.ConstraintIndex{typeof(func), typeof(set)}(func.value)
end

# Adding support for indicator constraints
function JuMP._build_indicator_constraint(
    error_fn::Function,
    lhs::JuMP.NonlinearExpr,
    constraint::JuMP.ScalarConstraint,
    ::Type{MOI.Indicator{A}},
) where {A}
    # Get the function and set
    fun, set = JuMP.jump_function(constraint), JuMP.moi_set(constraint)

    if lhs isa JuMP.NonlinearExpr && fun isa JuMP.NonlinearExpr
        return JuMP.build_constraint(
            error_fn,
            JuMP.NonlinearExpr(:(==), lhs, fun),
            MOI.EqualTo(0.0),
        )
    end

    if lhs isa JuMP.NonlinearExpr && fun isa JuMP.AffExpr
        rhs = _handle_affExpr(fun, set)
        return JuMP.build_constraint(error_fn, JuMP.NonlinearExpr(:(==), lhs, rhs), set)
    end

    return JuMP.VectorConstraint([lhs, fun], set)
end

# Adding support for transitions constraints
## Transition
struct TransitionAsModel
    model::JuMP.Model
    mode_variable::JuMP.VariableRef
    source_mode::Int
    destination_mode::Int
    switching_type::HybridSystems.AbstractSwitching
end

## Constructor for Transition
function TransitionAsModel(
    model::JuMP.GenericModel,
    mode_variable::JuMP.VariableRef,
    source_mode::Int,
    destination_mode::Int,
    switching_type::HybridSystems.AbstractSwitching = HybridSystems.ControlledSwitching(),
)
    return TransitionAsModel(
        model,
        mode_variable,
        source_mode,
        destination_mode,
        switching_type,
    )
end

JuMP._valid_model(::TransitionAsModel, ::Any) = nothing

function JuMP.model_convert(t::TransitionAsModel, con::Any)
    return JuMP.model_convert(t.model, con)
end

function _guard end
guard = JuMP.NonlinearOperator(_guard, :guard)

function _resetmaps end
resetmaps = JuMP.NonlinearOperator(_resetmaps, :resetmaps)

function _get_val_from_affExpr(set)
    if set isa MOI.LessThan
        return set.upper
    elseif set isa MOI.GreaterThan
        return set.lower
    elseif set isa MOI.EqualTo
        return set.value
    end
end

function _handle_affExpr(fun, set)
    if !(fun isa JuMP.AffExpr)
        error("Unsupported function type. `AffExpr` expected.")
    end

    if !(set isa MOI.LessThan || set isa MOI.GreaterThan || set isa MOI.EqualTo)
        error("Unsupported set type. `LessThan`, `GreaterThan` or `EqualTo` expected.")
    end

    return JuMP.NonlinearExpr(:-, fun, _get_val_from_affExpr(set))
end

function JuMP.add_constraint(t::TransitionAsModel, condition, ::String)
    if !(condition isa JuMP.ScalarConstraint)
        error("Unsupported constraint type. `ScalarConstraint` expected.")
    end

    if !(condition.func isa JuMP.NonlinearExpr || condition.func isa JuMP.AffExpr)
        error("Unsupported function type. `NonlinearExpr` or `AffExpr` expected.")
    end

    model = t.model
    mode_variable = t.mode_variable
    source_mode = t.source_mode
    destination_mode = t.destination_mode
    switching_type = t.switching_type

    # Build the lhs expression
    lhs_expr = if condition.func isa JuMP.NonlinearExpr
        JuMP.@expression(
            model,
            resetmaps(mode_variable == source_mode, mode_variable == destination_mode)
        )
    else
        JuMP.@expression(
            model,
            guard(mode_variable == source_mode, mode_variable == destination_mode)
        )
    end

    #=return JuMP._build_indicator_constraint(
        JuMP.error,
        lhs_expr,
        condition,
        MOI.Indicator{MOI.ACTIVATE_ON_ONE},
    )=#

    # Build the rhs expression
    set = MOI.EqualTo(0.0)
    rhs_expr = JuMP.jump_function(condition)
    if rhs_expr isa JuMP.AffExpr
        rhs_expr = _handle_affExpr(rhs_expr, condition.set)
        set = condition.set
    end

    # Combine the lhs and rhs expressions
    expr = JuMP.NonlinearExpr(:(==), lhs_expr, rhs_expr)

    # Build the constraint
    return JuMP.build_constraint(JuMP.error, expr, set)
end

function add_transition!(
    f::Function,
    model::JuMP.GenericModel,
    mode_variable::JuMP.VariableRef,
    from::Int,
    to::Int,
    switching_type::HybridSystems.AbstractSwitching = HybridSystems.ControlledSwitching(),
)
    t = TransitionAsModel(model, mode_variable, from, to, switching_type)
    f(t)
    return t
end

# Support incremental interface
MOI.supports_incremental_interface(::Optimizer) = true

function MOI.copy_to(model::Optimizer, src::MOI.ModelLike)
    return MOI.Utilities.default_copy_to(model, src)
end

_svec(vec, idx) = SVector([vec[i] for i in idx]...)
_svec_lower(vec, idx) = SVector([vec[i].lower for i in idx]...)
_svec_upper(vec, idx) = SVector([vec[i].upper for i in idx]...)

function _full(model, def, vars, vals)
    v = fill(def, length(model.variables))
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
        ) for o in model.obstacles
    ]
end

function dynamic!(model, mode = DEFAULT_MODE)
    variable_values = [v.value for v in model.variables]
    MOI.eval_constraint(
        model.modes[mode].evaluator,
        model.derivative_values,
        variable_values,
    )
    return
end

function symbolic(func::MathOptSymbolicAD._Function, xu)
    Symbolics.@variables(p[1:length(func.data)])
    e = MathOptSymbolicAD._expr_to_symbolics(
        func.model,
        func.expr,
        p,
        xu[func.ordered_variables],
    )
    return Symbolics.substitute(e, Dict(p[i] => func.data[i] for i in eachindex(p)))
end

function dynamic(model, x_idx, u_idx, mode = DEFAULT_MODE)
    # System eq x' = F_sys(x, u)
    Symbolics.@variables(x[1:length(x_idx)], u[1:length(u_idx)])
    xu = Vector{eltype(x)}(undef, length(model.variables))
    xu[x_idx] = x
    xu[u_idx] = u
    expr = [symbolic(func, xu) for func in model.modes[mode].evaluator.backend.constraints]
    # The second one is the inplace version that we don't want so we ignore it with `_`
    # `expression = Val{false}` is used to directly get a Julia function
    # `cse = true` is used to detect common sub-expressions (like `α` in the path planning example),
    # cfr. https://discourse.julialang.org/t/detecting-function-composition-in-symbolics-jl/115885/6
    # for a discussion on why it's not currently enabled by default
    F_sys, _ = Symbolics.build_function(
        expr,
        collect(x),
        collect(u);
        expression = Val{false},
        cse = true,
    )

    # To see the generated expression, use:
    #println(Symbolics.build_function(expr, collect(x), collect(u), cse = true))

    return F_sys
end

function _rect(lb, ub, lib, T)
    if lb == -Inf && ub == Inf
        r = FullSpace(1)
    elseif lb == -Inf
        r = hrep([HalfSpace([1], T(ub))])
    elseif ub == Inf
        r = hrep([HalfSpace([-1], -T(lb))])
    else
        r = HalfSpace([-1], -T(lb)) ∩ HalfSpace([1], T(ub))
    end

    return polyhedron(r, lib)
end

function _HyperRectangle_to_polyhedra(hrect, lib, T)
    n = length(hrect.lb) # number of dimensions
    lb_vector = hrect.lb
    ub_vector = hrect.ub

    result = Vector{Polyhedra.Rep}(undef, n)

    for i in 1:n
        result[i] = _rect(lb_vector[1], ub_vector[1], lib, T)
    end

    return reduce(∩, result)
end

function _hybrid_system(model, x_idx, u_idx, _X_, _U_)
    lib = CDDLib.Library()
    T::Type = Float64

    println(">>Assembling the hybrid system")

    # Create the automaton
    n_modes = length(model.modes)
    automaton = HybridSystems.GraphAutomaton(n_modes)

    # Add the transitions
    index = 1
    for (_, transition) in model.transitions
        src, dst = transition.source, transition.destination
        HybridSystems.add_transition!(automaton, src, dst, index)
        index += 1
    end

    # Assemble the continuous systems for each mode
    continuous_systems =
        Vector{MathematicalSystems.AbstractContinuousSystem}(undef, n_modes)

    funcs = [
        MathematicalSystems.ConstrainedAffineControlContinuousSystem(
            zeros(Float64, length(x_idx), 1),
            -0.5 .* ones(Float64, length(u_idx), 1),
            fill(0.0, 1),
            _X_,
            _U_,
        ),
        MathematicalSystems.ConstrainedAffineControlContinuousSystem(
            zeros(Float64, length(x_idx), 1),
            zeros(Float64, length(u_idx), 1),
            fill(1.0, 1),
            _X_,
            _U_,
        ),
    ]
    for mode in eachindex(model.modes)
        continuous_systems[mode] = funcs[mode]
    end

    # Assemble the guard for each transition
    #TODO: This is a bit of a hack. We should probably have a better way to handle this
    guard_constraints = [
        Dionysos.Utils.HyperRectangle([19.0], [Inf]),
        Dionysos.Utils.HyperRectangle([-Inf], [21.0]),
    ]

    guard_maps = Vector{MathematicalSystems.AbstractMap}(undef, length(model.transitions))
    index = 1
    for (_, transition) in model.transitions
        guard_maps[index] = MathematicalSystems.ConstrainedLinearControlMap(
            zeros(Float64, length(x_idx), 1),
            zeros(Float64, length(u_idx), 1),
            _HyperRectangle_to_polyhedra(guard_constraints[index], lib, T),
            _HyperRectangle_to_polyhedra(_U_, lib, T),
        )
        index += 1
    end

    # Assemble the switching conditions
    switchings = fill(HybridSystems.ControlledSwitching(), length(model.transitions))

    # Assemble the hybrid system
    system =
        HybridSystems.HybridSystem(automaton, continuous_systems, guard_maps, switchings)

    system.ext[:X] = _X_
    system.ext[:U] = _U_
    system.ext[:q_T] = n_modes
    system.ext[:ntransitions] = length(model.transitions)

    return system
end

function system(
    model,
    x_idx,
    u_idx;
    sysnoise = SVector(0.0, 0.0, 0.0),
    measnoise = SVector(0.0, 0.0, 0.0),
    tstep = 0.3,
    nsys = 5,
)
    _X_ = Dionysos.Utils.HyperRectangle(
        _svec_lower(model.variables, x_idx),
        _svec_upper(model.variables, x_idx),
    )
    _X_ = Dionysos.Utils.LazySetMinus(
        _X_,
        Dionysos.Utils.LazyUnionSetArray(obstacles(model, x_idx)),
    )
    _U_ = Dionysos.Utils.HyperRectangle(
        _svec_lower(model.variables, u_idx),
        _svec_upper(model.variables, u_idx),
    )
    if model.time_type == CONTINUOUS
        return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
            dynamic(model, x_idx, u_idx),
            Dionysos.Utils.get_dims(_X_),
            Dionysos.Utils.get_dims(_U_),
            _X_,
            _U_,
        )
    elseif model.time_type == DISCRETE
        # `dynamic(model, x_idx, u_idx)` is a function `(x_k, u) -> x_{k+1}`
        # However, `UniformGridAbstraction` will call it with `(x, u, t)` so with
        # an additional time argument. I would expect it to error because there is
        # too many arguments but it seems that `RuntimeGeneratedFunction` just
        # ignores any trailing arguments. To be safe, we still ignore the time
        # like so:
        dyn = dynamic(model, x_idx, u_idx)
        return MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem(
            (x, u) -> dyn(x, u),
            Dionysos.Utils.get_dims(_X_),
            Dionysos.Utils.get_dims(_U_),
            _X_,
            _U_,
        )
    elseif model.time_type == HYBRID
        return _hybrid_system(model, x_idx, u_idx, _X_, _U_)
    else
        error("The system type is unknown")
    end
end

function problem(model::Optimizer)
    x_idx = state_indices(model)
    u_idx = input_indices(model)
    _I_ = Dionysos.Utils.HyperRectangle(
        _svec([s.start.lower for s in model.variables], x_idx),
        _svec([s.start.upper for s in model.variables], x_idx),
    )
    _T_ = Dionysos.Utils.HyperRectangle(
        _svec([s.target.lower for s in model.variables], x_idx),
        _svec([s.target.upper for s in model.variables], x_idx),
    )

    sys = system(model, x_idx, u_idx)
    problem = Dionysos.Problem.OptimalControlProblem(
        sys,
        _I_,
        _T_,
        nothing,
        nothing,
        Dionysos.Problem.Infinity(),
    )

    return problem
end

function problem_hybrid(model::Optimizer, lib, T::Type, x_0 = [1.0], N = 11)
    x_idx = state_indices(model)
    u_idx = input_indices(model)
    _I_ = Dionysos.Utils.HyperRectangle(
        _svec([s.start.lower for s in model.variables], x_idx),
        _svec([s.start.upper for s in model.variables], x_idx),
    )
    _T_ = Dionysos.Utils.HyperRectangle(
        _svec([s.target.lower for s in model.variables], x_idx),
        _svec([s.target.upper for s in model.variables], x_idx),
    )

    sys = system(model, x_idx, u_idx)
    q_0 = length(x_idx)

    if q_0 != length(x_0)
        error(
            "The initial state must have the same dimension as the state variables {$q_0} != {$length(x_0)}",
        )
    end

    #TODO: This is a bit of a hack. We should probably have a better way to handle costs properly
    state_cost = Fill(Dionysos.Utils.ZeroFunction(), sys.ext[:q_T])
    transition_cost = Dionysos.Utils.QuadraticControlFunction(ones(T, 1, 1))

    problem = Dionysos.Problem.OptimalControlProblem(
        sys,
        (q_0, x_0),
        sys.ext[:q_T],
        Fill(state_cost, N),
        Fill(Fill(transition_cost, sys.ext[:ntransitions]), N),
        N,
    )
    return problem
end

function variable_type(model::Optimizer, vi::MOI.VariableIndex)
    if model.variables[vi.value].type == MODE
        return MODE
    elseif isnothing(model.modes[DEFAULT_MODE].dynamic[vi.value])
        return INPUT
    else
        return STATE
    end
end

function state_indices(model::Optimizer)
    return findall(eachindex(model.variables)) do i
        return variable_type(model, MOI.VariableIndex(i)) == STATE
    end
end

function input_indices(model::Optimizer)
    return findall(eachindex(model.variables)) do i
        return variable_type(model, MOI.VariableIndex(i)) == INPUT
    end
end

function mode_indices(model::Optimizer)
    return findall(eachindex(model.variables)) do i
        return variable_type(model, MOI.VariableIndex(i)) == MODE
    end
end

function _check_missing(model)
    for key in eachindex(model.modes)
        @show model.modes[key].dynamic
        if all(isnothing, model.modes[key].dynamic)
            error("Missing dynamics. i.e. ∂(x) = f(x, u) or Δ(x) = f(x, u)")
        end
    end

    transition_missing = isempty(model.transitions)
    guard_missing = false
    for key in eachindex(model.transitions)
        @show model.transitions[key].guard
        if all(isnothing, model.transitions[key].guard)
            guard_missing = true
            break
        end
    end

    resetmaps_missing = false
    for key in eachindex(model.transitions)
        @show model.transitions[key].resetmap
        if all(isnothing, model.transitions[key].resetmap)
            resetmaps_missing = true
            break
        end
    end

    if model.time_type == HYBRID && guard_missing && resetmaps_missing
        error("Missing guards and/or resetmaps. i.e. guard(x, u, t) and resetmaps(x, u)")
    end

    @show [v.type for v in model.variables]
end

function setup!(model::Optimizer)
    _check_missing(model)

    println(">>Setting up the model")
    input_index = 0
    state_index = 0
    mode_index = 0

    for i in eachindex(model.variables)
        var_type = variable_type(model, MOI.VariableIndex(i))
        if var_type == STATE
            # The set does not matter
            for mode in eachindex(model.modes)
                if isnothing(model.modes[mode].dynamic[i])
                    error("Missing dynamics. i.e. ∂(x) = f(x, u) or Δ(x) = f(x, u)")
                end

                MOI.Nonlinear.add_constraint(
                    model.modes[mode].nlp_model,
                    model.modes[mode].dynamic[i],
                    MOI.EqualTo(0.0),
                )
            end

            state_index += 1
            model.variables[i].index = state_index
        elseif var_type == INPUT
            input_index += 1
            model.variables[i].index = input_index
        elseif var_type == MODE
            mode_index += 1
            model.variables[i].index = mode_index
        end
    end
    model.derivative_values = fill(NaN, state_index)
    backend = MathOptSymbolicAD.DefaultBackend()
    vars = MOI.VariableIndex.(eachindex(model.variables))

    for mode in eachindex(model.modes)
        model.modes[mode].evaluator =
            MOI.Nonlinear.Evaluator(model.modes[mode].nlp_model, backend, vars)
        MOI.initialize(model.modes[mode].evaluator, Symbol[])
    end

    println(">>Model setup complete")
    return
end

function MOI.optimize!(model::Optimizer)
    setup!(model)

    if model.time_type == HYBRID
        println(">>Solving the hybrid system")
        qp_solver = JuMP.optimizer_with_attributes(
            OSQP.Optimizer,
            "eps_abs" => 1e-8,
            "eps_rel" => 1e-8,
            "max_iter" => 100000,
            MOI.Silent() => true,
        )

        mip_solver = JuMP.optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)

        cont_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, MOI.Silent() => true)

        miqp_solver = JuMP.optimizer_with_attributes(
            Pavito.Optimizer,
            "mip_solver" => mip_solver,
            "cont_solver" => cont_solver,
            MOI.Silent() => true,
        )

        algo = JuMP.optimizer_with_attributes(
            Dionysos.Optim.BemporadMorari.Optimizer{Float64},
            "continuous_solver" => qp_solver,
            "mixed_integer_solver" => miqp_solver,
            "indicator" => false,
            "log_level" => 0,
        )
        model.inner = MOI.instantiate(algo)
        MOI.set(
            model.inner,
            MOI.RawOptimizerAttribute("problem"),
            problem_hybrid(model, CDDLib.Library(), Float64),
        )
    else # CONTINUOUS or DISCRETE
        println(">>Solving the {model.time_type} system")
        MOI.set(model.inner, MOI.RawOptimizerAttribute("concrete_problem"), problem(model))
    end

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

function _diff end
∂ = JuMP.NonlinearOperator(_diff, :∂)

function _delta end
Δ = JuMP.NonlinearOperator(_delta, :Δ)

function _final end
final = JuMP.NonlinearOperator(_final, :final)

function _start end
start = JuMP.NonlinearOperator(_start, :start)

function rem end
rem_op = JuMP.NonlinearOperator(rem, :rem)

function _themode end
themode = JuMP.NonlinearOperator(_themode, :themode)

# Type piracy
function JuMP.parse_constraint_call(
    error_fn::Function,
    vectorized::Bool,
    ::Val{:∉},
    lhs,
    rhs,
)
    @assert !vectorized
    f, parse_code1 = JuMP._rewrite_expression(lhs)
    set, parse_code2 = JuMP._rewrite_expression(rhs)
    parse_code = quote
        $parse_code1
        $parse_code2
    end
    build_call = :(build_constraint($error_fn, $f, $(OuterSet)($set)))
    return parse_code, build_call
end

function Base.rem(x::JuMP.AbstractJuMPScalar, d)
    return rem_op(x, d)
end

function set_algorithm!(model::JuMP.GenericModel, algorithm::Any)
    model.moi_backend.optimizer.model.inner = MOI.instantiate(algorithm)
    return
end
