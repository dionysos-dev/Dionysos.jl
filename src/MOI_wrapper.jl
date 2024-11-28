import MathOptInterface as MOI
import StaticArrays as SA
import MathematicalSystems, HybridSystems
import JuMP
import MathOptSymbolicAD
import Symbolics

export ∂, Δ, final, start, rem, TransitionAsModel, add_transition!, guard, resetmaps, themode

@enum(VariableType, INPUT, STATE, MODE)

@enum(TimeType, UNKNOWN, CONTINUOUS, DISCRETE, HYBRID)

DEFAULT_MODE = 1

mutable struct Optimizer <: MOI.AbstractOptimizer
    inner::Any
    nlp_model::Any
    evaluator::Union{Nothing, MOI.Nonlinear.Evaluator}
    time_type::TimeType
    variable_values::Vector{Float64}
    derivative_values::Vector{Float64}
    variable_types::Vector{VariableType}
    variable_index::Vector{Int} # MOI.VariableIndex -> state/input/mode index
    lower::Vector{Float64}
    upper::Vector{Float64}
    start::Vector{MOI.Interval{Float64}}
    target::Vector{MOI.Interval{Float64}}
    dynamic::Dict{Int64, Vector{Union{Nothing, MOI.ScalarNonlinearFunction}}} # mode_name -> [dynamic_func for each state] (in continuous and discrete time there is only one mode_name)
    obstacles::Vector{Tuple{Vector{MOI.VariableIndex}, MOI.HyperRectangle}}
    objective_sense::MOI.OptimizationSense
    objective_function::MOI.AbstractScalarFunction
    nonlinear_index::Int
    mode_variable::Union{Nothing, MOI.VariableIndex}
    guards_table::Dict{
        String,
        Vector{Union{Nothing, Tuple{Int64, Int64, MOI.ScalarNonlinearFunction}}},
    }  # transition_index -> [(src, dst, guard_func) for each state]
    resetmaps_table::Dict{
        String,
        Vector{Union{Nothing, Tuple{Int64, Int64, MOI.ScalarNonlinearFunction}}},
    }  # transition_index -> [(src, dst, resetmap_func) for each state]
    integer_variables::Dict{MOI.VariableIndex, Union{MOI.Integer, MOI.ZeroOne}}
    indicator_constraints::Vector{
        Tuple{MOI.VariableIndex, Bool, MOI.ScalarNonlinearFunction},
    }
    function Optimizer()
        return new(
            MOI.instantiate(Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer),
            nothing,
            nothing,
            UNKNOWN,
            Float64[],
            Float64[],
            VariableType[],
            Int[],
            Float64[],
            Float64[],
            MOI.Interval{Float64}[],
            MOI.Interval{Float64}[],
            Dict{Int64, Vector{Union{Nothing, MOI.ScalarNonlinearFunction}}}(),
            Tuple{Vector{MOI.VariableIndex}, MOI.HyperRectangle}[],
            MOI.FEASIBILITY_SENSE,
            zero(MOI.ScalarAffineFunction{Float64}),
            0,
            nothing,
            Dict{
                String,
                Vector{Union{Nothing, Tuple{Int64, Int64, MOI.ScalarNonlinearFunction}}},
            }(),
            Dict{
                String,
                Vector{Union{Nothing, Tuple{Int64, Int64, MOI.ScalarNonlinearFunction}}},
            }(),
            Dict{MOI.VariableIndex, Union{MOI.Integer, MOI.ZeroOne}}(),
            Tuple{MOI.VariableIndex, Bool, MOI.ScalarNonlinearFunction}[],
        )
    end
end

MOI.is_empty(model::Optimizer) = isempty(model.variable_types)

function MOI.empty!(model::Optimizer)
    model.time_type = UNKNOWN
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
    model.mode_variable = nothing
    empty!(model.guards_table)
    empty!(model.resetmaps_table)
    empty!(model.integer_variables)
    empty!(model.indicator_constraints)
    return
end

function MOI.add_variable(model::Optimizer)
    println(">>Adding variable")
    push!(model.variable_values, NaN)
    push!(model.variable_types, INPUT)
    push!(model.variable_index, 0)
    push!(model.lower, -Inf)
    push!(model.upper, Inf)
    push!(model.start, MOI.Interval(-Inf, Inf))
    push!(model.target, MOI.Interval(-Inf, Inf))
    if !(haskey(model.dynamic, DEFAULT_MODE))
        model.dynamic[DEFAULT_MODE] = Vector{Union{Nothing, MOI.ScalarNonlinearFunction}}()
    end
    push!(model.dynamic[DEFAULT_MODE], nothing)
    return MOI.VariableIndex(length(model.upper))
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
    model.lower[func.value] = MOI.constant(set)
    return MOI.ConstraintIndex{typeof(func), typeof(set)}(func.value)
end

function MOI.add_constraint(model::Optimizer, func::MOI.VariableIndex, set::MOI.LessThan)
    model.upper[func.value] = MOI.constant(set)
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
            model.target[var.value] = set
            return _new_nonlinear_index!(model, func, set)
        end

        if var isa MOI.VariableIndex && func.head == :start
            model.start[var.value] = set
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
    if !haskey(model.dynamic, mode)
        model.dynamic[mode] = Vector{Union{Nothing, MOI.ScalarNonlinearFunction}}()
        # initialize the dynamic functions with nothing
        for _ in 1:length(model.variable_index)
            push!(model.dynamic[mode], nothing)
        end
    end

    model.time_type = type
    return model.dynamic[mode][val] = func
end

function _update_guard_table!(model, src, dst, func, val)
    if isempty(model.guards_table)
        model.guards_table = Dict{
            String,
            Vector{Union{Nothing, Tuple{Int64, Int64, MOI.ScalarNonlinearFunction}}},
        }()
    end

    transition = "$src->$dst"
    if !haskey(model.guards_table, transition)
        model.guards_table[transition] = []

        for _ in 1:length(model.variable_index)
            push!(model.guards_table[transition], nothing)
        end
    end

    return model.guards_table[transition][val] = (src, dst, func)
end

function _update_resetmaps_table!(model, src, dst, func, val)
    if isempty(model.resetmaps_table)
        model.resetmaps_table = Dict{
            String,
            Vector{Union{Nothing, Tuple{Int64, Int64, MOI.ScalarNonlinearFunction}}},
        }()
    end

    transition = "$src->$dst"
    if !haskey(model.resetmaps_table, transition)
        model.resetmaps_table[transition] = []

        for _ in 1:length(model.variable_index)
            push!(model.resetmaps_table[transition], nothing)
        end
    end

    return model.resetmaps_table[transition][val] = (src, dst, func)
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
        model.variable_types[lhs_var.value] = MODE

        _update_dynamic_dict!(model, HYBRID, rhs_rhs, rhs_var.value, lhs_val)
    end

    # handling guard
    if lhs.head == :guard
        lhs_var, lhs_src, lhs_dst =
            lhs.args[1].args[1], lhs.args[1].args[2], lhs.args[2].args[2]

        @show rhs
        @show rhs.args[1]
        @show set

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

# Adding support for primal start
function MOI.supports(::Optimizer, ::MOI.VariablePrimalStart, ::Type{MOI.VariableIndex})
    return true
end

function MOI.set(model::Optimizer, ::MOI.VariablePrimalStart, vi::MOI.VariableIndex, value)
    # create a MOI.Interval from the value if value is a scalar
    if !isa(value, MOI.Interval)
        value = MOI.Interval(value, value)
    end

    return model.start[vi.value] = value
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
    return TransitionAsModel(model, mode_variable, source_mode, destination_mode, switching_type)
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
        ) for o in model.obstacles
    ]
end

function dynamic!(model)
    MOI.eval_constraint(model.evaluator, model.derivative_values, model.variable_values)
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

function dynamic(model, x_idx, u_idx)
    # System eq x' = F_sys(x, u)
    Symbolics.@variables(x[1:length(x_idx)], u[1:length(u_idx)])
    xu = Vector{eltype(x)}(undef, length(model.variable_index))
    xu[x_idx] = x
    xu[u_idx] = u
    expr = [symbolic(func, xu) for func in model.evaluator.backend.constraints]
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

function _hybrid_system(model, x_idx, u_idx)
    error("Hybrid systems are not yet supported")
    return
end

function system(
    model,
    x_idx,
    u_idx;
    sysnoise = SA.SVector(0.0, 0.0, 0.0),
    measnoise = SA.SVector(0.0, 0.0, 0.0),
    tstep = 0.3,
    nsys = 5,
)
    _X_ =
        Dionysos.Utils.HyperRectangle(_svec(model.lower, x_idx), _svec(model.upper, x_idx))
    _X_ = Dionysos.Utils.LazySetMinus(
        _X_,
        Dionysos.Utils.LazyUnionSetArray(obstacles(model, x_idx)),
    )
    _U_ =
        Dionysos.Utils.HyperRectangle(_svec(model.lower, u_idx), _svec(model.upper, u_idx))
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
            (x, u, _) -> dyn(x, u),
            Dionysos.Utils.get_dims(_X_),
            Dionysos.Utils.get_dims(_U_),
            _X_,
            _U_,
        )
    elseif model.time_type == HYBRID
        return _hybrid_system(model, x_idx, u_idx)
    else
        error("The system type is unknown")
    end
end

function problem(model::Optimizer)
    x_idx = state_indices(model)
    u_idx = input_indices(model)
    _I_ = Dionysos.Utils.HyperRectangle(
        _svec([s.lower for s in model.start], x_idx),
        _svec([s.upper for s in model.start], x_idx),
    )
    _T_ = Dionysos.Utils.HyperRectangle(
        _svec([s.lower for s in model.target], x_idx),
        _svec([s.upper for s in model.target], x_idx),
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

function variable_type(model::Optimizer, vi::MOI.VariableIndex)
    if model.variable_types[vi.value] == MODE
        return MODE
    elseif isnothing(model.dynamic[DEFAULT_MODE][vi.value])
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

function mode_indices(model::Optimizer)
    return findall(eachindex(model.variable_index)) do i
        return variable_type(model, MOI.VariableIndex(i)) == MODE
    end
end

function _check_missing(model)
    for key in eachindex(model.dynamic)
        @show model.dynamic[key]
        if all(isnothing, model.dynamic[key])
            error("Missing dynamics. i.e. ∂(x) = f(x, u) or Δ(x) = f(x, u)")
        end
    end

    guard_missing = isempty(model.guards_table)
    for key in eachindex(model.guards_table)
        @show model.guards_table[key]
        if all(isnothing, model.guards_table[key])
            guard_missing = true
            break
        end
    end

    resetmaps_missing = isempty(model.resetmaps_table)
    for key in eachindex(model.resetmaps_table)
        @show model.resetmaps_table[key]
        if all(isnothing, model.resetmaps_table[key])
            resetmaps_missing = true
            break
        end
    end

    if model.time_type == HYBRID && guard_missing && resetmaps_missing
        error("Missing guards and/or resetmaps. i.e. guard(x, u, t) and resetmaps(x, u)")
    end

    @show model.variable_types
end

function setup!(model::Optimizer)
    _check_missing(model)

    println(">>Setting up the model")
    input_index = 0
    state_index = 0
    mode_index = 0
    model.nlp_model = MOI.Nonlinear.Model()
    for i in eachindex(model.variable_types)
        var_type = variable_type(model, MOI.VariableIndex(i))
        if var_type == STATE
            # The set does not matter
            MOI.Nonlinear.add_constraint(
                model.nlp_model,
                model.dynamic[DEFAULT_MODE][i],
                MOI.EqualTo(0.0),
            )
            state_index += 1
            model.variable_index[i] = state_index
        elseif var_type == INPUT
            input_index += 1
            model.variable_index[i] = input_index
        elseif var_type == MODE
            mode_index += 1
            model.variable_index[i] = mode_index
        end
    end
    model.derivative_values = fill(NaN, state_index)
    backend = MathOptSymbolicAD.DefaultBackend()
    vars = MOI.VariableIndex.(eachindex(model.lower))
    model.evaluator = MOI.Nonlinear.Evaluator(model.nlp_model, backend, vars)
    MOI.initialize(model.evaluator, Symbol[])
    println(">>Model setup complete")
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
