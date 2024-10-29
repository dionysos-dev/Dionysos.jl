# core.jl

module CoreControls

using LinearAlgebra  # For matrices and vector operations
using Plots          # For viz
using JuMP           # Mathematical modeling
using Ipopt  

using Dionysos
using StaticArrays

# Exports
## export macros
export @variable, @parameter, @constraint, @control 
export @statevar, @inputvar, @modevar
export @uniformgrid, @algorithm, @samplebased, @lazyellipsoid, @customalgorithm 
export  @objective, @minimize, @maximize

## export types
export ProblemType, Reachability, Safety, CustomProblem
export SystemType, Continuous, Discrete, Hybrid, Auto
export VariableType, StateVar, InputVar, ParameterVar, ModeVar
export Domain, Reals, PositiveReals, Integers
export Variable, Constraint, Objective, Control, Solution, Algorithm
export UniformGridAlgorithm, SampleBasedAlgorithm
export Minimize, Maximize

## export functions
export solve, Visualize, dot, diff


# Problem Types
abstract type ProblemType end
struct Reachability <: ProblemType end
struct Safety <: ProblemType end
struct CustomProblem <: ProblemType end

# System Types
abstract type SystemType end
struct Continuous <: SystemType end
struct Discrete <: SystemType end
struct Hybrid <: SystemType end
struct Auto <: SystemType end  # Automatically detect system type

# Variable Types
abstract type VariableType end
struct StateVar <: VariableType end
struct InputVar <: VariableType end
struct ParameterVar <: VariableType end
struct ModeVar <: VariableType end

# Domains
abstract type Domain end
struct Reals <: Domain end
struct PositiveReals <: Domain end
struct Integers <: Domain end

# Variable Struct
mutable struct Variable
    name::String
    var_type::VariableType
    domain::Domain
    bounds::Union{Tuple{Real, Real}, Nothing}
    value::Union{Nothing, Real}  # For parameters
    jump_vars::Vector{JuMP.VariableRef}  # For time-indexed variables
end

# Overload `getindex` to support time-indexed variable access
function Base.getindex(var::Variable, t::Int)
    return var.jump_vars[t + 1]  # Access the `jump_var` at time step `t`
end

# Constructor for Variable
function Variable(name::String, var_type::VariableType, domain::Domain; bounds::Union{Tuple{Real, Real}, Nothing} = nothing, value::Union{Nothing, Real} = nothing)
    return Variable(name, var_type, domain, bounds, value, Vector{JuMP.VariableRef}())
end

# Constraint Struct
mutable struct Constraint
    expr::Expr  # Use Julia's Expr type to store the expression

    Constraint(expr::Expr) = new(expr)
end

#=mutable struct Constraint
    lhs::Expr
    rhs::Expr
    is_implication::Bool

    function Constraint(expr::Expr)
        if expr.head == :call && expr.args[1] == :(==>) # Implication
            return new(expr.args[2], expr.args[3], true)
        else
            return new(expr, :(), false) # Normal constraint
        end
    end
end=#


macro constraint(model, expr)
    quote
        local model = $(esc(model))
        local constraint_expr = $(QuoteNode(expr))
        constraint = Constraint(constraint_expr)
        push!(model.constraints, constraint)
    end
end

# Objective Types
abstract type ObjectiveType end
struct Minimize <: ObjectiveType end
struct Maximize <: ObjectiveType end

# Objective Struct
mutable struct Objective
    objective_type::ObjectiveType
    expression::Expr

    Objective(objective_type::ObjectiveType, expression::Expr) = new(objective_type, expression)
end

# Control Model Struct
mutable struct Control
    name::String
    problem_type::ProblemType
    system_type::SystemType
    variables::Dict{String, Variable}
    constraints::Vector{Constraint}
    objective::Union{Objective, Nothing}
    model::JuMP.Model

    #function Control(;name::String = "Control", problem_type::ProblemType = CustomProblem(), system_type::SystemType = Auto())
    #    return new(name, problem_type, system_type, Dict{String, Variable}(), Vector{Constraint}(), nothing, JuMP.Model(Ipopt.Optimizer))
    #end
end

Control() = Control("Control", CustomProblem(), Auto(), Dict{String, Variable}(), Vector{Constraint}(), nothing, JuMP.Model(Ipopt.Optimizer))
Control(name::String, problem_type::ProblemType, system_type::SystemType) = Control(name, problem_type, system_type, Dict{String, Variable}(), Vector{Constraint}(), nothing, JuMP.Model(Ipopt.Optimizer))

#=macro control(expr)
    quote
        model = Control("Control", CustomProblem(), Auto(), Dict{String, Variable}(), Vector{Constraint}(), nothing, JuMP.Model(Ipopt.Optimizer))
        $(esc(expr))
        model
    end 
end=#

struct StateModel
    state_variable_index::Dict{Variable,Int}
    start_values::Vector{Float64}
    target_values::Vector{Float64}
    state::Dionysos.Utils.HyperRectangle
    obstacles::Vector{Dionysos.Utils.HyperRectangle}
    dynamic::Vector{Expr}
    function StateModel(state_variable_index::Dict{Variable,Int})
        nstates = length(state_variable_index)
        return new(
            state_variable_index,
            fill(NaN, nstates),
            fill(NaN, nstates),
            nothing,
            Hyperrectangle[],
            Vector{Expr}(undef, nstates),
        )
    end
end


macro control(args...)
    name = "Control"
    problem_type = :(CustomProblem())
    system_type = :(Auto())
    
    for arg in args
        if Meta.isexpr(arg, :(=))
            key, value = arg.args
            if key == :name
                name = value
            elseif key == :problem_type
                problem_type = :($(value))
            elseif key == :system_type
                system_type = :($(value))
            end
        end
    end
    
    quote
        Control($(esc(name)), $(esc(problem_type)), $(esc(system_type)), 
                Dict{String, Variable}(), Vector{Constraint}(), nothing, 
                JuMP.Model(Ipopt.Optimizer))
    end
end



# Helper function to detect system type
function detect_system_type(ctrl::Control)
    has_dot = any(expr_contains_function(constr.expr, :dot) for constr in ctrl.constraints)
    has_diff = any(expr_contains_function(constr.expr, :diff) for constr in ctrl.constraints)
    has_modes = any(var.var_type isa ModeVar for var in values(ctrl.variables))

    if has_modes
        return Hybrid()
    elseif has_dot
        return Continuous()
    elseif has_diff
        return Discrete()
    else
        return Continuous()  # Default to continuous if in doubt
    end
end

# Helper function to check if an expression contains a specific function
function expr_contains_function(expr::Expr, func_sym::Symbol)
    if expr.head == :call && expr.args[1] == func_sym
        return true
    elseif length(expr.args) > 0
        return any(expr_contains_function(arg, func_sym) for arg in expr.args if arg isa Expr)
    else
        return false
    end
end

# Solution Struct
mutable struct Solution
    format::String
    horizon::Int
    solution_data::Dict{String, Any}
end

# Constructor for Solution
function Solution(format::String, horizon::Int, solution_data::Dict{String, Any})
    return Solution(format, horizon, solution_data)
end

# Visualization Function
function Visualize(solution::Solution; plot::String, labels::Vector{String} = [])
    if plot == "state_trajectory"
        time = solution.solution_data["time"]
        for label in labels
            plot!(time, solution.solution_data[label], label=label)
        end
        xlabel!("Time")
        ylabel!("States")
        title!("State Trajectories")
        gui()
    elseif plot == "path"
        x = solution.solution_data[labels[1]]
        y = solution.solution_data[labels[2]]
        plot(x, y, label="Path")
        xlabel!(labels[1])
        ylabel!(labels[2])
        title!("Path Plot")
        gui()
    else
        error("Unknown plot type")
    end
end


# Algorithm Types
abstract type Algorithm end

struct UniformGridAlgorithm <: Algorithm
    state_grid::Dionysos.Domain.Grid
    input_grid::Vector{Float64}
    # Additional parameters can be added as needed
end

struct SampleBasedAlgorithm <: Algorithm
    num_samples::Int
    # Additional parameters can be added as needed
end

struct LazyEllipsoidAlgorithm <: Algorithm
    jacobian::Function
    # Additional parameters can be added as needed
end
struct CustomAlgorithm <: Algorithm
    fun::Function
end

#=
macro algorithm(name, args...)
    if name == :UniformGridAlgorithm
        grid_size = 100
        for arg in args
            if Meta.isexpr(arg, :(=))
                key, value = arg.args
                if key == :grid_size
                    grid_size = value
                end
            end
        end
        return :(UniformGridAlgorithm($grid_size))
    elseif name == :SampleBasedAlgorithm
        num_samples = 1000
        for arg in args
            if Meta.isexpr(arg, :(=))
                key, value = arg.args
                if key == :num_samples
                    num_samples = value
                end
            end
        end
        return :(SampleBasedAlgorithm($num_samples))
    else
        error("Unknown algorithm type")
    end
end=#

# macro for defining an algorithm
macro algorithm(algo_type, args...)
    if algo_type == "uniformgrid"
        return esc(uniformgrid(args...))
    elseif algo_type == :samplebased
        return esc(samplebased(args...))
    else
        error("Unknown algorithm type")
    end
end

macro uniformgrid(x, hx, input_list)
    quote
        state_grid = Dionysos.Domain.GridFree($(esc(x)), $(esc(hx)))
        UniformGridAlgorithm(state_grid, $(esc(input_list)))
    end
end

# get the number of state variables and the number of input variables or mode variables 
# and create the adequate grid and return the algorithm
macro uniformgrid(ctrl)
    quote
        local ctrl = $(esc(ctrl))
        nstates = count(v -> v.var_type isa StateVar, values(ctrl.variables))
        ninputs = count(v -> v.var_type isa InputVar || v.var_type isa ModeVar, values(ctrl.variables))
        
        x0 = SVector{nstates}(zeros(nstates))
        hx = SVector{nstates}(fill(2.0 / 4.0e3, nstates))
        
        state_grid = Dionysos.Domain.GridFree(x0, hx)
        UniformGridAlgorithm(state_grid, collect(1:ninputs))
    end
end
    
macro samplebased(num_samples)
    quote
        SampleBasedAlgorithm($num_samples)
    end
end

macro lazyellipsoid(jacobian)
    quote
        LazyEllipsoidAlgorithm($jacobian)
    end
end

# Macros for variable, parameter, constraint, and objective definition
macro variable(ctrl, name_str, var_type_expr, domain_expr, bounds_expr=:nothing)
    quote
        var_type = $(esc(var_type_expr))
        domain = $(esc(domain_expr))
        bounds = $(esc(bounds_expr)) === nothing ? nothing : $(esc(bounds_expr))
        var = CoreControls.Variable($(esc(name_str)), var_type, domain, bounds=bounds)
        $(esc(ctrl)).variables[$(esc(name_str))] = var
        $(esc(Symbol(name_str))) = var # Define variable in user's scope
    end
end

# function to define state variable
function statemeta(ctrl, name_str, domain_expr, bounds_expr=:nothing)
    @variable(ctrl, name_str, StateVar(), domain_expr, bounds_expr)
end

# Macro for state variable, which is a special case of variable
macro statevar(ctrl, name_str, domain_expr, bounds_expr=:nothing)
    quote
        var_type = StateVar()
        domain = $(esc(domain_expr))
        bounds = $(esc(bounds_expr)) === nothing ? nothing : $(esc(bounds_expr))
        var = CoreControls.Variable($(esc(name_str)), var_type, domain, bounds=bounds)
        $(esc(ctrl)).variables[$(esc(name_str))] = var
        $(esc(Symbol(name_str))) = var # Define variable in user's scope
    end
end

# Macro for input variable, which is a special case of variable
macro inputvar(ctrl, name_str, domain_expr, bounds_expr=:nothing)
    quote
        var_type = InputVar()
        domain = $(esc(domain_expr))
        bounds = $(esc(bounds_expr)) === nothing ? nothing : $(esc(bounds_expr))
        var = CoreControls.Variable($(esc(name_str)), var_type, domain, bounds=bounds)
        $(esc(ctrl)).variables[$(esc(name_str))] = var
        $(esc(Symbol(name_str))) = var # Define variable in user's scope
    end
end

# Macro for mode variable, which is a special case of variable
macro modevar(ctrl, name_str, n_modes)
    quote
        var_type = ModeVar()
        domain = Integers()
        bounds = (1, $(esc(n_modes)))
        var = CoreControls.Variable($(esc(name_str)), var_type, domain, bounds=bounds)
        $(esc(ctrl)).variables[$(esc(name_str))] = var
        $(esc(Symbol(name_str))) = var # Define variable in user's scope
    end
end

#=
macro parameter(ctrl, name_str, value_expr)
    quote
        
        value = $(esc(value_expr))
        var = CoreControls.Variable($(esc(name_str)), ParameterVar(), Reals(), value=value)
        $ctrl.variables[$(esc(name_str))] = var
        $(esc(Symbol(name_str))) = value # Define parameter value in user's scope
    end
end
=#

#=
macro parameter(model_expr, name_str, value_expr)
    quote
        let model = $(esc(model_expr))
            value = $(esc(value_expr))
            var = CoreControls.Variable($(esc(name_str)), ParameterVar(), Reals(), value=value)
            model.variables[$(esc(name_str))] = var
            $(esc(Symbol(name_str))) = value # Define parameter value in user's scope
        end
    end
end
=#

macro parameter(model_expr, name_str, value_expr)
    name_sym = Symbol(name_str)
    quote
        let model = $(esc(model_expr))
            value = $(esc(value_expr))
            var = CoreControls.Variable($(esc(name_str)), ParameterVar(), Reals(), value=value)
            model.variables[$(esc(name_str))] = var
            global $(esc(name_sym)) = value # Define parameter value as a global variable
        end
    end
end

#=macro constraint(ctrl, expr)
    quote
        constr = CoreControls.Constraint($expr)
        push!($ctrl.constraints, constr)
    end
end=#

macro minimize(ctrl, expr)
    quote
        local model = $(esc(ctrl))
        local expr = $(QuoteNode(expr))
        obj_type = Minimize()
        objective = Objective(obj_type, expr)
        model.objective = objective
    end
end

macro maximize(model, expr)
    quote
        local model = $(esc(ctrl))
        local expr = $(QuoteNode(expr))
        obj_type = Maximize()
        objective = Objective(obj_type, expr)
        model.objective = objective
    end
end

macro objective(ctrl, expr, obj_type)
    quote
        local model = $(esc(ctrl))
        local expr = $(QuoteNode(expr))
        objective = Objective($(esc(obj_type)), expr)
        model.objective = objective
    end
end

# Implementing dot() function for continuous-time dynamics
function dot(var::Variable, t::Int, Δt::Float64, ctrl::Control)
    return (var.jump_vars[t+1] - var.jump_vars[t]) / Δt
end

# Implementing diff() function for discrete-time dynamics
function diff(var::Variable, t::Int)
    return var.jump_vars[t+1] - var.jump_vars[t]
end

# Helper function to replace variables in expressions
function _replace_variables_time_indexed(expr::Expr, variables::Dict{String, Variable}, t::Int, Δt::Float64)
    args = [_replace_variables_time_indexed(arg, variables, t, Δt) for arg in expr.args]
    return Expr(expr.head, args...)
end

function replace_variables_time_indexed(state_model::StateModel, expr::Expr, variables::Dict{String, Variable}, t::Int, Δt::Float64)
    if expr.head == :call && expr.args[1] == :(==)
        lhs = expr.args[2]
        if Meta.isexpr(lhs, :call) && lhs.args[1] == :dot
            variable_name = lhs.args[2]
            var = variables[variable_name]
            state_model.dynamic[state_model.state_variable_index[var]] = _replace_variables_time_indexed(expr.args[2], variables, t, Δt)
        else
            dump(lhs)
            error("Not supported")
        end
    else
        dump(expr)
        error("Unsupported")
    end
    return
    #=if expr.head == :call && expr.args[1] == :(==>) # Implication
        lhs = replace_variables_time_indexed(expr.args[2], variables, t, Δt)
        rhs = replace_variables_time_indexed(expr.args[3], variables, t, Δt)
        return Expr(:call, :(==>), lhs, rhs)
    else=#if expr.head == :call
        func = expr.args[1]
        args = [_replace_variables_time_indexed(arg, variables, t, Δt) for arg in expr.args[2:end]]
        if func == :dot
            var_name = string(args[1])
            var = variables[var_name]
            return dot(var, t, Δt, var.model)
        elseif func == :diff
            error("WIP")
            var_name = string(args[1])
            var = variables[var_name]
            return diff(var, t)
        else
            return Expr(:call, func, args...)
        end
    elseif expr.head == :symbol
        var_name = string(expr)
        if haskey(variables, var_name)
            var = variables[var_name]
            if var.var_type isa ParameterVar
                return var.value
            else
                return var.jump_vars[t+1]
            end
        else
            return expr
        end
    else
        args = [replace_variables_time_indexed(arg, variables, t, Δt) for arg in expr.args]
        return Expr(expr.head, args...)
    end
end

# Function to apply constraints to the model
function apply_constraints_to_jump_model(ctrl::Control, jump_model::JuMP.Model, num_time_steps::Int, Δt::Float64)
    for t in 0:num_time_steps-1
        for constraint in ctrl.constraints
            if constraint.is_implication
                lhs = replace_variables_time_indexed(constraint.lhs, ctrl.variables, t, Δt)
                rhs = replace_variables_time_indexed(constraint.rhs, ctrl.variables, t, Δt)
                @constraint(jump_model, !lhs || rhs)
            else
                replaced_expr = replace_variables_time_indexed(constraint.expr, ctrl.variables, t, Δt)
                @constraint(jump_model, eval(replaced_expr))
            end
        end
    end
end

# Solver function
function solve(ctrl::Control, algorithm::Algorithm; horizon::Int = 1, Δt::Float64 = 1.0)

    println("Solving control model: $(ctrl.name)")

    # Detect system type if Auto
    if ctrl.system_type isa Auto
        ctrl.system_type = detect_system_type(ctrl)
        println("Detected system type: $(typeof(ctrl.system_type))")
    end

    # Initialize variables for each time step
    time_steps = 0:horizon
    for t in time_steps
        for var in values(ctrl.variables)
            if var.var_type isa StateVar || var.var_type isa InputVar || var.var_type isa ModeVar
                var_name = "$(var.name)_$t"
                if var.bounds === nothing
                    if var.domain isa Reals
                        jump_var = JuMP.@variable(ctrl.model, base_name=var_name)
                    elseif var.domain isa PositiveReals
                        # create a non-negative JuMP variable
                        jump_var = JuMP.@variable(ctrl.model, base_name=var_name)
                        JuMP.@constraint(ctrl.model, jump_var >= 0)
                    elseif var.domain isa Integers
                        jump_var = JuMP.@variable(ctrl.model, integer=true, base_name=var_name)
                    else
                        error("Unsupported domain for variable $(var.name)")
                    end
                else
                    lower, upper = var.bounds
                    if var.domain isa Reals || var.domain isa Integers
                        jump_var = JuMP.@variable(ctrl.model, base_name=var_name)
                        JuMP.@constraint(ctrl.model, lower <= jump_var <= upper)
                    else
                        error("Unsupported domain for variable $(var.name)")
                    end
                end
                push!(var.jump_vars, jump_var)
            end
            # Parameters are constants; no need to define in JuMP
        end
    end

    state_variable_indices = Dict{Variable,Int}()
    for var in ctrl.variables
        if v.var_type isa StateVar
            state_variable_indices[var] = length(state_variable_indices) + 1
        end
    end
    state_model = StateModel(state_variable_indices)

    # Add constraints to the model
    for constr in ctrl.constraints
        if ctrl.system_type isa Continuous || ctrl.system_type isa Hybrid
            for t in 0:horizon-1
                expr = replace_variables_time_indexed(state_model, constr.expr, ctrl.variables, t, Δt)
                @constraint(ctrl.model, expr)
            end
        elseif ctrl.system_type isa Discrete
            for t in 0:horizon-1
                expr = replace_variables_time_indexed(state_model, constr.expr, ctrl.variables, t, Δt)
                @constraint(ctrl.model, expr)
            end
        else
            error("Unknown system type")
        end
    end

    # Add objective to the model
    if ctrl.objective !== nothing
        total_cost = 0
        if ctrl.system_type isa Continuous || ctrl.system_type isa Hybrid
            for t in 0:horizon-1
                expr = replace_variables_time_indexed(ctrl.objective.expression, ctrl.variables, t, Δt)
                if ctrl.objective.objective_type isa Minimize
                    total_cost += expr * Δt
                elseif ctrl.objective.objective_type isa Maximize
                    total_cost -= expr * Δt
                else
                    error("Unsupported objective type")
                end
            end
        elseif ctrl.system_type isa Discrete
            for t in 0:horizon-1
                expr = replace_variables_time_indexed(ctrl.objective.expression, ctrl.variables, t, Δt)
                if ctrl.objective.objective_type isa Minimize
                    total_cost += expr
                elseif ctrl.objective.objective_type isa Maximize
                    total_cost -= expr
                else
                    error("Unsupported objective type")
                end
            end
        end
        @objective(ctrl.model, Min, total_cost)
    else
        error("Objective function not defined")
    end

    # Solve the model
    optimize!(ctrl.model)

    # Check solution status
    status = termination_status(ctrl.model)
    if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED
        println("Solution found.")
    else
        println("Solver terminated with status: $status")
    end

    # Extract solution data
    solution_data = Dict{String, Any}()
    time = [t * Δt for t in 0:horizon]
    solution_data["time"] = time
    for (name, var) in ctrl.variables
        if var.var_type isa StateVar || var.var_type isa InputVar || var.var_type isa ModeVar
            values = [value(var.jump_vars[t+1]) for t in 0:horizon]
            solution_data[name] = values
        end
    end

    return Solution("dict", horizon, solution_data)
end

# Function to print the Control model
function print(ctrl::Control)
    println("Control Model: $(ctrl.name)")
    println("Problem Type: $(typeof(ctrl.problem_type))")
    println("System Type: $(typeof(ctrl.system_type))")

    println("\nVariables:")
    for (name, var) in ctrl.variables
        println("  $(name): Type=$(typeof(var.var_type)), Domain=$(typeof(var.domain)), Bounds=$(var.bounds)")
    end

    println("\nConstraints:")
    for constr in ctrl.constraints
        println("  $(constr.expr)")
    end

    if ctrl.objective !== nothing
        println("\nObjective:")
        println("  Type=$(typeof(ctrl.objective.objective_type)), Expression=$(ctrl.objective.expression)")
    else
        println("\nObjective: Not defined")
    end
end

# exporting the macros
#export @variable, @parameter, @constraint, @objective
# Add this method to allow for dot notation
#FIXME: LoadError: syntax: invalid function name "ctrl::Control.solve"
###(ctrl::Control).solve(ctrl::Control, args...; kwargs...) = solve(ctrl, args...; kwargs...)

end  # module CoreControls
