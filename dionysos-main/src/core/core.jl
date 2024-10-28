# core.jl

module CoreControls

using LinearAlgebra  # For matrices and vector operations
using Plots          # For viz
using JuMP           # Mathematical modeling
#using Ipopt          # For optimization solver

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
end

# Constructor for Constraint
function Constraint(expr::Expr)
    return Constraint(expr)
end

# Objective Types
abstract type ObjectiveType end
struct Minimize <: ObjectiveType end
struct Maximize <: ObjectiveType end

# Objective Struct
mutable struct Objective
    objective_type::ObjectiveType
    expression::Expr
end

# Constructor for Objective
function Objective(objective_type::ObjectiveType, expression::Expr)
    return Objective(objective_type, expression)
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

    function Control(name::String = "Model", problem_type::ProblemType = CustomProblem(), system_type::SystemType = Auto())
        return new(name, problem_type, system_type, Dict{String, Variable}(), Vector{Constraint}(), nothing, JuMP.Model(Ipopt.Optimizer))
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

end  # module CoreControls
