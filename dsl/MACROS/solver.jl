# Define the Control class with grid type
mutable struct Control
    name::String
    grid_type::GridType
    variables::Vector{Variable}
    constraints::Vector{Constraint}
    objective::Union{Objective, Nothing}
    model::Model

    function Control(name::String = "Model", grid_type::GridType = UniformGrid())
        new(name, grid_type, Vector{Variable}(), Vector{Constraint}(), nothing, Model(Ipopt.Optimizer))
    end
end

# Macros for variable, constraint, and objective definition
macro variable(ctrl, name, domain, bounds=nothing)
    quote
        var = Variable($name, $domain, bounds=$bounds)
        push!($ctrl.variables, var)
        var
    end
end

macro constraint(ctrl, expr)
    quote
        constr = Constraint($expr)
        push!($ctrl.constraints, constr)
    end
end

macro objective(ctrl, type, expr)
    quote
        $ctrl.objective = Objective($type, $expr)
    end
end

# Solver method to handle different grid types (expand this as needed)
function (ctrl::Control).solve(; horizon::Int = 1)
    println("Solving with grid type: $(typeof(ctrl.grid_type))")

    # Initialize variables and constraints using JuMP and solve with Ipopt or other solver
    jump_vars = Dict{String, JuMP.VariableRef}()

    for var in ctrl.variables
        if var.bounds === nothing
            jump_var = @variable(ctrl.model, var.domain isa PositiveReals ? var >= 0 : var)
        else
            lower, upper = var.bounds
            jump_var = @variable(ctrl.model, lower <= var <= upper)
        end
        jump_vars[var.name] = jump_var
    end

    for constr in ctrl.constraints
        constraint_expr = constr.expr
        eval_expr = eval(Expr(:quote, constraint_expr))
        @constraint(ctrl.model, eval_expr)
    end

    if ctrl.objective !== nothing
        objective_expr = ctrl.objective.expression
        eval_obj_expr = eval(Expr(:quote, objective_expr))
        if ctrl.objective.objective_type isa Minimize
            @objective(ctrl.model, Min, eval_obj_expr)
        elseif ctrl.objective.objective_type isa Maximize
            @objective(ctrl.model, Max, eval_obj_expr)
        end
    end

    optimize!(ctrl.model)

    solution_data = Dict()
    for (name, var) in jump_vars
        solution_data[name] = value(var)
    end

    return Solution("dict", horizon, solution_data)
end
