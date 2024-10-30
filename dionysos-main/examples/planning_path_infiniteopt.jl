using InfiniteOpt, LinearAlgebra

# Create an InfiniteOpt model
model = InfiniteModel()

T = 10

# Define the time horizon (0 to T)
@infinite_parameter(model, t ∈ [0, T])

# Define the state variables: x1(t), x2(t), x3(t)
@variable(model, x[1:3], Infinite(t))

# Define the control variables: u1(t), u2(t)
@variable(model, -1 <= u[1:2] <= 1, Infinite(t))

# Set α(t) = arctan(tan(u2(t)) / 2)
@expression(model, α, atan(tan(u[2]) / 2))

@constraint(model, ∂(x[1], t) == u[1] * cos(α + x[3]) * sec(α))
@constraint(model, ∂(x[2], t) == u[1] * sin(α + x[3]) * sec(α))
@constraint(model, ∂(x[3], t) == u[1] * tan(u[2]))

# Initial conditions for the states
x1_initial = 0.0
x2_initial = 0.0
x3_initial = 0.0

@constraint(model, x[1](0) == x1_initial)
@constraint(model, x[2](0) == x2_initial)
@constraint(model, x[3](0) == x3_initial)

# Target position (goal state)
x_target = [10.0, 10, 0]

@constraint(model, x[1](T) == x_target[1])
@constraint(model, x[2](T) == x_target[2])
@constraint(model, x[3](T) == x_target[3])

struct OutsideHyperRectangle{T} <: MOI.AbstractVectorSet
    lower::Vector{T}
    upper::Vector{T}
    function OutsideHyperRectangle(lower::Vector{T}, upper::Vector{T}) where {T}
        l, u = length(lower), length(upper)
        if l != u
            throw(
                ArgumentError(
                    "length of lower (=$l) and upper (=$u) bounds must match.",
                ),
            )
        end
        return new{T}(lower, upper)
    end
end

MOI.dimension(set::OutsideHyperRectangle) = length(set.lower)

function Base.copy(set::OutsideHyperRectangle)
    return OutsideHyperRectangle(copy(set.lower), copy(set.upper))
end


# Obstacle boundaries (provided)
x1_lb = [1.0, 2.2, 2.2, 3.4, 4.6, 5.8, 5.8, 7.0, 8.2, 8.4, 9.3, 8.4, 9.3, 8.4, 9.3]
x1_ub = [1.2, 2.4, 2.4, 3.6, 4.8, 6.0, 6.0, 7.2, 8.4, 9.3, 10.0, 9.3, 10.0, 9.3, 10.0]
x2_lb = [0.0, 0.0, 6.0, 0.0, 1.0, 0.0, 7.0, 1.0, 0.0, 8.2, 7.0, 5.8, 4.6, 3.4, 2.2]
x2_ub = [9.0, 5.0, 10.0, 9.0, 10.0, 6.0, 10.0, 10.0, 8.5, 8.6, 7.4, 6.2, 5.0, 3.8, 2.6]

# Function to add rectangular obstacle avoidance constraints
for i in eachindex(x1_ub)
    rect = OutsideHyperRectangle(
        [x1_lb[i], x2_lb[i], -Inf],
        [x1_ub[i], x2_ub[i], Inf],
    )
    @constraint(model, x in rect)
end

# Objective: Minimize the distance to the target position at final time T
@objective(model, Min, (x[1](T) - x_target[1])^2 + (x[2](T) - x_target[2])^2)


import OrderedCollections

@enum(VariableType, INPUT, STATE, MODE)

variables = OrderedCollections.OrderedDict{GeneralVariableRef,VariableType}()

println(all_variables(model))
vars = all_variables(model)

for var in all_variables(model)
    if var.index_type == InfiniteOpt.InfiniteVariableIndex
        variables[var] = INPUT
    end
end

function detect_variable!(_, c, _, _)
    return
end

function detect_variable!(variables, _, v::GeneralVariableRef, ::MOI.Integer)
    if variables[v] != INPUT
        error("Cannot specify $v to be integer")
    end
    variables[v] = MODE
end

function _underiv(d::GeneralVariableRef)
    return GeneralVariableRef(d.model, d.raw_index, InfiniteOpt.InfiniteVariableIndex, d.param_index)
end

function detect_variable!(variables, _, f::NLPExpr, ::MOI.EqualTo)
    if f.tree_root.data == NodeData(:-) &&
        f.tree_root.child.data.value.index_type == InfiniteOpt.DerivativeIndex
        v = _underiv(f.tree_root.child.data.value)
        @show v
        variables[v] = STATE
    else
        error("Constraint $c not understood")
    end
end

for con_ref in all_constraints(model)
    con_obj = constraint_object(con_ref)
    detect_variable!(variables, con_ref, con_obj.func, con_obj.set)
end

variables

state_variable_index = Dict{InfiniteOpt.GeneralVariableRef,Int}()
input_variable_index = Dict{InfiniteOpt.GeneralVariableRef,Int}()
for (var, t) in variables
    if t == STATE
        state_variable_index[var] = length(state_variable_index) + 1
    elseif t == INPUT
        input_variable_index[var] = length(input_variable_index) + 1
    else
        error("TODO $t")
    end
end

import Dionysos

struct StateInputModel
    state_variable_index::Dict{InfiniteOpt.GeneralVariableRef,Int}
    input_variable_index::Dict{InfiniteOpt.GeneralVariableRef,Int}
    start_values::Vector{Float64}
    target_values::Vector{Float64}
    state::Union{Nothing,Dionysos.Utils.HyperRectangle}
    obstacles::Vector{Dionysos.Utils.HyperRectangle}
    dynamic::Vector{Expr}
    function StateInputModel(
        state_variable_index::Dict{InfiniteOpt.GeneralVariableRef,Int},
        input_variable_index::Dict{InfiniteOpt.GeneralVariableRef,Int},
    )
        nstates = length(state_variable_index)
        return new(
            state_variable_index,
            input_variable_index,
            fill(NaN, nstates),
            fill(NaN, nstates),
            nothing,
            Dionysos.Utils.HyperRectangle[],
            Vector{Expr}(undef, nstates),
        )
    end
end

simodel = StateInputModel(state_variable_index, input_variable_index)

function parse_constraint!(_, c, _, _)
    error(c)
end

function parse_constraint!(_, c, v::GeneralVariableRef, ::MOI.GreaterThan)
    if haskey(model.state_variable_index, v)
    elseif haskey(model.input_variable_index, v)
    else
        error("Unknown variable $v")
    end
end

function parse_constraint!(model, _, v::GeneralVariableRef, ::MOI.Integer)
    return
end

function parse_constraint!(model, _, f::NLPExpr, ::MOI.EqualTo)
    if f.tree_root.data == NodeData(:-) &&
        f.tree_root.child.data.value.index_type == InfiniteOpt.DerivativeIndex
        v = _underiv(f.tree_root.child.data.value)
        idx = model.state_variable_index[v]
        model.dynamic[idx] = NLPExpr(f.tree_root.child.sibling)
    else
        error("Constraint $c not understood")
    end
end


for con_ref in all_constraints(model)
    con_obj = constraint_object(con_ref)
    parse_constraint!(variables, con_ref, con_obj.func, con_obj.set)
end

simodel

# Solve the problem using a solver that supports nonlinear optimization
optimize!(model, optimizer_with_attributes(Ipopt.Optimizer))

include(joinpath(@__DIR__, "../src/core/core.jl"))


model = InfiniteModel()
@infinite_parameter(model, t ∈ [0, T])
@variable(model, x, Infinite(t))
e = sin(x) + 1

import MacroTools

macro test(expr)
    new_expr = MacroTools.postwalk(expr) do leaf
       if leaf isa Symbol
       return esc(leaf)
       else
       return leaf
       end
    end
    return expr
end

@test(x + 1)
