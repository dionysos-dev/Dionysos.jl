using JuMP

# Create an InfiniteOpt model
model = Model(Dionysos.UniformGrid)

# Define the state variables: x1(t), x2(t), x3(t)
@variable(model, x[1:3], start = 0.0)

# Define the control variables: u1(t), u2(t)
@variable(model, -1 <= u[1:2] <= 1)

# Set α(t) = arctan(tan(u2(t)) / 2)
@expression(model, α, atan(tan(u[2]) / 2))

function diff end

dot_diff = JuMP.NonlinearOperator(diff, :coucou)

@constraint(model, dot_diff(x[1]) == u[1] * cos(α + x[3]) * sec(α))
@constraint(model, dot_diff(x[2]) == u[1] * sin(α + x[3]) * sec(α))
@constraint(model, dot_diff(x[3]) == u[1] * tan(u[2]))

x_target = [10.0, 10, 0]

function _target end

target = JuMP.NonlinearOperator(_target, :target)

@constraint(model, target(x[1]) == x_target[1])
@constraint(model, target(x[2]) == x_target[2])
@constraint(model, target(x[3]) == x_target[3])

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

function obstacle(x_lb, y_lb, x_ub, y_ub)
    return MOI.HyperRectangle(
        [x_lb, y_lb, -Inf],
        [x_ub, y_ub, Inf],
    )
end

# Function to add rectangular obstacle avoidance constraints

#rect = JuMP.NonlinearOperator(obstacle, :rect)

function _not_in end

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
    build_call = :(build_constraint($error_fn, $f, $OuterSet($set)))
    return parse_code, build_call
end

struct OuterSet{S<:MOI.AbstractVectorSet} <: MOI.AbstractVectorSet
    inner::S
end

MOI.dimension(set::OuterSet) = MOI.dimension(set.inner)
Base.copy(set::OuterSet) = OuterSet(copy(set.inner))

for i in eachindex(x1_ub)
    @constraint(model, x ∉ obstacle(x1_lb[i], x2_lb[i], x1_ub[i], x2_ub[i]))
end

# Objective: Minimize the distance to the target position at final time T
@objective(model, Min, (target(x[1]) - x_target[1])^2 + (target(x[2]) - x_target[2])^2)

optimize!(model)
