using JuMP
include(joinpath(@__DIR__, "../src/core/opt.jl"))

function _diff end
∂ = JuMP.NonlinearOperator(_diff, :diff)

function _final end
final = JuMP.NonlinearOperator(_final, :final)

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
    build_call = :(build_constraint($error_fn, $f, $(Opt.OuterSet)($set)))
    return parse_code, build_call
end

function obstacle(x_lb, y_lb, x_ub, y_ub)
    return MOI.HyperRectangle(
        [x_lb, y_lb, -Inf],
        [x_ub, y_ub, Inf],
    )
end

# Create an InfiniteOpt model
model = Model(Opt.Optimizer)

# Define the state variables: x1(t), x2(t), x3(t)
#x_low, x_upp = [-0.1, -0.1, -pi - 0.4], [10.1, 10.1, pi + 0.4]
x_low, x_upp = [0.0, 0.0, -pi - 0.4], [4.0, 10.0, pi + 0.4]
x_start = [0.4, 0.4, 0.0]
@variable(model, x_low[i] <= x[i=1:3] <= x_upp[i], start = x_start[i])
#@variable(model, 0.0 <= x1 <= 4.0, start = 0.0)
#@variable(model, 0.0 <= x2 <= 4.0, start = 0.0)

# Define the control variables: u1(t), u2(t)
@variable(model, -1 <= u[1:2] <= 1)

# Set α(t) = arctan(tan(u2(t)) / 2)
@expression(model, α, atan(tan(u[2]) / 2))

@constraint(model, ∂(x[1]) == u[1] * cos(α + x[3]) * sec(α))
@constraint(model, ∂(x[2]) == u[1] * sin(α + x[3]) * sec(α))
@constraint(model, ∂(x[3]) == u[1] * tan(u[2]))

#x_target = [10.0, 10, 0]
x_target = [3.3, 0.5, 0]

@constraint(model, final(x[1]) in MOI.Interval(3.0, 3.6))
@constraint(model, final(x[2]) in MOI.Interval(0.3, 0.8))
@constraint(model, final(x[3]) in MOI.Interval(-100.0, 100.0))

# Obstacle boun daries (provided)
x1_lb = [1.0, 2.2, 2.2] #, 3.4, 4.6, 5.8, 5.8, 7.0, 8.2, 8.4, 9.3, 8.4, 9.3, 8.4, 9.3]
x1_ub = [1.2, 2.4, 2.4] #, 3.6, 4.8, 6.0, 6.0, 7.2, 8.4, 9.3, 10.0, 9.3, 10.0, 9.3, 10.0]
x2_lb = [0.0, 0.0, 6.0] #, 0.0, 1.0, 0.0, 7.0, 1.0, 0.0, 8.2, 7.0, 5.8, 4.6, 3.4, 2.2]
x2_ub = [9.0, 5.0, 10.0] #, 9.0, 10.0, 6.0, 10.0, 10.0, 8.5, 8.6, 7.4, 6.2, 5.0, 3.8, 2.6]

# Function to add rectangular obstacle avoidance constraints

for i in eachindex(x1_ub)
    @constraint(model, x[1:2] ∉ MOI.HyperRectangle([x1_lb[i], x2_lb[i]], [x1_ub[i], x2_ub[i]]))
end

#set_attribute(model, "L_growthbound", ...)

x0 = SVector(0.0, 0.0, 0.0);
h = SVector(0.2, 0.2, 0.2);
set_attribute(model, "state_grid", Dionysos.Domain.GridFree(x0, h))

# Definition of the grid of the input-space on which the abstraction is based (origin `u0` and input-space discretization `h`):
u0 = SVector(0.0, 0.0);
h = SVector(0.3, 0.3);
set_attribute(model, "input_grid", Dionysos.Domain.GridFree(u0, h))

# Objective: Minimize the distance to the target position at final time T
@objective(model, Min, (target(x[1]) - x_target[1])^2 + (target(x[2]) - x_target[2])^2)

#set_optimizer(model, Opt.Optimizer)

optimize!(model)

# Get the results
abstract_system = get_attribute(model, "abstract_system");
abstract_problem = get_attribute(model, "abstract_problem");
abstract_controller = get_attribute(model, "abstract_controller");
concrete_controller = get_attribute(model, "concrete_controller")
concrete_problem = get_attribute(model, "concrete_problem");
concrete_system = concrete_problem.system

using Test
@test length(abstract_controller.data) == 19400 #src

# ### Trajectory display
# We choose a stopping criterion `reached` and the maximal number of steps `nsteps` for the sampled system, i.e. the total elapsed time: `nstep`*`tstep`
# as well as the true initial state `x0` which is contained in the initial state-space `_I_` defined previously.
nstep = 100
function reached(x)
    if x ∈ concrete_problem.target_set
        return true
    else
        return false
    end
end
using StaticArrays
x0 = SVector(0.4, 0.4, 0.0)
import Dionysos
control_trajectory = Dionysos.System.get_closed_loop_trajectory(
    concrete_system.f,
    concrete_controller,
    x0,
    nstep;
    stopping = reached,
)

using Plots

# Here we display the coordinate projection on the two first components of the state space along the trajectory.
fig = plot(; aspect_ratio = :equal);
# We display the concrete domain
plot!(concrete_system.X; color = :yellow, opacity = 0.5);

# We display the abstract domain
plot!(abstract_system.Xdom; color = :blue, opacity = 0.5);

# We display the concrete specifications
plot!(concrete_problem.initial_set; color = :green, opacity = 0.2);
plot!(concrete_problem.target_set; dims = [1, 2], color = :red, opacity = 0.2);

# We display the abstract specifications
plot!(
    Dionysos.Symbolic.get_domain_from_symbols(abstract_system, abstract_problem.initial_set);
    color = :green,
);
plot!(
    Dionysos.Symbolic.get_domain_from_symbols(abstract_system, abstract_problem.target_set);
    color = :red,
);

# We display the concrete trajectory
plot!(control_trajectory; ms = 0.5)
