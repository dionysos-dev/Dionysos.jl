using Test     #src
# # Example: Path planning problem solved by [Uniform grid abstraction](https://github.com/dionysos-dev/Dionysos.jl/blob/master/docs/src/manual/manual.md#solvers).
#
#md # [![Binder](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Path planning.ipynb)
#md # [![nbviewer](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Path planning.ipynb)
#
# This example was borrowed from [1, IX. Examples, A] whose dynamics comes from the model given in [2, Ch. 2.4].
# This is a **reachability problem** for a **continuous system**.
#
# Let us consider the 3-dimensional state space control system of the form
# ```math
# \dot{x} = f(x, u)
# ```
# with $f: \mathbb{R}^3 × U ↦ \mathbb{R}^3$ given by
# ```math
# f(x,(u_1,u_2)) = \begin{bmatrix} u_1 \cos(α+x_3)\cos(α)^{-1} \\ u_1 \sin(α+x_3)\cos(α)^{-1} \\ u_1 \tan(u_2)  \end{bmatrix}
# ```
# and with $U = [−1, 1] \times [−1, 1]$ and $α = \arctan(\tan(u_2)/2)$. Here, $(x_1, x_2)$ is the position and $x_3$ is the
# orientation of the vehicle in the 2-dimensional plane. The control inputs $u_1$ and $u_2$ are the rear
# wheel velocity and the steering angle.
# The control objective is to drive the vehicle which is situated in a maze made of obstacles from an initial position to a target position.
#
#
# In order to study the concrete system and its symbolic abstraction in a unified framework, we will solve the problem
# for the sampled system with a sampling time $\tau$.
# For the construction of the relations in the abstraction, it is necessary to over-approximate attainable sets of
# a particular cell. In this example, we consider the used of a growth bound function  [1, VIII.2, VIII.5] which is one of the possible methods to over-approximate
# attainable sets of a particular cell based on the state reach by its center.
#
# For this reachability problem, the abstraction controller is built by solving a fixed-point equation which consists in computing the pre-image
# of the target set.

# First, let us import [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl) and [Plots](https://github.com/JuliaPlots/Plots.jl).
using StaticArrays, Plots, Revise

# At this point, we import Dionysos and JuMP.
using Dionysos, JuMP

# ### Definition of the problem
#include(joinpath(@__DIR__, "../../../../src/MOi_wrapper.jl"))

# We define the problem using JuMP as follows.
# We first create a JuMP model:
model = Model(Dionysos.Optimizer)

# Define the prediction horizon
N = 5
discretization_step = 0.1

# Define the state variables: x1(t), x2(t), x3(t) for t = 1, ..., N
x_low = [-3.5, -2.6, -pi]
x_upp = -x_low
x_initial = [0.0, 0.0, 0.0]  # Initial state
#x_initial = [1.0, -1.7, 0.0]
@variable(model, x_low[i] <= x[i = 1:3] <= x_upp[i])#, start = x_initial[i])

# Define the control variables: u1(t), u2(t) for t = 1, ..., N-1
@variable(model, -1 <= u[1:2] <= 1)

# Define the dynamics

@constraint(model, ∂(x[1]) == x[1] + u[1] * cos(x[3]))
@constraint(model, ∂(x[2]) == x[2] + u[1] * sin(x[3]))
#@constraint(model, ∂(x[3) == mod(x[3] + u[2, t], 2 * pi))
@constraint(model, ∂(x[3]) == x[3] + u[2])


# Define the initial and target sets

x_target = [0.5, 0.5, -pi]  # Target state
#x_target = [sqrt(32)/3, sqrt(20)/3, -pi]

## constraint on the initial state+discretization_step
#@constraint(model, start(x[1]) in MOI.Interval(x_initial[1] - discretization_step, x_initial[1] + discretization_step))
#@constraint(model, x[:, 1] in MOI.Interval.(x_initial .- discretization_step, x_initial .+ discretization_step))

## constraint on the target state+0.2
#@constraint(model, x[:, N] == x_target)
#FIXME: final should be able to recognized that we are walking with an horizon N and set the final corectely.
# writing final(x[?, N]) in Interval is redundant. We should be able to write final(x[?]) in Interval

#FIXME: if a julia variable is provided we should be able to handle it or to give more clear feedback
#@constraint(model, final(x[1,N]) in MOI.Interval(x_target[1] - discretization_step, x_target[1] + discretization_step))
#@constraint(model, final(x[2,N]) in MOI.Interval(x_target[2] - discretization_step, x_target[2] + discretization_step))
#@constraint(model, final(x[3,N]) in MOI.Interval(x_target[3], x_target[3]))

#@constraint(model, start(x[1]) in MOI.Interval(x_initial[1] - discretization_step, x_initial[1] + discretization_step))
#@constraint(model, start(x[2]) in MOI.Interval(x_initial[2] - discretization_step, x_initial[2] + discretization_step))
#@constraint(model, start(x[3]) in MOI.Interval(x_initial[3], x_initial[3]))

@constraint(model, start(x[1]) in MOI.Interval(-0.2, 0.2))
@constraint(model, start(x[2]) in MOI.Interval(-0.2, 0.2))
@constraint(model, start(x[3]) in MOI.Interval(-0.2, 0.2))

#@constraint(model, start(x[1]) in MOI.Interval(0.8, 1.2))
#@constraint(model, start(x[2]) in MOI.Interval(-1.9, -1.5))
#@constraint(model, start(x[3]) in MOI.Interval(-0.2, 0.2))

@constraint(model, final(x[1]) in MOI.Interval(0.3, 0.7))
@constraint(model, final(x[2]) in MOI.Interval(0.3, 0.7))
@constraint(model, final(x[3]) in MOI.Interval(-3.14, 3.14))

# Obstacle boundaries (provided)
function extract_rectangles(matrix)
    if isempty(matrix)
        return []
    end

    n, m = size(matrix)
    tlx, tly, brx, bry = Int[], Int[], Int[], Int[]

    # Build histogram heights
    for i in 1:n
        j = 1
        while j <= m
            if matrix[i, j] == 1
                j += 1
                continue
            end
            push!(tlx, j)
            push!(tly, i)
            while j <= m && matrix[i, j] == 0
                j += 1
            end
            push!(brx, j - 1)
            push!(bry, i)
        end
    end

    return zip(tlx, tly, brx, bry)
end

function get_obstacles(lb,ub, h)
    #lb_x1 = -3.5, ub_x1 = 3.5, lb_x2 = -2.6, ub_x2 = 2.6, h = 0.1
    lb_x1, lb_x2, lb_x3 = lb
    ub_x1, ub_x2, ub_x3 = ub
    

    # Define the obstacles
    x1 = range(lb_x1; stop = ub_x1, step = h)
    x2 = range(lb_x2; stop = ub_x2, step = h)
    steps1, steps2 = length(x1), length(x2)

    X1 = x1' .* ones(steps2)
    X2 = ones(steps1)' .* x2

    Z1 = (X1 .^ 2 .- X2 .^ 2) .<= 4
    Z2 = (4 .* X2 .^ 2 .- X1 .^ 2) .<= 16

    # Find the upper and lower bounds of X1 and X2 for the obstacle 
    grid = Z1 .& Z2

    return [
        MOI.HyperRectangle(
            [x1[x1lb], x2[x2lb]],
            [x1[x1ub], x2[x2ub]]
        ) for (x1lb, x2lb, x1ub, x2ub) in extract_rectangles(grid)
    ]
end

obstacles = get_obstacles(x_low, x_upp, discretization_step)
# Function to add rectangular obstacle avoidance constraints

for obstacle in obstacles
    @constraint(model, x[1:2] ∉ obstacle)
end

# ### Definition of the abstraction

# Definition of the grid of the state-space on which the abstraction is based (origin `x0` and state-space discretization `h`):

# We define the growth bound function of $f$:
function jacobian_bound(u)
    α = abs(-u[1] * sin(u[2]))
    β = abs(u[1] * cos(u[2]))
    #@show β
    return StaticArrays.SMatrix{3, 3}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, α, β, 0.0)
end
set_attribute(model, "jacobian_bound", jacobian_bound)

function growth_bound(r, u, _)
    β = u[1] * r[3]
    return StaticArrays.SVector{3}(β, β, 0.0)
end
set_attribute(model, "growthbound_map", growth_bound)

function sys_inv(x, u, _)
    return StaticArrays.SVector{3}(
        x[1] - u[1] * cos((x[3] - u[2]) % (2 * π)),
        x[2] - u[1] * sin((x[3] - u[2]) % (2 * π)),
        (x[3] - u[2]) % (2 * π), # to check
    )
end
set_attribute(model, "sys_inv_map", sys_inv)

set_attribute(model, "time_step", 1.0)

x0 = SVector(0.0, 0.0, 0.0);
h = SVector(0.1, 0.1, 0.2);
#h = SVector(discretization_step, discretization_step, discretization_step);
set_attribute(model, "state_grid", Dionysos.Domain.GridFree(x0, h))

# Definition of the grid of the input-space on which the abstraction is based (origin `u0` and input-space discretization `h`):
u0 = SVector(1.1, 0.0);
h = SVector(0.3, 0.3);
set_attribute(model, "input_grid", Dionysos.Domain.GridFree(u0, h))

optimize!(model)

# Get the results
abstract_system = get_attribute(model, "abstract_system");
abstract_problem = get_attribute(model, "abstract_problem");
abstract_controller = get_attribute(model, "abstract_controller");
concrete_controller = get_attribute(model, "concrete_controller")
concrete_problem = get_attribute(model, "concrete_problem");
concrete_system = concrete_problem.system

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

control_trajectory = Dionysos.System.get_closed_loop_trajectory(
    get_attribute(model, "discretized_system"),
    concrete_controller,
    x_initial,
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
    Dionysos.Symbolic.get_domain_from_symbols(
        abstract_system,
        abstract_problem.initial_set,
    );
    color = :green,
);
plot!(
    Dionysos.Symbolic.get_domain_from_symbols(abstract_system, abstract_problem.target_set);
    color = :red,
);

# We display the concrete trajectory
plot!(control_trajectory; ms = 0.5)

# ### References
# 1. G. Reissig, A. Weber and M. Rungger, "Feedback Refinement Relations for the Synthesis of Symbolic Controllers," in IEEE Transactions on Automatic Control, vol. 62, no. 4, pp. 1781-1796.
# 2. K. J. Aström and R. M. Murray, Feedback systems. Princeton University Press, Princeton, NJ, 2008.
