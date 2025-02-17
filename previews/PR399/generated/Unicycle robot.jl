using StaticArrays, Plots

using Dionysos, JuMP

model = Model(Dionysos.Optimizer)

hx = 0.1

x_low = [-3.5, -2.6, -pi]
x_upp = -x_low
@variable(model, x_low[i] <= x[i = 1:3] <= x_upp[i])

@variable(model, -1 <= u[1:2] <= 1)

@constraint(model, Δ(x[1]) == x[1] + u[1] * cos(x[3]))
@constraint(model, Δ(x[2]) == x[2] + u[1] * sin(x[3]))
@constraint(model, Δ(x[3]) == x[3] + u[2])

x_initial = [1.0, -1.7, 0.0]
x_target = [sqrt(32) / 3, sqrt(20) / 3, -pi]

@constraint(model, start(x[1]) in MOI.Interval(x_initial[1] - hx, x_initial[1] + hx))
@constraint(model, start(x[2]) in MOI.Interval(x_initial[2] - hx, x_initial[2] + hx))
@constraint(model, start(x[3]) in MOI.Interval(x_initial[3] - hx, x_initial[3] + hx))

@constraint(model, final(x[1]) in MOI.Interval(x_target[1] - hx, x_target[1] + hx))
@constraint(model, final(x[2]) in MOI.Interval(x_target[2] - hx, x_target[2] + hx))
@constraint(model, final(x[3]) in MOI.Interval{Float64}(-pi, pi))

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

function get_obstacles(lb, ub, h)
    # lb_x1 = -3.5, ub_x1 = 3.5, lb_x2 = -2.6, ub_x2 = 2.6, h = 0.1
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
        MOI.HyperRectangle([x1[x1lb], x2[x2lb]], [x1[x1ub], x2[x2ub]]) for
        (x1lb, x2lb, x1ub, x2ub) in extract_rectangles(grid)
    ]
end

obstacles = get_obstacles(x_low, x_upp, hx)

for obstacle in obstacles
    @constraint(model, x[1:2] ∉ obstacle)
end

function growth_bound(r, u)
    β = u[1] * r[3]
    return StaticArrays.SVector{3}(β, β, 0.0)
end
set_attribute(model, "growthbound_map", growth_bound)

x0 = SVector(0.0, 0.0, 0.0);
h = SVector(hx, hx, 0.2);
set_attribute(model, "state_grid", Dionysos.Domain.GridFree(x0, h))

u0 = SVector(1.1, 0.0);
h = SVector(0.3, 0.3);
set_attribute(model, "input_grid", Dionysos.Domain.GridFree(u0, h))

optimize!(model)

abstract_system = get_attribute(model, "abstract_system");
abstract_problem = get_attribute(model, "abstract_problem");
abstract_controller = get_attribute(model, "abstract_controller");
concrete_controller = get_attribute(model, "concrete_controller")
concrete_problem = get_attribute(model, "concrete_problem");
concrete_system = concrete_problem.system

nstep = 100
function reached(x)
    if x ∈ concrete_problem.target_set
        return true
    else
        return false
    end
end

control_trajectory = Dionysos.System.get_closed_loop_trajectory(
    get_attribute(model, "discrete_time_system"),
    concrete_controller,
    x_initial,
    nstep;
    stopping = reached,
)

using Plots

fig = plot(; aspect_ratio = :equal);

plot!(concrete_system.X; color = :yellow, opacity = 0.5);

plot!(abstract_system.Xdom; color = :blue, opacity = 0.5);

plot!(concrete_problem.initial_set; color = :green, opacity = 0.2);
plot!(concrete_problem.target_set; dims = [1, 2], color = :red, opacity = 0.2);

plot!(
    Dionysos.Symbolic.get_domain_from_states(abstract_system, abstract_problem.initial_set);
    color = :green,
);
plot!(
    Dionysos.Symbolic.get_domain_from_states(abstract_system, abstract_problem.target_set);
    color = :red,
);

plot!(control_trajectory; ms = 0.5)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
