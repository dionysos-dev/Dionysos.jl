using StaticArrays, Plots

using Dionysos, JuMP

model = Model(Dionysos.Optimizer)

hx = 0.05
l = 1.0
g = 9.81

x_low, x_upp = [-π, -10.0], [π + pi, 10.0]
@variable(model, x_low[i] <= x[i = 1:2] <= x_upp[i])
nothing #hide

@variable(model, -3.0 <= u <= 3.0)

@constraint(model, ∂(x[1]) == x[2])
@constraint(model, ∂(x[2]) == -(g / l) * sin(x[1]) + u)

x1_initial, x2_initial = (5.0 * pi / 180.0) .* [-1, 1], 0.5 .* [-1, 1]
x1_target, x2_target = pi .+ (5.0 * pi / 180.0) .* [-1, 1], 1.0 .* [-1, 1]

@constraint(model, start(x[1]) in MOI.Interval(x1_initial...))
@constraint(model, start(x[2]) in MOI.Interval(x2_initial...))

@constraint(model, final(x[1]) in MOI.Interval(x1_target...))
@constraint(model, final(x[2]) in MOI.Interval(x2_target...))

function jacobian_bound(u)
    return SMatrix{2, 2}(0.0, 1.0, (g / l), 0)
end
set_attribute(model, "jacobian_bound", jacobian_bound)

set_attribute(model, "time_step", 0.1)

x0 = SVector(0.0, 0.0);
h = SVector(hx, hx);
set_attribute(model, "state_grid", Dionysos.Domain.GridFree(x0, h))

u0 = SVector(0.0);
h = SVector(0.3);
set_attribute(model, "input_grid", Dionysos.Domain.GridFree(u0, h))

optimize!(model)

abstract_system = get_attribute(model, "abstract_system");
abstract_problem = get_attribute(model, "abstract_problem");
abstract_controller = get_attribute(model, "abstract_controller");
concrete_controller = get_attribute(model, "concrete_controller")
concrete_problem = get_attribute(model, "concrete_problem");
concrete_system = concrete_problem.system
nothing #hide

nstep = 100
function reached(x)
    if x ∈ concrete_problem.target_set
        return true
    else
        return false
    end
end
x0 = SVector(Dionysos.Utils.sample(concrete_problem.initial_set)...)
control_trajectory = Dionysos.System.get_closed_loop_trajectory(
    get_attribute(model, "discrete_time_system"),
    concrete_controller,
    x0,
    nstep;
    stopping = reached,
)

using Plots

fig = plot(; aspect_ratio = :equal);

plot!(concrete_system.X; color = :grey, label = "");

plot!(
    Dionysos.Symbolic.get_domain_from_states(abstract_system, abstract_problem.initial_set);
    color = :green,
    opacity = 1.0,
);
plot!(
    Dionysos.Symbolic.get_domain_from_states(abstract_system, abstract_problem.target_set);
    color = :red,
    opacity = 1.0,
);
plot!(control_trajectory; markersize = 1, arrows = false)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
