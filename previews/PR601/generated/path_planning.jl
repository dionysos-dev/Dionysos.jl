using StaticArrays, JuMP, Plots
using Symbolics, MathOptSymbolicAD

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const AB = DI.Optim.Abstraction;

model = Model(Dionysos.Optimizer);

x_low, x_upp = [0.0, 0.0, -pi - 0.4], [4.0, 10.0, pi + 0.4]
@variable(model, x_low[i] <= x[i = 1:3] <= x_upp[i], start = [0.4, 0.4, 0.0][i])
@variable(model, -1 <= u[1:2] <= 1);

@expression(model, α, atan(tan(u[2]) / 2))
@constraint(model, ∂(x[1]) == u[1] * cos(α + x[3]) * sec(α))
@constraint(model, ∂(x[2]) == u[1] * sin(α + x[3]) * sec(α))
@constraint(model, ∂(x[3]) == u[1] * tan(u[2]));

@constraint(model, final(x[1]) in MOI.Interval(3.0, 3.6))
@constraint(model, final(x[2]) in MOI.Interval(0.3, 0.8));

walls = [([1.0, 0.0], [1.2, 9.0]), ([2.2, 0.0], [2.4, 5.0]), ([2.2, 6.0], [2.4, 10.0])]
for (lo, hi) in walls
    @constraint(model, x[1:2] ∉ MOI.HyperRectangle(lo, hi))
end

function jacobian_bound(u)
    β = abs(u[1] / cos(atan(tan(u[2]) / 2)))
    return SMatrix{3, 3}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, β, β, 0.0)
end

set_attribute(model, "jacobian_bound", jacobian_bound)
set_attribute(model, "approx_mode", AB.UniformGridAbstraction.GROWTH)
set_attribute(model, "time_step", 0.3)
set_attribute(
    model,
    "state_grid",
    MP.GridFree(SVector(0.0, 0.0, 0.0), SVector(0.2, 0.2, 0.2)),
)
set_attribute(model, "input_grid", MP.GridFree(SVector(0.0, 0.0), SVector(0.3, 0.3)))
set_attribute(model, "print_level", 0);

optimize!(model);

termination_status(model)

abstract_system = get_attribute(model, "abstract_system");
abstract_value_function = get_attribute(model, "abstract_value_function");
concrete_problem = get_attribute(model, "concrete_problem");
concrete_controller = get_attribute(model, "concrete_controller");

trajectory = Dionysos.simulate(model, SVector(0.4, 0.4, 0.0); nsteps = 100);

fig = plot(; aspect_ratio = :equal)
plot!(concrete_problem.system.X; color = :grey, opacity = 1.0, label = "")
plot!(abstract_system; value_function = abstract_value_function)
plot!(
    UT.project_set(concrete_problem.initial_set, [1, 2]);
    color = :green,
    opacity = 0.4,
    label = "Initial set",
)
plot!(
    UT.project_set(concrete_problem.target_set, [1, 2]);
    color = :red,
    opacity = 0.5,
    label = "Target set",
)
plot!(trajectory; ms = 2.0, color = :blue)

include(
    joinpath(
        dirname(dirname(pathof(Dionysos))),
        "problems",
        "PathPlanning",
        "path_planning.jl",
    ),
);

obstacles = [UT.box(SVector(lo...), SVector(hi...)) for (lo, hi) in walls];

anim = Dionysos.animate_trajectory_dashboard(
    PathPlanning.system_plot!(; obstacles = obstacles, xlims = (0, 4), ylims = (0, 10)),
    trajectory;
    xdims = (1, 2),
    udims = (1, 2),
    Δt = 0.3,
    frame_step = 2,
    xlims_state = (0, 4),
    ylims_state = (0, 10),
    xlabel_state = "x₁",
    ylabel_state = "x₂",
);
gif(anim; fps = 5)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
