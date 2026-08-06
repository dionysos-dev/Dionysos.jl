using StaticArrays, JuMP, Plots
using Symbolics, MathOptSymbolicAD

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping

l, g = 1.0, 9.81

model = Model(Dionysos.Optimizer);

@variable(model, -π <= x1 <= π)
@variable(model, -10.0 <= x2 <= 10.0)
@variable(model, -3.0 <= u <= 3.0);

@constraint(model, ∂(x1) == x2)
@constraint(model, ∂(x2) == -(g / l) * sin(x1) + u)

@constraint(model, start(x1) in MOI.Interval(-0.09, 0.09))
@constraint(model, start(x2) in MOI.Interval(-0.5, 0.5))

@constraint(model, final(x1) in MOI.Interval(π - 15π / 180, π + 15π / 180))
@constraint(model, final(x2) in MOI.Interval(-1.0, 1.0))

@constraint(model, [x1] ∉ MOI.HyperRectangle([-π + 16π / 180], [-π + 38π / 180]));

set_attribute(model, "jacobian_bound", u -> SMatrix{2, 2}(0.0, 1.0, g / l, 0.0))
set_attribute(model, "time_step", 0.1)

h = SVector(3π / 180, 0.05)
set_attribute(model, "state_grid", MP.GridFree(SVector(-π + h[1] / 2, 0.0), h))
set_attribute(model, "input_grid", MP.GridFree(SVector(0.0), SVector(0.3)))

set_attribute(model, "use_periodic_mapping", true)
set_attribute(model, "periodic_dims", SVector(1))
set_attribute(model, "periodic_periods", SVector(2π))
set_attribute(model, "periodic_start", SVector(-π))
set_attribute(model, "print_level", 0);

optimize!(model);

termination_status(model)

concrete_problem = get_attribute(model, "concrete_problem");
concrete_controller = get_attribute(model, "concrete_controller");
abstract_system = get_attribute(model, "abstract_system");
abstract_value_function = get_attribute(model, "abstract_value_function");

trajectory = Dionysos.simulate(model, SVector(0.0, 0.0); nsteps = 200)

fig = plot(; aspect_ratio = :equal)
plot!(concrete_problem)
plot!(trajectory; ms = 1.0, arrows = false, color = :blue)

include(
    joinpath(
        dirname(dirname(pathof(Dionysos))),
        "problems",
        "Pendulum",
        "simple_pendulum.jl",
    ),
)

anim = Dionysos.animate_trajectory_dashboard(
    SimplePendulum.system_plot!(),
    trajectory;
    xdims = (1, 2),
    udims = (1,),
    Δt = 0.1,
    frame_step = 2,
    xlabel_state = "θ [rad]",
    ylabel_state = "ω [rad/s]",
    ylabel_input = "τ [N·m]",
)

gif(anim; fps = 4)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
