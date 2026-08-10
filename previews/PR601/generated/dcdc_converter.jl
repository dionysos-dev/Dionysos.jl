using StaticArrays, JuMP, Plots

using Symbolics, MathOptSymbolicAD

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const AB = DI.Optim.Abstraction;

vs, xL, xC, r0, rL, rC = 1.0, 3.0, 70.0, 1.0, 0.05, 0.005

A1 = SMatrix{2, 2}(-rL / xL, 0.0, 0.0, -1.0 / (xC * (r0 + rC)))
A2 = SMatrix{2, 2}(
    -(rL + r0 * rC / (r0 + rC)) / xL,
    5.0 * r0 / ((r0 + rC) * xC),
    -r0 / ((r0 + rC) * xL * 5.0),
    -1.0 / (xC * (r0 + rC)),
)
b = SVector(vs / xL, 0.0)

dynamic = (x, u) -> u[1] == 1 ? A1 * x + b : A2 * x + b;

A2_abs = SMatrix{2, 2}(
    -(rL + r0 * rC / (r0 + rC)) / xL,
    5.0 * r0 / ((r0 + rC) * xC),
    r0 / ((r0 + rC) * xL * 5.0),
    -1.0 / (xC * (r0 + rC)),
)

jacobian_bound = u -> u[1] == 1 ? A1 : A2_abs;

x_low, x_upp = [1.15, 5.45], [1.55, 5.85]
safe = UT.box(SVector(x_low...), SVector(x_upp...))
initial = UT.box(SVector(1.19, 5.59), SVector(1.21, 5.61));

model = direct_model(Dionysos.Optimizer());
@variable(model, x_low[i] <= x[i = 1:2] <= x_upp[i])
@variable(model, 1 <= u <= 2);

set_role!(x, Dionysos.STATE)
set_attribute(model, "dynamics", dynamic)
set_attribute(model, "time_domain", Dionysos.CONTINUOUS);

@constraint(model, x in Start(initial))
@constraint(model, x in Always(safe));

hx = SVector(2.0 / 4.0e3, 2.0 / 4.0e3)
set_attribute(model, "state_grid", MP.GridFree(SVector(0.0, 0.0), hx))
set_attribute(model, "input_grid", MP.GridFree(SVector(1), SVector(1)))
set_attribute(model, "jacobian_bound", jacobian_bound)
set_attribute(model, "time_step", 0.5)
set_attribute(model, "approx_mode", AB.UniformGridAbstraction.GROWTH)
set_attribute(model, "print_level", 0);

optimize!(model);

termination_status(model)

concrete_problem = get_attribute(model, "concrete_problem");
concrete_controller = get_attribute(model, "concrete_controller");
invariant_set = get_attribute(model, "invariant_set");

trajectory = Dionysos.simulate(model, SVector(1.2, 5.6); nsteps = 150);

fig = plot(; aspect_ratio = :equal)
plot!(concrete_problem)
plot!(trajectory; ms = 2.0, color = :blue)

include(
    joinpath(dirname(dirname(pathof(Dionysos))), "problems", "DCDC", "dcdc_converter.jl"),
);

anim = Dionysos.animate_trajectory_dashboard(
    DCDC.system_plot!(),
    trajectory;
    xdims = (1, 2),
    udims = (1,),
    Δt = 0.5,
    frame_step = 2,
    xlabel_state = "iL [A]",
    ylabel_state = "vC [V]",
    ylabel_input = "switch position",
);
gif(anim; fps = 8)

origin = SVector(0.0, 0.0)
η = (2 / 4.0) * 10^(-3)
ϵ = 0.1 * 0.01
P = SMatrix{2, 2}(1.0224, 0.0084, 0.0084, 1.0031)
state_grid = MP.GridEllipsoidalRectangular(origin, SVector(η, η), P / ϵ)

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("input_grid"),
    MP.GridFree(SVector(1), SVector(1)),
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), jacobian_bound)
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.5)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.CENTER_SIMULATION,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
MOI.optimize!(optimizer);

MOI.get(optimizer, MOI.TerminationStatus())

abstract_system2 = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"));
controller2 = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"));

traj2 = ST.get_closed_loop_trajectory(
    MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system")),
    controller2,
    SVector(1.2, 5.6),
    300,
);


fig2 = plot(; aspect_ratio = :equal)
plot!(concrete_problem)
plot!(traj2; ms = 2.0, color = :blue)

hybrid = Model(Dionysos.Optimizer);
@variable(hybrid, 1.15 <= xh[i = 1:2] <= 1.55)
set_lower_bound(xh[2], 5.45)
set_upper_bound(xh[2], 5.85)

@mode(hybrid, closed)
@mode(hybrid, opened);

@constraint(closed, ∂(xh[1]) == A1[1, 1] * xh[1] + A1[1, 2] * xh[2] + b[1])
@constraint(closed, ∂(xh[2]) == A1[2, 1] * xh[1] + A1[2, 2] * xh[2] + b[2])
@constraint(opened, ∂(xh[1]) == A2[1, 1] * xh[1] + A2[1, 2] * xh[2] + b[1])
@constraint(opened, ∂(xh[2]) == A2[2, 1] * xh[1] + A2[2, 2] * xh[2] + b[2]);

add_transition!(hybrid, closed => opened) do t
    return @constraint(t, xh in Guard(safe))
end
add_transition!(hybrid, opened => closed) do t
    return @constraint(t, xh in Guard(safe))
end

@constraint(closed, xh in Always(safe))
@constraint(opened, xh in Always(safe));

for (md, Ab) in ((closed, A1), (opened, A2_abs))
    set_attribute(md, "state_grid", MP.GridFree(SVector(0.0, 0.0), 4 .* hx))
    set_attribute(md, "time_step", 0.5)
    set_attribute(md, "approx_mode", AB.UniformGridAbstraction.GROWTH)
    set_attribute(md, "jacobian_bound", u -> Ab)
    set_attribute(md, "print_level", 0)
end

set_attribute(hybrid, "print_level", 0)

optimize!(hybrid);

termination_status(hybrid)

hybrid_traj = Dionysos.simulate(hybrid, (SVector(1.2, 5.6), 1); nsteps = 200);

fig3 = plot(; aspect_ratio = :equal)
plot!(concrete_problem)
plot!(hybrid_traj; ms = 2.0, color = :blue)

hybrid_switch = ST.Trajectory(
    ST.states(hybrid_traj);
    inputs = [SVector(Float64(m)) for m in ST.modes(hybrid_traj)[1:(end - 1)]],
);

anim_hybrid = Dionysos.animate_trajectory_dashboard(
    DCDC.system_plot!(),
    hybrid_switch;
    xdims = (1, 2),
    udims = (1,),
    Δt = 0.5,
    frame_step = 2,
    xlabel_state = "iL [A]",
    ylabel_state = "vC [V]",
    ylabel_input = "switch position",
);
gif(anim_hybrid; fps = 8)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
