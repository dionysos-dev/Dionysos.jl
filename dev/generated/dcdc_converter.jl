using StaticArrays, JuMP, Plots

using Symbolics, MathOptSymbolicAD

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const AB = DI.Optim.Abstraction

include(
    joinpath(dirname(dirname(pathof(Dionysos))), "problems", "DCDC", "dcdc_converter.jl"),
)

x_low, x_upp = [1.15, 5.45], [1.55, 5.85]
safe = UT.box(SVector(x_low...), SVector(x_upp...))
initial = UT.box(SVector(1.19, 5.59), SVector(1.21, 5.61))

model = direct_model(Dionysos.Optimizer());
@variable(model, x_low[i] <= x[i = 1:2] <= x_upp[i])
@variable(model, 1 <= u <= 2);

set_role!(x, Dionysos.STATE)
set_attribute(model, "dynamics", DCDC.dynamic())
set_attribute(model, "time_domain", Dionysos.CONTINUOUS);

@constraint(model, x in Start(initial))
@constraint(model, x in Always(safe));

hx = SVector(2.0 / 4.0e3, 2.0 / 4.0e3)
set_attribute(model, "state_grid", MP.GridFree(SVector(0.0, 0.0), hx))
set_attribute(model, "input_grid", MP.GridFree(SVector(1), SVector(1)))
set_attribute(model, "jacobian_bound", DCDC.jacobian_bound())
set_attribute(model, "time_step", 0.5)
set_attribute(model, "approx_mode", AB.UniformGridAbstraction.GROWTH)
set_attribute(model, "print_level", 0);

optimize!(model);

termination_status(model)

concrete_problem = get_attribute(model, "concrete_problem");
concrete_controller = get_attribute(model, "concrete_controller");
invariant_set = get_attribute(model, "invariant_set");

trajectory = Dionysos.simulate(model, SVector(1.2, 5.6); nsteps = 300)

fig = plot(; aspect_ratio = :equal)
plot!(concrete_problem)
plot!(trajectory; with_arrows = false, ms = 2.0, color = :blue)

anim = Dionysos.animate_trajectory_dashboard(
    DCDC.system_plot!(),
    trajectory;
    xdims = (1, 2),
    udims = (1,),
    Δt = 0.5,
    frame_step = 6,
    xlabel_state = "iL [A]",
    ylabel_state = "vC [V]",
    ylabel_input = "switch position",
)
gif(anim; fps = 5)

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
MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), DCDC.jacobian_bound())
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
plot!(traj2; with_arrows = false, ms = 2.0, color = :blue)

hybrid = Model(Dionysos.Optimizer);
@variable(hybrid, 1.15 <= xh[i = 1:2] <= 1.55)
set_lower_bound(xh[2], 5.45)
set_upper_bound(xh[2], 5.85)

@mode(hybrid, closed)
@mode(hybrid, opened);

A1, A2 = DCDC.A1(), DCDC.A2()
bvec = SVector(DCDC.Params().vs / DCDC.Params().xL, 0.0)

@constraint(closed, ∂(xh[1]) == A1[1, 1] * xh[1] + A1[1, 2] * xh[2] + bvec[1])
@constraint(closed, ∂(xh[2]) == A1[2, 1] * xh[1] + A1[2, 2] * xh[2] + bvec[2])
@constraint(opened, ∂(xh[1]) == A2[1, 1] * xh[1] + A2[1, 2] * xh[2] + bvec[1])
@constraint(opened, ∂(xh[2]) == A2[2, 1] * xh[1] + A2[2, 2] * xh[2] + bvec[2]);

add_transition!(hybrid, closed => opened) do t
    return @constraint(t, xh in Guard(safe))
end
add_transition!(hybrid, opened => closed) do t
    return @constraint(t, xh in Guard(safe))
end

@constraint(closed, xh in Always(safe))
@constraint(opened, xh in Always(safe));

for md in (closed, opened)
    set_attribute(md, "state_grid", MP.GridFree(SVector(0.0, 0.0), 4 .* hx))
    set_attribute(md, "time_step", 0.5)
    set_attribute(md, "approx_mode", AB.UniformGridAbstraction.GROWTH)
    set_attribute(md, "jacobian_bound", u -> DCDC.jacobian_bound()(SVector(1)))
    set_attribute(md, "print_level", 0)
end

optimize!(hybrid);

termination_status(hybrid)

input_map = get_attribute(hybrid, "abstract_system").input_mapping;
(input_map.continuous_inputs, input_map.switching_inputs)

hybrid_traj = Dionysos.simulate(hybrid, (SVector(1.2, 5.6), 1); nsteps = 200)


sort(unique(ST.modes(hybrid_traj)))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
