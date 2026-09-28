using StaticArrays, JuMP, Plots
import LazySets

using Symbolics, MathOptSymbolicAD

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const AB = DI.Optim.Abstraction;

m, f0, f1, f2 = 1650.0, 0.1, 5.0, 0.25
v_lead, v_desired, τ_h = 13.89, 24.0, 1.8
a_max = 0.3 * 9.81;

z_low, z_upp = 0.0, 100.0
v_low, v_upp = 0.0, 30.0

safe = LazySets.HPolytope([
    LazySets.HalfSpace(SVector(1.0, 0.0), z_upp),
    LazySets.HalfSpace(SVector(-1.0, 0.0), -z_low),
    LazySets.HalfSpace(SVector(0.0, 1.0), v_upp),
    LazySets.HalfSpace(SVector(0.0, -1.0), -v_low),
    LazySets.HalfSpace(SVector(-1.0, τ_h), 0.0),   # z ≥ τ_h v
]);

model = Model(Dionysos.Optimizer)
@variable(model, z_low <= z <= z_upp)
@variable(model, v_low <= v <= v_upp)
@variable(model, -a_max <= a <= a_max)

@constraint(model, ∂(z) == v_lead - v)
@constraint(model, ∂(v) == a - (f0 + f1 * v + f2 * v^2) / m)

initial = LazySets.Hyperrectangle(; low = SVector(59.0, 19.8), high = SVector(61.0, 20.2))

@constraint(model, [z, v] in Start(initial))
@constraint(model, [z, v] in Always(safe));

jacobian_bound = u -> SMatrix{2, 2}(0.0, 0.0, 1.0, -(f1 + 2 * f2 * v_low) / m);

set_attribute(model, "state_grid", MP.GridFree(SVector(0.0, 0.0), SVector(1.0, 0.2)))
set_attribute(model, "input_grid", MP.GridFree(SVector(0.0), SVector(0.4)))
set_attribute(model, "jacobian_bound", jacobian_bound)
set_attribute(model, "time_step", 0.5)
set_attribute(model, "approx_mode", AB.UniformGridAbstraction.GROWTH)
set_attribute(model, "print_level", 0)

optimize!(model);

termination_status(model)

invariant_set = get_attribute(model, "invariant_set")
abstract_system = get_attribute(model, "abstract_system")
XMapping = DI.Symbolic.get_state_mapping(abstract_system);

function min_safe_gap(v; a_brake = a_max)
    v_star = v_lead + τ_h * a_brake
    return τ_h * v + max(v - v_star, 0.0)^2 / (2 * a_brake)
end

a_outer = a_max + (f0 + f1 * v_upp + f2 * v_upp^2) / m
speeds = range(v_low, v_upp; length = 200);

fig = plot(; xlabel = "gap [m]", ylabel = "ego speed [m/s]", legend = :bottomright)
plot!(fig, get_attribute(model, "concrete_problem"))
plot!(
    fig,
    (invariant_set, XMapping);
    color = :blue,
    linecolor = :blue,
    label = "Invariant set",
)
plot!(
    fig,
    [min_safe_gap(s) for s in speeds],
    speeds;
    lw = 2,
    color = :black,
    label = "inner",
)
plot!(
    fig,
    [min_safe_gap(s; a_brake = a_outer) for s in speeds],
    speeds;
    lw = 2,
    ls = :dash,
    color = :black,
    label = "outer",
)

band(v_c; ε = 0.5, margin = 5.0, gap_high = z_upp) = LazySets.Hyperrectangle(;
    low = SVector(τ_h * (v_c + ε) + margin, v_c - ε),
    high = SVector(gap_high, v_c + ε),
)

target = UT.set_union([band(v_desired), band(v_lead; gap_high = 45.0)])

acc = Model(Dionysos.Optimizer)
@variable(acc, z_low <= za <= z_upp)
@variable(acc, v_low <= va <= v_upp)
@variable(acc, -a_max <= aa <= a_max)

@constraint(acc, ∂(za) == v_lead - va)
@constraint(acc, ∂(va) == aa - (f0 + f1 * va + f2 * va^2) / m)

far_behind =
    LazySets.Hyperrectangle(; low = SVector(89.0, 12.8), high = SVector(91.0, 13.2))

@constraint(acc, [za, va] in Start(far_behind))
@constraint(acc, [za, va] in Always(safe))
@constraint(acc, [za, va] in EventuallyAlways(target; stay_on_first_entry = true))

for (k, val) in (
    ("state_grid", MP.GridFree(SVector(0.0, 0.0), SVector(1.0, 0.2))),
    ("input_grid", MP.GridFree(SVector(0.0), SVector(0.4))),
    ("jacobian_bound", jacobian_bound),
    ("time_step", 0.5),
    ("approx_mode", AB.UniformGridAbstraction.GROWTH),
    ("print_level", 0),
)
    set_attribute(acc, k, val)
end

optimize!(acc);

termination_status(acc)

winning_set = get_attribute(acc, "winning_set")
acc_mapping = DI.Symbolic.get_state_mapping(get_attribute(acc, "abstract_system"));

trajectory =
    Dionysos.simulate(acc, SVector(90.0, 13.0); nsteps = 160, stopping = _ -> false);

let xs = collect(ST.states(trajectory))
    (gap_start = xs[1][1], gap_end = xs[end][1], top_speed = maximum(x[2] for x in xs))
end

fig = plot(; xlabel = "gap [m]", ylabel = "ego speed [m/s]", legend = :bottomright)
plot!(fig, get_attribute(acc, "concrete_problem"))
plot!(
    fig,
    (winning_set, acc_mapping);
    color = :yellow,
    linecolor = :yellow,
    label = "Winning set",
)
plot!(fig, trajectory; ms = 2.0, color = :blue)

include(
    joinpath(
        dirname(dirname(pathof(Dionysos))),
        "problems",
        "AdaptiveCruiseControl",
        "adaptive_cruise_control.jl",
    ),
);

anim = Dionysos.animate_trajectory_dashboard(
    AdaptiveCruiseControl.system_plot!(),
    trajectory;
    xdims = (1, 2),
    udims = (1,),
    Δt = 0.5,
    xlabel_state = "gap [m]",
    ylabel_state = "ego speed [m/s]",
    xlabel_input = "time [s]",
    ylabel_input = "acceleration [m/s²]",
);
gif(anim; fps = 4)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
