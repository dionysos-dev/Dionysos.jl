using StaticArrays, JuMP, Plots
using Symbolics, MathOptSymbolicAD

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const AB = DI.Optim.Abstraction;

Ta, α, β = 18.0, 0.1, 2.0

model = Model(Dionysos.Optimizer);

@variable(model, 17.0 <= T <= 25.0, start = 18.0)
@variable(model, 0.0 <= u <= 1.0);

@mode(model, off)
@mode(model, on);

@constraint(off, ∂(T) == -α * (T - Ta))
@constraint(on, ∂(T) == -α * (T - Ta) + β * u);

@constraint(off, u == 0.0)
@constraint(on, 0.2 <= u <= 1.0);

add_transition!(model, off => on) do t
    return @constraint(t, T <= 19.0)
end

add_transition!(model, on => off) do t
    return @constraint(t, T >= 21.0)
end;

comfortable = UT.box(SVector(20.5), SVector(23.0))

@constraint(off, [T] in Final(comfortable))
@constraint(on, [T] in Final(comfortable));

for m in (off, on)
    set_attribute(m, "state_grid", MP.GridFree(SVector(0.0), SVector(0.2)))
    set_attribute(m, "input_grid", MP.GridFree(SVector(0.0), SVector(0.2)))
    set_attribute(m, "time_step", 0.5)
    set_attribute(m, "approx_mode", AB.UniformGridAbstraction.GROWTH)
    set_attribute(m, "jacobian_bound", u -> SMatrix{1, 1}(-α))
    set_attribute(m, "print_level", 0)
end

optimize!(model);

termination_status(model)

concrete_problem = get_attribute(model, "concrete_problem");

concrete_problem.system

trajectory = Dionysos.simulate(model, (SVector(18.0), 1); nsteps = 60);

include(
    joinpath(
        dirname(dirname(pathof(Dionysos))),
        "problems",
        "Thermostat",
        "thermostat_hybrid_system.jl",
    ),
);

anim = Dionysos.animate_trajectory_dashboard(
    ThermostatHybridSystem.system_plot!(; problem = concrete_problem),
    trajectory;
    xdims = (1,),
    udims = (1,),
    Δt = 0.5,
    xlabel_state = "time [s]",
    ylabel_state = "T [°C]",
    ylabel_input = "heating power",
);
gif(anim; fps = 6)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
