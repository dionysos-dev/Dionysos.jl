using Test

using StaticArrays, Plots

using JuMP, Dionysos

k = 0.5 # Maximum heating rate (degrees per unit time)
α = 0.1 # Ambient cooling rate (degrees per unit time)
T_r = 20.0 # Reference temperature
δ = 1.0 # Tolerance

model = Model(Dionysos.Optimizer)

T_low, T_upp = 0.0, 30.0
T_start = 10.0
@variable(model, T_low <= T <= T_upp, start = T_start)

@variable(model, 0 <= u <= 1)

@variable(model, 1 <= mode <= 2, integer = true)

@constraint(model, mode == 2 => {∂(T) == k * u})
@constraint(model, mode == 1 => {∂(T) == -α})

@constraint(model, guard(mode == 2, mode == 1) => {T >= T_r - δ})

@constraint(model, guard(mode == 1, mode == 2) => {T <= T_r + δ})

@constraint(model, start(T) in MOI.Interval(9.5, 10.5))

@constraint(model, final(T) in MOI.Interval(19.0, 21.0))

set_attribute(model, "time_step", 0.1)

x0 = SVector(0.0);
h = SVector(0.1);
set_attribute(model, "state_grid", Dionysos.Domain.GridFree(x0, h))

u0 = SVector(0.0);
h = SVector(0.3);
set_attribute(model, "input_grid", Dionysos.Domain.GridFree(u0, h))

optimize!(model)

optimizer = get_attribute(model, "inner")

MOI.get(optimizer, MOI.SolveTimeSec())

termination = MOI.get(optimizer, MOI.TerminationStatus())

objective_value = MOI.get(optimizer, MOI.ObjectiveValue())

xu = MOI.get(optimizer, Dionysos.System.ContinuousTrajectoryAttribute());

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
