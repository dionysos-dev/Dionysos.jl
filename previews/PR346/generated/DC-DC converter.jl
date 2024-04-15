using StaticArrays, Plots

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "dc_dc.jl"))

concrete_problem = DCDC.problem(; approx_mode = "growth")
concrete_system = concrete_problem.system

x0 = SVector(0.0, 0.0)
hx = SVector(2.0 / 4.0e3, 2.0 / 4.0e3)
state_grid = DO.GridFree(x0, hx)
u0 = SVector(1)
hu = SVector(1)
input_grid = DO.GridFree(u0, hu)

using JuMP
optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.optimize!(optimizer)

abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

nstep = 300
x0 = SVector(1.2, 5.6)
control_trajectory =
    ST.get_closed_loop_trajectory(concrete_system.f, concrete_controller, x0, nstep)

fig = plot(; aspect_ratio = :equal);
plot!(concrete_system.X);
plot!(control_trajectory)

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "dc_dc.jl"))

concrete_problem = DCDC.problem(; approx_mode = "δ-GAS Lyapunov")
concrete_system = concrete_problem.system

origin = SVector(0.0, 0.0)
η = (2 / 4.0) * 10^(-3)

ϵ = 0.1 * 0.01
P = SMatrix{2, 2}(1.0224, 0.0084, 0.0084, 1.0031)
state_grid = DO.GridEllipsoidalRectangular(origin, SVector(η, η), P / ϵ, concrete_system.X)

u0 = SVector(1)
hu = SVector(1)
input_grid = DO.GridFree(u0, hu)

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("δGAS"), true)
MOI.optimize!(optimizer)

abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

nstep = 300
x0 = SVector(1.2, 5.6)
control_trajectory =
    ST.get_closed_loop_trajectory(concrete_system.f, concrete_controller, x0, nstep)

fig = plot(; aspect_ratio = :equal);
plot!(concrete_system.X);
plot!(control_trajectory)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
