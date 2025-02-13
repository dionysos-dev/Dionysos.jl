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

concrete_problem = DCDC.problem()
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
MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), DCDC.jacobian_bound())
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.5)

MOI.optimize!(optimizer)

abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
abstraction_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstraction_construction_time_sec"))
println("Time to construct the abstraction: $(abstraction_time)")
abstract_problem_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem_time_sec"))
println("Time to solve the abstract problem: $(abstract_problem_time)")

invariant_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("invariant_set"))
uninvariant_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("uninvariant_set"))

nstep = 300
x0 = SVector(1.2, 5.6)
control_trajectory = ST.get_closed_loop_trajectory(
    MOI.get(optimizer, MOI.RawOptimizerAttribute("discretized_system")),
    concrete_controller,
    x0,
    nstep,
);

fig = plot(; aspect_ratio = :equal);
plot!(concrete_system.X; label = "");
plot!(control_trajectory)

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "dc_dc.jl"))

concrete_problem = DCDC.problem()
concrete_system = concrete_problem.system

origin = SVector(0.0, 0.0)
η = (2 / 4.0) * 10^(-3)

ϵ = 0.1 * 0.01
P = SMatrix{2, 2}(1.0224, 0.0084, 0.0084, 1.0031)
state_grid = DO.GridEllipsoidalRectangular(origin, SVector(η, η), P / ϵ)

u0 = SVector(1)
hu = SVector(1)
input_grid = DO.GridFree(u0, hu)

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), DCDC.jacobian_bound())
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    Dionysos.Optim.Abstraction.UniformGridAbstraction.DELTA_GAS,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("δGAS"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.5)
MOI.optimize!(optimizer)

abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

nstep = 300
x0 = SVector(1.2, 5.6)
control_trajectory = ST.get_closed_loop_trajectory(
    MOI.get(optimizer, MOI.RawOptimizerAttribute("discretized_system")),
    concrete_controller,
    x0,
    nstep,
)

fig = plot(; aspect_ratio = :equal);
plot!(concrete_system.X; label = "");
plot!(concrete_problem.initial_set; color = :red, label = "");
plot!(control_trajectory)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
