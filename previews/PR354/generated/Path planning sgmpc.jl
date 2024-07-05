using StaticArrays, Plots

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const PR = DI.Problem
const OP = DI.Optim
const AB = OP.Abstraction

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "path_planning_sgmpc.jl"))

initial = SVector(1.0, -1.7, 0.0)
#initial = SVector(0.0, 0.0, 0.0)
#target = SVector(0.5, 0.5, -pi)
target = SVector(-0.5, 0.5, -pi)
#target = SVector(2.6, 2.0, -pi)
#target = SVector(-2.6, 2.0, -pi)
#target = SVector(sqrt(32.0 / 3.0), sqrt(20.0 / 3.0), -pi)

concrete_problem = PathPlanningSgMPC.problem(; sgmpc = true, initial = initial, target = target)
concrete_system = concrete_problem.system;

x0 = SVector(0.0, 0.0, 0.0);
h = SVector(0.1, 0.1, 0.2);
state_grid = DO.GridFree(x0, h);

u0 = SVector(1.1, 0.0);
h = SVector(0.3, 0.3);
input_grid = DO.GridFree(u0, h);

using JuMP
optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.optimize!(optimizer)

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

nstep = 100
function reached(x)
    if x ∈ concrete_problem.target_set
        return true
    else
        return false
    end
end
x0 = initial #SVector(1.1, -1.6, 0.0)
control_trajectory = ST.get_closed_loop_trajectory(
    concrete_system.f,
    concrete_controller,
    x0,
    nstep;
    stopping = reached,
)

fig = plot(; aspect_ratio = :equal);

plot!(concrete_system.X; color = :yellow, opacity = 0.5);

plot!(abstract_system.Xdom; color = :blue, opacity = 0.5);

plot!(concrete_problem.initial_set; color = :green, opacity = 0.2);
plot!(concrete_problem.target_set; dims = [1, 2], color = :red, opacity = 0.2);

plot!(
    SY.get_domain_from_symbols(abstract_system, abstract_problem.initial_set);
    color = :green,
);
plot!(
    SY.get_domain_from_symbols(abstract_system, abstract_problem.target_set);
    color = :red,
);

display(control_trajectory)

plot!(control_trajectory; ms = 0.5)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
