using StaticArrays, Plots

using Dionysos
const DI = Dionysos
const DO = DI.Domain
const CO = DI.Control
const UT = DI.Utils
const OP = DI.Optim
const AB = OP.Abstraction

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "PathPlanning.jl"))

problem = PathPlanning.problem(simple=true, approx_mode="growth");

F_sys = problem.system.f;
_X_ = problem.system.X;
_U_ = problem.system.U;

x0 = SVector(0.0, 0.0, 0.0);
h = SVector(0.2, 0.2, 0.2);
state_grid = DO.GridFree(x0, h);

u0 = SVector(0.0, 0.0);
h = SVector(0.3, 0.3);
input_grid = DO.GridFree(u0, h);

using JuMP
optimizer = MOI.instantiate(AB.SCOTS.Optimizer)
MOI.set(optimizer, MOI.RawOptimizerAttribute("problem"), problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.optimize!(optimizer)

abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("controller"))

nstep = 100
function reached(x)
    if x∈problem.target_set
        return true
    else
        return false
    end
end
x0 = SVector(0.4, 0.4, 0.0)
x_traj, u_traj = CO.get_closed_loop_trajectory(problem.system.f, controller, x0, nstep; stopping=reached)

fig = plot(aspect_ratio=:equal)

plot!(problem.system.X, color=:yellow, opacity=0.5)

abstract_system = AB.SCOTS.get_abstract_system(optimizer)
plot!(abstract_system.Xdom, color=:blue, opacity=0.5)

plot!(problem.initial_set, color=:green, opacity=0.2)
plot!(problem.target_set; dims=[1,2], color=:red, opacity=0.2)

abstract_problem = AB.SCOTS.get_abstract_problem(optimizer)
plot!(abstract_problem.initial_set, color=:green)
plot!(abstract_problem.target_set, color=:red)

plot!(fig, UT.DrawTrajectory(x_traj), ms=0.5)

display(fig)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

