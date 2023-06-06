using StaticArrays, Plots

using Dionysos
const DI = Dionysos
const DO = DI.Domain
const OP = DI.Optim
const CO = DI.Control
const UT = DI.Utils

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "PathPlanning.jl"))

problem = PathPlanning.problem();

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
optimizer = MOI.instantiate(OP.Abstraction.Optimizer)
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
Plots.plot!(problem.system.X; dims=[1,2], color=:yellow, opacity=0.5)
Plots.plot!(problem.initial_set; dims=[1,2], color=:green)
Plots.plot!(problem.target_set; dims=[1,2], color=:red)
Plots.plot!(fig, UT.DrawTrajectory(x_traj))
display(fig)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

