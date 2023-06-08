import CDDLib
import GLPK
import OSQP
using JuMP
import Pavito
import HiGHS
import Ipopt

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const CO = DI.Control
const OP = DI.Optim

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "GolLazarBelta.jl"))

problem = GolLazarBelta.problem(CDDLib.Library(), Float64);

qp_solver = optimizer_with_attributes(
    OSQP.Optimizer,
    "eps_abs" => 1e-8,
    "eps_rel" => 1e-8,
    "max_iter" => 100000,
    MOI.Silent() => true
);

mip_solver = optimizer_with_attributes(
    HiGHS.Optimizer,
    MOI.Silent() => true
);

cont_solver = optimizer_with_attributes(
    Ipopt.Optimizer,
    MOI.Silent() => true
);

miqp_solver = optimizer_with_attributes(
    Pavito.Optimizer,
    "mip_solver" => mip_solver,
    "cont_solver" => cont_solver,
    MOI.Silent() => true
);


algo = optimizer_with_attributes(
    OP.BemporadMorari.Optimizer{Float64},
    "continuous_solver" => qp_solver,
    "mixed_integer_solver" => miqp_solver,
    "indicator" => false,
    "log_level" => 0
);

optimizer = MOI.instantiate(algo)
MOI.set(optimizer, MOI.RawOptimizerAttribute("problem"), problem)
MOI.optimize!(optimizer)

MOI.get(optimizer, MOI.SolveTimeSec())

termination = MOI.get(optimizer, MOI.TerminationStatus())

objective_value = MOI.get(optimizer, MOI.ObjectiveValue())

xu = MOI.get(optimizer, CO.ContinuousTrajectoryAttribute());

using Plots
using Polyhedra
using HybridSystems
using Suppressor

##Initialize our canvas
fig = plot(aspect_ratio=:equal, xtickfontsize=10, ytickfontsize=10, guidefontsize=16, titlefontsize=14)
xlims!(-10.5, 3.0)
ylims!(-10.5, 3.0)

##Plot the discrete modes
for mode in states(problem.system)
    t = (problem.system.ext[:q_T] in [mode, mode + 11]) ? "XT" : (mode == problem.system.ext[:q_A] ? "A" : (mode == problem.system.ext[:q_B] ? "B" :
            mode <= 11 ? string(mode) : string(mode - 11)))
    set = stateset(problem.system, mode)
    plot!(set, color=:white)
    UT.text_in_set_plot!(fig, set, t)
end

##Plot obstacles
for i in eachindex(problem.system.ext[:obstacles])
    set = problem.system.ext[:obstacles][i]
    plot!(set, color=:black, opacity=0.5)
    println(typeof(set))
    UT.text_in_set_plot!(fig, set, "O$i")
end

##Plot trajectory
x0 = problem.initial_set[2]
x_traj = [x0, xu.x...]
plot!(fig, UT.DrawTrajectory(x_traj))

##Plot initial point
plot!(fig, UT.DrawPoint(x0), color=:blue)
annotate!(fig, x0[1], x0[2]-0.5, "x0")
display(fig)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

