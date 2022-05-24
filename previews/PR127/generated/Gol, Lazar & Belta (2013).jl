import CDDLib
import GLPK
import OSQP
using JuMP
import Pavito
import Cbc
import Ipopt

using Dionysos
using Dionysos.Control
using Dionysos.Problem

include(joinpath(dirname(dirname(pathof(Dionysos))), "examples", "gol_lazar_belta.jl"))

system = gol_lazar_belta(CDDLib.Library(), Float64);

x0 = [1.0, -6.0];
q0 = 3;

N = 11;

state_cost = Fill(ZeroFunction(), nmodes(system))
transition_cost = QuadraticControlFunction(ones(1, 1))

problem = OptimalControlProblem(
    system,
    q0, x0,
    Fill(state_cost, N),
    Fill(Fill(transition_cost, ntransitions(system)), N),
    system.ext[:q_T],
    N
);

qp_solver = optimizer_with_attributes(
    OSQP.Optimizer,
    "eps_abs" => 1e-8,
    "eps_rel" => 1e-8,
    "max_iter" => 100000,
    MOI.Silent() => true
);

mip_solver = optimizer_with_attributes(
    Cbc.Optimizer,
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


algo = optimizer_with_attributes(BemporadMorari.Optimizer{Float64},
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

xu = MOI.get(optimizer, ContinuousTrajectoryAttribute());

using PyPlot
using Colors

##Auxiliary function for annotating
function text_in_set_plot!(ax, po, t;  fillcolor = :white, linecolor = :black, fillalpha = 1)
    ##solve finding center (other solvers? https://jump.dev/JuMP.jl/dev/installation/#Supported-solvers)
    solver = optimizer_with_attributes(GLPK.Optimizer, "presolve" => GLPK.GLP_ON)
    poly = matplotlib.patches.Polygon(get_ordered_vertices(po))
    poly.set_facecolor(fillcolor)
    poly.set_edgecolor(linecolor)
    poly.set_alpha(fillalpha)
    ax.add_patch(poly)

    if t !== nothing
        c, r = hchebyshevcenter(hrep(po), solver, verbose=0)
        ax.annotate(t, c, ha="center", va="center")
    end
end

##Initialize our canvas
PyPlot.pygui(true)
fig = PyPlot.figure()

ax = PyPlot.axes(aspect = "equal")
ax.set_xlim(-10.5,3)
ax.set_ylim(-10.5,3)

##Show the discrete modes
for mode in states(system)
    t = (system.ext[:q_T] in [mode, mode + 11]) ? "XT" : (mode == system.ext[:q_A] ? "A" : (mode == system.ext[:q_B] ? "B" :
            mode <= 11 ? string(mode) : string(mode - 11)))
    text_in_set_plot!(ax, stateset(system, mode), t, fillcolor = "none", linecolor = :black)
end

##Plot obstacles
for i in eachindex(system.ext[:obstacles])
    text_in_set_plot!(ax, system.ext[:obstacles][i], "O$i", fillcolor = :black, fillalpha = 0.1)
end


##Initial state
ax.scatter([x0[1]], [x0[2]])
ax.annotate("x0", [x0[1], x0[2]-0.5], ha="center", va="center")

##Split the vector into x1 and x2
x1 = [xu.x[j][1] for j in eachindex(xu.x)]
x2 = [xu.x[j][2] for j in eachindex(xu.x)]

##Plot the trajectory
ax.scatter(x1, x2)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

