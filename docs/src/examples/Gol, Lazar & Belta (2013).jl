using Test     #src
# # Example: Gol, Lazar and Belta (2013)
#
#md # [![Binder](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Gol%2C Lazar %26 Belta (2013).ipynb)
#md # [![nbviewer](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Gol%2C Lazar %26 Belta (2013).ipynb)
# This example reproduces parts of the numerical results of [3]. A similar example reproducing all results of [3] is available as a codeocean capsule in [4].
#
# This example was borrowed from [1, Example VIII.A] and tackles
# an optimal control for the hybrid system with state evolution governed by
# ```math
# x(k+1) = \begin{bmatrix} 1 & 1 \\ 0 & 1 \end{bmatrix}x(k) + \begin{bmatrix} 0.5 \\ 1.0 \end{bmatrix} u(k)
# ```

# The goal is to take the state vector toward a target set **XT** by visiting one of the squares
# **A** or **B** and avoiding the obstacles **O1** and **O2**


# First, let us import [CDDLib](https://github.com/JuliaPolyhedra/CDDLib.jl),
# [GLPK](https://github.com/jump-dev/GLPK.jl), [OSQP](https://github.com/oxfordcontrol/OSQP.jl),
# [JuMP](https://github.com/jump-dev/JuMP.jl), [Pavito](https://github.com/jump-dev/Pavito.jl)
# and [Ipopt](https://github.com/jump-dev/Ipopt.jl)
import CDDLib
import GLPK
import OSQP
using JuMP
import Pavito
import HiGHS
import Ipopt

# At this point we import Dionysos
using Dionysos
using Dionysos.Control
using Dionysos.Problem

# And the file defining the hybrid system for this problem
include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "GolLazarBelta.jl"))

# Now we instantiate our optimal control problem using the function provided by [GolLazarBelta.jl](@__REPO_ROOT_URL__/problems/GolLazarBelta.jl)
problem = GolLazarBelta.problem(CDDLib.Library(), Float64);

# Finally, we select the method presented in [2] as our optimizer

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


algo = optimizer_with_attributes(BemporadMorari.Optimizer{Float64},
    "continuous_solver" => qp_solver,
    "mixed_integer_solver" => miqp_solver,
    "indicator" => false,
    "log_level" => 0
);

# and use it to solve the given problem, with the help of the abstraction layer
# MathOptInterface provided by [JuMP](https://github.com/jump-dev/JuMP.jl)
optimizer = MOI.instantiate(algo)
MOI.set(optimizer, MOI.RawOptimizerAttribute("problem"), problem)
MOI.optimize!(optimizer)

# We check the solver time
MOI.get(optimizer, MOI.SolveTimeSec())

# the termination status
termination = MOI.get(optimizer, MOI.TerminationStatus())

# the objective value
objective_value = MOI.get(optimizer, MOI.ObjectiveValue())

@test objective_value ≈ 11.38 atol=1e-2     #src

# and recover the corresponding continuous trajectory
xu = MOI.get(optimizer, ContinuousTrajectoryAttribute());

# A little bit of data visualization now:

using PyPlot
using Colors
using Polyhedra
using HybridSystems

##Auxiliary function for annotating
function text_in_set_plot!(ax, po, t;  fillcolor = :white, linecolor = :black, fillalpha = 1)
    ##solve finding center (other solvers? https://jump.dev/JuMP.jl/dev/installation/#Supported-solvers)
    solver = optimizer_with_attributes(GLPK.Optimizer, "presolve" => GLPK.GLP_ON)
    poly = matplotlib.patches.Polygon(GolLazarBelta.get_ordered_vertices(po))
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
PyPlot.pygui(true) #jl
fig = PyPlot.figure()

ax = PyPlot.axes(aspect = "equal")
ax.set_xlim(-10.5,3)
ax.set_ylim(-10.5,3)

##Show the discrete modes
for mode in states(problem.system)
    t = (problem.system.ext[:q_T] in [mode, mode + 11]) ? "XT" : (mode == problem.system.ext[:q_A] ? "A" : (mode == problem.system.ext[:q_B] ? "B" :
            mode <= 11 ? string(mode) : string(mode - 11)))
    text_in_set_plot!(ax, stateset(problem.system, mode), t, fillcolor = "none", linecolor = :black)
end

##Plot obstacles
for i in eachindex(problem.system.ext[:obstacles])
    text_in_set_plot!(ax, problem.system.ext[:obstacles][i], "O$i", fillcolor = :black, fillalpha = 0.1)
end


##Initial state
ax.scatter([problem.initial_set[2][1]], [problem.initial_set[2][2]])
ax.annotate("x0", [problem.initial_set[2][1], problem.initial_set[2][2]-0.5], ha="center", va="center")

##Split the vector into x1 and x2
x1 = [xu.x[j][1] for j in eachindex(xu.x)]
x2 = [xu.x[j][2] for j in eachindex(xu.x)]

##Plot the trajectory
ax.scatter(x1, x2)
gcf() #md
gcf() #nb

# ### References
#
# 1. Gol, E. A., Lazar, M., & Belta, C. (2013). Language-guided controller synthesis for linear systems. IEEE Transactions on Automatic Control, 59(5), 1163-1176.
# 1. Bemporad, A., & Morari, M. (1999). Control of systems integrating logic, dynamics, and constraints. Automatica, 35(3), 407-427.
# 1. Legat B., Bouchat J., Jungers R. M. (2021). Abstraction-based branch and bound approach to Q-learning for hybrid optimal control. 3rd Annual Learning for Dynamics & Control Conference, 2021.
# 1. Legat B., Bouchat J., Jungers R. M. (2021). Abstraction-based branch and bound approach to Q-learning for hybrid optimal control. https://www.codeocean.com/. https://doi.org/10.24433/CO.6650697.v1.
