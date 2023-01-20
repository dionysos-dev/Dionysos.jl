using Test     #src
# # Example: State-Feedback Abstractions - Grid 
#
#md # [![Binder](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/State-Feedback Abstractions - Grid.ipynb)
#md # [![nbviewer](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/State-Feedback Abstractions - Grid.ipynb)
#
# This example was presented in [1, Example V.B] and tackles
# an optimal control for a piecewise-affine system with state evolution governed by
# ```math
# x(k+1) \in  F\big(x(k),u(k)\big),~~x(0)=x_0
# ```
# with
# ```math
# F(x,u) = \{A_{\psi(x)}x+B_{\psi(x)}u+g_{\psi(x)}\}\oplus\Omega_{\psi(x)},~~x(0)=x_0
# ```
# where $\psi:\mathcal{X}\rightarrow\{1,\dots,N_p\}$ selects one of the $N_p=3$ subsystems depending on
# the part of the state-space that the state $x$ is. We define these subsystems as 
# ```math
# 	A_1=\begin{bmatrix}
# 1.01 & 0.3\\
# -0.1 & 1.01
# \end{bmatrix},~B_1=\begin{bmatrix}
# 1&0\\ 0 & 1
# \end{bmatrix},~g_1=\begin{bmatrix}
# -0.1\\-0.1
# \end{bmatrix}
# ```
# $A_2=A_1^\top,~ A_3=A_1,~B_2=B_3=B_1,~g_2=0$ and $g_3=-g_1$. These systems are three spiral sources 
# with unstable equilibria at $x_{e1}=[-0.9635~~0.3654]^\top,~x_{e2}=0,$ and $x_{e3}=-x_{e1}$. We also 
# define the additive-noise sets $\Omega_1=\Omega_2=\Omega_3=[-0.05,0.05]^2$, the control-input set 
# $\mathcal{U}=[-0.5,0.5]^2$ and the state space $\mathcal{X}=[-2,2]^2$. The $N_p=3$ partitions of 
# $\mathcal{X}$ are $\mathcal{X}_1= \{x\in\mathcal{X}~:~x_1\leq-1 \},~\mathcal{X}_3= \{x\in\mathcal{X}~:~x_1>1 \},$ 
# and $\mathcal{X}_2=\mathcal{X}\setminus(\mathcal{X}_1\cup\mathcal{X}_3)$. The goal is to bring the state $x$ from 
# the initial set $\mathcal{X}_0$ to a final set $\mathcal{X}_*$, while avoiding the obstacle $\mathcal{O}$. 
# The associated stage-cost function is $\J(x,u)=[x~u~1]Q[x~u~1]^\top$  with $Q=\diag(10^{-2}I,0)$ which evenly 
# penalize states and inputs far away from the origin. This define the optimal control problem instantiate as follows

include("../problems/PWAsys.jl")

dt = 0.01;
const problem = PWAsys.problem(lib, dt, Usz)
const system = problem.system
#
#
# To build a deterministic state-feedback abstraction in alternating simulation relation  with the system as described 
# in [1], a set of balls of radius 0.2 covering the state space is adopted as cells $\xi\in\X_d$. We assume that inside cells intersecting the boundary of partitions of $\X$ the selected piecewise-affine mode is the same all over its interior and given by the mode defined at its center. An alternative to this assumption is to split these cells and use the S-Procedure to incorporate the respective cuts into the design problem, but we do not proceed in this way to favor a clear illustration of the results in this example. 

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
import Cbc
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
