using Test     #src
# # Example: DC-DC converter solved by [Naive abstraction](https://github.com/dionysos-dev/Dionysos.jl/blob/master/docs/src/manual/manual.md#solvers).
#
#md # [![Binder](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/DC-DC converter.ipynb)
#md # [![nbviewer](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/DC-DC converter.ipynb)
#
# We consider a boost DC-DC converter which has been widely studied from the point of view of hybrid control, see for example in  [1, V.A],[2],[3].
# This is a **safety problem** for a **switching system**.
#
# ![Boost DC-DC converter.](https://github.com/dionysos-dev/Dionysos.jl/blob/master/docs/assets/dcdcboost.jpg?raw=true)
#
# The state of the system is given by $x(t) = \begin{bmatrix} i_l(t) & v_c(t) \end{bmatrix}^\top$.
# The switching system has two modes consisting in two-dimensional affine dynamics:
# ```math
# \dot{x} = f_p(x) = A_p x + b_p,\quad p=1,2
# ```
# with
# ```math
# A_1 = \begin{bmatrix} -\frac{r_l}{x_l} &0 \\ 0 & -\frac{1}{x_c}\frac{1}{r_0+r_c}  \end{bmatrix}, A_2= \begin{bmatrix} -\frac{1}{x_l}\left(r_l+\frac{r_0r_c}{r_0+r_c}\right) & -\frac{1}{x_l}\frac{r_0}{r_0+r_c}  \\ \frac{1}{x_c}\frac{r_0}{r_0+r_c}   & -\frac{1}{x_c}\frac{1}{r_0+r_c}  \end{bmatrix}, b_1 = b_2 = \begin{bmatrix} \frac{v_s}{x_l}\\0\end{bmatrix}.
# ```
# The goal is to design a controller to keep the state of the system in a safety region around the reference desired value, using as input only the switching
# signal. In order to study the concrete system and its symbolic abstraction in a unified framework, we will solve the problem
# for the sampled system with a sampling time $\tau$. For the construction of the relations in the abstraction, it is necessary to over-approximate attainable sets of
# a particular cell. In this example, we consider the use of a growth bound function  [4, VIII.2, VIII.5] which is one of the possible methods to over-approximate
# attainable sets of a particular cell based on the state reach by its center.
#

# First, let us import [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl) and [Plots](https://github.com/JuliaPlots/Plots.jl).

import CDDLib
import GLPK
import OSQP
using JuMP
import Pavito
import HiGHS
import Ipopt

# At this point, we import the useful Dionysos sub-modules.
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

# ### Definition of the system
# we can import the module containing the DCDC problem like this 
#include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "pc-osqp.jl"))
include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "pc-osqp-rmaps.jl"))

problem = PCOSQPRM.problem(CDDLib.Library(), Float64);

# Finally, we select the method presented in [2] as our optimizer

qp_solver = optimizer_with_attributes(
    OSQP.Optimizer,
    "eps_abs" => 1e-8,
    "eps_rel" => 1e-8,
    "max_iter" => 100000,
    MOI.Silent() => true,
);

mip_solver = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true);

cont_solver = optimizer_with_attributes(Ipopt.Optimizer, MOI.Silent() => true);

miqp_solver = optimizer_with_attributes(
    Pavito.Optimizer,
    "mip_solver" => mip_solver,
    "cont_solver" => cont_solver,
    MOI.Silent() => true,
);

algo = optimizer_with_attributes(
    OP.BemporadMorari.Optimizer{Float64},
    "continuous_solver" => qp_solver,
    "mixed_integer_solver" => miqp_solver,
    "indicator" => false,
    "log_level" => 0,
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

##@test objective_value ≈ 11.38 atol = 1e-2     #src

# and recover the corresponding continuous trajectory
xu = MOI.get(optimizer, ST.ContinuousTrajectoryAttribute());

# ## A little bit of data visualization now:

using Plots
using Polyhedra
using HybridSystems
using Suppressor

#Initialize our canvas
fig = plot(;
    aspect_ratio = :equal,
    xtickfontsize = 10,
    ytickfontsize = 10,
    guidefontsize = 16,
    titlefontsize = 14,
);
xlims!(-10.5, 3.0)
ylims!(-10.5, 3.0)

#Plot the discrete modes
for mode in states(problem.system)
    t =
        (problem.system.ext[:modes] in [mode, mode + 11]) ? "XT" :
        (
            mode == problem.system.ext[:q_A] ? "A" :
            (
                mode == problem.system.ext[:q_B] ? "B" :
                mode <= 11 ? string(mode) : string(mode - 11)
            )
        )
    set = stateset(problem.system, mode)
    plot!(set; color = :white)
    UT.text_in_set_plot!(fig, set, t)
end

#Plot obstacles
for i in eachindex(problem.system.ext[:obstacles])
    set = problem.system.ext[:obstacles][i]
    plot!(set; color = :black, opacity = 0.5)
    UT.text_in_set_plot!(fig, set, "O$i")
end

#Plot trajectory
x0 = problem.initial_set[2]
x_traj = [x0, xu.x...]
plot!(fig, UT.DrawTrajectory(x_traj));

#Plot initial point
plot!(fig, UT.DrawPoint(x0); color = :blue)
annotate!(fig, x0[1], x0[2] - 0.5, "x0")

# ### References
#
# 1. Gol, E. A., Lazar, M., & Belta, C. (2013). Language-guided controller synthesis for linear systems. IEEE Transactions on Automatic Control, 59(5), 1163-1176.
# 1. Bemporad, A., & Morari, M. (1999). Control of systems integrating logic, dynamics, and constraints. Automatica, 35(3), 407-427.
# 1. Legat B., Bouchat J., Jungers R. M. (2021). Abstraction-based branch and bound approach to Q-learning for hybrid optimal control. 3rd Annual Learning for Dynamics & Control Conference, 2021.
# 1. Legat B., Bouchat J., Jungers R. M. (2021). Abstraction-based branch and bound approach to Q-learning for hybrid optimal control. https://www.codeocean.com/. https://doi.org/10.24433/CO.6650697.v1.
