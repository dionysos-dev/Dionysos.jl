using Test     #src
# # Example: DC-DC converter solved by [Uniform grid abstraction](https://github.com/dionysos-dev/Dionysos.jl/blob/master/docs/src/manual/manual.md#solvers).
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

using StaticArrays, Plots

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
include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "dc_dc.jl"))

# and we can instantiate the DC system with the provided system
problem = DCDC.problem(; approx_mode = DCDC.GROWTH)

controller = problem.solve(;
    method = :uniform_abstraction,
    initial_state = [0.0, 0.0],
    initial_input = [1.0],
)
# Didn't work because the default step sizes are too coarse / too small
problem.set_state_stepsize([2.0 / 4.0e3, 2.0 / 4.0e3])
problem.set_input_stepsize([1.0])
controller = problem.solve()
# Now it worked

# ### Trajectory display
# We choose the number of steps `nsteps` for the sampled system, i.e. the total elapsed time: `nstep`*`tstep`
# as well as the true initial state `x0` which is contained in the initial state-space defined previously.
nstep = 300
x0 = [1.2, 5.6]
control_trajectory = ST.get_closed_loop_trajectory(problem.system.f, controller, x0, nstep)

fig = plot(; aspect_ratio = :equal);
plot!(problem.system.X);
plot!(control_trajectory)

# # Example: DC-DC converter solved by [Uniform grid abstraction] (https://github.com/dionysos-dev/Dionysos.jl/blob/master/docs/src/manual/manual.md#solvers) by exploiting the incremental stability of the system.
# ### Definition of the system
# we can import the module containing the DCDC problem like this 

# and we can instantiate the DC system with the provided system
problem = DCDC.problem(; approx_mode = DCDC.DELTA_GAS)

origin = SVector(0.0, 0.0)
η = (2 / 4.0) * 10^(-3)

# Note: In the following, `P` and `ϵ` are computed by hand, but their computation is not crucial since they only affect the visualization of the abstraction. See https://github.com/dionysos-dev/Dionysos.jl/issues/345
ϵ = 0.1 * 0.01
P = SMatrix{2, 2}(1.0224, 0.0084, 0.0084, 1.0031)
state_grid = DO.GridEllipsoidalRectangular(origin, SVector(η, η), P / ϵ, concrete_system.X)

# Here we give a more complex grid to the solver
controller = problem.solve(;
    method = :uniform_abstraction,
    δ_gas = true,
    state_grid = DO.GridEllipsoidalRectangular(
        origin,
        SVector(η, η),
        P / ϵ,
        concrete_system.X,
    ),
    initial_input = [1.0],
    input_stepsize = [1.0],
)

# ### Trajectory display
# We choose the number of steps `nsteps` for the sampled system, i.e. the total elapsed time: `nstep`*`tstep`
# as well as the true initial state `x0` which is contained in the initial state-space defined previously.
nstep = 300
x0 = [1.2, 5.6]
control_trajectory = ST.get_closed_loop_trajectory(problem.system.f, controller, x0, nstep)

fig = plot(; aspect_ratio = :equal);
plot!(problem.system.X);
plot!(control_trajectory)

# ### References
# 1. A. Girard, G. Pola and P. Tabuada, "Approximately Bisimilar Symbolic Models for Incrementally Stable Switched Systems," in IEEE Transactions on Automatic Control, vol. 55, no. 1, pp. 116-126, Jan. 2010.
# 2. S. Mouelhi, A. Girard, and G. Gössler. “CoSyMA: a tool for controller synthesis using multi-scale abstractions”. In: HSCC. ACM. 2013, pp. 83–88.
# 3. A. Girard. “Controller synthesis for safety and reachability via approximate bisimulation”. In: Automatica 48.5 (2012), pp. 947–953.
# 4. G. Reissig, A. Weber and M. Rungger, "Feedback Refinement Relations for the Synthesis of Symbolic Controllers," in IEEE Transactions on Automatic Control, vol. 62, no. 4, pp. 1781-1796.
