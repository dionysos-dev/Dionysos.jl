using Test     #src
# # Example: DC-DC converter
#
#md # [![Binder](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/DC-DC converter.ipynb)
#md # [![nbviewer](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/DC-DC converter.ipynb)
#
# We consider a boost DC-DC converter which has been widely studied from the point of view of hybrid control, see for example in  [1, V.A],[2],[3].
# This is a **safety problem** for a **switching system**.
#
# ![Boost DC-DC converter.](https://github.com/JulienCalbert/Dionysos.jl/blob/master/docs/assets/dcdcboost.jpg?raw=true)
# The state of the system is given by $x(t) = \begin{bmatrix} i_l(t) & v_c(t) \end{bmatrix}^\top$.
# The switching system has two modes consisting in two-dimensional affine dynamics:
# ```math
# \dot{x} = f_p(x) = A_p x + b_p,\quad p=1,2
# ```
# with
# ```math
# A_1 = \begin{bmatrix} -\frac{r_l}{x_l} &0 \\ 0 & -\frac{1}{x_c}\frac{1}{r_0+r_c}  \end{bmatrix}, A_2= \begin{bmatrix} -\frac{1}{x_l}\left(r_l+\frac{r_0r_c}{r_0+r_c}\right) & -\frac{1}{x_l}\frac{r_0}{r_0+r_c}  \\ \frac{1}{x_c}\frac{r_0}{r_0+r_c}   & -\frac{1}{x_c}\frac{1}{r_0+r_c}  \end{bmatrix}, b = \begin{bmatrix} \frac{v_s}{x_l}\\0\end{bmatrix}.
# ```
# The goal is to design a controller to keep the state of the system in a safety region around the reference desired value, using as input only the switching
# signal.
#
#
# In order to study the concrete system and its symbolic abstraction in a unified framework, we will solve the problem
# for the sampled system with a sampling time $\tau$.
#
# The abstraction is based on a feedback refinment relation [4,V.2 Definition].
# Basically, this is equivalent to an alternating simulation relationship with the additional constraint that the input of the
# concrete and symbolic system preserving the relation must be identical.
# This allows to easily determine the controller of the concrete system from the abstraction controller by simply adding a quantization step.
#
# For the construction of the relations in the abstraction, it is necessary to over-approximate attainable sets of
# a particular cell. In this example, we consider the used of a growth bound function  [4, VIII.2, VIII.5] which is one of the possible methods to over-approximate
# attainable sets of a particular cell based on the state reach by its center. Therefore, it is used
# to compute the relations in the abstraction based on the feedback refinement relation.
#

# First, let us import [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl).

using StaticArrays

# At this point, we import the useful Dionysos sub-module for this problem: [Abstraction](@__REPO_ROOT_URL__/src/Abstraction/abstraction.jl).
using Dionysos
using Dionysos.Abstraction
AB = Dionysos.Abstraction;

# ### Definition of the system

# Definition of the parameters of the system:
vs = 1.0; rL = 0.05; xL = 3.0; rC = 0.005; xC = 70.0; r0 = 1.0;

# Definition of the dynamics functions $f_p$ of the system:
b = SVector(vs/xL, 0.0);
A1 = SMatrix{2,2}(-rL/xL, 0.0, 0.0, -1.0/xC/(r0+rC));
A2 = SMatrix{2,2}(-(rL+r0*rC/(r0+rC))/xL, 5.0*r0/(r0+rC)/xC,
    -r0/(r0+rC)/xL/5.0, -1.0/xC/(r0+rC));
F_sys = let b = b, A1 = A1, A2 = A2
    (x, u) -> u[1] == 1 ? A1*x + b : A2*x + b
end;
# Definition of the growth bound functions of $f_p$:
ngrowthbound = 5;
A2_abs = SMatrix{2,2}(-(rL+r0*rC/(r0+rC))/xL, 5.0*r0/(r0+rC)/xC,
                      r0/(r0+rC)/xL/5.0, -1.0/xC/(r0+rC));
L_growthbound = let A1 = A1, A2_abs = A2_abs
    u -> u[1] == 1 ? A1 : A2_abs
end;
# Here it is considered that there is no system and measurement noise:
sysnoise = SVector(0.0, 0.0);
measnoise = SVector(0.0, 0.0);

# Definition of the discretization time step parameters: `tstep` and `nsys`:
tstep = 0.5;
nsys = 5;

# Finally, we build the control system:
contsys = AB.NewControlSystemGrowthRK4(tstep, F_sys, L_growthbound, sysnoise,
                                       measnoise, nsys, ngrowthbound);

# ### Definition of the control problem
# Definition of the state-space (limited to be rectangle):
_X_ = AB.HyperRectangle(SVector(1.15, 5.45), SVector(1.55, 5.85));

# Definition of the input-space, the later discretization of the input ensures that it can only take the values $1$ or $2$:
_U_ = AB.HyperRectangle(SVector(1), SVector(2));

# ### Definition of the abstraction

# Definition of the grid of the state-space on which the abstraction is based (origin `x0` and state-space discretization `h`):
x0 = SVector(0.0, 0.0);
h = SVector(2.0/4.0e3, 2.0/4.0e3);
Xgrid = AB.GridFree(x0, h);
# Construction of the struct `DomainList` containing the feasible cells of the state-space.
# Note, we used `AB.INNER` to make sure to add cells entirely contained in the domain because we are working with a safety problem.
Xfull = AB.DomainList(Xgrid);
AB.add_set!(Xfull, _X_, AB.INNER)

# Definition of the grid of the input-space on which the abstraction is based (origin `u0` and input-space discretization `h`):
u0 = SVector(1);
h = SVector(1);
Ugrid = AB.GridFree(u0, h);
# Construction of the struct `DomainList` containing the quantized inputs:
Ufull = AB.DomainList(Ugrid);
AB.add_set!(Ufull, _U_, AB.OUTER);

# Construction of the abstraction:
symmodel = AB.NewSymbolicModelListList(Xfull, Ufull);
@time AB.compute_symmodel_from_controlsystem!(symmodel, contsys)

# ### Construction of the controller
# In this problem, we consider both: the initial state-space and the safety state-space are equal to the entire state-space.
#
# Computation of the initial symbolic states:
Xinit = AB.DomainList(Xgrid);
union!(Xinit, Xfull)
initlist = Int[]
for pos in AB.enum_pos(Xinit)
    push!(initlist, AB.get_state_by_xpos(symmodel, pos))
end
# Computation of the safety symbolic states:
Xsafe = AB.DomainList(Xgrid)
union!(Xsafe, Xfull)
safelist = Int[]
for pos in AB.enum_pos(Xsafe)
    push!(safelist, AB.get_state_by_xpos(symmodel, pos))
end
# Construction of the controller:
contr = AB.NewControllerList();
@time AB.compute_controller_safe!(contr, symmodel.autom, initlist, safelist)

# ### Trajectory display
# We choose the number of steps `nsteps` for the sampled system, i.e. the total elapsed time: `nstep`*`tstep`
# as well as the true initial state `x0` which is contained in the initial state-space defined previously.
nstep = 300;
x0 = SVector(1.2, 5.6);
# To complete

# ### References
# 1. A. Girard, G. Pola and P. Tabuada, "Approximately Bisimilar Symbolic Models for Incrementally Stable Switched Systems," in IEEE Transactions on Automatic Control, vol. 55, no. 1, pp. 116-126, Jan. 2010.
# 2. S. Mouelhi, A. Girard, and G. Gössler. “CoSyMA: a tool for controller synthesis using multi-scale abstractions”. In: HSCC. ACM. 2013, pp. 83–88.
# 3. A. Girard. “Controller synthesis for safety and reachability via approximate bisimulation”. In: Automatica 48.5 (2012), pp. 947–953.
# 4. G. Reissig, A. Weber and M. Rungger, "Feedback Refinement Relations for the Synthesis of Symbolic Controllers," in IEEE Transactions on Automatic Control, vol. 62, no. 4, pp. 1781-1796.
