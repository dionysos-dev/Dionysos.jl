using Test     #src
# # Example: Path planning problem
#
#md # [![Binder](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Path planning.ipynb)
#md # [![nbviewer](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Path planning.ipynb)
#
# This example was borrowed from [1, IX. Examples, A] whose dynamics comes from the model given in [2, Ch. 2.4].
# This is a **safety problem** for a **continuous system*.
#
# Let us consider the 3-dimensional state space control system of the form
# ```math
# \dot{x} = f(x, u)
# ```
# with $f: \mathbb{R}^3 × U ↦ \mathbb{R}^3$ given by
# ```math
# f(x,(u_1,u_2)) = \begin{bmatrix} u_1 \cos(α+x_3)\cos(α^{-1}) \\ u_1 \sin(α+x_3)\cos(α^{-1}) \\ u_1 \tan(u_2)  \end{bmatrix}
# ```
# and with $U = [−1, 1] \times [−1, 1]$ and $α = \arctan(\tan(u_2)/2)$. Here, $(x_1, x_2)$ is the position and $x_3$ is the
# orientation of the vehicle in the 2-dimensional plane. The control inputs $u_1$ and $u_2$ are the rear
# wheel velocity and the steering angle.
# The control objective is to drive the vehicle which is situated in a maze made of obstacles from an initial position to a target position.
#
#
# In order to study the concrete system and its symbolic abstraction in a unified framework, we will solve the problem
# for the sampled system with a sampling time $\tau$.
#
# The abstraction is based on a feedback refinment relation [1,V.2 Definition].
# Basically, this is equivalent to an alternating simulation relationship with the additional constraint that the input of the
# concrete and symbolic system preserving the relation must be identical.
# This allows to easily determine the controller of the concrete system from the abstraction controller by simply adding a quantization step.
#
# For the construction of the relations in the abstraction, it is necessary to over-approximate attainable sets of
# a particular cell. In this example, we consider the used of a growth bound function  [1, VIII.2, VIII.5] which is one of the possible methods to over-approximate
# attainable sets of a particular cell based on the state reach by its center. Therefore, it is used
# to compute the relations in the abstraction based on the feedback refinement relation.
#
# For this reachability problem, the abstraction controller is built by solving a fixed-point equation which consists in computing the the pre-image
# of the target set.

# First, let us import [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl).

using StaticArrays

# At this point, we import the useful Dionysos sub-module for this problem: [Abstraction](@__REPO_ROOT_URL__/src/Abstraction/abstraction.jl).
using Dionysos
using Dionysos.Abstraction
AB = Dionysos.Abstraction;

# ### Definition of the control problem
# Definition of the state-space (limited to be rectangle):
_X_ = AB.HyperRectangle(SVector(0.0, 0.0, -pi-0.4), SVector(4.0, 10.0, pi+0.4));

# Definition of the obstacles (limited to be rectangle):
X1_lb = [1.0, 2.2,  2.2, 3.4,  4.6, 5.8,  5.8,  7.0, 8.2, 8.4,  9.3, 8.4,  9.3, 8.4,  9.3];
X1_ub = [1.2, 2.4,  2.4, 3.6,  4.8, 6.0,  6.0,  7.2, 8.4, 9.3, 10.0, 9.3, 10.0, 9.3, 10.0];
X2_lb = [0.0, 0.0,  6.0, 0.0,  1.0, 0.0,  7.0,  1.0, 0.0, 8.2,  7.0, 5.8,  4.6, 3.4,  2.2];
X2_ub = [9.0, 5.0, 10.0, 9.0, 10.0, 6.0, 10.0, 10.0, 8.5, 8.6,  7.4, 6.2,  5.0, 3.8,  2.6];

# Definition of the input-space (limited to be rectangle):
_U_ = AB.HyperRectangle(SVector(-1.0, -1.0), SVector(1.0, 1.0));

# Definition of the initial state-space (here it consists in a single point):
_I_ = AB.HyperRectangle(SVector(0.4, 0.4, 0.0), SVector(0.4, 0.4, 0.0));

# Definition of the target state-space (limited to be hyper-rectangle):
_T_ = AB.HyperRectangle(SVector(3.0, 0.5, -100.0), SVector(3.6, 0.8, 100.0));

# ### Definition of the system
# We define the discretization time step `tstep` and the dynamics function $f$ of the system:
tstep = 0.3;
nsys = 5;
ngrowthbound = 5;
function F_sys(x, u)
    α = atan(tan(u[2])/2)
    return SVector{3}(
        u[1]*cos(α + x[3])/cos(α),
        u[1]*sin(α + x[3])/cos(α),
        u[1]*tan(u[2]))
end;
# Definition of the growth bound function of $f$:
function L_growthbound(u)
    β = abs(u[1]/cos(atan(tan(u[2])/2)))
    return SMatrix{3,3}(
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        β, β, 0.0)
end;
# Here it is considered that there is no system and measurement noise:
sysnoise = SVector(0.0, 0.0, 0.0);
measnoise = SVector(0.0, 0.0, 0.0);

# Finally, we build the control system:
contsys = AB.NewControlSystemGrowthRK4(tstep, F_sys, L_growthbound, sysnoise,
                                       measnoise, nsys, ngrowthbound);

# ### Definition of the abstraction

# Definition of the grid of the state-space on which the abstraction is based (origin `x0` and state-space discretization `h`):
x0 = SVector(0.0, 0.0, 0.0);
h = SVector(0.2, 0.2, 0.2);
Xgrid = AB.GridFree(x0, h);
# Construction of the struct `DomainList` containing the feasible cells of the state-space:
Xfull = AB.DomainList(Xgrid);
AB.add_set!(Xfull, _X_, AB.OUTER)
for (x1lb, x2lb, x1ub, x2ub) in zip(X1_lb, X2_lb, X1_ub, X2_ub)
    box = AB.HyperRectangle(SVector(x1lb, x2lb, _X_.lb[3]), SVector(x1ub, x2ub, _X_.ub[3]))
    if box ⊆ _X_ && isempty(box ∩ _I_) && isempty(box ∩ _T_)
        AB.remove_set!(Xfull, box, AB.OUTER)
    end
end

# Definition of the grid of the input-space on which the abstraction is based (origin `u0` and input-space discretization `h`):
u0 = SVector(0.0, 0.0);
h = SVector(0.3, 0.3);
Ugrid = AB.GridFree(u0, h);
# Construction of the struct `DomainList` containing the quantized inputs:
Ufull = AB.DomainList(Ugrid);
AB.add_set!(Ufull, _U_, AB.OUTER)

# Construction of the abstraction:
symmodel = AB.NewSymbolicModelListList(Xfull, Ufull);
@time AB.compute_symmodel_from_controlsystem!(symmodel, contsys)

# ### Construction of the controller
# Computation of the initial symbolic states:
Xinit = AB.DomainList(Xgrid);
AB.add_subset!(Xinit, Xfull, _I_, AB.OUTER)
initlist = Int[]
for pos in AB.enum_pos(Xinit)
    push!(initlist, AB.get_state_by_xpos(symmodel, pos))
end
# Computation of the target symbolic states:
Xtarget = AB.DomainList(Xgrid)
AB.add_subset!(Xtarget, Xfull, _T_, AB.OUTER)
targetlist = Int[]
for pos in AB.enum_pos(Xtarget)
    push!(targetlist, AB.get_state_by_xpos(symmodel, pos))
end
# Construction of the controller:
contr = AB.NewControllerList();
@time AB.compute_controller_reach!(contr, symmodel.autom, initlist, targetlist)

# ### Trajectory display
# We choose the number of steps `nsteps` for the sampled system, i.e. the total elapsed time: `nstep`*`tstep`
# as well as the true initial state `x0` which is contained in the initial state-space `_I_` defined previously.
nstep = 100;
x0 = SVector(0.4, 0.4, 0.0);
# Here we display the coordinate projection on the two first components of the state space along the trajectory.
#
# to complete

# ### References
# 1. G. Reissig, A. Weber and M. Rungger, "Feedback Refinement Relations for the Synthesis of Symbolic Controllers," in IEEE Transactions on Automatic Control, vol. 62, no. 4, pp. 1781-1796.
# 2. K. J. Aström and R. M. Murray, Feedback systems. Princeton University Press, Princeton, NJ, 2008.
