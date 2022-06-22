```@meta
EditURL = "<unknown>/docs/src/examples/Path planning 2.jl"
```

# Example: Path planning problem

[![Binder](https://mybinder.org/badge_logo.svg)](<unknown>/generated/Path planning.ipynb)
[![nbviewer](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](<unknown>/generated/Path planning.ipynb)

This example was borrowed from [1, IX. Examples, A] whose dynamics comes from the model given in [2, Ch. 2.4].
This is a **reachability problem** for a **continuous system**.

Let us consider the 3-dimensional state space control system of the form
```math
\dot{x} = f(x, u)
```
with $f: \mathbb{R}^3 × U ↦ \mathbb{R}^3$ given by
```math
f(x,(u_1,u_2)) = \begin{bmatrix} u_1 \cos(α+x_3)\cos(α^{-1}) \\ u_1 \sin(α+x_3)\cos(α^{-1}) \\ u_1 \tan(u_2)  \end{bmatrix}
```
and with $U = [−1, 1] \times [−1, 1]$ and $α = \arctan(\tan(u_2)/2)$. Here, $(x_1, x_2)$ is the position and $x_3$ is the
orientation of the vehicle in the 2-dimensional plane. The control inputs $u_1$ and $u_2$ are the rear
wheel velocity and the steering angle.
The control objective is to drive the vehicle which is situated in a maze made of obstacles from an initial position to a target position.


In order to study the concrete system and its symbolic abstraction in a unified framework, we will solve the problem
for the sampled system with a sampling time $\tau$.

The abstraction is based on a feedback refinment relation [1,V.2 Definition].
Basically, this is equivalent to an alternating simulation relationship with the additional constraint that the input of the
concrete and symbolic system preserving the relation must be identical.
This allows to easily determine the controller of the concrete system from the abstraction controller by simply adding a quantization step.

For the construction of the relations in the abstraction, it is necessary to over-approximate attainable sets of
a particular cell. In this example, we consider the used of a growth bound function  [1, VIII.2, VIII.5] which is one of the possible methods to over-approximate
attainable sets of a particular cell based on the state reach by its center. Therefore, it is used
to compute the relations in the abstraction based on the feedback refinement relation.

For this reachability problem, the abstraction controller is built by solving a fixed-point equation which consists in computing the the pre-image
of the target set.

First, let us import [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl).

```@example Path Planning 2
using StaticArrays
```

At this point, we import the useful Dionysos sub-module for this problem: [Abstraction](<unknown>/src/Abstraction/abstraction.jl).

```@example Path Planning 2
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const CO = DI.Control
const SY = DI.Symbolic
```

And the file defining the hybrid system for this problem

```@example Path Planning 2
include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "PathPlanning.jl"))
```

### Definition of the problem

Now we instantiate the problem using the function provided by [PathPlanning.jl](<unknown>/problems/PathPlanning.jl)

```@example Path Planning 2
problem = PathPlanning.problem();
nothing #hide
```

`F_sys` is the function, `_X_` the state domain and `_U_` the input domain

```@example Path Planning 2
F_sys = problem.system.f;
_X_ = problem.system.X;
_U_ = problem.system.U;
nothing #hide
```

We define the growth bound function of $f$:

```@example Path Planning 2
ngrowthbound = 5;
function L_growthbound(u)
    β = abs(u[1]/cos(atan(tan(u[2])/2)))
    return SMatrix{3,3}(
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        β, β, 0.0)
end;
nothing #hide
```

Here it is considered that there is no system and measurement noise:

```@example Path Planning 2
sysnoise = SVector(0.0, 0.0, 0.0);
measnoise = SVector(0.0, 0.0, 0.0);
nothing #hide
```

We define the discretization time step parameters: `tstep` and `nsys`:

```@example Path Planning 2
tstep = 0.3;
nsys = 5;
nothing #hide
```

Finally, we build the control system:

```@example Path Planning 2
contsys = ST.NewControlSystemGrowthRK4(tstep, F_sys, L_growthbound, sysnoise,
                                       measnoise, nsys, ngrowthbound);
nothing #hide
```

### Definition of the abstraction

Definition of the grid of the state-space on which the abstraction is based (origin `x0` and state-space discretization `h`):

```@example Path Planning 2
x0 = SVector(0.0, 0.0, 0.0);
h = SVector(0.2, 0.2, 0.2);
Xgrid = DO.GridFree(x0, h);
nothing #hide
```

Construction of the struct `DomainList` containing the feasible cells of the state-space:

```@example Path Planning 2
Xfull = DO.DomainList(Xgrid);
DO.add_set!(Xfull, _X_, DO.OUTER);
nothing #hide
```

Definition of the grid of the input-space on which the abstraction is based (origin `u0` and input-space discretization `h`):

```@example Path Planning 2
u0 = SVector(0.0, 0.0);
h = SVector(0.3, 0.3);
Ugrid = DO.GridFree(u0, h);
nothing #hide
```

Construction of the struct `DomainList` containing the quantized inputs:

```@example Path Planning 2
Ufull = DO.DomainList(Ugrid);
DO.add_set!(Ufull, _U_, DO.OUTER);
nothing #hide
```

Construction of the abstraction:

```@example Path Planning 2
symmodel = SY.NewSymbolicModelListList(Xfull, Ufull);
@time SY.compute_symmodel_from_controlsystem!(symmodel, contsys)
```

### Construction of the controller

`_I_` is the initial state domain and `_T_` is the target state domain

```@example Path Planning 2
_I_ = problem.initial_set;
_T_ = problem.target_set;
nothing #hide
```

Computation of the initial symbolic states:

```@example Path Planning 2
Xinit = DO.DomainList(Xgrid);
DO.add_subset!(Xinit, Xfull, _I_, DO.OUTER)
initlist = [SY.get_state_by_xpos(symmodel, pos) for pos in DO.enum_pos(Xinit)];
nothing #hide
```

Computation of the target symbolic states:

```@example Path Planning 2
Xtarget = DO.DomainList(Xgrid)
DO.add_subset!(Xtarget, Xfull, _T_, DO.OUTER)
targetlist = [SY.get_state_by_xpos(symmodel, pos) for pos in DO.enum_pos(Xtarget)];
nothing #hide
```

Construction of the controller:

```@example Path Planning 2
contr = CO.NewControllerList();
@time CO.compute_controller_reach!(contr, symmodel.autom, initlist, targetlist)
```

### Trajectory display
We choose the number of steps `nsteps` for the sampled system, i.e. the total elapsed time: `nstep`*`tstep`
as well as the true initial state `x0` which is contained in the initial state-space `_I_` defined previously.

```@example Path Planning 2
nstep = 100;
x0 = SVector(0.4, 0.4, 0.0);
nothing #hide
```

Here we display the coordinate projection on the two first components of the state space along the trajectory.

To complete

### References
1. G. Reissig, A. Weber and M. Rungger, "Feedback Refinement Relations for the Synthesis of Symbolic Controllers," in IEEE Transactions on Automatic Control, vol. 62, no. 4, pp. 1781-1796.
2. K. J. Aström and R. M. Murray, Feedback systems. Princeton University Press, Princeton, NJ, 2008.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

