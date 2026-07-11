# # Getting Started
#
#md # [![Binder](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Getting Started.ipynb)
#md # [![nbviewer](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Getting Started.ipynb)
#
# In this file we will visit the basic functionalities provided by Dionysos for the optimal control of complex systems. In summary, the topics covered are
# - Grids and discretizations
# - Dynamical system declaration
# - Continuous and discrete state image mapping
# - Plotting
# 
# First, let us import a few packages that are necessary to run this example.
using Dionysos
using StaticArrays, MathematicalSystems
using LinearAlgebra
using Plots

# The main package [Dionysos](https://github.com/dionysos-dev/Dionysos.jl) provides most important data structures that we will need.
# Additionally  [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl) provides faster implementation of Arrays (which have static memory allocation),
# [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) allows us to perform some additional operations.

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic

# Additionally, we will short the submodules accondingly 
#
# We use `UT.box` to represent the boundary of the state space `rectX` and the input space `rectU`.
rectX = UT.box(SVector(-2, -2), SVector(2, 2));
rectU = UT.box(SVector(-5), SVector(5));

# A discretization of the state space is declared using the `GridFree` structure, which requires the definition of a center `x0` and 
# a vector `hx` of discretization steps in each direction.
x0 = SVector(0.0, 0.0)
hx = SVector(1.0/5, 1.0/5)
Xgrid = MP.GridFree(x0, hx);

u0 = SVector(0.0)
hu = SVector(1.0/5)
Ugrid = MP.GridFree(u0, hu);

# `Xgrid` represents the state space grid and holds information of `x0` and `hx`, but is not a collection of cells. Indeed, a cell can be efficiently represented by a tuple of `Int`, for instance 'pos', with which the corresponding cartesian position can be computed by `x0 + hx .* pos` or using functions to be shown. In Dionysos, a set of cells is called a `Mapping`
# and the `ExplicitGridMapping` structure is the simplest used to represent this set. 
Xmap = MP.ExplicitGridMapping(Xgrid)
Umap = MP.ExplicitGridMapping(Ugrid)

MP.cover!(Xmap, rectX, MP.INNER)
MP.cover!(Umap, rectU, MP.INNER)

Xset = MP.FullStateSet{2}()  # default "all states"
Uset = MP.FullStateSet{1}()

# Now we have to define our dynamical system. For the sake of simplicity, note that we consider a linear time-invariant dynamical system but the functions 
# defining it allow the definition of a generic nonlinear and time-dependent system. We also define a step time `tstep` for discretizing the continuous-time dynamic.
# The parameters
tstep = 0.1

A = SMatrix{2, 2}(0.0, 1.0, -3.0, 1.0)
B = SMatrix{2, 1}(0.0, 1.0)

F_sys = (x, u) -> A*x + B*u

jacobian_bound = u -> abs.(A)

concrete_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
    F_sys,
    2,
    1,
    nothing,
    nothing,
)

continuous_approx =
    ST.ContinuousTimeGrowthBound(concrete_system; jacobian_bound = jacobian_bound)
discrete_approx = ST.discretize(continuous_approx, tstep)

concrete_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
    F_sys,
    2,
    1,
    nothing,
    nothing,
)
continuous_approx =
    ST.ContinuousTimeGrowthBound(concrete_system; jacobian_bound = jacobian_bound)
discrete_approx = ST.discretize(continuous_approx, tstep)

abstract_system = SY.SymbolicModelList(
    Xmap,
    Umap;
    Xset = Xset,
    Rset = Xset,  # allowed targets
    Uset = Uset,
)
SY.compute_abstract_system_from_concrete_system!(abstract_system, discrete_approx)

xpos = MP.get_pos_by_coord(Xgrid, SVector(1.1, 1.3))
x_center = MP.get_coord_by_pos(Xgrid, xpos)

q = MP.get_state_by_pos(Xmap, xpos)
abstract_input = 1
u = SY.get_concrete_input(abstract_system, abstract_input)

post = Int[]
SY.compute_post!(post, SY.get_automaton(abstract_system), q, abstract_input)

# Build a "post set" by adding post state
post_set = SY.get_state_set_from_states(abstract_system, post)

# Let us visualize this
fig = plot(; aspect_ratio = :equal)
dims = [1, 2]
xlims!(-2, 2);
ylims!(-2, 2)

plot!((Xset, Xmap); fc = "grey", dims = dims, label = "X", efficient = false)
plot!(
    (MP.stateset_from_states(Xmap, [q]), Xmap);
    fc = "blue",
    dims = dims,
    label = "x cell",
)
plot!((post_set, Xmap); fc = "green", dims = dims, label = "Post", efficient = false)

# In the previous picture, we have the state space lattice in grey, the chosen cell in blue and 
# the corresponding Post image in green.
