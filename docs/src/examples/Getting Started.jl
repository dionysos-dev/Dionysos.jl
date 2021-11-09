# # Getting Started
#
#md # [![Binder](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Getting Started.ipynb)
#md # [![nbviewer](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Getting Started.ipynb)
#
# In this file we will visit the basic functionalities provided by Dionysos for the optimal control of complex systems.
# 
# First, let us import a few packages that are necessary to run this example.
using Dionysos
using StaticArrays
using LinearAlgebra
using PyPlot

# The main package [Dionysos](https://github.com/dionysos-dev/Dionysos.jl) provides most important data structures that we will need.
# Additionally  [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl) provides faster implementation of Arrays (which have static memory allocation),
# [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) allows us to perform some additional operations and [PyPlot](https://github.com/JuliaPy/PyPlot.jl) is
# important for data visualization.

include("../../../src/plotting.jl")

# The submodule [plotting.jl](@__REPO_ROOT_URL__/src/plotting.jl) has functions that will be useful for 2D-visualization of the functions that we are implementing.

AB = Dionysos.Abstraction;

# Additionally, we will short the [Abstraction](@__REPO_ROOT_URL__/src/Abstraction/abstraction.jl) submodule as `AB` 
#
# We use Hyper
_X_ = AB.HyperRectangle(SVector(-2, -2), SVector(2, 2));
_U_ = AB.HyperRectangle(SVector(-5), SVector(5));

x0 = SVector(0.0, 0.0);
h = SVector(1.0/5, 1.0/5);
Xgrid = AB.GridFree(x0, h);
# Construction of the struct `DomainList` containing the feasible cells of the state-space.
# Note, we used `AB.INNER` to make sure to add cells entirely contained in the domain because we are working with a safety problem.
Xfull = AB.DomainList(Xgrid);

AB.add_set!(Xfull, _X_, AB.INNER)


# Definition of the grid of the input-space on which the abstraction is based (origin `u0` and input-space discretization `h`):
u0 = SVector(0.0);
h = SVector(1.0/5);
Ugrid = AB.GridFree(u0, h);
# Construction of the struct `DomainList` containing the quantized inputs:
Ufull = AB.DomainList(Ugrid);
AB.add_set!(Ufull, _U_, AB.INNER);

tstep = 0.1;
nsys=10;
ngrowthbound=10;
A = SMatrix{2,2}(0.0, 1.0,
                    -3.0, 1.0);
B = SMatrix{2,1}(0.0, 1.0);

F_sys = let A = A
    (x,u) -> A*x + B*u
end;
L_growthbound = x -> abs.(A)

DF_sys =  (x,u) -> A;
bound_DF = x -> 0;
bound_DDF = x -> 0;

measnoise = SVector(0.0, 0.0);
sysnoise = SVector(0.0, 0.0);
contsys = AB.NewControlSystemGrowthRK4(tstep, F_sys, L_growthbound, sysnoise,
                                       measnoise, nsys, ngrowthbound);

# Construction of the abstraction:
symmodel = AB.NewSymbolicModelListList(Xfull, Ufull);
AB.compute_symmodel_from_controlsystem!(symmodel, contsys)

Xinit = AB.DomainList(Xgrid);
union!(Xinit, Xfull)
initlist = [AB.get_state_by_xpos(symmodel, pos) for pos in AB.enum_pos(Xinit)];



xpos = AB.get_pos_by_coord(Xgrid, SVector(1.1, 1.3))
upos = AB.get_pos_by_coord(Ugrid, SVector(-1))

x = AB.get_coord_by_pos(Xgrid, xpos)
u = AB.get_coord_by_pos(Ugrid, upos)


Xspecific = AB.DomainList(Xgrid)
AB.add_pos!(Xspecific, xpos)

Uspecific = AB.DomainList(Ugrid)
AB.add_pos!(Uspecific, upos)

PyPlot.pygui(true)
fig = PyPlot.figure()

ax = PyPlot.axes(aspect = "equal")

Plot.domain!(ax, 1:2, Xinit, fc = "green")
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)

Plot.domain!(ax, 1:2, Xspecific, fc = "blue")


r = Xgrid.h/2.0
Fr = contsys.growthbound_map(r, u, -contsys.tstep) + contsys.measnoise
Fx = contsys.sys_inv_map(x, u, contsys.tstep)
rectI = AB.get_pos_lims_outer(Xgrid, AB.HyperRectangle(Fx - Fr, Fx + Fr))


Xtarget = AB.DomainList(Xgrid)
for y in Iterators.product(AB._ranges(rectI)...)
    target = AB.get_state_by_xpos(symmodel, y)
    AB.add_pos!(Xtarget, y)
end

Plot.domain!(ax, 1:2, Xtarget, fc = "red")
Plot.cell_pre_image!(ax, 1:2, Xspecific, Uspecific, contsys)

Plot.cell_image!(ax, 1:2, Xspecific, Uspecific, contsys)

#Plot.cell_approx!(ax, 1:2, Xsimple, Usimple, contsys)

#Plot.cell_approx!(ax, 1:2, Xsimple, Usimple, contsys)
