using StaticArrays
using Dionysos
using PyPlot

include("../../../src/plotting.jl")
AB = Dionysos.Abstraction;



# State-space defined 
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

tstep = 0.1
nsys=1;
A = SMatrix{2,2}(0.0, 1.0,
                    -3.0, 1.0);
B = SMatrix{2,1}(0.0, 1.0);

F_sys = let A = A
    (x,u) -> A*x + B*u
end;
DF_sys =  (x,u) -> A;
bound_DF = x -> 0;
bound_DDF = x -> 0;

measnoise = SVector(0.0, 0.0);

contsys = AB.NewControlSystemLinearizedRK4(tstep, F_sys, DF_sys, bound_DF, bound_DDF,
                    measnoise, nsys)

# Construction of the abstraction:
symmodel = AB.NewSymbolicModelListList(Xfull, Ufull);
AB.compute_symmodel_from_controlsystem!(symmodel, contsys)

Xinit = AB.DomainList(Xgrid);
union!(Xinit, Xfull)
initlist = [AB.get_state_by_xpos(symmodel, pos) for pos in AB.enum_pos(Xinit)];



xpos = AB.get_pos_by_coord(Xgrid, SVector(1.1, 1.3))
upos = AB.get_pos_by_coord(Ugrid, SVector(1))
x = AB.get_coord_by_pos(Xgrid, xpos)
u = AB.get_coord_by_pos(Ugrid, upos)
Xsimple = AB.DomainList(Xgrid)
AB.add_pos!(Xsimple, xpos)
Usimple = AB.DomainList(Ugrid)
AB.add_pos!(Usimple, upos)

PyPlot.pygui(true)
fig = PyPlot.figure()

ax = PyPlot.axes(aspect = "equal")

Plot.domain!(ax, 1:2, Xinit, fc = "green")
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)

Plot.domain!(ax, 1:2, Xsimple, fc = "red")
Plot.cell_image!(ax, 1:2, Xsimple, Usimple, contsys)
#Plot.cell_approx!(ax, 1:2, Xsimple, Usimple, contsys)

#Plot.cell_approx!(ax, 1:2, Xsimple, Usimple, contsys)
