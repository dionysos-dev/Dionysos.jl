include("../src/abstraction.jl")

module TestMain

using Test
using Main.Abstraction
AB = Main.Abstraction

sleep(0.1) # used for good printing
println("Started test")

@testset "FromControlSystem" begin
lb = (0.0, 0.0)
ub = (10.0, 11.0)
x0 = (0.0, 0.0)
h = (1.0, 2.0)
Xgrid = AB.NewGridSpaceList(x0, h)
AB.add_set!(Xgrid, AB.HyperRectangle(lb, ub), AB.OUTER)
Xfull = AB.NewSubSet(Xgrid)
AB.add_all!(Xfull)

lb = (-1.0,)
ub = (1.0,)
u0 = (0.0,)
h = (0.5,)
Ugrid = AB.NewGridSpaceList(u0, h)
AB.add_set!(Ugrid, AB.HyperRectangle(lb, ub), AB.OUTER)

tstep = 5.0
nsys = 3
nbound = 3
# F_sys(x, u) = [1.0-cos(x[2]), -x[1] + u[1]]
# L_bound(u) = [0.0 1.0; 1.0 0.0]
F_sys(x, u) = (u[1], -cos(x[1]))
L_bound(r, u) = (0.0, r[1])
sysnoise = (1.0, 1.0).*0.1
measnoise = (1.0, 1.0).*0.0

contsys = AB.NewControlSystemRK4(
    tstep, F_sys, L_bound, sysnoise, measnoise, nsys, nbound)
symmodel = AB.NewSymbolicModelListList(Xgrid, Ugrid)
AB.compute_symmodel_from_controlsystem!(symmodel, contsys)
@test AB.get_ntrans(symmodel.autom) == 1145

xpos = (1, 2)
upos = (1,)
x = AB.get_coord_by_pos(Xgrid, xpos)
u = AB.get_coord_by_pos(Ugrid, upos)
source = AB.get_state_by_xpos(symmodel, xpos)
symbol = AB.get_symbol_by_upos(symmodel, upos)

Xsimple = AB.NewSubSet(Xgrid)
AB.add_pos!(Xsimple, xpos)
Usimple = AB.NewSubSet(Ugrid)
AB.add_pos!(Usimple, upos)
Ysimple = AB.NewSubSet(Xgrid)
targetlist = Int[]
AB.compute_post!(targetlist, symmodel.autom, source, symbol)
for target in targetlist
    AB.add_pos!(Ysimple, AB.get_xpos_by_state(symmodel, target))
end

@static if get(ENV, "TRAVIS", "false") == "false"
    include("../src/plotting.jl")
    using PyPlot
    fig = PyPlot.figure()
    ax = fig.gca()
    ax.set_xlim([-1.0, 11.0])
    ax.set_ylim([-2.0, 14.0])
    Plot.subset!(ax, 1:2, Xfull, fa = 0.1)
    Plot.subset!(ax, 1:2, Xsimple)
    Plot.subset!(ax, 1:2, Ysimple)
    Plot.trajectory_open_loop!(ax, 1:2, contsys, x, u, 50)
    Plot.cell_image!(ax, 1:2, Xsimple, Usimple, contsys)
    Plot.cell_approx!(ax, 1:2, Xsimple, Usimple, contsys)
end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
