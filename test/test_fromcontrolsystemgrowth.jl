include("../src/abstraction.jl")

module TestMain

using Test
using StaticArrays
using Main.Abstraction
AB = Main.Abstraction

sleep(0.1) # used for good printing
println("Started test")

@testset "FromControlSystem" begin
lb = SVector(0.0, 0.0)
ub = SVector(10.0, 11.0)
x0 = SVector(0.0, 0.0)
h = SVector(1.0, 2.0)
Xgrid = AB.GridFree(x0, h)
Xfull = AB.DomainList(Xgrid)
AB.add_set!(Xfull, AB.HyperRectangle(lb, ub), AB.OUTER)

lb = SVector(-1.0)
ub = SVector(1.0)
u0 = SVector(0.0)
h = SVector(0.5)
Ugrid = AB.GridFree(u0, h)
Ufull = AB.DomainList(Ugrid)
AB.add_set!(Ufull, AB.HyperRectangle(lb, ub), AB.OUTER)

tstep = 5.0
nsys = 3
ngrowthbound = 3
# F_sys(x, u) = [1.0-cos(x[2]), -x[1] + u[1]]
# L_growthbound(u) = [0.0 1.0; 1.0 0.0]
F_sys(x, u) = SVector(u[1], -cos(x[1]))
L_growthbound(u) = SMatrix{2,2}(0.0, 1.0, 0.0, 0.0)
sysnoise = SVector(1.0, 1.0)*0.1
measnoise = SVector(1.0, 1.0)*0.0

contsys = AB.NewControlSystemGrowthRK4(
    tstep, F_sys, L_growthbound, sysnoise, measnoise, nsys, ngrowthbound)
symmodel = AB.NewSymbolicModelListList(Xfull, Ufull)
AB.compute_symmodel_from_controlsystem!(symmodel, contsys)
@test AB.get_ntrans(symmodel.autom) == 1145

xpos = (1, 2)
upos = (1,)
x = AB.get_coord_by_pos(Xgrid, xpos)
u = AB.get_coord_by_pos(Ugrid, upos)
source = AB.get_state_by_xpos(symmodel, xpos)
symbol = AB.get_symbol_by_upos(symmodel, upos)

Xsimple = AB.DomainList(Xgrid)
AB.add_pos!(Xsimple, xpos)
Usimple = AB.DomainList(Ugrid)
AB.add_pos!(Usimple, upos)
Ysimple = AB.DomainList(Xgrid)
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
    ax.set_xlim((-1.0, 11.0))
    ax.set_ylim((-2.0, 14.0))
    Plot.domain!(ax, 1:2, Xfull, fa = 0.1)
    Plot.domain!(ax, 1:2, Xsimple)
    Plot.domain!(ax, 1:2, Ysimple)
    Plot.trajectory_open_loop!(ax, 1:2, contsys, x, u, 50)
    Plot.cell_image!(ax, 1:2, Xsimple, Usimple, contsys)
    Plot.cell_approx!(ax, 1:2, Xsimple, Usimple, contsys)
end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
