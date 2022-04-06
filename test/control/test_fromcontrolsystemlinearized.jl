module TestMain

using Test
using StaticArrays
using ..Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const CO = DI.Control
const SY = DI.Symbolic

sleep(0.1) # used for good printing
println("Started test")

@testset "FromControlSystem" begin
lb = SVector(0.0, 0.0)
ub = SVector(10.0, 11.0)
x0 = SVector(0.0, 0.0)
h = SVector(1.0, 2.0)
Xgrid = DO.GridFree(x0, h)
Xfull = DO.DomainList(Xgrid)
DO.add_set!(Xfull, UT.HyperRectangle(lb, ub), DO.OUTER)

lb = SVector(-1.0)
ub = SVector(1.0)
u0 = SVector(0.0)
h = SVector(0.5)
Ugrid = DO.GridFree(u0, h)
Ufull = DO.DomainList(Ugrid)
DO.add_set!(Ufull, UT.HyperRectangle(lb, ub), DO.OUTER)

tstep = 0.5
nsys = 3
# F_sys(x, u) = (1.0-cos(x(2)), -x(1) + u(1)
# L_growthbound(u) = (0.0 1.0; 1.0 0.0)
F_sys(x, u) = SVector(u[1], -cos(x[1]))
DF_sys(x, u) = SMatrix{2,2}(0.0, sin(x[1]), 0.0, 0.0)
bound_DF(u) = 1.0
# DDF_1 = [0.0 0.0; 0.0 0.0]
# DDF_2 = [cos(x[1]) 0.0; 0.0 0.0]
bound_DDF(u) = 1.0
measnoise = SVector(1.0, 1.0)*0.0

contsys = ST.NewControlSystemLinearizedRK4(
    tstep, F_sys, DF_sys, bound_DF, bound_DDF, measnoise, nsys)
symmodel = SY.NewSymbolicModelListList(Xfull, Ufull)
SY.compute_symmodel_from_controlsystem!(symmodel, contsys)
@test SY.ntransitions(symmodel.autom) == 2155

xpos = (1, 2)
upos = (1,)
x = DO.get_coord_by_pos(Xgrid, xpos)
u = DO.get_coord_by_pos(Ugrid, upos)
source = SY.get_state_by_xpos(symmodel, xpos)
symbol = SY.get_symbol_by_upos(symmodel, upos)

Xsimple = DO.DomainList(Xgrid)
DO.add_pos!(Xsimple, xpos)
Usimple = DO.DomainList(Ugrid)
DO.add_pos!(Usimple, upos)
Ysimple = DO.DomainList(Xgrid)
targetlist = Int[]
SY.compute_post!(targetlist, symmodel.autom, source, symbol)
for target in targetlist
    DO.add_pos!(Ysimple, SY.get_xpos_by_state(symmodel, target))
end

# @static if get(ENV, "CI", "false") == "false"
#     include("../../src/Abstraction/plotting.jl")
#     using PyPlot
#     fig = PyPlot.figure()
#     ax = fig.gca()
#     ax.set_xlim((-1.0, 11.0))
#     ax.set_ylim((-2.0, 14.0))
#     Plot.domain!(ax, 1:2, Xfull, fa = 0.1)
#     Plot.domain!(ax, 1:2, Xsimple)
#     Plot.domain!(ax, 1:2, Ysimple)
#     Plot.trajectory_open_loop!(ax, 1:2, contsys, x, u, 50)
#     Plot.cell_image!(ax, 1:2, Xsimple, Usimple, contsys)
#     Plot.cell_approx!(ax, 1:2, Xsimple, Usimple, contsys)
# end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
