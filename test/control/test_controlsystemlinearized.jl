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

@testset "ControlSystemLinearized" begin
lb = SVector(0.0, 0.0)
ub = SVector(10.0, 11.0)
x0 = SVector(0.0, 0.0)
h = SVector(1.0, 2.0)
Xgrid = DO.GridFree(x0, h)
Xfull = DO.DomainList(Xgrid)
@test DO.get_ncells(Xfull) == 0
DO.add_set!(Xfull, UT.HyperRectangle(lb, ub), DO.OUTER)
@test DO.get_ncells(Xfull) == 77

lb = SVector(-1.0)
ub = SVector(1.0)
u0 = SVector(0.0)
h = SVector(0.5)
Ugrid = DO.GridFree(u0, h)
Ufull = DO.DomainList(Ugrid)
DO.add_set!(Ufull, UT.HyperRectangle(lb, ub), DO.OUTER)
@test DO.get_ncells(Ufull) == 5

tstep = 1.0
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

xpos = (1, 1)
upos = (1,)
x = DO.get_coord_by_pos(Xgrid, xpos)
u = DO.get_coord_by_pos(Ugrid, upos)
Xsimple = DO.DomainList(Xgrid)
DO.add_pos!(Xsimple, xpos)
@test DO.get_ncells(Xsimple) == 1
Usimple = DO.DomainList(Ugrid)
DO.add_pos!(Usimple, upos)
@test DO.get_ncells(Usimple) == 1

# @static if get(ENV, "CI", "false") == "false"
#     include("../../src/Abstraction/plotting.jl")
#     using PyPlot
#     fig = PyPlot.figure()
#     ax = fig.gca()
#     ax.set_xlim((-10.0, 10.0))
#     ax.set_ylim((-10.0, 10.0))
#     Plot.domain!(ax, 1:2, Xsimple)
#     Plot.trajectory_open_loop!(ax, 1:2, contsys, x, u, 50)
#     Plot.cell_image!(ax, 1:2, Xsimple, Usimple, contsys)
#     Plot.cell_approx!(ax, 1:2, Xsimple, Usimple, contsys)
# end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
