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

@testset "ControllerSafe" begin
lb = SVector(-5.0, -5.0)
ub = SVector(5.0, 5.0)
x0 = SVector(0.0, 0.0)
h = SVector(0.47, 0.23)
Xgrid = DO.GridFree(x0, h)
Xfull = DO.DomainList(Xgrid)
DO.add_set!(Xfull, UT.HyperRectangle(lb, ub), DO.OUTER)

lb = SVector(-4.0)
ub = SVector(4.0)
u0 = SVector(0.0)
h = SVector(0.5)
Ugrid = DO.GridFree(u0, h)
Ufull = DO.DomainList(Ugrid)
DO.add_set!(Ufull, UT.HyperRectangle(lb, ub), DO.OUTER)

tstep = 0.2
nsys = 3
ngrowthbound = 3
F_sys(x, u) = SVector(u[1], -x[2] + u[1])
sysnoise = SVector(1.0, 1.0)*0.001
measnoise = SVector(1.0, 1.0)*0.001
L_growthbound(u) = SMatrix{2,2}(0.0, 0.0, 0.0, -1.0)

contsys = ST.NewControlSystemGrowthRK4(
    tstep, F_sys, L_growthbound, sysnoise, measnoise, nsys, ngrowthbound)
symmodel = SY.NewSymbolicModelListList(Xfull, Ufull)
SY.compute_symmodel_from_controlsystem!(symmodel, contsys)
@test SY.ntransitions(symmodel.autom) == 62077

Xinit = DO.DomainList(Xgrid)
DO.add_set!(Xinit, UT.HyperRectangle(SVector(-3.0, -3.0), SVector(-2.9, -2.9)), DO.OUTER)
initlist = Int[]
for pos in DO.enum_pos(Xinit)
    push!(initlist, SY.get_state_by_xpos(symmodel, pos))
end
Xsafe = DO.DomainList(Xgrid)
union!(Xsafe, Xfull)
DO.remove_set!(Xsafe, UT.HyperRectangle(SVector(-1.0, -2.0), SVector(-1.1, 4.0)), DO.OUTER)
safelist = Int[]
for pos in DO.enum_pos(Xsafe)
    push!(safelist, SY.get_state_by_xpos(symmodel, pos))
end

contr = CO.NewControllerList()
CO.compute_controller_safe!(contr, symmodel.autom, initlist, safelist)
@test length(contr) == 15043

invlist = Int[]
for source in 1:symmodel.autom.nstates
    if !isempty(UT.fix_and_eliminate_first(contr, source))
        push!(invlist, source)
    end
end
Xinv = DO.DomainList(Xgrid)
Yinv = DO.DomainList(Xgrid)
correct = true
for source in invlist
    DO.add_pos!(Xinv, SY.get_xpos_by_state(symmodel, source))
    if !correct
        break
    end
    targetlist = Int[]
    for symbol in UT.fix_and_eliminate_first(contr, source)
        SY.compute_post!(targetlist, symmodel.autom, source, symbol)
    end
    for target in targetlist
        DO.add_pos!(Yinv, SY.get_xpos_by_state(symmodel, target))
    end
    correct = correct && targetlist ⊆ safelist
end
@test correct

xpos = DO.get_somepos(Xinit)
x0 = DO.get_coord_by_pos(Xgrid, xpos)

# @static if get(ENV, "CI", "false") == "false"
#     include("../../src/Abstraction/plotting.jl")
#     using PyPlot
#     fig = PyPlot.figure()
#     ax = fig.gca()
#     ax.set_xlim((-5.5, 5.5))
#     ax.set_ylim((-5.3, 5.3))
#     Plot.domain!(ax, 1:2, Xfull, fa = 0.0)
#     Plot.domain!(ax, 1:2, Xinit)
#     Plot.domain!(ax, 1:2, Xsafe, fa = 0.1)
#     Plot.domain!(ax, 1:2, Xinv, fa = 0.1, fc = "yellow")
#     Plot.domain!(ax, 1:2, Yinv; fc = "blue")
#     Plot.trajectory_closed_loop!(
#         ax, 1:2, contsys, symmodel, contr, x0, 100, randchoose = true)
# end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
