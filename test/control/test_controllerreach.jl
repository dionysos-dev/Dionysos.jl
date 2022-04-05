include("../../src/Dionysos.jl")

module Test


# using TestAbstraction
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

@testset "ControllerReach" begin
lb = SVector(-5.0, -5.0)
ub = SVector(5.0, 5.0)
x0 = SVector(0.0, 0.0)
h = SVector(0.47, 0.23)
Xgrid = DO.GridFree(x0, h)
Xfull = DO.DomainList(Xgrid)
DO.add_set!(Xfull, UT.HyperRectangle(lb, ub), DO.OUTER)
DO.remove_set!(Xfull, UT.HyperRectangle(SVector(-1.0, -2.0), SVector(-1.1, 4.0)), DO.OUTER)

lb = SVector(-2.0)
ub = SVector(2.0)
u0 = SVector(0.0)
h = SVector(1.0)
Ugrid = DO.GridFree(u0, h)
Ufull = DO.DomainList(Ugrid)
DO.add_set!(Ufull, UT.HyperRectangle(lb, ub), DO.OUTER)

tstep = 1.0
nsys = 3
ngrowthbound = 3
F_sys(x, u) = SVector(1.0, u[1])
sysnoise = SVector(1.0, 1.0)*0.001
measnoise = SVector(1.0, 1.0)*0.001
L_growthbound(u) = SMatrix{2,2}(0.0, 0.0, 0.0, 0.0)

contsys = ST.NewControlSystemGrowthRK4(
    tstep, F_sys, L_growthbound, sysnoise, measnoise, nsys, ngrowthbound)
symmodel = SY.NewSymbolicModelListList(Xfull, Ufull)
SY.compute_symmodel_from_controlsystem!(symmodel, contsys)
SY.compute_symmodel_from_controlsystem!(symmodel, contsys)
display(symmodel)

Xinit = DO.DomainList(Xgrid)
DO.add_set!(Xinit, UT.HyperRectangle(SVector(-3.0, -3.0), SVector(-2.9, -2.9)), DO.OUTER)
initlist = Int[]
for pos in DO.enum_pos(Xinit)
    push!(initlist, SY.get_state_by_xpos(symmodel, pos))
end
Xtarget = DO.DomainList(Xgrid)
DO.add_set!(Xtarget, UT.HyperRectangle(SVector(0.0, 0.0), SVector(4.0, 4.0)), DO.OUTER)
targetlist = Int[]
for pos in DO.enum_pos(Xtarget)
    push!(targetlist, SY.get_state_by_xpos(symmodel, pos))
end

contr = CO.NewControllerList()
CO.compute_controller_reach!(contr, symmodel.autom, initlist, targetlist)
@test length(contr) == 412
if VERSION >= v"1.5"
    function f(autom, initlist, targetlist)
        contr = CO.NewControllerList()
        initset, targetset, num_targets_unreachable, current_targets, next_targets = CO._data(contr, autom, initlist, targetlist)
        # Preallocates to make sure `_compute_controller_reach` does not need to allocate
        sizehint!(contr.data, 500)
        sizehint!(current_targets, 50)
        sizehint!(next_targets, 200)
        @allocated CO._compute_controller_reach!(contr, autom, initset, targetset, num_targets_unreachable, current_targets, next_targets)
    end
    f(symmodel.autom, initlist, targetlist)
    @test f(symmodel.autom, initlist, targetlist) == 3224896 # check this
end

xpos = DO.get_somepos(Xinit)
x0 = DO.get_coord_by_pos(Xgrid, xpos)
Xsimple = DO.DomainList(Xgrid)
XUYsimple_ = Any[]

for i in 1:6
    source = SY.get_state_by_xpos(symmodel, xpos)
    Xs = DO.DomainList(Xgrid)
    Ys = DO.DomainList(Xgrid)
    Us = DO.DomainList(Ugrid)
    DO.add_pos!(Xs, xpos)
    targetlist = Int[]
    for (symbol,) in UT.fix_and_eliminate_first(contr, source)
        SY.compute_post!(targetlist, symmodel.autom, source, symbol)
        DO.add_pos!(Us, SY.get_upos_by_symbol(symmodel, symbol))
    end
    for target in targetlist
        DO.add_pos!(Ys, SY.get_xpos_by_state(symmodel, target))
    end
    push!(XUYsimple_, (Xs, Us, Ys))
    if !isempty(Ys)
        xpos = DO.get_somepos(Ys)
    end
end

# @static if get(ENV, "CI", "false") == "false"
#     include("../../src/utils/plotting/plotting.jl")
#     using PyPlot
#     fig = PyPlot.figure()
#     ax = fig.gca()
#     ax.set_xlim((-5.5, 5.5))
#     ax.set_ylim((-5.3, 5.3))
#     Plot.domain!(ax, 1:2, Xfull, fa = 0.1)
#     Plot.domain!(ax, 1:2, Xinit)
#     Plot.domain!(ax, 1:2, Xtarget)
#     for (Xs, Us, Ys) in XUYsimple_
#         Plot.domain!(ax, 1:2, Xs, fc = "green")
#         Plot.domain!(ax, 1:2, Ys, fc = "blue")
#         Plot.cell_image!(ax, 1:2, Xs, Us, contsys)
#         Plot.cell_approx!(ax, 1:2, Xs, Us, contsys)
#     end
#     Plot.trajectory_closed_loop!(ax, 1:2, contsys, symmodel, contr, x0, 10)
# end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
