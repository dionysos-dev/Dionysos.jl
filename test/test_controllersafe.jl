include("../src/abstraction.jl")

module TestMain

using Test
using StaticArrays
using Main.Abstraction
AB = Main.Abstraction

sleep(0.1) # used for good printing
println("Started test")

@testset "ControllerSafe" begin
lb = SVector(-5.0, -5.0)
ub = SVector(5.0, 5.0)
x0 = SVector(0.0, 0.0)
h = SVector(0.47, 0.23)
Xgrid = AB.GridFree(x0, h)
Xfull = AB.DomainList(Xgrid)
AB.add_set!(Xfull, AB.HyperRectangle(lb, ub), AB.OUTER)

lb = SVector(-4.0)
ub = SVector(4.0)
u0 = SVector(0.0)
h = SVector(0.5)
Ugrid = AB.GridFree(u0, h)
Ufull = AB.DomainList(Ugrid)
AB.add_set!(Ufull, AB.HyperRectangle(lb, ub), AB.OUTER)

tstep = 0.2
nsys = 3
ngrowthbound = 3
F_sys(x, u) = SVector(u[1], -x[2] + u[1])
sysnoise = SVector(1.0, 1.0)*0.001
measnoise = SVector(1.0, 1.0)*0.001
L_growthbound(u) = SMatrix{2,2}(0.0, 0.0, 0.0, -1.0)

contsys = AB.NewControlSystemGrowthRK4(
    tstep, F_sys, L_growthbound, sysnoise, measnoise, nsys, ngrowthbound)
symmodel = AB.NewSymbolicModelListList(Xfull, Ufull)
AB.compute_symmodel_from_controlsystem!(symmodel, contsys)
@test AB.get_ntrans(symmodel.autom) == 62077

Xinit = AB.DomainList(Xgrid)
AB.add_set!(Xinit, AB.HyperRectangle(SVector(-3.0, -3.0), SVector(-2.9, -2.9)), AB.OUTER)
initlist = Int[]
for pos in AB.enum_pos(Xinit)
    push!(initlist, AB.get_state_by_xpos(symmodel, pos))
end
Xsafe = AB.DomainList(Xgrid)
union!(Xsafe, Xfull)
AB.remove_set!(Xsafe, AB.HyperRectangle(SVector(-1.0, -2.0), SVector(-1.1, 4.0)), AB.OUTER)
safelist = Int[]
for pos in AB.enum_pos(Xsafe)
    push!(safelist, AB.get_state_by_xpos(symmodel, pos))
end

contr = AB.NewControllerList()
AB.compute_controller_safe!(contr, symmodel.autom, initlist, safelist)
@test AB.get_npairs(contr) == 15043

invlist = Int[]
for source in 1:symmodel.autom.nstates
    symbollist = Int[]
    AB.compute_enabled_symbols!(symbollist, contr, source)
    if !isempty(symbollist)
        push!(invlist, source)
    end
end
Xinv = AB.DomainList(Xgrid)
Yinv = AB.DomainList(Xgrid)
correct = true
for source in invlist
    AB.add_pos!(Xinv, AB.get_xpos_by_state(symmodel, source))
    if !correct
        break
    end
    symbollist = Int[]
    targetlist = Int[]
    AB.compute_enabled_symbols!(symbollist, contr, source)
    for symbol in symbollist
        AB.compute_post!(targetlist, symmodel.autom, source, symbol)
    end
    for target in targetlist
        AB.add_pos!(Yinv, AB.get_xpos_by_state(symmodel, target))
    end
    correct = correct && targetlist ⊆ safelist
end
@test correct

xpos = AB.get_somepos(Xinit)
x0 = AB.get_coord_by_pos(Xgrid, xpos)

@static if get(ENV, "TRAVIS", "false") == "false"
    include("../src/plotting.jl")
    using PyPlot
    fig = PyPlot.figure()
    ax = fig.gca()
    ax.set_xlim((-5.5, 5.5))
    ax.set_ylim((-5.3, 5.3))
    Plot.domain!(ax, 1:2, Xfull, fa = 0.0)
    Plot.domain!(ax, 1:2, Xinit)
    Plot.domain!(ax, 1:2, Xsafe, fa = 0.1)
    Plot.domain!(ax, 1:2, Xinv, fa = 0.1, fc = "yellow")
    Plot.domain!(ax, 1:2, Yinv; fc = "blue")
    Plot.trajectory_closed_loop!(
        ax, 1:2, contsys, symmodel, contr, x0, 100, randchoose = true)
end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
