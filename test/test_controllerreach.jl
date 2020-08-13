include("../src/abstraction.jl")

module TestMain

using Test
using StaticArrays
using Main.Abstraction
AB = Main.Abstraction

sleep(0.1) # used for good printing
println("Started test")

@testset "ControllerReach" begin
lb = SVector(-5.0, -5.0)
ub = SVector(5.0, 5.0)
x0 = SVector(0.0, 0.0)
h = SVector(0.47, 0.23)
Xgrid = AB.NewGrid(x0, h)
Xfull = AB.NewDomainList(Xgrid)
AB.add_set!(Xfull, AB.HyperRectangle(lb, ub), AB.OUTER)
AB.remove_set!(Xfull, AB.HyperRectangle((-1.0, -2.0), (-1.1, 4.0)), AB.OUTER)

lb = SVector(-2.0)
ub = SVector(2.0)
u0 = SVector(0.0)
h = SVector(1.0)
Ugrid = AB.NewGrid(u0, h)
Ufull = AB.NewDomainList(Ugrid)
AB.add_set!(Ufull, AB.HyperRectangle(lb, ub), AB.OUTER)

tstep = 1.0
nsys = 3
ngrowthbound = 3
F_sys(x, u) = SVector(1.0, u[1])
sysnoise = SVector(1.0, 1.0)*0.001
measnoise = SVector(1.0, 1.0)*0.001
L_growthbound(u) = SMatrix{2,2}(0.0, 0.0, 0.0, 0.0)

contsys = AB.NewControlSystemGrowthRK4(
    tstep, F_sys, L_growthbound, sysnoise, measnoise, nsys, ngrowthbound)
symmodel = AB.NewSymbolicModelListList(Xfull, Ufull)
AB.compute_symmodel_from_controlsystem!(symmodel, contsys)
display(symmodel)

Xinit = AB.NewDomainList(Xgrid)
AB.add_set!(Xinit, AB.HyperRectangle(SVector(-3.0, -3.0), SVector(-2.9, -2.9)), AB.OUTER)
initlist = Int[]
for pos in AB.enum_pos(Xinit)
    push!(initlist, AB.get_state_by_xpos(symmodel, pos))
end
Xtarget = AB.NewDomainList(Xgrid)
AB.add_set!(Xtarget, AB.HyperRectangle(SVector(0.0, 0.0), SVector(4.0, 4.0)), AB.OUTER)
targetlist = Int[]
for pos in AB.enum_pos(Xtarget)
    push!(targetlist, AB.get_state_by_xpos(symmodel, pos))
end

contr = AB.NewControllerList()
AB.compute_controller_reach!(contr, symmodel.autom, initlist, targetlist)
@test AB.get_npairs(contr) == 836

xpos = AB.get_somepos(Xinit)
x0 = AB.get_coord_by_pos(Xgrid, xpos)
Xsimple = AB.NewDomainList(Xgrid)
XUYsimple_ = Any[]

for i = 1:6
    source = AB.get_state_by_xpos(symmodel, xpos)
    Xs = AB.NewDomainList(Xgrid)
    Ys = AB.NewDomainList(Xgrid)
    Us = AB.NewDomainList(Ugrid)
    AB.add_pos!(Xs, xpos)
    symbollist = Int[]
    targetlist = Int[]
    AB.compute_enabled_symbols!(symbollist, contr, source)
    for symbol in symbollist
        AB.compute_post!(targetlist, symmodel.autom, source, symbol)
        AB.add_pos!(Us, AB.get_upos_by_symbol(symmodel, symbol))
    end
    for target in targetlist
        AB.add_pos!(Ys, AB.get_xpos_by_state(symmodel, target))
    end
    push!(XUYsimple_, (Xs, Us, Ys))
    if !isempty(Ys)
        xpos = AB.get_somepos(Ys)
    end
end

@static if get(ENV, "TRAVIS", "false") == "false"
    include("../src/plotting.jl")
    using PyPlot
    fig = PyPlot.figure()
    ax = fig.gca()
    ax.set_xlim((-5.5, 5.5))
    ax.set_ylim((-5.3, 5.3))
    Plot.domain!(ax, 1:2, Xfull, fa = 0.1)
    Plot.domain!(ax, 1:2, Xinit)
    Plot.domain!(ax, 1:2, Xtarget)
    for (Xs, Us, Ys) in XUYsimple_
        Plot.domain!(ax, 1:2, Xs, fc = "green")
        Plot.domain!(ax, 1:2, Ys, fc = "blue")
        Plot.cell_image!(ax, 1:2, Xs, Us, contsys)
        Plot.cell_approx!(ax, 1:2, Xs, Us, contsys)
    end
    Plot.trajectory_closed_loop!(ax, 1:2, contsys, symmodel, contr, x0, 10)
end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
