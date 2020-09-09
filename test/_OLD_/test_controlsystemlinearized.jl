include("../src/abstraction.jl")

module TestMain

using Test
using StaticArrays
using Main.Abstraction
AB = Main.Abstraction

sleep(0.1) # used for good printing
println("Started test")

@testset "ControlSystemLinearized" begin
lb = SVector(0.0, 0.0)
ub = SVector(10.0, 11.0)
x0 = SVector(0.0, 0.0)
h = SVector(1.0, 2.0)
Xgrid = AB.GridFree(x0, h)
Xfull = AB.DomainList(Xgrid)
@test AB.get_ncells(Xfull) == 0
AB.add_set!(Xfull, AB.HyperRectangle(lb, ub), AB.OUTER)
@test AB.get_ncells(Xfull) == 77

lb = SVector(-1.0)
ub = SVector(1.0)
u0 = SVector(0.0)
h = SVector(0.5)
Ugrid = AB.GridFree(u0, h)
Ufull = AB.DomainList(Ugrid)
AB.add_set!(Ufull, AB.HyperRectangle(lb, ub), AB.OUTER)
@test AB.get_ncells(Ufull) == 5

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

contsys = AB.NewControlSystemLinearizedRK4(
    tstep, F_sys, DF_sys, bound_DF, bound_DDF, measnoise, nsys)

xpos = (1, 1)
upos = (1,)
x = AB.get_coord_by_pos(Xgrid, xpos)
u = AB.get_coord_by_pos(Ugrid, upos)
Xsimple = AB.DomainList(Xgrid)
AB.add_pos!(Xsimple, xpos)
@test AB.get_ncells(Xsimple) == 1
Usimple = AB.DomainList(Ugrid)
AB.add_pos!(Usimple, upos)
@test AB.get_ncells(Usimple) == 1

@static if get(ENV, "TRAVIS", "false") == "false"
    include("../src/plotting.jl")
    using PyPlot
    fig = PyPlot.figure()
    ax = fig.gca()
    ax.set_xlim((-10.0, 10.0))
    ax.set_ylim((-10.0, 10.0))
    Plot.domain!(ax, 1:2, Xsimple)
    Plot.trajectory_open_loop!(ax, 1:2, contsys, x, u, 50)
    Plot.cell_image!(ax, 1:2, Xsimple, Usimple, contsys)
    Plot.cell_approx!(ax, 1:2, Xsimple, Usimple, contsys)
end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
