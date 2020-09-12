include(joinpath(@__DIR__, "../../src/Abstraction/abstraction.jl"))

module TestMain

using Test
using Main.Abstraction
ABS = Main.Abstraction
using StaticArrays

println("")

@testset "Abstraction/automaton-growth" begin
xorig = SVector(0.0, 0.0)
xh = SVector(1.0, 2.0)
uorig = SVector(0.0)
uh = SVector(0.5)
mng = ABS.ListGridManager(xorig, xh, uorig, uh)

XDom = ABS.AddXDomain!(mng)
rect = ABS.HyperRectangle(SVector(0.0, 0.0), SVector(10.0, 11.0))
ABS.add_cells!(mng, XDom, rect, ABS.OUTER)
ABS.create_symbols!(mng, XDom)
stateset = ABS.AddStateSet!(mng)
ABS.add_symbols!(mng, stateset, XDom)

UDom = ABS.AddUDomain!(mng)
rect = rect = ABS.HyperRectangle(SVector(-1.0), SVector(1.0))
ABS.add_cells!(mng, UDom, rect, ABS.OUTER)
ABS.create_symbols!(mng, UDom)
labelset = ABS.AddLabelSet!(mng)
ABS.add_symbols!(mng, labelset, UDom)

tstep = 5.0
nsys = 3
ngrowthbound = 3
F_sys(x, u) = SVector(u[1], -cos(x[1]))
L_growthbound(u) = SMatrix{2,2}(0.0, 1.0, 0.0, 0.0)
sysnoise = SVector(1.0, 1.0)*0.1
measnoise = SVector(1.0, 1.0)*0.0

contsys = ABS.ControlSystemGrowthRK4(
    tstep,
    F_sys, L_growthbound,
    sysnoise, measnoise,
    nsys, ngrowthbound)

TT = ABS.transitiontype(mng)
autom = ABS.AddAutomaton!(mng)
ABS.add_transitions!(mng, autom, stateset, labelset, stateset, contsys)
@test ABS.get_ntransitions(mng, autom) === 1145

targetset = ABS.AddStateSet!(mng)
# XCT = ABS.xcelltype(mng)
# UCT = ABS.ucelltype(mng)
get_xcell = ABS.coord2cell(mng.XDisc)
get_ucell = ABS.coord2cell(mng.UDisc)
get_state = ABS.cell2symbol(mng.XMap)
xcell = get_xcell(SVector(1.0, 4.0))
ucell = get_ucell(SVector(0.5))
XDom = ABS.AddXDomain!(mng)
UDom = ABS.AddUDomain!(mng)
ABS.add_cell!(mng, XDom, xcell)
ABS.add_cell!(mng, UDom, ucell)
sourceset = ABS.AddStateSet!(mng)
labelset = ABS.AddLabelSet!(mng)
ABS.add_symbols!(mng, sourceset, XDom)
ABS.add_symbols!(mng, labelset, UDom)
ABS.add_symbols!(mng, targetset, autom, sourceset, labelset)
coordlist = (SVector(Float64.(xi)) for xi in Iterators.product(3:5, 2:12))
celllist = Set(get_xcell(x) for x in coordlist)
symbollist = Set(get_state(cell) for cell in celllist)
@test Set(ABS.enum_symbols(mng, targetset)) == symbollist
print("")
end

end  # module TestMain
