include(joinpath(@__DIR__, "../../src/Abstraction/abstraction.jl"))

module TestMain

using Test
using Main.Abstraction
ABS = Main.Abstraction
using StaticArrays

println("")

@testset "Abstraction/symbolic" begin
xorig = SVector(0.0, 0.0)
xh = SVector(1.0, 2.0)
uorig = SVector(0)
uh = SVector(1)
mng = ABS.ListGridManager(xorig, xh, uorig, uh)

XDom1 = ABS.AddXDomain!(mng)
rect = ABS.HyperRectangle(SVector(-2.0, -4.0), SVector(2.0, 4.0))
ABS.add_cells!(mng, XDom1, rect, ABS.INNER)
XCT = ABS.xcelltype(mng)

ABS.create_symbols!(mng, XDom1)
celllist = Set(XCT((a, b)) for a = -1:1, b = -1:1)
@test Set(ABS.enum_states(mng)) == celllist
XDom2 = ABS.AddXDomain!(mng)
rect = ABS.HyperRectangle(SVector(0.0, 0.0), SVector(40.0, 40.0))
ABS.add_cells!(mng, XDom2, rect, ABS.INNER)
stateset = ABS.AddStateSet!(mng)
symb = ABS.get_some_symbol(mng, mng.XMap)
@test symb !== nothing
ABS.add_symbol!(mng, stateset, symb)
@test ABS.get_nsymbols(mng, stateset) === 1
@test empty!(mng, stateset) === stateset
@test ABS.get_nsymbols(mng, stateset) === 0
ABS.add_symbols!(mng, stateset, XDom2)
@test ABS.get_nsymbols(mng, stateset) === 1
XDom3 = ABS.AddXDomain!(mng)
ABS.add_cells!(mng, XDom3, stateset)
@test Set(ABS.enum_cells(mng, XDom3)) == Set([XCT((1, 1))])

UDom1 = ABS.AddUDomain!(mng)
rect = ABS.HyperRectangle(SVector(-5), SVector(2))
ABS.add_cells!(mng, UDom1, rect, ABS.INNER)
UCT = ABS.ucelltype(mng)

ABS.create_symbols!(mng, UDom1)
celllist = Set(UCT((a,)) for a = -4:1)
@test Set(ABS.enum_labels(mng)) == celllist
UDom2 = ABS.AddUDomain!(mng)
rect = ABS.HyperRectangle(SVector(0), SVector(0))
ABS.add_cells!(mng, UDom2, rect, ABS.OUTER)
labelset = ABS.AddLabelSet!(mng)
ABS.add_symbols!(mng, labelset, UDom2)
@test ABS.get_nsymbols(mng, labelset) === 1
UDom3 = ABS.AddUDomain!(mng)
ABS.add_cells!(mng, UDom3, labelset)
@test Set(ABS.enum_cells(mng, UDom3)) == Set([UCT((0,))])
print("")
end

end  # module TestMain
