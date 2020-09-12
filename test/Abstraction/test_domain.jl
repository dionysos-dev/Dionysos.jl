include(joinpath(@__DIR__, "../../src/Abstraction/abstraction.jl"))

module TestMain

using Test
using Main.Abstraction
ABS = Main.Abstraction
using StaticArrays

println("")

@testset "Abstraction/domain" begin
xorig = SVector(0.0, 0.0)
xh = SVector(1.0, 2.0)
uorig = SVector(0)
uh = SVector(1)
mng = ABS.ListGridManager(xorig, xh, uorig, uh)

INT_MAX = typemax(Int)
INT_MIN = typemin(Int)
@test !ABS._compare_size(ABS.HyperRectangle(SVector(0, 0), SVector(100, 100)), 1000)
@test ABS._compare_size(ABS.HyperRectangle(SVector(0, 0), SVector(100, 100)), 11_000)
@test !ABS._compare_size(ABS.HyperRectangle(SVector(0, INT_MIN + 1), SVector(100, 100)), 11_000)

XDom1 = ABS.AddXDomain!(mng)
XCT = ABS.xcelltype(mng)
xcell = XCT((1, 2))
ABS.add_cell!(mng, XDom1, xcell)
@test ABS.get_ncells(mng, XDom1) === 1
@test ABS.get_some_cell(mng, XDom1) === xcell
@test empty!(mng, XDom1) === XDom1
@test ABS.get_ncells(mng, XDom1) === 0
@test ABS.get_some_cell(mng, XDom1) === nothing
rect = ABS.HyperRectangle(SVector(-2.0, -4.0), SVector(2.0, 4.0))
ABS.add_cells!(mng, XDom1, rect, ABS.INNER)
ABS.add_cells!(mng, XDom1, rect, ABS.INNER)
@test ABS.get_ncells(mng, XDom1) === 9
XDom2 = ABS.AddXDomain!(mng)
rect = ABS.HyperRectangle(SVector(-2000.0, -4000.0), SVector(10.0, 10.0))
ABS.add_cells!(mng, XDom2, XDom1, rect, ABS.INNER)
@test ABS.get_ncells(mng, XDom1) === 9
rect = ABS.HyperRectangle(SVector(-20.0, -Inf), SVector(Inf, 10.0))
ABS.add_cells!(mng, XDom2, XDom1, rect, ABS.INNER)
@test ABS.get_ncells(mng, XDom1) === 9
celllist = Set(XCT((a, b)) for a = -1:1, b = -1:1)
@test Set(ABS.enum_cells(mng, XDom2)) == celllist
rect = ABS.HyperRectangle(SVector(-Inf, -Inf), SVector(0.0, 10.0))
ABS.remove_cells!(mng, XDom2, rect, ABS.INNER)
@test ABS.get_ncells(mng, XDom2) === 6

UDom1 = ABS.AddUDomain!(mng)
UCT = ABS.ucelltype(mng)
ucell = UCT((1,))
rect = ABS.HyperRectangle(SVector(-5), SVector(2))
ABS.add_cells!(mng, UDom1, rect, ABS.INNER)
celllist = Set(UCT((a,)) for a = -4:1)
@test Set(ABS.enum_cells(mng, UDom1)) == celllist
print("")
end

end  # module TestMain
