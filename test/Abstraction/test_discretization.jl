include(joinpath(@__DIR__, "../../src/Abstraction/abstraction.jl"))

module TestMain

using Test
using Main.Abstraction
ABS = Main.Abstraction
using StaticArrays

ListGridDisc = (orig, h, posmin) -> ABS.ListGridDisc(orig, h)
BDDGridDisc = (orig, h, posmin) -> ABS.BDDGridDisc(orig, h, posmin)

println("")

for GridDisc in [ListGridDisc, BDDGridDisc]
@testset "Abstraction/discretization" begin
orig = SVector(0.0, 0.0)
h = SVector(1.0, 2.0)
posmin = SVector(-5, -5)
disc = GridDisc(orig, h, posmin)
x = SVector(3.14, 15.9)
coord2cell = ABS.coord2cell(disc)
cell2pos = ABS.cell2pos(disc)
pos2cell = ABS.pos2cell(disc)
cell2coord = ABS.cell2coord(disc)
coord2pos = ABS.coord2pos(disc)
cell = coord2cell(x)
# Main.@code_warntype hash(cell, UInt(0))
@test coord2pos(x) === SVector(3, 8)
@test cell === pos2cell(SVector(3, 8))
@test cell2pos(cell) === SVector(3, 8)
@test pos2cell((3, 8)) === cell
@test cell2coord(cell) === SVector(3.0, 16.0)

coord2pos_set = ABS.coord2pos_set(disc)
INT_MAX = typemax(Int)
INT_MIN = typemin(Int)
rect = ABS.HyperRectangle(SVector(0.1, 0.1), SVector(1.1, 1.1))
@test coord2pos_set(rect, ABS.INNER) ===
    ABS.HyperRectangle(SVector(1, 1), SVector(0, 0))
@test coord2pos_set(rect, ABS.OUTER) ===
    ABS.HyperRectangle(SVector(0, 0), SVector(1, 1))
rect = ABS.HyperRectangle(SVector(-1e100, 2.71), SVector(Inf, Inf))
@test coord2pos_set(rect, ABS.INNER) ===
    ABS.HyperRectangle(SVector(INT_MIN, 2), SVector(INT_MAX, INT_MAX))
@test coord2pos_set(rect, ABS.OUTER) ===
    ABS.HyperRectangle(SVector(INT_MIN, 1), SVector(INT_MAX, INT_MAX))
print("")
end
end

end  # module TestMain
