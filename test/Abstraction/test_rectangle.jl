include(joinpath(@__DIR__, "../../src/Abstraction/abstraction.jl"))

module TestMain

using Test
using Main.Abstraction
ABS = Main.Abstraction
using StaticArrays

println("")

@testset "Abstraction/rectangle" begin
rectI1 = ABS.HyperRectangle(SVector(0, 1), SVector(1, 2))
rectI2 = ABS.HyperRectangle(SVector(10, 1), SVector(5, 2))
@test !isempty(rectI1)
@test intersect(rectI1, rectI1) === rectI1
@test isempty(rectI2)
@test isempty(intersect(rectI1, rectI2))
@test rectI1 ⊆ rectI1
@test rectI2 ⊆ rectI1
print("")
end

end  # module TestMain
