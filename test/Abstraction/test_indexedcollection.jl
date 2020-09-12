include(joinpath(@__DIR__, "../../src/Abstraction/abstraction.jl"))

module TestMain

using Test
using Main.Abstraction
ABS = Main.Abstraction
using StaticArrays

println("")

@testset "Abstraction/indexedcollection" begin
coll = ABS.IndColl{String}()
@test push!(coll, "a") === coll
@test coll.e2i == Dict("a" => UInt(1))
@test coll.i2e == ["a"]
@test push!(coll, "b") === coll
@test coll.e2i == Dict(["a" => UInt(1), "b" => UInt(2)])
@test coll.i2e == ["a", "b"]
@test push!(coll, "a") === coll
@test coll.e2i == Dict(["a" => UInt(1), "b" => UInt(2)])
@test coll.i2e == ["a", "b"]
@test collect(ABS.get_elems(coll)) == ["a", "b"]
@test collect(ABS.get_indexes(coll)) == [1, 2]
@test length(coll) === 2
@test append!(coll, ["e", "f"]) === coll
@test length(coll) === 4
print("")
end

end  # module TestMain
