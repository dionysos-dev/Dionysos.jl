include(joinpath(@__DIR__, "../../src/BDD/BDD.jl"))

module TestMain

using Test
using Main.BDD
using CUDD

println("")

@testset "IntTupleSet" begin
    set = BDD.IntTupleSet{3}()
    @test isempty(set)
    @test sprint(show, MIME"text/plain"(), set) == "Main.BDD.IntTupleSet{3,Int64} with 0×0×0 bits"
    list = [(1, 2, 3), (1, 5, 6), (1, 9, 8), (1, 6, 8), (4, 5, 6), (7, 6, 4), (1, 2, 3)]
    m = maximum(maximum(x) for x in list)
    for x in Iterators.product((1:m for i in 1:3)...)
        @test x ∉ set
    end
    @test sprint(show, MIME"text/plain"(), set) == "Main.BDD.IntTupleSet{3,Int64} with 0×0×0 bits"
    @test collect(set) isa Vector{NTuple{3,Int64}}
    @test isempty(collect(set))
    set1 = BDD.IntTupleSet{3,Int16}()
    set2 = BDD.IntTupleSet{3,UInt16}()
    for x in list
        push!(set1, Int16.(x))
        push!(set2, UInt16.(x))
        @test push!(set, x) === set
        for x in Iterators.product((1:m for i in 1:3)...)
            @test (x ∈ set) === (Int16.(x) ∈ set1) === (UInt16.(x) ∈ set2)
            @test (x ∈ set1) === (x ∈ set2) === false
        end
    end
    for x in Iterators.product((1:m for i in 1:3)...)
        @test (x ∈ set) === (x ∈ list)
        @test (UInt16.(x) ∈ set2) === (x ∈ list)
    end
    for x in set
        @test x ∈ list
    end
    for x in set2
        @test eltype(x) === UInt16
        @test Int.(x) ∈ list
    end
    @test collect(set) == sort!(unique!(list), by = x -> (x[3], x[2], x[1]))
    rem_list = [(1, 2, 3), (4, 5, 6), (1, 3, 2), (1, 2)]
    for x in rem_list
        @test delete!(set, x) === set
    end
    setdiff!(list, rem_list)
    @test collect(set) == sort!(unique!(list), by = x -> (x[3], x[2], x[1]))
    for x in rem_list
        @test delete!(set, x) === set
    end
    @test delete!(set, "hello") === set
    @test collect(set) == list
    @test isempty(empty(set))
    @test !isempty(set)
    @test isempty(Base.emptymutable(set))
    @test !isempty(set)
    @test isempty(empty!(set))
    @test isempty(set)
    @test sprint(show, MIME"text/plain"(), set) == "Main.BDD.IntTupleSet{3,Int64} with 3×4×4 bits"
    @test push!(set, (40, 1, 1)) === set
    @test sprint(show, MIME"text/plain"(), set) in [
        "Main.BDD.IntTupleSet{3,Int64} with 6×4×4 bits", # Julia v1.0
        "Main.BDD.IntTupleSet{3,Int64} with 6×4×4 bits:\n  (40, 1, 1)" # Julia v1.5
    ]
    print("")
end

end  # module TestMain
