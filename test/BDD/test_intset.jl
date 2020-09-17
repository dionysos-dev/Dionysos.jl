include(joinpath(@__DIR__, "../../src/BDD/BDD.jl"))
using Test
using CUDD

println("")

@testset "BDD/IntSet" begin
    set = BDD.IntSet()
    @test isempty(set)
    @test sprint(show, MIME"text/plain"(), set) == "Main.BDD.IntSet{Int64} with 0 bits"
    list = [1, 10, 5, 1, 4, 8, 3, 4, 2, 14, 28, 13]
    m = maximum(list)
    for x in 0:m
        @test x ∉ set
    end
    @test sprint(show, MIME"text/plain"(), set) == "Main.BDD.IntSet{Int64} with 0 bits"
    @test collect(set) isa Vector{Int}
    @test isempty(collect(set))
    set1 = BDD.IntSet{UInt64}()
    for x in list
        push!(set1, UInt64(x))
        @test push!(set, x) === set
        for x in 0:m
            @test (x ∈ set) === (UInt64(x) ∈ set1)
            @test x ∉ set1
        end
        @test collect(set) == collect(set1)
    end
    @test isempty(empty(set))
    @test !isempty(set)
    @test isempty(Base.emptymutable(set))
    @test !isempty(set)
    @test isempty(empty!(set))
    @test isempty(set)
    @test sprint(show, MIME"text/plain"(), set) == "Main.BDD.IntSet{Int64} with 5 bits"
    push!(set, 40)
    @test sprint(show, MIME"text/plain"(), set) in [
        "Main.BDD.IntSet{Int64} with 6 bits",       # Julia v1.0
        "Main.BDD.IntSet{Int64} with 6 bits:\n  40" # Julia v1.5
    ]
    push!(set, 42)
    @test Set(set) == Set([40, 42])
    delete!(set, 41)
    @test Set(set) == Set([40, 42])
    delete!(set, 42)
    @test Set(set) == Set([40])
    delete!(set, 40)
    @test isempty(set)
    set = 6
    print("")
end
