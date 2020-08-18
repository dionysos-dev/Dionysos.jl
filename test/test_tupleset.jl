include("../src/tupleset.jl")

module TestMain

using Test
using Main.BDD

@testset "TupleUIntSet" begin
    dim = 3
    set = BDD.TupleUIntSet(dim)
    @test isempty(set)
    _list = [(1, 2, 3), (1, 5, 6), (1, 9, 8), (1, 6, 8), (4, 5, 6), (7, 6, 4)]
    list = [UInt.(x) for x in _list]
    m = maximum(maximum(x) for x in list)
    for i in Iterators.product((1:m for i in 1:dim)...)
        @test !(i in set)
    end
    for x in list
        push!(set, x)
        # @test push!(set, l) === set
        # for i in 0:m
        #     @test (i in set) == (i in expected)
        # end
        # @test collect(set) == collect(expected)
    end
    @test isempty(empty(set))
    @test !isempty(set)
    @test isempty(Base.emptymutable(set))
    @test !isempty(set)
    @test isempty(empty!(set))
    @test isempty(set)
    @test sprint(show, MIME"text/plain"(), set) == "Main.BDD.BitSet with 5 bits"
    push!(set, 40)
    @test sprint(show, MIME"text/plain"(), set) in [
        "Main.BDD.BitSet with 6 bits",       # Julia v1.0
        "Main.BDD.BitSet with 6 bits:\n  40" # Julia v1.5
    ]
end

end  # module TestMain
