using Test
using Dionysos
BDD = Dionysos.Utils.BDD
@testset "BDD.BitSet" begin
    set = BDD.BitSet()
    @test isempty(set)
    @test sprint(show, MIME"text/plain"(), set) == "Dionysos.Utils.BDD.BitSet with 0 bits"
    list = [1, 10, 5, 1, 4, 8, 3, 4, 2, 14, 28, 13]
    m = maximum(list)
    for i in 0:m
        @test !(i in set)
    end
    @test sprint(show, MIME"text/plain"(), set) == "Dionysos.Utils.BDD.BitSet with 0 bits"
    @test collect(set) isa Vector{Int}
    @test isempty(collect(set))
    expected = BitSet()
    for l in list
        push!(expected, l)
        @test push!(set, l) === set
        for i in 0:m
            @test (i in set) == (i in expected)
        end
        @test collect(set) == collect(expected)
    end
    @test isempty(empty(set))
    @test !isempty(set)
    @test isempty(Base.emptymutable(set))
    @test !isempty(set)
    @test isempty(empty!(set))
    @test isempty(set)
    @test sprint(show, MIME"text/plain"(), set) == "Dionysos.Utils.BDD.BitSet with 5 bits"
    push!(set, 40)
    @test sprint(show, MIME"text/plain"(), set) in [
        "Dionysos.Utils.BDD.BitSet with 6 bits",       # Julia v1.0
        "Dionysos.Utils.BDD.BitSet with 6 bits:\n  40", # Julia v1.5
    ]
end
