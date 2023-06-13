module TestMain

using Test
using Dionysos
using CUDD

BDD = Dionysos.Utils.BDD
@testset "IntTupleSet" begin
    set = BDD.IntTupleSet{3, Int16}()
    @test isempty(set)
    @test sprint(show, MIME"text/plain"(), set) ==
          "Dionysos.Utils.BDD.IntTupleSet{3, Int16} with 0×0×0 bits"
    _list = [(1, 2, 3), (1, 5, 6), (1, 9, 8), (1, 6, 8), (4, 5, 6), (7, 6, 4), (1, 2, 3)]
    list = [Int16.(x) for x in _list]
    m = maximum(maximum(x) for x in list)
    for x in Iterators.product((1:m for i in 1:3)...)
        @test !(x ∈ set)
    end
    @test sprint(show, MIME"text/plain"(), set) ==
          "Dionysos.Utils.BDD.IntTupleSet{3, Int16} with 0×0×0 bits"
    @test collect(set) isa Vector{NTuple{3, Int16}}
    @test isempty(collect(set))
    set1 = BDD.IntTupleSet{3, Int16}()
    set2 = BDD.IntTupleSet{3, UInt16}()
    for x in list
        push!(set1, x)
        push!(set2, UInt16.(x))
        @test push!(set, x) === set
        for x in Iterators.product((1:m for i in 1:3)...)
            @test (x ∈ set) == (x ∈ set1) == (UInt16.(x) ∈ set2)
        end
    end
    for x in Iterators.product((1:m for i in 1:3)...)
        @test (x ∈ set) == (x ∈ list)
        @test (UInt16.(x) ∈ set2) == (x ∈ list)
    end
    for x in set
        @test (x ∈ list)
    end
    for x in set2
        @test eltype(x) == UInt16
        @test (Int16.(x) ∈ list)
    end
    @test collect(set) == sort!(unique!(list); by = x -> (x[3], x[2], x[1]))
    _rem_list = [(1, 2, 3), (4, 5, 6), (1, 3, 2), (1, 2)]
    rem_list = [Int16.(x) for x in _rem_list]
    for x in rem_list
        @test delete!(set, x) === set
    end
    setdiff!(list, rem_list)
    @test collect(set) == sort!(unique!(list); by = x -> (x[3], x[2], x[1]))
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
    @test sprint(show, MIME"text/plain"(), set) ==
          "Dionysos.Utils.BDD.IntTupleSet{3, Int16} with 3×4×4 bits"
    @test push!(set, Int16.((40, 1, 1))) === set
    @test sprint(show, MIME"text/plain"(), set) in [
        "Dionysos.Utils.BDD.IntTupleSet{3, Int16} with 6×4×4 bits", # Julia v1.0
        "Dionysos.Utils.BDD.IntTupleSet{3, Int16} with 6×4×4 bits:\n  (40, 1, 1)", # Julia v1.5
    ]
end

end  # module TestMain
