module TestTreeAndSortedTupleSet

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

@testset "SortedTupleSet" begin
    @testset "constructor and basic push" begin
        s = UT.SortedTupleSet{2, Tuple{Int, Int}}()

        @test length(s) == 0
        @test s.is_sorted == true
        @test UT.get_data(s) == Tuple{Int, Int}[]

        UT.push_new!(s, (2, 20))
        @test length(s) == 1
        @test s.is_sorted == false
        @test UT.get_data(s) == [(2, 20)]

        UT.push_new!(s, (1, 10))
        @test length(s) == 2
        @test s.is_sorted == false
    end

    @testset "ensure_sorted!" begin
        s = UT.SortedTupleSet{2, Tuple{Int, Int}}()
        UT.append_new!(s, [(3, 30), (1, 10), (2, 20)])

        @test s.is_sorted == false

        UT.ensure_sorted!(s)

        @test s.is_sorted == true
        @test UT.get_data(s) == [(1, 10), (2, 20), (3, 30)]
    end

    @testset "delete! with custom comparison" begin
        s = UT.SortedTupleSet{2, Tuple{Int, Int}}()
        UT.append_new!(s, [(1, 10), (2, 20), (2, 99), (3, 30)])

        Base.delete!(s, (2, 0), (a, b) -> a[1] == b[1])

        @test length(s) == 2
        @test Set(UT.get_data(s)) == Set([(1, 10), (3, 30)])
    end

    @testset "empty!" begin
        s = UT.SortedTupleSet{2, Tuple{Int, Int}}()
        UT.append_new!(s, [(1, 10), (2, 20)])

        out = empty!(s)

        @test out == true
        @test length(s) == 0
        @test s.is_sorted == true
        @test isempty(UT.get_data(s))
    end

    @testset "fix_and_eliminate_first on pairs" begin
        s = UT.SortedTupleSet{2, Tuple{Int, Int}}()
        UT.append_new!(s, [(2, 20), (1, 10), (1, 11), (3, 30)])

        vals = collect(UT.fix_and_eliminate_first(s, 1))

        @test s.is_sorted == true
        @test vals == [(10,), (11,)]
    end

    @testset "fix_and_eliminate_first on triples" begin
        s = UT.SortedTupleSet{3, Tuple{Int, Int, Int}}()
        UT.append_new!(s, [(2, 20, 200), (1, 10, 100), (1, 11, 101), (3, 30, 300)])

        vals = collect(UT.fix_and_eliminate_first(s, 1))

        @test vals == [(10, 100), (11, 101)]
    end

    @testset "fix_and_eliminate_tail!" begin
        s = UT.SortedTupleSet{3, Tuple{Int, Int, Int}}()
        UT.append_new!(s, [(1, 4, 5), (2, 4, 5), (3, 7, 8), (9, 4, 5)])

        out = Int[]
        UT.fix_and_eliminate_tail!(out, s, (4, 5))

        @test out == [1, 2, 9]
    end

    @testset "append_new! on Base.Set" begin
        s = Set([1, 2])
        UT.append_new!(s, [2, 3, 4])

        @test s == Set([1, 2, 3, 4])
    end
end

end # module
