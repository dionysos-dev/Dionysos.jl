module BDD
module TestMain

    using Test
    using Dionysos
    using CUDD

    BDD = Dionysos.Utils.BDD
    using Test

    @testset "BitSet" begin
        # create an empty set
        set = BDD.BitSet()
        @test isempty(set)
        @test length(set) == 0

        # add 1 and 2 to the list
        @test push!(set, 1) === set
        @test !(isempty(set))
        @test length(set) == 1
        @test 1 in set

        @test push!(set, 2) === set
        @test length(set) == 2
        @test 2 in set

        # test that 3 and 4 are not in the list
        @test !(3 in set)
        @test !(4 in set)

        # delete 1 from list
        @test delete!(set, 1) === set
        @test 1 ∉ set
        @test length(set) == 1

        # set list to be empty again
        empty!(set)
        @test isempty(set)
        @test length(set) == 0

        # put elements 1 and 2 in the list again
        set = BDD.BitSet()
        push!(set, 1)
        push!(set, 2)
        elements = [x for x in set]
        @test elements == [1, 2]

        # test inclusion
        @test 1 in set
        @test 2 in set
        @test 3 ∉ set
        @test 4 ∉ set
    end

end #end module BDD
end #end module TestMain
