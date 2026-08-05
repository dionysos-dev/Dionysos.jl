module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

@testset "INCL_MODE enum" begin
    @test UT.INNER isa UT.INCL_MODE
    @test UT.OUTER isa UT.INCL_MODE
    @test UT.CENTER isa UT.INCL_MODE
    @test length(instances(UT.INCL_MODE)) == 3
    # the three modes are distinct
    @test allunique((UT.INNER, UT.OUTER, UT.CENTER))
end

@testset "invert_incl_mode (soundness inverts the hole mode)" begin
    # Cutting a hole B out of A: to stay sound, an INNER cover of A needs an
    # OUTER cover of the hole, and vice versa; CENTER has no dual so it maps to
    # OUTER as well.
    @test UT.invert_incl_mode(UT.INNER) === UT.OUTER
    @test UT.invert_incl_mode(UT.OUTER) === UT.INNER
    @test UT.invert_incl_mode(UT.CENTER) === UT.OUTER

    # INNER/OUTER form an involution; CENTER does not.
    @test UT.invert_incl_mode(UT.invert_incl_mode(UT.INNER)) === UT.INNER
    @test UT.invert_incl_mode(UT.invert_incl_mode(UT.OUTER)) === UT.OUTER
end

end # module
