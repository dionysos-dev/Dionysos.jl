module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

@testset "periodic wrap_coord" begin
    periodic_dims = SVector(1)
    periods = SVector(2.0)              # wrap dim 1 into [start, start + 2)

    # default start = 0  =>  domain [0, 2)
    @test UT.wrap_coord(SVector(2.5), periodic_dims, periods)[1] ≈ 0.5
    @test UT.wrap_coord(SVector(-0.5), periodic_dims, periods)[1] ≈ 1.5
    @test UT.wrap_coord(SVector(1.0), periodic_dims, periods)[1] ≈ 1.0   # already inside

    # offset start = -1  =>  domain [-1, 1)
    w = UT.wrap_coord(SVector(3.5), periodic_dims, periods; start = SVector(-1.0))
    @test w[1] ≈ -0.5
    @test -1.0 <= w[1] < 1.0

    # only the periodic dimension is wrapped; other dims pass through
    x2 = SVector(2.5, 5.0)
    w2 = UT.wrap_coord(x2, periodic_dims, periods)
    @test w2[1] ≈ 0.5
    @test w2[2] ≈ 5.0

    # the closure form agrees with the direct call
    f = UT.get_periodic_wrapper(periodic_dims, periods)
    @test f(SVector(2.5))[1] ≈ 0.5
end

end # module TestMain
