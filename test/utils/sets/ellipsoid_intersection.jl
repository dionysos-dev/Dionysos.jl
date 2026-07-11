module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import LazySets

@testset "ellipsoid intersection / disjointness" begin
    # Unit balls (shape matrix Q = I) of radius 1.
    Q = SMatrix{2, 2}(1.0, 0.0, 0.0, 1.0)
    E0 = LazySets.Ellipsoid(SVector(0.0, 0.0), Q)
    E_far = LazySets.Ellipsoid(SVector(3.0, 0.0), Q)   # centers 3 apart > 1+1 -> disjoint
    E_near = LazySets.Ellipsoid(SVector(1.5, 0.0), Q)  # centers 1.5 apart < 2 -> overlap
    E_touch = LazySets.Ellipsoid(SVector(0.2, 0.0), Q) # concentric-ish -> overlap

    @test UT.is_disjoint(E0, E_far)
    @test !UT.is_disjoint(E0, E_near)
    @test !UT.is_disjoint(E0, E_touch)

    # Symmetry of the predicate.
    @test UT.is_disjoint(E_far, E0) == UT.is_disjoint(E0, E_far)

    # A witness scaling from get_ℓ_ast_intersect exists for overlapping ellipsoids.
    gstar, λstar = UT.get_ℓ_ast_intersect(E0, E_near)
    @test isfinite(λstar)
end

end # module TestMain
