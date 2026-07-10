module TestMain

using Test
using Dionysos
import LinearAlgebra as LA
import LazySets
const DI = Dionysos
const UT = DI.Utils

# The golden values below are in the quadratic-form convention
# {x : (x−c)ᵀP(x−c) ≤ 1}; LazySets stores the shape matrix Q = P⁻¹.
ell(P, c) = LazySets.Ellipsoid(c, UT._symmetrize(LA.inv(P)))

@testset "EllipsoidBasics" begin
    c1 = [1.0; 1.0]
    P1 = [
        0.4 -0.1
        -0.1 0.5
    ]
    E1 = ell(P1, c1)
    c2 = [4.0; 1.0]
    P2 = [
        1.0 0.0
        0.0 1.0
    ]
    E2 = ell(P2, c2)
    c3 = [4.0; 1.0]
    P3 = [
        0.25 0.0
        0.0 1.0
    ]
    E3 = ell(P3, c3)
    # sublevel scalings compose: (P/2)·2 = P
    E4 = UT.get_sublevel_set(UT.get_sublevel_set(E1, 2.0), 0.5)

    @test UT.get_center(E1) == [1.0; 1.0]
    @test isapprox(UT.get_quadratic_form(E1), P1, atol = 1e-12)
    @test LA.isposdef(UT.get_quadratic_form(E1)) &&
          LA.isposdef(UT.get_quadratic_form(E2)) &&
          LA.isposdef(UT.get_quadratic_form(E3))
    @test UT.get_dim(E1) == 2
    @test UT.center_distance(E1, E2) == 3.0
    @test UT.center_distance(E1, [5.0; 4.0]) == 5.0
    @test abs(UT.get_volume(E2) - π) <= 10e-6
    @test isapprox(UT.get_farthest_point(E3, [0.0; 1.0]), [0.0; 1.0], atol = 1e-12)
    @test isapprox(UT.get_farthest_point(E3, [1.0; 0.0]), [2.0; 0.0], atol = 1e-12)
    B2 = UT.get_min_bounding_box(E2)
    @test isapprox(LazySets.low(B2), UT.get_center(E2) .- [1.0; 1.0], atol = 1e-9)
    @test isapprox(LazySets.high(B2), UT.get_center(E2) .+ [1.0; 1.0], atol = 1e-9)
    B3 = UT.get_min_bounding_box(E3)
    @test isapprox(LazySets.low(B3), UT.get_center(E3) .- [2.0; 1.0], atol = 1e-9)
    @test isapprox(LazySets.high(B3), UT.get_center(E3) .+ [2.0; 1.0], atol = 1e-9)
    @test isapprox(UT.get_quadratic_form(E4), P1, atol = 1e-12)
    # E3 has distinct semi-axes (2 along x, 1 along y), so the i-th longest
    # axis is well defined; the endpoint order depends on the (arbitrary)
    # eigenvector sign, hence the sort.
    q1, q2 = sort(collect(UT.get_axis_points(E3, 1)))
    @test isapprox(q1, [2.0, 1.0], atol = 1e-9) && isapprox(q2, [6.0, 1.0], atol = 1e-9)
    q1, q2 = sort(collect(UT.get_axis_points(E3, 2)))
    @test isapprox(q1, [4.0, 0.0], atol = 1e-9) && isapprox(q2, [4.0, 2.0], atol = 1e-9)
    @test isapprox(UT.get_length_semiaxis(E3), [2.0, 1.0], atol = 1e-9)
end

@testset "EllipsoidOperations" begin
    c0 = [1.6; 1.4]
    P0 = [
        0.4 -0.1
        -0.1 0.5
    ]

    c = [1.5; 1.5]
    P = [
        4.0 0.5
        0.5 6.0
    ]

    E0 = ell(P0, c0)
    E = ell(P, c)

    ############################################
    a = 0.8
    E1 = ell(P0, c0 + [a; a])

    @test !UT.is_disjoint(E1, E)
    @test !UT.is_disjoint(E, E1)
    E1scaled = UT.scale_for_noninclusion_contact_point(E1, E)
    err =
        LA.norm(UT.get_center(E1scaled) - [2.4; 2.2]) +
        LA.norm(UT.get_quadratic_form(E1scaled) - [2.388 -0.597; -0.597 2.985])
    @test err <= 10e-4
    @test UT.is_included(E1, E) == false
    @test UT.is_included(E, E1) == true
    @test UT.is_included(E, E1scaled) == false
    @test UT.is_included(E1scaled, E) == false
    E1scaled2 = UT.scale_for_inclusion_contact_point(E1, E)
    err =
        LA.norm(UT.get_center(E1scaled2) - [2.4; 2.2]) +
        LA.norm(UT.get_quadratic_form(E1scaled2) - [0.466 -0.116; -0.116 0.582])
    @test err <= 10e-2
    @test UT.is_included(
        E,
        ell(UT.get_quadratic_form(E1scaled2) * 0.99, UT.get_center(E1scaled2)),
    ) == true
    ############################################
    a = 1.2
    E1 = ell(P0, c0 + [a; a])

    @test !UT.is_disjoint(E1, E)
    @test !UT.is_disjoint(E, E1)
    E1scaled = UT.scale_for_noninclusion_contact_point(E1, E)
    err =
        LA.norm(UT.get_center(E1scaled) - [2.8; 2.599]) +
        LA.norm(UT.get_quadratic_form(E1scaled) - [0.725 -0.181; -0.181 0.906])
    @test err <= 1e-2
    @test UT.is_included(E1, E) == false
    @test UT.is_included(E, E1) == false
    @test UT.is_included(E, E1scaled) == false
    @test UT.is_included(E1scaled, E) == false
    E1scaled2 = UT.scale_for_inclusion_contact_point(E1, E)
    err =
        LA.norm(UT.get_center(E1scaled2) - [2.8; 2.599]) +
        LA.norm(UT.get_quadratic_form(E1scaled2) - [0.254 -0.063; -0.063 0.318])
    @test err <= 10e-2
    @test UT.is_included(
        E,
        ell(UT.get_quadratic_form(E1scaled2) * 0.99, UT.get_center(E1scaled2)),
    ) == true
    ############################################
    a = 2.5
    E1 = ell(P0, c0 + [a; a])

    @test UT.is_disjoint(E1, E)
    @test UT.is_disjoint(E, E1)
    E1scaled = UT.scale_for_noninclusion_contact_point(E1, E)
    err =
        LA.norm(UT.get_center(E1scaled) - [4.1; 3.9]) +
        LA.norm(UT.get_quadratic_form(E1scaled) - [0.1195 -0.02989; -0.029897 0.14948])
    @test err <= 10e-4
    @test UT.is_included(E1, E) == false
    @test UT.is_included(E, E1) == false
    @test UT.is_included(E, E1scaled) == false
    @test UT.is_included(E1scaled, E) == false
    E1scaled2 = UT.scale_for_inclusion_contact_point(E1, E)
    err =
        LA.norm(UT.get_center(E1scaled2) - [4.1; 3.9]) + LA.norm(
            UT.get_quadratic_form(E1scaled2) - [0.073295 -0.018323; -0.01832382 0.0916196],
        )
    @test err <= 10e-4
    @test UT.is_included(
        E,
        ell(UT.get_quadratic_form(E1scaled2) * 0.99, UT.get_center(E1scaled2)),
    ) == true

    #############################################
    Id = [
        1.0 0.0
        0.0 1.0
    ]
    E5 = LazySets.affine_map(Id, E0, c)
    @test isapprox(UT.get_quadratic_form(E5), P0, atol = 1e-9)
    @test isapprox(UT.get_center(E5), [3.1; 2.9], atol = 1e-12)

    E6 = ell(Id, [-2.0; -2.0])
    E7 = ell(Id, [2.0; 2.0])
    @test UT.compress_if_intersection(E6, E7) == E6

    ############################################
    P1 = [
        0.17 -0.09
        -0.09 0.26
    ]
    E8 = ell(P1, c0)
    @test isapprox(sqrt(UT.get_quadratic_form(E8)), P0, atol = 1e-6)
end

println("End test")
end
