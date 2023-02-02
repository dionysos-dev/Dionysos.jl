module TestMain

using Test
using Dionysos
using LinearAlgebra
using IntervalArithmetic
const DI = Dionysos
const UT = DI.Utils

@testset "EllipsoidBasics" begin
    c1 = [1.0; 1.0]
    P1 = [0.4 -0.1;
          -0.1 0.5]
    E1 = UT.Ellipsoid(P1, c1)
    c2 = [4.0; 1.0]
    P2 = [1.0 0.0;
          0.0 1.0]
    E2 = UT.Ellipsoid(P2, c2)
    c3 = [4.0; 1.0]
    P3 = [0.25 0.0;
          0.0 1.0]
    E3 = UT.Ellipsoid(P3, c3)

    @test UT.get_center(E1) == [1.0; 1.0]
    @test UT.get_shape(E1) == [0.4 -0.1;-0.1 0.5]
    @test isposdef(UT.get_shape(E1)) && isposdef(UT.get_shape(E2)) && isposdef(UT.get_shape(E3)) == true
    @test UT.get_dims(E1) == 2   
    @test UT.centerDistance(E1, E2) == 3.0
    @test UT.pointCenterDistance(E1, [5.0;4.0]) == 5.0
    @test abs(UT.get_volume(E2)-π) <= 10e-6
    @test UT.get_farthest_point(E3, [0.0;1.0]) == [0.0;1.0]
    @test UT.get_farthest_point(E3, [1.0;0.0]) == [2.0;0.0]
    @test UT.get_min_bounding_box(E2) == IntervalBox(UT.get_center(E2) .- [1.0;1.0], UT.get_center(E2).+ [1.0;1.0])
    @test UT.get_min_bounding_box(E3) == IntervalBox(UT.get_center(E3) .- [2.0;1.0], UT.get_center(E3).+ [2.0;1.0])
end

@testset "EllipsoidOperations" begin
    c0 = [1.6; 1.4]
    P0 = [0.4 -0.1;
          -0.1 0.5]

    c = [1.5; 1.5]
    P = [4.0 0.5;       
         0.5 6.0]

    E0 = UT.Ellipsoid(P0, c0)
    E = UT.Ellipsoid(P, c)

    ############################################
    a = 0.8
    E1 = UT.Ellipsoid(P0, c0+[a;a])

    @test UT.intersect(E1, E) == true
    @test UT.intersect(E, E1) == true
    E1scaled = UT.scale_for_noninclusion_contact_point(E1, E) 
    err = norm(E1scaled.c-[2.4; 2.2]) + norm(E1scaled.P-[2.388 -0.597; -0.597 2.985])
    @test err <= 10e-4
    @test Base.in(E1, E) == false
    @test Base.in(E, E1) == true
    @test Base.in(E, E1scaled) == false
    @test Base.in(E1scaled, E) == false
    E1scaled2 = UT.scale_for_inclusion_contact_point(E1, E) 
    err = norm(E1scaled2.c-[2.4; 2.2]) + norm(E1scaled2.P-[0.466 -0.116; -0.116 0.582])
    @test err <= 10e-2
    @test Base.in(E, UT.Ellipsoid(E1scaled2.P*0.99, E1scaled2.c)) == true
    ############################################
    a = 1.2 
    E1 = UT.Ellipsoid(P0, c0+[a;a])

    @test UT.intersect(E1, E) == true
    @test UT.intersect(E, E1) == true
    E1scaled = UT.scale_for_noninclusion_contact_point(E1, E) 
    err = norm(E1scaled.c-[2.8; 2.599]) + norm(E1scaled.P-[0.725 -0.181; -0.181 0.906])
    @test err <= 1e-2
    @test Base.in(E1, E) == false
    @test Base.in(E, E1) == false
    @test Base.in(E, E1scaled) == false
    @test Base.in(E1scaled, E) == false
    E1scaled2 = UT.scale_for_inclusion_contact_point(E1, E) 
    err = norm(E1scaled2.c-[2.8; 2.599]) + norm(E1scaled2.P-[0.254 -0.063; -0.063 0.318])
    @test err <= 10e-2
    @test Base.in(E, UT.Ellipsoid(E1scaled2.P*0.99, E1scaled2.c)) == true
    ############################################
    a = 2.5 
    E1 = UT.Ellipsoid(P0, c0+[a;a])

    @test UT.intersect(E1, E) == false
    @test UT.intersect(E, E1) == false
    E1scaled = UT.scale_for_noninclusion_contact_point(E1, E) 
    err = norm(E1scaled.c-[4.1; 3.9]) + norm(E1scaled.P-[0.1195 -0.02989; -0.029897 0.14948])
    @test err <= 10e-4
    @test Base.in(E1, E) == false
    @test Base.in(E, E1) == false
    @test Base.in(E, E1scaled) == false
    @test Base.in(E1scaled, E) == false
    E1scaled2 = UT.scale_for_inclusion_contact_point(E1, E) 
    err = norm(E1scaled2.c-[4.1; 3.9]) + norm(E1scaled2.P-[0.073295 -0.018323; -0.01832382 0.0916196])
    @test err <= 10e-4
    @test Base.in(E, UT.Ellipsoid(E1scaled2.P*0.99, E1scaled2.c)) == true

end

println("End test")
end