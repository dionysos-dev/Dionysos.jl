module TestContinuousDomain
using Test
using StaticArrays
using Dionysos
const DI = Dionysos
const DO = DI.Domain
const UT = DI.Utils

@testset "ContinuousDomain" begin
    orig  = SVector(0.0, 0.0)
    bound = 10.0
    dom = DO.ContinuousBoundedEllipsoidDomain(orig, bound)
    println(dom)
    @test dom.orig === orig
    @test dom.bound == bound
    @test dom.ellips isa Set{UT.Ellipsoid}
    @test isempty(dom.ellips)

    P = SMatrix{2, 2}(2.0, 1.0, 1.0, 3.0)
    ell = UT.Ellipsoid(P, orig)
    dom2 = DO.ContinuousBoundedEllipsoidDomain(orig, bound, Set([ell]))
    println(dom2)
    @test dom2.ellips isa Set{typeof(ell)}
    @test length(dom2.ellips) == 1
    @test ell in dom2.ellips
end
end 
