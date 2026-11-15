module TestMain

using Test
using Dionysos
using Clarabel, Ipopt, JuMP
using HybridSystems

const DI = Dionysos
const UT = DI.Utils
const SY = DI.Symbolic

println("Started test")

@testset "Get min bounding box" begin
    opt_qp = optimizer_with_attributes(Ipopt.Optimizer, MOI.Silent() => true)
    P = [
        1.0 0.0
        0.0 0.25
    ]
    box = UT.get_min_bounding_box(UT.Ellipsoid(P, zeros(2)); optimizer = opt_qp)
    x = [interval.hi for interval in box]
    @test x ≈ [1.0, 2.0] atol = 1e-2
    P = [1.5 -0.5; -0.5 2.5]
    box = UT.get_min_bounding_box(UT.Ellipsoid(P, zeros(2)); optimizer = opt_qp)
    x = [interval.hi for interval in box]
    @test x ≈ [0.84515, 0.65465] atol = 1e-2
end

@testset "Provide P" begin
    opt_sdp = optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => true)
    A = [
        0.0 1.0 0.0
        0.0 0.0 1.0
        2.0 1.0 5.0
    ]
    B = [0.0; 0.0; 1.0][:, :]
    c = zeros(3)

    sys = HybridSystems.ConstrainedAffineControlDiscreteSystem(A, B, c, Nothing, Nothing)
    ans, K, P, gamma = SY._provide_P(sys, opt_sdp)
    @test K ≈ [-2.0 -1.0 -5.0] atol = 1e-2
end

println("End test")

end  # module TestMain
