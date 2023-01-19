module TestMain

using Test
using Dionysos
using SDPA, Ipopt, JuMP
using HybridSystems

const DI = Dionysos
const SY = DI.Symbolic

sleep(0.1) # used for good printing
println("Started test")

@testset "Get min bounding box" begin
    opt_qp = optimizer_with_attributes(Ipopt.Optimizer, MOI.Silent() => true)
    P = [1.0 0.0;
         0.0 0.25]
    x = SY._get_min_bounding_box(P, opt_qp)
    @test x ≈ [1.0, 2.0] atol=1e-2
    P = [1.5 -0.5; -0.5 2.5]
    x = SY._get_min_bounding_box(P, opt_qp)
    @test x ≈ [0.84515, 0.65465] atol=1e-2
end

@testset "Provide P" begin
    opt_sdp = optimizer_with_attributes(SDPA.Optimizer, MOI.Silent() => true)
    A = [0.0 1.0 0.0; 
         0.0 0.0 1.0; 
         2.0 1.0 5.0]
    B = [0.0; 0.0; 1.0][:,:]
    c = zeros(3)

    sys = HybridSystems.ConstrainedAffineControlDiscreteSystem(A,B,c,Nothing, Nothing)
    ans, K, P, gamma = SY._provide_P(sys, opt_sdp)
    @test K ≈ [-1.97742  -1.0  -5.0] atol=1e-2
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
