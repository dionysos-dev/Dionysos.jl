module TestMain

using Test
using Dionysos
using Clarabel, JuMP
using HybridSystems

const DI = Dionysos
const ST = DI.System

println("Started test")

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
    ans, K, P, gamma = ST._provide_P(sys, opt_sdp)
    @test K ≈ [-2.0 -1.0 -5.0] atol = 1e-1
end

println("End test")

end  # module TestMain
