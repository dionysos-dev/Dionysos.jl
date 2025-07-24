module TestMain

using Test
using StaticArrays, MathematicalSystems
using Dionysos
const SY = Dionysos.Symbolic

sleep(0.1) # used for good printing
println("Started test")

@testset "VectorContinuousSystem - basic usage" begin
    A1 = [1.0;;]
    B1 = [1.0;;]
    c1 = [0.0]
    X1 = Dionysos.Utils.HyperRectangle([0.0], [1.0])
    U1 = Dionysos.Utils.HyperRectangle([-1.0], [1.0])
    sys1 = MathematicalSystems.ConstrainedAffineControlContinuousSystem(A1, B1, c1, X1, U1)

    A2 = [2.0;;]
    B2 = [2.0;;]
    c2 = [1.0]
    X2 = Dionysos.Utils.HyperRectangle([1.0], [2.0])
    U2 = Dionysos.Utils.HyperRectangle([0.0], [2.0])
    sys2 = MathematicalSystems.ConstrainedAffineControlContinuousSystem(A2, B2, c2, X2, U2)

    v_sys = SY.VectorContinuousSystem([sys1, sys2])

    @test SY.statedim(v_sys) == 2
    @test SY.inputdim(v_sys) == 2

    @test SY.stateset(v_sys) == (
        Dionysos.Utils.HyperRectangle([0.0], [1.0]),
        Dionysos.Utils.HyperRectangle([1.0], [2.0]),
    )
    @test SY.inputset(v_sys) == (
        Dionysos.Utils.HyperRectangle([-1.0], [1.0]),
        Dionysos.Utils.HyperRectangle([0.0], [2.0]),
    )
end

sleep(0.1) # used for good printing
println("End test")
end # module
