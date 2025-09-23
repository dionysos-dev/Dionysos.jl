module TestMain

using Test
using StaticArrays, MathematicalSystems
using Dionysos
const SY = Dionysos.Symbolic
const UT = Dionysos.Utils

sleep(0.1) # used for good printing
println("Started tests")

@testset "VectorContinuousSystem Construction" begin
    # Test basic construction with two affine control systems
    A1 = [1.0;;]
    B1 = [1.0;;]
    c1 = [0.0]
    X1 = UT.HyperRectangle([0.0], [1.0])
    U1 = UT.HyperRectangle([-1.0], [1.0])
    sys1 = MathematicalSystems.ConstrainedAffineControlContinuousSystem(A1, B1, c1, X1, U1)

    A2 = [2.0;;]
    B2 = [2.0;;]
    c2 = [1.0]
    X2 = UT.HyperRectangle([1.0], [2.0])
    U2 = UT.HyperRectangle([0.0], [2.0])
    sys2 = MathematicalSystems.ConstrainedAffineControlContinuousSystem(A2, B2, c2, X2, U2)

    v_sys = SY.VectorContinuousSystem([sys1, sys2])

    @test v_sys isa SY.VectorContinuousSystem
    @test length(v_sys.systems) == 2
    @test v_sys.systems[1] === sys1
    @test v_sys.systems[2] === sys2
end

@testset "VectorContinuousSystem Dimensions" begin
    # Create systems with different dimensions
    A1 = [1.0;;]
    B1 = [1.0;;]
    c1 = [0.0]
    X1 = UT.HyperRectangle([0.0], [1.0])
    U1 = UT.HyperRectangle([-1.0], [1.0])
    sys1 = MathematicalSystems.ConstrainedAffineControlContinuousSystem(A1, B1, c1, X1, U1)

    A2 = [2.0;;]
    B2 = [2.0;;]
    c2 = [1.0]
    X2 = UT.HyperRectangle([1.0], [2.0])
    U2 = UT.HyperRectangle([0.0], [2.0])
    sys2 = MathematicalSystems.ConstrainedAffineControlContinuousSystem(A2, B2, c2, X2, U2)

    v_sys = SY.VectorContinuousSystem([sys1, sys2])

    @test SY.statedim(v_sys) == 2  # 1 + 1
    @test SY.inputdim(v_sys) == 2  # 1 + 1
end

@testset "VectorContinuousSystem State and Input Sets" begin
    A1 = [1.0;;]
    B1 = [1.0;;]
    c1 = [0.0]
    X1 = UT.HyperRectangle([0.0], [1.0])
    U1 = UT.HyperRectangle([-1.0], [1.0])
    sys1 = MathematicalSystems.ConstrainedAffineControlContinuousSystem(A1, B1, c1, X1, U1)

    A2 = [2.0;;]
    B2 = [2.0;;]
    c2 = [1.0]
    X2 = UT.HyperRectangle([1.0], [2.0])
    U2 = UT.HyperRectangle([0.0], [2.0])
    sys2 = MathematicalSystems.ConstrainedAffineControlContinuousSystem(A2, B2, c2, X2, U2)

    v_sys = SY.VectorContinuousSystem([sys1, sys2])

    @test SY.stateset(v_sys) ==
          (UT.HyperRectangle([0.0], [1.0]), UT.HyperRectangle([1.0], [2.0]))
    @test SY.inputset(v_sys) ==
          (UT.HyperRectangle([-1.0], [1.0]), UT.HyperRectangle([0.0], [2.0]))
end

@testset "VectorContinuousSystem Mixed System Types" begin
    # Test with different system types
    A1 = [1.0;;]
    X1 = UT.HyperRectangle([0.0], [1.0])
    linear_sys = MathematicalSystems.ConstrainedLinearContinuousSystem(A1, X1)

    A2 = [2.0;;]
    B2 = [2.0;;]
    c2 = [1.0]
    X2 = UT.HyperRectangle([1.0], [2.0])
    U2 = UT.HyperRectangle([0.0], [2.0])
    affine_sys =
        MathematicalSystems.ConstrainedAffineControlContinuousSystem(A2, B2, c2, X2, U2)

    v_sys = SY.VectorContinuousSystem([linear_sys, affine_sys])

    @test SY.statedim(v_sys) == 2  # 1 + 1
    @test SY.inputdim(v_sys) == 1  # 0 + 1 (linear system has no input)
end

@testset "VectorContinuousSystem Single System" begin
    # Test with single system
    A = [1.0;;]
    B = [1.0;;]
    c = [0.0]
    X = UT.HyperRectangle([0.0], [1.0])
    U = UT.HyperRectangle([-1.0], [1.0])
    sys = MathematicalSystems.ConstrainedAffineControlContinuousSystem(A, B, c, X, U)

    v_sys = SY.VectorContinuousSystem([sys])

    @test SY.statedim(v_sys) == 1
    @test SY.inputdim(v_sys) == 1
    @test length(SY.stateset(v_sys)) == 1
    @test length(SY.inputset(v_sys)) == 1
end

@testset "VectorContinuousSystem Higher Dimensions" begin
    # Test with higher dimensional systems
    A1 = [1.0 0.0; 0.0 1.0]  # 2D system
    B1 = [1.0; 0.0;;]        # 2D state, 1D input (matrix format)
    c1 = [0.0, 0.0]
    X1 = UT.HyperRectangle([0.0, 0.0], [1.0, 1.0])
    U1 = UT.HyperRectangle([-1.0], [1.0])
    sys1 = MathematicalSystems.ConstrainedAffineControlContinuousSystem(A1, B1, c1, X1, U1)

    A2 = [2.0;;]  # 1D system
    B2 = [1.0;;]
    c2 = [1.0]
    X2 = UT.HyperRectangle([1.0], [2.0])
    U2 = UT.HyperRectangle([0.0], [2.0])
    sys2 = MathematicalSystems.ConstrainedAffineControlContinuousSystem(A2, B2, c2, X2, U2)

    v_sys = SY.VectorContinuousSystem([sys1, sys2])

    @test SY.statedim(v_sys) == 3  # 2 + 1
    @test SY.inputdim(v_sys) == 2  # 1 + 1
    @test length(SY.stateset(v_sys)) == 2
    @test length(SY.inputset(v_sys)) == 2
end

sleep(0.1) # used for good printing
println("End test")

end # module
