module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import MathematicalSystems as MS
import LazySets

@testset "VectorContinuousSystem Construction" begin
    # Test basic construction with two affine control systems
    A1 = [1.0;;]
    B1 = [1.0;;]
    c1 = [0.0]
    X1 = LazySets.Hyperrectangle(; low = [0.0], high = [1.0])
    U1 = LazySets.Hyperrectangle(; low = [-1.0], high = [1.0])
    sys1 = MS.ConstrainedAffineControlContinuousSystem(A1, B1, c1, X1, U1)

    A2 = [2.0;;]
    B2 = [2.0;;]
    c2 = [1.0]
    X2 = LazySets.Hyperrectangle(; low = [1.0], high = [2.0])
    U2 = LazySets.Hyperrectangle(; low = [0.0], high = [2.0])
    sys2 = MS.ConstrainedAffineControlContinuousSystem(A2, B2, c2, X2, U2)

    v_sys = ST.VectorContinuousSystem([sys1, sys2])

    @test v_sys isa ST.VectorContinuousSystem
    @test length(v_sys.systems) == 2
    @test v_sys.systems[1] === sys1
    @test v_sys.systems[2] === sys2
end

@testset "VectorContinuousSystem Dimensions" begin
    # Create systems with different dimensions
    A1 = [1.0;;]
    B1 = [1.0;;]
    c1 = [0.0]
    X1 = LazySets.Hyperrectangle(; low = [0.0], high = [1.0])
    U1 = LazySets.Hyperrectangle(; low = [-1.0], high = [1.0])
    sys1 = MS.ConstrainedAffineControlContinuousSystem(A1, B1, c1, X1, U1)

    A2 = [2.0;;]
    B2 = [2.0;;]
    c2 = [1.0]
    X2 = LazySets.Hyperrectangle(; low = [1.0], high = [2.0])
    U2 = LazySets.Hyperrectangle(; low = [0.0], high = [2.0])
    sys2 = MS.ConstrainedAffineControlContinuousSystem(A2, B2, c2, X2, U2)

    v_sys = ST.VectorContinuousSystem([sys1, sys2])

    @test MS.statedim(v_sys) == 2  # 1 + 1
    @test MS.inputdim(v_sys) == 2  # 1 + 1
end

@testset "VectorContinuousSystem State and Input Sets" begin
    A1 = [1.0;;]
    B1 = [1.0;;]
    c1 = [0.0]
    X1 = LazySets.Hyperrectangle(; low = [0.0], high = [1.0])
    U1 = LazySets.Hyperrectangle(; low = [-1.0], high = [1.0])
    sys1 = MS.ConstrainedAffineControlContinuousSystem(A1, B1, c1, X1, U1)

    A2 = [2.0;;]
    B2 = [2.0;;]
    c2 = [1.0]
    X2 = LazySets.Hyperrectangle(; low = [1.0], high = [2.0])
    U2 = LazySets.Hyperrectangle(; low = [0.0], high = [2.0])
    sys2 = MS.ConstrainedAffineControlContinuousSystem(A2, B2, c2, X2, U2)

    v_sys = ST.VectorContinuousSystem([sys1, sys2])

    @test MS.stateset(v_sys) == (
        LazySets.Hyperrectangle(; low = [0.0], high = [1.0]),
        LazySets.Hyperrectangle(; low = [1.0], high = [2.0]),
    )
    @test MS.inputset(v_sys) == (
        LazySets.Hyperrectangle(; low = [-1.0], high = [1.0]),
        LazySets.Hyperrectangle(; low = [0.0], high = [2.0]),
    )
end

@testset "VectorContinuousSystem Mixed System Types" begin
    # Test with different system types
    A1 = [1.0;;]
    X1 = LazySets.Hyperrectangle(; low = [0.0], high = [1.0])
    linear_sys = MS.ConstrainedLinearContinuousSystem(A1, X1)

    A2 = [2.0;;]
    B2 = [2.0;;]
    c2 = [1.0]
    X2 = LazySets.Hyperrectangle(; low = [1.0], high = [2.0])
    U2 = LazySets.Hyperrectangle(; low = [0.0], high = [2.0])
    affine_sys = MS.ConstrainedAffineControlContinuousSystem(A2, B2, c2, X2, U2)

    v_sys = ST.VectorContinuousSystem([linear_sys, affine_sys])

    @test MS.statedim(v_sys) == 2  # 1 + 1
    @test MS.inputdim(v_sys) == 1  # 0 + 1 (linear system has no input)
end

@testset "VectorContinuousSystem Single System" begin
    # Test with single system
    A = [1.0;;]
    B = [1.0;;]
    c = [0.0]
    X = LazySets.Hyperrectangle(; low = [0.0], high = [1.0])
    U = LazySets.Hyperrectangle(; low = [-1.0], high = [1.0])
    sys = MS.ConstrainedAffineControlContinuousSystem(A, B, c, X, U)

    v_sys = ST.VectorContinuousSystem([sys])

    @test MS.statedim(v_sys) == 1
    @test MS.inputdim(v_sys) == 1
    @test length(MS.stateset(v_sys)) == 1
    @test length(MS.inputset(v_sys)) == 1
end

@testset "VectorContinuousSystem Higher Dimensions" begin
    # Test with higher dimensional systems
    A1 = [1.0 0.0; 0.0 1.0]  # 2D system
    B1 = [1.0; 0.0;;]        # 2D state, 1D input (matrix format)
    c1 = [0.0, 0.0]
    X1 = LazySets.Hyperrectangle(; low = [0.0, 0.0], high = [1.0, 1.0])
    U1 = LazySets.Hyperrectangle(; low = [-1.0], high = [1.0])
    sys1 = MS.ConstrainedAffineControlContinuousSystem(A1, B1, c1, X1, U1)

    A2 = [2.0;;]  # 1D system
    B2 = [1.0;;]
    c2 = [1.0]
    X2 = LazySets.Hyperrectangle(; low = [1.0], high = [2.0])
    U2 = LazySets.Hyperrectangle(; low = [0.0], high = [2.0])
    sys2 = MS.ConstrainedAffineControlContinuousSystem(A2, B2, c2, X2, U2)

    v_sys = ST.VectorContinuousSystem([sys1, sys2])

    @test MS.statedim(v_sys) == 3  # 2 + 1
    @test MS.inputdim(v_sys) == 2  # 1 + 1
    @test length(MS.stateset(v_sys)) == 2
    @test length(MS.inputset(v_sys)) == 2
end

end # module
