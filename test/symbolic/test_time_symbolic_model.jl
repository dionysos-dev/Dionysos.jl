module TestMain

using Test
using StaticArrays, MathematicalSystems
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic

sleep(0.1) # used for good printing
println("Started tests")

@testset "TimeSymbolicModel Time Active (Identity Matrix)" begin
    # Test case where time evolves (A ≈ I)
    A = reshape([1.0], 1, 1)
    X = UT.HyperRectangle([0.0], [2.0])
    sys_time = MathematicalSystems.ConstrainedLinearContinuousSystem(A, X)
    tstep = 0.5
    tm = SY.TimeSymbolicModel(sys_time, tstep)

    # Test basic properties
    @test tm.tsteps == [0.0, 0.5, 1.0, 1.5, 2.0]
    @test tm.is_time_active == true
    @test tm.time_domain === X

    # Test time conversion functions
    @test SY.floor_time2int(tm, 0.75) == 2
    @test SY.int2time(tm, 3) == 1.0
    @test SY.ceil_time2int(tm, 1.25) == 4
    @test SY.ceil_time2int(tm, 2.5) == 5  # Beyond domain
    @test SY.floor_time2int(tm, 1.1) == 3
    @test SY.floor_time2int(tm, 1.3) == 3
    @test SY.floor_time2int(tm, 1.9) == 4
end

@testset "TimeSymbolicModel Time Frozen (Zero Matrix)" begin
    # Test case where time is frozen (A ≈ 0)
    A = reshape([0.0], 1, 1)
    X = UT.HyperRectangle([0.0], [2.0])
    sys_time = MathematicalSystems.ConstrainedLinearContinuousSystem(A, X)
    tstep = 0.5
    tm = SY.TimeSymbolicModel(sys_time, tstep)

    # Test basic properties
    @test tm.tsteps == [0.0]
    @test tm.is_time_active == false
    @test tm.time_domain === X

    # Test time conversion functions (all should return 1 for frozen time)
    @test SY.floor_time2int(tm, 1.5) == 1
    @test SY.int2time(tm, 1) == 0.0
    @test SY.ceil_time2int(tm, 1.5) == 1
    @test SY.floor_time2int(tm, 1.5) == 1
end

@testset "TimeSymbolicModel Boundary Cases Active Time" begin
    A = reshape([1.0], 1, 1)
    X = UT.HyperRectangle([0.0], [2.0])
    sys_time = MathematicalSystems.ConstrainedLinearContinuousSystem(A, X)
    tstep = 0.5
    tm = SY.TimeSymbolicModel(sys_time, tstep)

    # Test exact time values (floor behavior for floor_time2int)
    @test SY.floor_time2int(tm, 0.0) == 1    # Exact start
    @test SY.floor_time2int(tm, 0.5) == 2    # Exact step (0.5 <= 0.5, so index 2)
    @test SY.floor_time2int(tm, 1.0) == 3    # Exact step (1.0 <= 1.0, so index 3)
    @test SY.floor_time2int(tm, 2.0) == 5    # Exact end (2.0 <= 2.0, so index 5)

    # Test ceil function at boundaries
    @test SY.ceil_time2int(tm, 0.0) == 1
    @test SY.ceil_time2int(tm, 0.5) == 2
    @test SY.ceil_time2int(tm, 1.0) == 3
    @test SY.ceil_time2int(tm, 2.0) == 5

    # Test floor function at boundaries (should be same as floor_time2int)
    @test SY.floor_time2int(tm, 0.0) == 1
    @test SY.floor_time2int(tm, 0.5) == 2    # floor(0.5) means largest index with value <= 0.5
    @test SY.floor_time2int(tm, 1.0) == 3    # floor(1.0) means largest index with value <= 1.0
    @test SY.floor_time2int(tm, 2.0) == 5    # floor(2.0) means largest index with value <= 2.0
end

@testset "TimeSymbolicModel Edge Cases" begin
    A = reshape([1.0], 1, 1)
    X = UT.HyperRectangle([0.0], [2.0])
    sys_time = MathematicalSystems.ConstrainedLinearContinuousSystem(A, X)
    tstep = 0.5
    tm = SY.TimeSymbolicModel(sys_time, tstep)

    # Test negative time (should clamp to first index)
    @test SY.floor_time2int(tm, -1.0) == 1
    @test SY.ceil_time2int(tm, -1.0) == 1

    # Test time beyond domain  
    @test SY.floor_time2int(tm, 3.0) == 5  # Last valid index (2.0 is at index 5)
    @test SY.ceil_time2int(tm, 3.0) == 5   # Beyond domain, clamp to last

    # Test int2time boundary cases
    @test SY.int2time(tm, 1) == 0.0
    @test SY.int2time(tm, 5) == 2.0
end

@testset "TimeSymbolicModel Different Step Sizes" begin
    # Test with smaller step size
    A = reshape([1.0], 1, 1)
    X = UT.HyperRectangle([0.0], [1.0])
    sys_time = MathematicalSystems.ConstrainedLinearContinuousSystem(A, X)
    tstep = 0.25
    tm = SY.TimeSymbolicModel(sys_time, tstep)

    @test tm.tsteps == [0.0, 0.25, 0.5, 0.75, 1.0]
    @test length(tm.tsteps) == 5
    @test tm.is_time_active == true

    # Test with larger step size
    tstep_large = 1.0
    tm_large = SY.TimeSymbolicModel(sys_time, tstep_large)
    @test tm_large.tsteps == [0.0, 1.0]
    @test length(tm_large.tsteps) == 2
end

@testset "TimeSymbolicModel Different Domain Sizes" begin
    # Test with different domain bounds
    A = reshape([1.0], 1, 1)
    X = UT.HyperRectangle([-1.0], [1.0])  # Domain from -1 to 1
    sys_time = MathematicalSystems.ConstrainedLinearContinuousSystem(A, X)
    tstep = 0.5
    tm = SY.TimeSymbolicModel(sys_time, tstep)

    @test tm.tsteps == [-1.0, -0.5, 0.0, 0.5, 1.0]
    @test tm.is_time_active == true

    # Test conversion functions with negative times
    @test SY.int2time(tm, 1) == -1.0
    @test SY.int2time(tm, 3) == 0.0
    @test SY.floor_time2int(tm, -0.25) == 2  # Between -0.5 and 0.0
end

@testset "TimeSymbolicModel Error Cases" begin
    # Test with invalid matrix (neither identity nor zero)
    A = reshape([0.5], 1, 1)  # Not identity, not zero
    X = UT.HyperRectangle([0.0], [2.0])
    sys_time = MathematicalSystems.ConstrainedLinearContinuousSystem(A, X)
    tstep = 0.5

    @test_throws ErrorException SY.TimeSymbolicModel(sys_time, tstep)
end

@testset "TimeSymbolicModel Type Stability" begin
    A = reshape([1.0], 1, 1)
    X = UT.HyperRectangle([0.0], [2.0])
    sys_time = MathematicalSystems.ConstrainedLinearContinuousSystem(A, X)
    tstep = 0.5
    tm = SY.TimeSymbolicModel(sys_time, tstep)

    # Test return types
    @test SY.floor_time2int(tm, 1.0) isa Int
    @test SY.int2time(tm, 2) isa Float64
    @test SY.ceil_time2int(tm, 1.0) isa Int
    @test SY.floor_time2int(tm, 1.0) isa Int
end

sleep(0.1) # used for good printing
println("End test")

end
