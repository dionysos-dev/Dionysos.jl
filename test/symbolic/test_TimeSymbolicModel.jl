module TestMain

using Test
using StaticArrays, Plots, MathematicalSystems
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
  # doit exister
sleep(0.1) # used for good printing
println("Started test")

@testset "SymbolicModel_active" begin
    A = reshape([1.0], 1, 1)  # Matrice 1x1
    X = UT.HyperRectangle([0.0], [2.0])  # domaine du temps de 0 à 2 avec Dionysos
    sys_time = MathematicalSystems.ConstrainedLinearContinuousSystem(A, X)
    tstep = 0.5
    tm = SY.BuildTimeSymbolicModel(sys_time, tstep)
    @test tm.tsteps == [0.0, 0.5, 1.0, 1.5, 2.0]
    @test tm.is_active == true
    @test SY.time2int(tm, 0.75) == 2  
    @test SY.int2time(tm, 3) == 1.0  
    @test SY.ceil_time2int(tm, 1.25) == 4  
    @test SY.ceil_time2int(tm, 2.5) == 5  
    @test SY.floor_time2int(tm, 1.1) == 3  
    @test SY.floor_time2int(tm, 1.3) == 3  
    @test SY.floor_time2int(tm, 1.9) == 4  
end
@testset "SymbolicModel_frozen" begin
    A = reshape([0.0], 1, 1)  # Matrice 1x1
    X = UT.HyperRectangle([0.0], [2.0])  # domaine du temps de 0 à 2 avec Dionysos
    sys_time = MathematicalSystems.ConstrainedLinearContinuousSystem(A, X)
    tstep = 0.5
    tm = SY.BuildTimeSymbolicModel(sys_time, tstep)
    @test tm.tsteps == [0.0]
    @test tm.is_active == false
    @test SY.time2int(tm, 1.5) == 1  
    @test SY.int2time(tm, 1) == 0.0  
    @test SY.ceil_time2int(tm, 1.5) == 1  
    @test SY.floor_time2int(tm, 1.5) == 1  
end
sleep(0.1) # used for good printing
println("End test")
end