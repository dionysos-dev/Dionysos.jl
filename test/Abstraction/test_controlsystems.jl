include(joinpath(@__DIR__, "../../src/Abstraction/abstraction.jl"))

module TestMain

using Test
using LinearAlgebra
using Main.Abstraction
ABS = Main.Abstraction
using StaticArrays

println("")

@testset "Abstraction/controlsystems" begin
tstep = 1.0
nsys = 1_000
ngrowthbound = 1_000
F_sys(x, u) = SVector(u[1], -x[1])
L_growthbound(u) = SMatrix{2,2}(0.0, 1.0, 0.0, 0.0)
sysnoise = SVector(1.0, 1.0)
measnoise = SVector(1.0, 1.0)

contsys = ABS.ControlSystemGrowthRK4(
    tstep,
    F_sys, L_growthbound,
    sysnoise, measnoise,
    nsys, ngrowthbound)

@test contsys.tstep === tstep
@test contsys.sysnoise === sysnoise
@test contsys.measnoise === measnoise
@test ABS.xspacetype(contsys) === (2, Float64)

@test contsys.sys_map(SVector(1.0, 1.0), SVector(1.0), 2.0) ≈ SVector(3.0, -3.0)
@test contsys.growthbound_map(SVector(1.0, 1.0), SVector(1.0), 2.0) ≈ SVector(3.0, 7.0)

DF_sys(x, u) = SMatrix{2,2}(0.0, -1.0, 0.0, 0.0)
bound_DF(u) = 1.0
bound_DDF(u) = 1.0

contsys = ABS.ControlSystemLinearizedRK4(
    tstep,
    F_sys, DF_sys, bound_DF, bound_DDF,
    measnoise,
    nsys)

@test contsys.tstep === tstep
@test contsys.measnoise === measnoise

@test contsys.sys_map(SVector(1.0, 1.0), SVector(1.0), 2.0) ≈ SVector(3.0, -3.0)
@test all(contsys.linsys_map(SVector(1.0, 1.0), SMatrix{2,2}(1.0I), SVector(1.0), 2.0) .≈
    (SVector(3.0, -3.0), SMatrix{2,2}(1.0, -2.0, 0.0, 1.0)))
@test contsys.error_map(1.0, SVector(1.0), 2.0) ≈ exp(2)*(exp(2)-1)/2
print("")
end

end  # module TestMain
