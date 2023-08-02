module Test_BR

include(joinpath(@__DIR__, "..", "..", "src", "ZMPBipedRobot.jl"))
import .ZMPBipedRobot
const ZMProbot = ZMPBipedRobot
using Test

sleep(0.1) # used for good printing
println("Started test")

@testset "Biped Robot" begin
    URDFfileName = "ZMP_2DBipedRobot.urdf"
    br = ZMProbot.BipedRobot(;
        readFile = true,
        URDFfileName = URDFfileName,
        paramFileName = "param_test.jl",
    )

    @test br.L_leg == 0.172
    @test br.L_thigh == 0.20125
    @test br.offset_hip_to_motor == 0.04025
    @test br.offset_ankle_to_foot == 0.009
    @test br.d == 0.052
end

sleep(0.1) # used for good printing
println("End test")

end # End Main Module 
