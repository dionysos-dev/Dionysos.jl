module TestPC

include(joinpath(@__DIR__,"..", "..", "src", "ZMPBipedRobot.jl"))
import .ZMPBipedRobot as ZMProbot
using Test 

sleep(0.1) # used for good printing
println("Started test")

br = ZMProbot.BipedRobot(;
    readFile = true,
    paramFileName = "param_test.jl",
)

fp = ZMProbot.FootPlanner(br = br)
zt = ZMProbot.ZMPTrajectory(br = br, fp = fp, check = true)
@testset "ZMP Trajectory Generator" begin 
    @test length(zt.ZMP) == length(fp.center)
end 

sleep(0.1) # used for good printing
println("End test")

end # End Main Module 
