module Test_SF

include(joinpath(@__DIR__,"..", "..", "src", "ZMPBipedRobot.jl"))
import .ZMPBipedRobot as ZMProbot
using Test 

sleep(0.1) # used for good printing
println("Started test")


@testset "Swing Foot Trajectory Generator" begin 
    
    URDFfileName = "ZMP_2DBipedRobot.urdf"

    br = ZMProbot.BipedRobot(;
        readFile = true,
        URDFfileName = URDFfileName,
        paramFileName = "param_test.jl",
    )
    fp = ZMProbot.FootPlanner(br = br)
    zt = ZMProbot.ZMPTrajectory(br = br, fp = fp, check = false)
    sf = ZMProbot.SwingFootTrajectory(br = br, fp = fp, zt = zt)

    tref = reduce(vcat, zt.timeVec)
    stepL_plot = reduce(hcat, sf.stepL)
    stepR_plot = reduce(hcat, sf.stepR)

    @test length(sf.stepL) == length(fp.center) 
    @test length(sf.stepR) == length(fp.center) 
    @test length(stepR_plot[1, :]) == length(tref)
    @test length(stepL_plot[1, :]) == length(tref)
end 

sleep(0.1) # used for good printing
println("End test")

end # End Main Module 
