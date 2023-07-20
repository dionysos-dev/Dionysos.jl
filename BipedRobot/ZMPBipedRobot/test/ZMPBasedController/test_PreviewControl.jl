module Test_PC

include(joinpath(@__DIR__,"..", "..", "src", "ZMPBipedRobot.jl"))
import .ZMPBipedRobot as ZMProbot
using Test 

sleep(0.1) # used for good printing
println("Started test")

@testset "Preview Control" begin 
        
    br = ZMProbot.BipedRobot(;
    readFile = true,
    paramFileName = "param_test.jl",
    )

    pc = ZMProbot.PreviewController(br = br)
    @test length(pc.Gx) == 3
    @test length(pc.Gi) == 1
    @test length(pc.Gd) == (br.previewTime/br.Ts + 1)
end 

sleep(0.1) # used for good printing
println("End test")

end # End Main Module 
