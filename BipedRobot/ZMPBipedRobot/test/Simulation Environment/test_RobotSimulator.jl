module Test_RS

include(joinpath(@__DIR__,"..", "..", "src", "ZMPBipedRobot.jl"))
import .ZMPBipedRobot as ZMProbot
using Test 

sleep(0.1) # used for good printing
println("Started test")

@testset "Robot Simulator" begin 
    
end 

sleep(0.1) # used for good printing
println("End test")

end # End Main Module 
