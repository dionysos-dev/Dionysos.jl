module TestFP

include(joinpath(@__DIR__,"..", "..", "src", "ZMPBipedRobot.jl"))
import .ZMPBipedRobot as ZMProbot
using Test 

sleep(0.1) # used for good printing
println("Started test")

br = ZMProbot.BipedRobot(;
    readFile = true,
    paramFileName = "param_test.jl",
)
# Straight path 
t = vec(0 : 100)
br.yPath = 1.18 .+ 0.0 .* t
br.xPath =  0.01*t;
br.initial_position = [br.xPath[1]; br.yPath[1]; 0.0]
fp = ZMProbot.FootPlanner(br = br)
@testset "Foot Planner Diagonal Straight Path" begin 
    @test fp.center[end][1:2] ≈ [br.xPath[end]; br.yPath[end]] atol = 0.01
    @test length(fp.center) > 0  
end 

# Diagonal Straight path 
t = vec(0 : 100);            
br.xPath =  0.01*t
br.yPath = 1.18 .+  0.01*t
br.initial_position = [br.xPath[1], br.yPath[1], pi/4]
fp = ZMProbot.FootPlanner(br = br, check = false)
@testset "Foot Planner Diagonal Straight Path" begin 
    @test fp.center[end][1:2] ≈ [br.xPath[end]; br.yPath[end]] atol = 0.01
    @test length(fp.center) > 0  
end 

# Circular path test 
t = vec(100 : -1: 75); 
br.xPath = -0 .- 1.18 * sin.(2*pi/100 .*t);
br.yPath = 0 .+ 1.18 * cos.(2*pi/100 .*t);
br.initial_position = [br.xPath[1], br.yPath[1], 0]
fp = ZMProbot.FootPlanner(br = br)
@testset "Foot Planner Circular Path" begin 
    @test fp.center[end][1:2] ≈ [br.xPath[end]; br.yPath[end]] atol = 0.01
    @test length(fp.center) > 0  
end 

# Sinusoidale path test 
t = collect(100 : - 1: 0)
br.yPath = 1.18 * cos.(4*pi/100 .* t)
br.xPath = 0.01 .* (100 .- t)
br.initial_position = [br.xPath[1], br.yPath[1], 0]
fp = ZMProbot.FootPlanner(br = br, check = false)

@testset "Foot Planner Sinusoidale Path" begin 
    @test fp.center[end][1:2] ≈ [br.xPath[end]; br.yPath[end]] atol = 0.01
    @test length(fp.center) > 0  
end 

sleep(0.1) # used for good printing
println("End test")

end # End Main Module 
