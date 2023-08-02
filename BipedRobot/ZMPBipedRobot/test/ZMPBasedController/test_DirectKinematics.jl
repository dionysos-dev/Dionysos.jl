module Test_DK

include(joinpath(@__DIR__, "..", "..", "src", "ZMPBipedRobot.jl"))
import .ZMPBipedRobot
const ZMProbot = ZMPBipedRobot
using Test

sleep(0.1) # used for good printing
println("Started test")

@testset "Direct Kinematics" begin
    URDFfileName = "ZMP_2DBipedRobot.urdf"

    br = ZMProbot.BipedRobot(;
        readFile = true,
        URDFfileName = URDFfileName,
        paramFileName = "param_test.jl",
    )
    # Simple test with a given value of q1 and q2 
    q1 = 0:(pi / 100):pi
    q2 = (-2 * pi):(pi / 50):0

    # Foot frame location in hip frame 
    p_foot2hip = ZMProbot.foot2hip(br)
    # Hip frame location in foot frame 
    p_hip2foot = ZMProbot.hip2foot(br)

    # Expected Direct Kinematics  
    dk = ZMProbot.twoLinksDirectKinematics.(q1, q2, br.L_thigh, br.L_leg)
    xz = reduce(vcat, dk)

    # Computed Direct Kinematics using D-H convention 
    xhip_foot = p_foot2hip[1].(q1, q2)
    zhip_foot = p_foot2hip[3].(q1, q2) .+ br.offset_hip_to_motor # remove the offset to compare with expected value 
    xfoot_hip = p_hip2foot[1].(q1, q2)
    zfoot_hip = p_hip2foot[2].(q1, q2)

    @test [isapprox(xz[:, 1], xhip_foot; atol = 1e-6) for i in 1:length(xhip_foot)] == ones(length(xhip_foot))
    @test [isapprox(xz[:, 2], zhip_foot; atol = 1e-6) for i in 1:length(xhip_foot)] == ones(length(xhip_foot))
end

sleep(0.1) # used for good printing
println("End test")

end # End Main Module 
