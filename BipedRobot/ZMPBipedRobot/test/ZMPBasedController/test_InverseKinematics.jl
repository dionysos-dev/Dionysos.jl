module Test_IK

include(joinpath(@__DIR__, "..", "..", "src", "ZMPBipedRobot.jl"))
import .ZMPBipedRobot
const ZMProbot = ZMPBipedRobot
using Test

sleep(0.1) # used for good printing
println("Started test")

@testset "Inverse Kinematics" begin
    URDFfileName = "ZMP_2DBipedRobot.urdf"
    br = ZMProbot.BipedRobot(;
        readFile = true,
        URDFfileName = URDFfileName,
        paramFileName = "param_test.jl",
    )

    # Construct the Preview Controller
    pc = ZMProbot.PreviewController(; br = br)

    # Run the Foot Planer Algorithm and get the foot position 
    fp = ZMProbot.FootPlanner(; br = br)

    # Get the ZMP reference trajectory 
    zt = ZMProbot.ZMPTrajectory(; br = br, fp = fp)

    # Convert the ZMP reference trajectory into CoM trajectory
    ct = ZMProbot.CoMTrajectory(; br = br, pc = pc, zt = zt)

    # Get the Swing Foot trajectory 
    sf = ZMProbot.SwingFootTrajectory(; br = br, fp = fp, zt = zt)

    # Get the joint trajectory from the all path 
    ik = ZMProbot.InverseKinematics(; br = br, fp = fp, ct = ct, sf = sf)
    ik.q_r[:, 1]
    # Get local trajectories in the sagittal plane 
    stepLocalR, stepLocalL = ZMProbot.computeGlobal2Local(br, fp, ct, sf)
    steplocalR_plot = reduce(hcat, stepLocalR)
    steplocalL_plot = reduce(hcat, stepLocalL)

    # Hip frame location in foot frame 
    p_foot2hip = ZMProbot.foot2hip(br)

    # Computed Direct Kinematics using D-H convention 
    # Here, all angles are negative because the simulator have another convention (opposite sign)
    xhip_Rfoot = p_foot2hip[1].(-ik.q_r[:, 1], -ik.q_r[:, 2])
    zhip_Rfoot = p_foot2hip[3].(-ik.q_r[:, 1], -ik.q_r[:, 2]) .+ br.offset_hip_to_motor
    xhip_Lfoot = p_foot2hip[1].(-ik.q_l[:, 1], -ik.q_l[:, 2])
    zhip_Lfoot = p_foot2hip[3].(-ik.q_l[:, 1], -ik.q_l[:, 2]) .+ br.offset_hip_to_motor

    @test [
        isapprox(xhip_Lfoot, steplocalL_plot[1, :]; atol = 1e-6) for
        i in 1:length(xhip_Rfoot)
    ] == ones(length(xhip_Rfoot))
    @test [
        isapprox(zhip_Lfoot, steplocalL_plot[3, :]; atol = 1e-6) for
        i in 1:length(xhip_Rfoot)
    ] == ones(length(xhip_Rfoot))
    @test [
        isapprox(xhip_Rfoot, steplocalR_plot[1, :]; atol = 1e-6) for
        i in 1:length(xhip_Rfoot)
    ] == ones(length(xhip_Rfoot))
    @test [
        isapprox(zhip_Rfoot, steplocalR_plot[3, :]; atol = 1e-6) for
        i in 1:length(xhip_Rfoot)
    ] == ones(length(xhip_Rfoot))
end

sleep(0.1) # used for good printing
println("End test")

end # End Main Module 
