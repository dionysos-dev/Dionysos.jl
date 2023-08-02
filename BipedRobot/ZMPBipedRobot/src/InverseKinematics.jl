"""

A Inverse Kinematics perform a analytical computation of the inverse kinematics of a double pendulum.

A `InverseKinematics` has 2 differents outputs :
    1.  `q_r` which is an matrix which containt all references joints angles for the right foot
    2.  `q_l` which is an matrix which containt all references joints angles for the left foot

The input of of the `CoMTrajectory` require : 
    1.  `br` is a structure `BipedRobot` which have all the description of the robot
    2.  `fp` is a structure `FootPlanner`which have the position of the right, left and center point of the foot during the path
    3.  `ct` is a structure `CoMTrajectory` which have the computed CoM
    4.  `sf` is a structure `SwingFootTrajectory` which have the global position of both feet during the step 

"""
struct InverseKinematics
    q_r::Matrix
    q_l::Matrix
end

"""

Direct kinematics of a double pendulum
"""
function twoLinksDirectKinematics(
    θ1::T,
    θ2::T,
    L1::T,
    L2::T,
) where {T <: Union{Int64, Float64}}
    x = L1 * sin(θ1) + L2 * sin(θ1 + θ2)
    z = -L1 * cos(θ1) - L2 * cos(θ1 + θ2)
    return [x z]
end

"""

Inverse kinematics of a double pendulum
"""
function twoLinksInverseKinematics(
    x::T,
    y::T,
    L1::T,
    L2::T,
) where {T <: Union{Int64, Float64}}
    D = (x^2 + y^2 - L1^2 - L2^2) / (2 * L1 * L2) # = cos θ1 
    if (D^2 > 1)
        D = 1
        println("Singular configuration or not reachable for ($(x), $(y)) !")
    end
    E = -sqrt(1 - D^2) # = +/- sin θ1 
    γ = atan(L2 * E / (L1 + L2 * D))
    β = atan(-x, -y)

    # Solution 
    θ1 = γ + β
    θ2 = acos(D)
    return [θ1 θ2]
end

"""
Constructor
"""
function InverseKinematics(;
    br::BipedRobot,
    fp::FootPlanner,
    ct::CoMTrajectory,
    sf::SwingFootTrajectory,
    check::Bool = false,
)
    stepLocalR, stepLocalL = computeGlobal2Local(br, fp, ct, sf)

    steplocalR_plot = reduce(hcat, stepLocalR)
    steplocalL_plot = reduce(hcat, stepLocalL)

    if (check)
        plt = plot(; xlabel = "X [m]", ylabel = "Z [m]", layout = (1, 2), dpi = 600)
        scatter!(
            plt[1],
            steplocalL_plot[1, :],
            steplocalL_plot[3, :];
            ms = 2,
            msw = 0.1,
            label = false,
            title = "Local Left Foot",
        )
        scatter!(
            plt[2],
            steplocalR_plot[1, :],
            steplocalR_plot[3, :];
            ms = 2,
            msw = 0.1,
            label = false,
            title = "Local Right Foot",
        )
        display(plt)
    end
    L1 = br.L_thigh
    L2 = br.L_leg
    qr = twoLinksInverseKinematics.(steplocalR_plot[1, :], steplocalR_plot[3, :], L1, L2)
    ql = twoLinksInverseKinematics.(steplocalL_plot[1, :], steplocalL_plot[3, :], L1, L2)
    qr = reduce(vcat, qr)
    ql = reduce(vcat, ql)

    if (check)
        plt_θ = plot(;
            title = "Inverse Kinematrics",
            xlabel = L"$t$ [s]",
            legend = false,
            dpi = 600,
            layout = (2, 2),
        )
        zt = ZMPTrajectory(; br = br, fp = fp, check = false)
        tplot = reduce(vcat, zt.timeVec)
        qref = [ql[:, 1] qr[:, 1] ql[:, 2] qr[:, 2]]
        for (side_idx, side) in enumerate(["Left", "Right"])
            for (joint, name) in enumerate(["Leg", "Knee"])
                plot!(
                    plt_θ[joint, side_idx],
                    tplot,
                    qref[:, (2 * joint - 1) + (side_idx - 1)];
                    label = false,
                    title = "$(side) $(name)",
                )
                str = latexstring("q_$(2 + (2*joint - 1) + (side_idx - 1))") * " [rad]"
                ylabel!(plt_θ[joint, side_idx], str)
            end
        end
        display(plt_θ)
    end

    return InverseKinematics(qr, ql)
end

"""

Here we convert a global trajectory into a local trajectory trought homogenous transformation 
"""
function computeGlobal2Local(
    br::BipedRobot,
    fp::FootPlanner,
    ct::CoMTrajectory,
    sf::SwingFootTrajectory,
)
    center = fp.center
    stepR = sf.stepR
    stepL = sf.stepL
    hip_height = ct.hip_height
    d1 = br.offset_hip_to_motor
    d4 = br.offset_ankle_to_foot

    stepLocalR = Array{Float64}[]
    stepLocalL = Array{Float64}[]

    p_l_global =
        p_r_global = for stepNum in 1:length(stepL)

            # Angle of rotation between the reference frame and the local frame 
            θ = center[stepNum][3]

            # Change global coordinate into hip frame  
            # Rotation matrix 
            R_ref2foot = [
                cos(θ) -sin(θ) 0
                sin(θ) cos(θ) 0
                0 0 1
            ]

            p_l_global = zeros(4, length(hip_height[stepNum][1, :]))
            p_r_global = zeros(4, length(hip_height[stepNum][1, :]))
            p_l_local = zeros(4, length(hip_height[stepNum][1, :]))
            p_r_local = zeros(4, length(hip_height[stepNum][1, :]))

            for k in 1:length(hip_height[stepNum][1, :])
                d_ref2foot = hip_height[stepNum][1:3, k] .- [0; 0; d1]
                H = [R_ref2foot d_ref2foot; 0 0 0 1]
                H = inv(H)

                p_l_global[:, k] = [stepL[stepNum][:, k] .+ [0, 0, d4]; 1]
                p_l_local[:, k] = *(H, p_l_global[:, k])

                p_r_global[:, k] = [stepR[stepNum][:, k] .+ [0, 0, d4]; 1]
                p_r_local[:, k] = *(H, p_r_global[:, k])
            end
            push!(stepLocalR, p_r_local)
            push!(stepLocalL, p_l_local)
        end
    return stepLocalR, stepLocalL
end
