"""

A Swing Foot Trajectory as defined in sec. 3.1 of : 
    1.  Khusainov R, Sagitov A, Klimchik A, Magid E. Arbitrary Trajectory Foot Planner for Bipedal Walking: 
        In: Proceedings of the 14th International Conference on Informatics in Control, Automation and Robotics [Internet]. 
        Madrid, Spain: SCITEPRESS - Science and Technology Publications; 2017 [cité 20 déc 2022]. p. 417‑24. 
        Disponible sur: http://www.scitepress.org/DigitalLibrary/Link.aspx?doi=10.5220/0006442504170424


A `SwingFootTrajectory` has 2 differents outputs :
    1.  `stepR` which is an array and which the eleemnt `i` in the array is a coordinate matrix of the global position of the right foot at
        the step `i`
    2.  `stepL` which is an array and which the eleemnt `i` in the array is a coordinate matrix of the global position of the left foot at
        the step `i`

The input of of the `CoMTrajectory` require : 
    1.  `br` is a structure `BipedRobot` which have all the description of the robot
    2.  `fp` is a structure `FootPlanner`which have the position of the right, left and center point of the foot during the path
    3.  `zt` is a structure `ZMPTrajectory` which have the reference ZMP and the time reference

"""
struct SwingFootTrajectory
    stepR::Array
    stepL::Array
end

"""

Function that describe the X or Y motion of the foot during the step as describe in the paper 
"""
function xyFoot(
    br::BipedRobot,
    t::Union{Vector, StepRangeLen},
    xs::Union{Int64, Float64},
    xf::Union{Int64, Float64},
)
    δ = br.δ
    Tstep = br.Tstep
    Tver = br.Tver
    tDSP = δ * Tstep
    return 0.5 * (xf - xs) .* (1 .- cos.(pi * (t .- tDSP) / (Tstep - tDSP - Tver))) .+ xs
end

"""

Function that describe the Z motion of the foot during the step as describe in the paper 
"""
function zFoot(br::BipedRobot, t::Union{Vector, StepRangeLen})
    δ = br.δ
    Tstep = br.Tstep
    tDSP = δ * Tstep
    hstep = br.hstep
    return 0.5 * hstep .* (1 .- cos.(2 .* pi .* (t .- tDSP) ./ (Tstep - tDSP)))
end

"""

Constructor 
"""
function SwingFootTrajectory(;
    br::BipedRobot,
    fp::FootPlanner,
    zt::ZMPTrajectory,
    check::Bool = false,
)
    stepR, stepL = computeSwingFootTrajectory(br, fp, zt)
    if (check)
        stepR_plot = reduce(hcat, stepR)
        stepL_plot = reduce(hcat, stepL)
        step_plot = reduce(vcat, [stepL_plot, stepR_plot])

        plt = plot(;
            title = "Swing Foot Trajectory",
            xlabel = "X[m]",
            ylabel = "Y[m]",
            zlabel = "Z[m]",
            dpi = 600,
        )
        i = 1
        for foot in ["left", "right"]
            plot!(
                step_plot[i + 2 * (i - 1), :],
                step_plot[i + 1 + 2 * (i - 1), :],
                step_plot[i + 2 + 2 * (i - 1), :];
                label = foot,
            )
            i = i + 1
        end
        display(plt)
    end
    return SwingFootTrajectory(stepR, stepL)
end

"""

Compute all trajectories (X, Y, Z) of each step (and both feet) of the motion
"""
function computeSwingFootTrajectory(br::BipedRobot, fp::FootPlanner, zt::ZMPTrajectory)
    isLeftSupport = br.isLeftSupport
    tstart = br.Tdelay
    Twait = br.Twait
    Ts = br.Ts
    Tstep = br.Tstep
    tver = br.Tver
    δ = br.δ

    offset = br.d
    right = fp.right
    left = fp.left
    center = fp.center
    timeVec = zt.timeVec
    ZMP = zt.ZMP

    stepR = Array{Float64}[]
    stepL = Array{Float64}[]

    left_flag = ~isLeftSupport
    tDSP = δ * Tstep

    for stepNum in 1:length(ZMP)
        tstep = timeVec[stepNum]
        if (stepNum < length(ZMP))
            if (left_flag == true)
                lastStep = left[stepNum]        # Position of the last step left
                nextStep = left[stepNum + 1]    # Position of the next step left
            else
                lastStep = right[stepNum]        # Position of the last step right
                nextStep = right[stepNum + 1]    # Position of the next step right
            end
        elseif (stepNum == length(ZMP))
            temp = center[end]                  # Get last position 
            lastFootPos = footPosition(temp[1], temp[2], temp[3], offset, left_flag)
            if (left_flag == true)
                lastStep = left[stepNum]        # Position of the last step left
                nextStep = lastFootPos         # Position of the next step left (beside the right foot)
            else
                lastStep = right[stepNum]        # Position of the last step right
                nextStep = lastFootPos         # Position of the next step right (beside the left foot)
            end
        end

        xf = nextStep[1]       # x-position of the next foot 
        xs = lastStep[1]       # x-position of the last foot 
        yf = nextStep[2]       # y-position of the next foot 
        ys = lastStep[2]       # y-position of the last foot 

        # Taking account of the initial delay 
        if (stepNum == 1)
            delay = tstart + Twait
        else
            delay = tstep[1] - Ts # - Ts because in the way that timeVec are constructed
        end

        # Reset the temp variables 
        stepX = []
        stepY = []
        stepZ = []
        i = 1

        # At the beginning (Double support phase), no mouvement 
        while (tstep[i] <= delay + tDSP)
            stepX = reduce(vcat, [stepX, xs])
            stepY = reduce(vcat, [stepY, ys])
            stepZ = reduce(vcat, [stepZ, 0])
            i += 1
        end

        # Evaluate the position of the foot in local time (only between tDSP < t < t0)
        stepZ = reduce(vcat, [stepZ, zFoot(br, tstep[i:end] .- delay)])

        # Compute the vertical time index
        idxVer = Int(round(tver / Ts))

        # Evaluate the position of the foot in local time (only between tDSP < t <= t0 - tver)
        stepX = reduce(vcat, [stepX, xyFoot(br, tstep[i:(end - idxVer)] .- delay, xs, xf)])
        stepY = reduce(vcat, [stepY, xyFoot(br, tstep[i:(end - idxVer)] .- delay, ys, yf)])

        # Evaluate the position of the foot in local time (only between t0 - tver < t < t0)
        stepX = reduce(vcat, [stepX, xf * ones(idxVer)])
        stepY = reduce(vcat, [stepY, yf * ones(idxVer)])

        # Group up all computated values in one 
        swingFoot = hcat(stepX, stepY, stepZ)'

        if (left_flag == true)  # Left foot = swing foot 
            # define the stand foot as the coordinate of the foot position 
            standFoot =
                vcat(right[stepNum][1:2] .* ones(2, length(tstep)), zeros(1, length(tstep)))
            push!(stepL, swingFoot)
            push!(stepR, standFoot)
        else                    # Right foot = swing foot 
            # define the stand foot as the coordinate of the foot position 
            standFoot =
                vcat(left[stepNum][1:2] .* ones(2, length(tstep)), zeros(1, length(tstep)))
            push!(stepR, swingFoot)
            push!(stepL, standFoot)
        end
        left_flag = ~left_flag             # Change foot swing 
    end
    return stepR, stepL
end
