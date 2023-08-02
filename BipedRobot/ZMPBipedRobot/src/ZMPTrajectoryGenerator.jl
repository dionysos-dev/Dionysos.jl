"""

A ZMP Trajectory Generator as conform to the definition of the ZMP stated by : 
    1. Vukobratović M, Stepanenko J. On the stability of anthropomorphic systems. Mathematical Biosciences. oct 1972;15(1‑2):1‑37. 

A `ZMPTrajectory` has 2 differents outputs :
    1.  `ZMP` which is an array and which the eleemnt `i` in the array is a reference ZMP vector for the step `i`
    2.  `timeVec` which is an array and which the eleemnt `i` in the array is a time vector for the step `i`

The input of of the `ZMPTrajectory` require : 
    1.  `br` is a structure `BipedRobot` which have all the description of the robot
    2.  `fp` is a structure `FootPlanner` which have the position of the right, left and center point of the foot during the path
"""
struct ZMPTrajectory
    ZMP::Array
    timeVec::Array
end

"""

Constructor 
"""
function ZMPTrajectory(; br::BipedRobot, fp::FootPlanner, check::Bool = false)
    ZMP, timeVec = computeZMPTrajectory(br::BipedRobot, fp::FootPlanner)
    if (check)
        ZMPplot = reduce(hcat, ZMP)
        plt =
            plot(; title = "ZMP Trajectory", xlabel = "X [m]", ylabel = "Y [m]", dpi = 600)
        xpath = br.xPath
        ypath = br.yPath
        right_plot = reduce(hcat, fp.right)
        left_plot = reduce(hcat, fp.left)
        center = fp.center
        plot!(xpath, ypath; label = "Reference path", lw = 2)
        scatter!(left_plot[1, :], left_plot[2, :]; shape = :rect, label = "Left")
        scatter!(right_plot[1, :], right_plot[2, :]; shape = :rect, label = "Right")
        scatter!(
            getindex.(center, 1),
            getindex.(center, 2);
            label = "Center",
            mc = :black,
            markershape = :xcross,
        )
        plot!(ZMPplot[1, :], ZMPplot[2, :]; label = "Reference ZMP", lw = 2)
        display(plt)
    end
    return ZMPTrajectory(ZMP, timeVec)
end

"""

Contruct the reference ZMP trajectory by fusing 3 part :
    1. Initial delay at the beginning 
    2. Double support phase (both fooots are on the ground) 
    3. Swing phase or Single support phase (one foot on the ground)

First, the initial delay is taken into account.
Then, we have to loop trought all the steps and for each step, we fuse a double support phase and a single support phase. 

The double support phase is describe by a cubic spline with zero derivative at the beginning and the end.

The single support phase is describe by a constant value (define by the position of the support foot) 
during the remmaining time step. 

"""
function computeZMPTrajectory(br::BipedRobot, fp::FootPlanner)
    isLeftSupport = br.isLeftSupport
    δ = br.δ
    Tstep = br.Tstep
    Ts = br.Ts
    Tdelay = br.Tdelay
    Twait = br.Twait
    right = fp.right
    left = fp.left
    center = fp.center

    DSPtime = (1:(δ * (Tstep / Ts))) * Ts
    SSPtime = (1:((1 - δ) * (Tstep / Ts))) * Ts
    stepNum = 1                                             # step index
    left_flag = ~isLeftSupport                                # Initiate the move foot 

    timeVec = Array{Float64}[]
    ZMP = Array{Float64}[]

    # Delay vector 
    nDelay = Int(round((Tdelay + Twait) / Ts))
    t = (0:nDelay) * Ts
    zmp = [center[1][1]; center[1][2]] .* ones(2, nDelay + 1)

    # Initiate 
    lastZMP = [center[1][1]; center[1][2]]
    lastTime = t[end]

    while (stepNum <= length(center))
        if (left_flag == true)
            nextZMP = right[stepNum][1:2] # ZMP ref of the next step 
        else
            nextZMP = left[stepNum][1:2] # ZMP of the next step       
        end
        left_flag = ~left_flag    # Change support foot 

        # Creating a cubic spline in the X direction for connecting 2 differents ZMPx 
        xSplineX = DSPtime
        coeff = getSplineCoeff(xSplineX[1], xSplineX[end], lastZMP[1], nextZMP[1])                                     # X position of the point for the spline (time data) 
        ySplineX = spline(xSplineX, coeff)

        # Same of y-coordinate
        xSplineY = DSPtime
        coeff = getSplineCoeff(xSplineX[1], xSplineX[end], lastZMP[2], nextZMP[2])                                     # X position of the point for the spline (time data) 
        ySplineY = spline(xSplineY, coeff)

        # Contruct the first part with a spline for the ZMP value 
        t = vcat(t, lastTime .+ DSPtime)  # Construct the time vector 
        lastTime = t[end]                  # get the last value of t
        zmp = vcat(zmp', hcat(ySplineX, ySplineY)) # Construct the ZMP vector 
        zmp = zmp'

        # Construct the last part for the time vector 
        t = vcat(t, lastTime .+ SSPtime)

        # For the last part, the ZMP stay constant 
        zmp = vcat(zmp', nextZMP' .* ones(length(SSPtime), 2))# Construct the ZMP vector 
        zmp = zmp'

        push!(timeVec, t)      # add to the variable 
        push!(ZMP, zmp)        # add to the variable 

        lastZMP = zmp[1:2, end]     # Get last value 
        lastTime = t[end]               # Get last value 
        t = []                 # Reset the vector 
        zmp = reshape([], 2, 0) # Reset the vector 

        stepNum = stepNum + 1  # next step 
    end
    return ZMP, timeVec
end
