"""

A CoM Trajectory Generator perform a classical solver to a discrete-time system :

    `x(k+1) = Ad x(k) + Bd u(k)`
    `y(k+1) = Cd x(k+1)`

Where `u(k)` is given by a preview controller. 

A `CoMTrajectory` has 2 differents outputs :
    1.  `CoM` which is an array and which the eleemnt `i` in the array is a computed CoM vector for the step `i`
    2.  `p` which is an array and which the eleemnt `i` in the array is a computed ZMP vector for the step `i`
    3.  `robot_height` which is an array and which the eleemnt `i` in the array is a height shiffted CoM vector for the step `i`
        This will be usefull to the inverse kinematics block because we don't neglect the mass of the legs 

The input of of the `CoMTrajectory` require : 
    1.  `br` is a structure `BipedRobot` which have all the description of the robot
    2.  `pc` is a structure `PreviewController` which have the gains of the preview controller
    3.  `fp` is a structure `FootPlanner`which have the position of the right, left and center point of the foot during the path


"""

struct CoMTrajectory
    CoM::Array
    p::Array
    hip_height::Array
end

"""

Constructor 
"""
function CoMTrajectory(;
    br::BipedRobot,
    pc::PreviewController,
    zt::ZMPTrajectory,
    check::Bool = false,
)
    CoM, ZMP, hip_height = computeCoMTrajectory(br, pc, zt)
    if (check)
        CoMplot = reduce(hcat, CoM)
        ZMPplot = reduce(hcat, ZMP)
        timeVec = zt.timeVec
        ZMP_ref = zt.ZMP
        ZMP_refplot = reduce(hcat, ZMP_ref)
        tplot = reduce(vcat, timeVec)
        plt = plot(;
            title = "CoM Trajectory",
            xlabel = "t [s]",
            ylabel = "Y[m]",
            layout = (2, 1),
            dpi = 600,
        )
        for (i, dir) in enumerate(["X[m]", "Y[m]"])
            plot!(
                plt[i],
                tplot,
                ZMP_refplot[i, :];
                ylabel = dir,
                label = "Reference ZMP",
                lw = 2,
            )
            plot!(plt[i], tplot, ZMPplot[i, :]; ylabel = dir, label = "ZMP", lw = 2)
            plot!(plt[i], tplot, CoMplot[i, :]; ylabel = dir, label = "CoM", lw = 2)
        end
        display(plt)
    end
    return CoMTrajectory(CoM, ZMP, hip_height)
end

"""

Solve the discrete time system in 2D : only on the XY plane
"""
function compute2DCoMTrajectory(br::BipedRobot, pc::PreviewController, zt::ZMPTrajectory)
    Ad = br.Ad
    Bd = br.Bd
    Cd = br.Cd
    Ts = br.Ts
    previewTime = br.previewTime

    ZMP = zt.ZMP
    timeVec = zt.timeVec
    longTime = reduce(vcat, timeVec)
    longZMP = reduce(hcat, ZMP)
    kmax = length(longTime)
    h = Int(round(previewTime / Ts) + 1)
    Gx = pc.Gx
    Gi = pc.Gi
    Gd = pc.Gd

    xx = br.xinit
    xy = br.yinit
    intEi_k_1 = [0; 0]

    p = [
        xx[1, 1] - zc / g * xx[3, 1]     # Initial ZMP on x component 
        xy[1, 1] - zc / g * xy[3, 1]     # Initial ZMP on y component    
    ]

    # Simulation of the pattern generator with preview control 
    # This is done using x(k+1) = Ad * x(k) + Bd * u(k)
    # with u(k) = -Gx * x(k) - Gi * e(k) - \sum_{i = 1}^{i = N_l} Gd(i) * r(i)
    for k in 1:(kmax - 1)

        # Compute the error 
        e_k = p[:, k] - longZMP[:, k]

        # Compute the integral action 
        intEi_k = intEi_k_1 - Gi .* e_k #.* Ts;
        intEi_k_1 = intEi_k

        # Compute the preview action 
        preview_k = [0.0; 0.0]
        for l in 1:h
            if (k + l < kmax)
                preview_k = preview_k - Gd[l] .* longZMP[:, k + l]
            else    # if we exceed the vector of reference, we took the last value 
                preview_k = preview_k - Gd[l] .* longZMP[:, kmax]
            end
        end

        # Compute the control input for both x-direction and y-direction 
        u_k = [
            intEi_k[1] .- Gx * xx[:, k] .+ preview_k[1]
            intEi_k[2] .- Gx * xy[:, k] .+ preview_k[2]
        ]
        # Compute the next state vector 
        xx_k_1 = Ad * xx[:, k] + Bd * u_k[1]
        xy_k_1 = Ad * xy[:, k] + Bd * u_k[2]

        xx = reduce(hcat, [xx, xx_k_1])
        xy = reduce(hcat, [xy, xy_k_1])

        # Compute the next controlled output
        p_k_1 = [Cd * xx[:, k + 1]; Cd * xy[:, k + 1]]
        p = reduce(hcat, [p, p_k_1])
    end
    CoM = Matrix(transpose(reduce(hcat, [xx[1, :], xy[1, :]])))
    return CoM, p
end

"""

Construct a 3D version of the CoM. And store the CoM into array to get the vector for each step. 

The z position of the CoM is computed by a constant. At the beginning (step 1), we just impose a third order polynomial behavior of the z component :

    `CoMz(t) = a_0 + a_1 * t + a_2 * t^2 + a_3 * t^3 `

Where `t  ∈ [0; Tdelay + Tstep]` and a_i  are constants 
"""
function computeCoMTrajectory(br::BipedRobot, pc::PreviewController, zt::ZMPTrajectory)
    zc = br.zc
    Δz = br.Δz
    L1 = br.L_thigh
    L2 = br.L_leg
    d1 = br.offset_hip_to_motor
    d4 = br.offset_ankle_to_foot
    ΔCoMz = br.ΔCoMz
    max_height = (L1 + L2 + d1 + d4)

    Tdelay = br.Tdelay
    longCoM, longp = compute2DCoMTrajectory(br, pc, zt)

    p = Array{Float64}[]
    CoM = Array{Float64}[]
    hip_height = Array{Float64}[]

    init = 1
    stepNum = 1
    ZMP = zt.ZMP
    timeVec = zt.timeVec

    for stepNum in 1:length(ZMP)
        tstep = timeVec[stepNum]
        len = length(tstep)
        if (stepNum == 1)
            CoMz = zeros(1, length(tstep))
            new_height = zeros(1, length(tstep))
            coeffCoM = getSplineCoeff(0.0, Float64(Tdelay), zc + ΔCoMz, zc)
            coeffHeight = getSplineCoeff(0.0, Float64(Tdelay), max_height, max_height - Δz)
            for idx in eachindex(tstep)
                if (tstep[idx] <= Tdelay)
                    CoMz[idx] = spline(tstep[idx], coeffCoM)[1]
                    new_height[idx] = spline(tstep[idx], coeffHeight)[1]
                else
                    CoMz[idx] = zc
                    new_height[idx] = -Δz .+ max_height
                end
            end
        else
            CoMz = zc * ones(1, len)
            new_height = (-Δz .+ max_height) * ones(1, len)
        end
        tempCoM = longCoM[:, init:(init + len - 1)]
        temph = vcat(tempCoM, new_height)
        tempCoM = vcat(tempCoM, CoMz)
        tempZMP = longp[:, init:(init + len - 1)]

        push!(CoM, tempCoM)
        push!(hip_height, temph)
        push!(p, tempZMP)
        init = init + len
    end
    return CoM, p, hip_height
end
