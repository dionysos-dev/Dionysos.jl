"""

A Walking Optimization perform solve an optimisation problem based on : 
    `x = argmin/argmax f(x)`
    `g(x) <= 0`
    `h(x) = 0`

where `g(x)` and `h(x)` are constraints functions and `x` is an candidate (ZMP based controller)
The methodology used here is to simulate each candidate from a discrete input-space. 

A `WalkingOptimization` has 2 differents outputs :
    1.  `x` which is the value of a candidate 
    2.  `f` which is the objective value of the candidate `x`

"""
mutable struct WalkingOptimization
    # Input space
    Δz::Real
    Tstep::Real
    Lmax::Real
    δ::Real
    Tver::Real
    hstep::Real

    # Auto defined input space
    zc::Real
    previewTime::Real

    # Output space 
    energy::Real
    walkingSpeed::Real
end

"""

Compute the energy from a torque and joint speed vector 
"""
function computeEnergy(torque::Vector{T}, omega::Vector{T}, dt::T) where {T <: Real}
    energy = zero(T)
    n = length(torque)
    power = zeros(T, n)
    torque = abs.(torque)
    omega = abs.(omega)
    for i in 1:(n - 1)
        avg_torque = (torque[i] + torque[i + 1]) / 2
        avg_omega = (omega[i] + omega[i + 1]) / 2
        energy += avg_torque * avg_omega * dt
        power[i] = avg_torque * avg_omega
    end
    return energy, power
end

"""

Return a `BipedRobot` data structure based on a `WalkingOptimization` data structure 
"""
function defineBipedRobot(wo::WalkingOptimization)
    include(joinpath(packagepath(), "param_test.jl"))
    urdfpath() = joinpath(packagepath(), "ZMP_2DBipedRobot.urdf")
    Δz = wo.Δz
    Tstep = wo.Tstep
    Lmax = wo.Lmax
    δ = wo.δ
    Tver = wo.Tver
    hstep = wo.hstep
    zc = wo.zc
    previewTime = wo.previewTime

    A, B, C = cartTableModel(zc, g)
    Ad, Bd, Cd = continuous2discrete(A, B, C, Ts)

    p = rank(Cd)       # rank of Cd
    n = rank(Ad)       # rank of Ad

    # Define cost function weights
    Qe = q_e * eye(p)     # tracking error weight matrix
    Qx = 0.0 * eye(n)       # incremental state weight matrix
    R = r_u * eye(p)       # control vector weight matrix

    visual = URDFVisuals(urdfpath())
    e = root(visual.xdoc)

    L1 = 0.0
    L2 = 0.0
    d = 0.0
    offset_hip_to_motor = 0.0
    offset_ankle_to_foot = 0.0

    for (i, joint) in enumerate(e["joint"])
        if (attribute(joint, "name")) == "l_hip_to_motor"
            s = attribute(joint["origin"]..., "xyz")
            numbers = split(s, " ")
            offset_hip_to_motor = abs(parse(Float64, numbers[3]))
            d = 2 * abs(parse(Float64, numbers[2]))
        end
        if (attribute(joint, "name")) == "l_thigh_link_to_motor"
            s = attribute(joint["origin"]..., "xyz")
            numbers = split(s, " ")
            L1 = abs(parse(Float64, numbers[3]))
        end
        if (attribute(joint, "name")) == "l_ankle"
            s = attribute(joint["origin"]..., "xyz")
            numbers = split(s, " ")
            L2 = abs(parse(Float64, numbers[3]))
        end
    end

    for (i, link) in enumerate(e["link"])
        if (attribute(link, "name") == "l_foot_link")
            s = attribute(link["visual"][1]["geometry"][1]["box"]..., "size")
            numbers = split(s, " ")
            offset_ankle_to_foot = abs(parse(Float64, numbers[3]))
        end
    end

    return BipedRobot(
        Ts,
        zc,
        g,
        Ad,
        Bd,
        Cd,
        zmax,
        Δz,
        L1,
        L2,
        offset_hip_to_motor,
        offset_ankle_to_foot,
        Qx,
        Qe,
        R,
        previewTime,
        Lmax,
        θ_max,
        d,
        initial_position,
        isLeftSupport,
        xPath,
        yPath,
        Tstep,
        Tdelay,
        Twait,
        δ,
        xinit,
        yinit,
        ΔCoMz,
        hstep,
        Tver,
    )
end

"""

Compute `zc` and `previewTime` based on other parameters 
"""
function computeAutoDefineParameters!(wo::WalkingOptimization)
    rs = RobotSimulator(;
        fileName = "ZMP_2DBipedRobot.urdf",
        symbolic = false,
        add_contact_points = true,
        add_gravity = true,
        add_flat_ground = true,
    )

    br = defineBipedRobot(wo)
    tstep = 0:(br.Ts):(br.Tdelay)
    max_height = (br.L_thigh + br.L_leg)
    CoMz = zeros(length(tstep))
    zc_init = 0.21902787835056647
    zc = zc_init - br.Δz

    new_height = zeros(length(tstep))
    coeffCoM = getSplineCoeff(0.0, Float64(Tdelay), zc_init, zc)
    coeffHeight = getSplineCoeff(0.0, Float64(Tdelay), max_height, max_height - br.Δz)
    for idx in eachindex(tstep)
        if (tstep[idx] <= Tdelay)
            new_height[idx] = spline(tstep[idx], coeffHeight)[1]
            CoMz[idx] = spline(tstep[idx], coeffCoM)[1]
        end
    end
    L1 = br.L_thigh
    L2 = br.L_leg
    q = twoLinksInverseKinematics.(zeros(length(tstep)), -new_height, L1, L2)
    qr = reduce(vcat, q)
    ql = reduce(vcat, q)
    qref = [ql[:, 1] qr[:, 1] ql[:, 2] qr[:, 2]]

    ## TODO : change to a adjustable position controller 
    tend = tstep[end]
    Kp = 10000.0
    Ki = 100.0
    Kd = 100.0
    ctrl = true
    Δt = 1e-3
    controller! = trajectory_controller!(rs, tstep, qref, Δt, Kp, Ki, Kd, ctrl)
    RigidBodyDynamics.simulate(rs.state, tend; Δt = Δt, controller!)
    CoMsim = reduce(hcat, rs.CoM)

    wo.zc = CoMsim[3, end]

    br = defineBipedRobot(wo)
    _, _, Gd = computeGains(br)
    Gd = vec(reduce(vcat, Gd))
    indexArray = findall(x -> x >= -1, Gd)
    wo.previewTime = (indexArray[1] + 1) * br.Ts
    return br.Tdelay = wo.previewTime
end

"""

Return the Energy objective function 
"""
function EnergyObjectiveFunction(
    torque::VecOrMat{T},
    omega::VecOrMat{T},
    dt::T,
) where {T <: Real}
    energy = zeros(4)
    power = zeros(4, length(torque[1, :]))
    for actuator in range(1, 4)
        energy[actuator], power[actuator, :] =
            computeEnergy(torque[actuator, :], omega[:, actuator], dt)
        # println("Energy $(actuator - dof_offset) = $(E[actuator - dof_offset]) J")
    end
    total_energy = sum(energy)
    total_power = sum(power[actuator, :] for actuator in 1:4)
    mean_total_power = sum(total_power) / length(power[1, :])
    return total_energy, total_power, mean_total_power
end

"""

Return the torque objective function 
"""
function TorquesObjectiveFunction(tau::VecOrMat{T}, Δt::T) where {T <: Real}
    tau_sqrt = tau .^ 2
    temp = zeros(4)
    Δt = 1e-3
    for i in 1:4
        temp[i] = sum((tau_sqrt[i, 1:(end - 1)] + tau_sqrt[i, 2:end]) / 2 * Δt)
    end
    return sum(temp)
end

"""

Return the speed objective function 
"""
function SpeedObjectiveFunction(speed::Vector{T}) where {T <: Real}
    return sum(speed) / length(speed)
end

"""

Simulate the candidate `x` for 3 step 
"""
function simulate(
    rs::RobotSimulator,
    qref::Matrix,
    tref::Union{Vector, StepRangeLen},
    tend::Float64,
    Kp::Float64,
    Ki::Float64,
    Kd::Float64,
    Δt::Float64,
)
    ctrl = true
    boom = [0, 0]
    foot = [0, 0]
    actuators = [0, 0, 0, 0]
    config = vec(reduce(hcat, [boom; actuators; foot]))
    set_configuration!(rs.state, config)
    zero_velocity!(rs.state)
    controller! = trajectory_controller!(rs, tref, qref, Δt, Kp, Ki, Kd, ctrl)
    ts, qs, vs = RigidBodyDynamics.simulate(rs.state, tend; Δt = Δt, controller!)
    qsim = reduce(hcat, qs)'
    vsim = reduce(hcat, vs)'
    tsim = reduce(vcat, ts)
    torque_sim = reduce(hcat, rs.torques)
    CoMsim = reduce(hcat, rs.CoM)

    return tsim, qsim, vsim, torque_sim, CoMsim
end

"""

Compute the objective value of a candidate `x`
"""
function ZMPbipedObjective(x)
    # TODO : Adapt the controller 
    Kp = 10000.0
    Ki = 100.0
    Kd = 100.0

    # println("x = $(x)")
    Δz = x[1]
    Tstep = x[2]
    Lmax = x[3]
    δ = x[4]
    Tver = x[5]
    hstep = x[6]
    zc = 0.4275
    previewTime = 5
    rs = RobotSimulator(;
        fileName = "ZMP_2DBipedRobot.urdf",
        symbolic = false,
        add_contact_points = true,
        add_gravity = true,
        add_flat_ground = true,
    )

    wo = WalkingOptimization(Δz, Tstep, Lmax, δ, Tver, hstep, zc, previewTime, NaN, NaN)
    computeAutoDefineParameters!(wo)
    plot_result = false
    br = defineBipedRobot(wo)
    pc = PreviewController(; br = br, check = plot_result)
    fp = FootPlanner(; br = br, check = plot_result)
    zt = ZMPTrajectory(; br = br, fp = fp, check = plot_result)
    ct = CoMTrajectory(; br = br, pc = pc, zt = zt, check = plot_result)
    sf = SwingFootTrajectory(; br = br, fp = fp, zt = zt, check = plot_result)
    ik = InverseKinematics(; br = br, fp = fp, ct = ct, sf = sf, check = plot_result)

    qr = ik.q_r
    ql = ik.q_l
    timeVec = zt.timeVec
    tplot = reduce(vcat, timeVec)
    qref = [ql[:, 1] qr[:, 1] ql[:, 2] qr[:, 2]]
    tend_sim = br.Tdelay + 3.0 * br.Tstep + br.Twait
    # tend_sim = tplot[end]
    tsim, qsim, vsim, torque_sim, CoMsim =
        simulate(rs, qref, tplot, tend_sim, Kp, Ki, Kd, 1e-3)

    tau_hip_max = 3.4
    tau_knee_max = 2.5

    omega_hip_max = 3.14
    omega_knee_max = 5.23

    tau_constraints =
        (maximum(abs.(torque_sim[3, :])) < tau_hip_max) &&
        (maximum(abs.(torque_sim[4, :])) < tau_hip_max) &&
        (maximum(abs.(torque_sim[5, :])) < tau_knee_max) &&
        (maximum(abs.(torque_sim[6, :])) < tau_knee_max)
    omega_constraints =
        (maximum(abs.(vsim[:, 3])) < omega_hip_max) &&
        (maximum(abs.(vsim[:, 4])) < omega_hip_max) &&
        (maximum(abs.(vsim[:, 5])) < omega_knee_max) &&
        (maximum(abs.(vsim[:, 6])) < omega_knee_max)

    Δt = 1e-3
    tstart = br.Tdelay + br.Tstep + br.Twait
    tend = tstart + br.Tstep
    i_start = Int(round(tstart / Δt + 1))
    i_end = Int(round(tend / Δt + 1))
    L2 = br.L_leg
    L1 = br.L_thigh

    if (tau_constraints && omega_constraints)
        walkSpeed = SpeedObjectiveFunction(vsim[i_start:i_end, 1])
        energyPerStep, _, _ = EnergyObjectiveFunction(
            torque_sim[3:6, i_start:i_end],
            vsim[i_start:i_end, 3:6],
            Δt,
        )
        int_tau = TorquesObjectiveFunction(torque_sim[3:6, i_start:i_end], Δt) / Lmax
        energyPerStep = energyPerStep / Lmax
    else
        walkSpeed = -Inf
        energyPerStep = Inf
        int_tau = Inf
    end
    g = [
        x[3] - sqrt(
            2 * L2^2 * (1 - cos(acos(((L1 + L2 - x[1])^2 - L1^2 - L2^2) / (2 * L1 * L2)))),
        ),
    ]
    # h = [x[2]-1; x[4]-0.2; x[5]; x[6]-0.01]
    h = [0.0]
    return [walkSpeed, energyPerStep, int_tau], g, h
end
