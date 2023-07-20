mutable struct CLcontrol
    ΔZMPx::Deque{Float64}
    joint_state::Deque{Vector}
    ∫e_xdt_prev::Float64
    previous_state::Vector
end

function oneStepPC!(
    br::BipedRobot,
    pc::PreviewController,
    zt::ZMPTrajectory,
    CLc::CLcontrol,
)
    Ad = br.Ad
    Bd = br.Bd
    Cd = br.Cd
    Ts = br.Ts
    previewTime = br.previewTime
    ZMP = zt.ZMP
    zc = br.zc
    g = br.g

    timeVec = zt.timeVec
    longTime = reduce(vcat, timeVec)
    longZMP = reduce(hcat, ZMP)
    kmax = length(longTime)
    h = Int(round(previewTime / Ts) + 1)
    Gx = pc.Gx
    Gi = pc.Gi
    Gd = pc.Gd

    xx_prev = CLc.previous_state[:, 1]
    px_prev = 0#Cd * xx_prev

    # Preview action 
    ux_prev = 0
    h = Int(round(previewTime / Ts) + 1)
    # println("preview length : $(h)")
    for (id, dzmp) in enumerate(CLc.ΔZMPx)
        # println("Gd[id] : $(Gd[id])")
        ux_prev = ux_prev - Gd[id] .* dzmp
    end
    # println("ux_prev = $(ux_prev)")
    # State feedback action 
    ux_sf = -(Gx * xx_prev)[1]
    # println("ux_sf = $(ux_sf)")
    # Integral action
    e_x = px_prev[1] - popfirst!(CLc.ΔZMPx)
    ∫e_xdt = CLc.∫e_xdt_prev - e_x
    CLc.∫e_xdt_prev = ∫e_xdt
    # println("∫e_xdt = $(∫e_xdt)")

    ux_pc = Gi[1] * ∫e_xdt + ux_prev + ux_sf

    xx_next = Ad * xx_prev + Bd * ux_pc

    return CLc.previous_state = xx_next
end

function ClosedLoopZMPcontroller!(
    br::BipedRobot,
    rs::RobotSimulator,
    Δt::Float64,
    Kp::Float64,
    Ki::Float64,
    Kd::Float64,
    debugMode::Bool = false,
)
    time_idx = 1
    stepNum = 1
    left_support = br.isLeftSupport
    sim_index = 0
    mechanism = rs.mechanism
    state = rs.state
    dynamics_results = DynamicsResult(rs.mechanism)
    d1 = br.offset_hip_to_motor
    d4 = br.offset_ankle_to_foot
    ddl = 2 # For 2 non-actuated foot 
    qref = [0; 0; 0; 0]
    qactual = [0; 0; 0; 0]
    desired_q = [0; 0; 0; 0]

    # Init measure data
    nddl = num_velocities(mechanism)
    tau = zeros(nddl)
    com = center_of_mass(state).v
    push!(rs.torques, tau)
    push!(rs.CoM, com)
    push!(rs.ZMP, com[1:2])

    # Init the controller 
    plot_result = false
    pc = PreviewController(; br = br, check = plot_result)
    fp = FootPlanner(; br = br, check = plot_result)
    zt = ZMPTrajectory(; br = br, fp = fp, check = plot_result)
    sf = SwingFootTrajectory(; br = br, fp = fp, zt = zt, check = plot_result)
    ct = CoMTrajectory(; br = br, pc = pc, zt = zt, check = plot_result)
    dk = hip2foot(br)
    pid = PID(Kp, Ki, Kd)
    ∫edt = zeros(4)

    # Get the actual reference vectors
    tref = reduce(vcat, zt.timeVec)
    ZMPref = reduce(hcat, zt.ZMP)
    center = fp.center
    stepR = reduce(hcat, sf.stepR)
    stepL = reduce(hcat, sf.stepL)
    hip_height = reduce(hcat, ct.hip_height)

    # Declare the FIFO Data Structures 
    fifo_ΔZMPx = Deque{Float64}()
    fifo_states = Deque{Vector}()
    fifo_ref_states = Deque{Vector}()

    # Set initial values for the FIFO Data Structures
    h = Int(round(br.previewTime / br.Ts) + 1)
    com = center_of_mass(state).v
    for i in 1:(h - 1)
        push!(fifo_ΔZMPx, com[1])
        push!(fifo_states, [0; 0; 0; 0])
        push!(fifo_ref_states, [0; 0; 0; 0])
    end
    previous_state = zeros(3)

    # New struct for the closed loop control 
    clc = CLcontrol(fifo_ΔZMPx, fifo_states, 0.0, previous_state)

    function controller!(τ, t, state)
        if (t >= tref[time_idx])

            # Update the ZMP measured and reference
            p = measureZMP(rs, dynamics_results, state)
            px_ref = ZMPref[1, 1]

            # Put the new data in a FIFO DataStructures
            ΔZMPx = p[1] - px_ref
            qactual = configuration(state)[(end - 2 - 3):(end - 2)]
            push!(fifo_ΔZMPx, ΔZMPx)
            push!(fifo_states, qactual)
            push!(fifo_ref_states, qref)

            # Look for support foot 
            if t >= zt.timeVec[stepNum][end]
                stepNum = stepNum + 1
            end

            # Compute the variation of the Pelvis structure in x-direction
            oneStepPC!(br, pc, zt, clc)
            xx = clc.previous_state[:]
            ΔCoMx = xx[1]

            # Get the actual configuration 
            qactual = popfirst!(fifo_states)

            # Seek for the reference foot (support foot)
            if left_support * mod(stepNum, 2) == true # Left foot as reference
                q_dk = [qactual[1] qactual[3]] # -> Left leg joints 
            else # Right foot as reference
                q_dk = [qactual[2] qactual[4]] # -> Right leg joints 
            end

            # Forward Kinematics to compute the actual Pelvis configuration using the support foot as reference
            xhip = dk[1].(q_dk[1], q_dk[2])
            yhip = dk[2].(q_dk[1], q_dk[2]) .+ d4

            # Compute the corrected postion of the Pelvis configuration
            Δhip_coord = [ΔCoMx; 0; yhip - hip_height[3, time_idx]]
            hip_coord = [xhip; 0; yhip - d1] .- Δhip_coord  # We don't need to correct the ZMP in y-direction (plannar robot)

            # Get the local Pelvis configuration 
            θ = center[stepNum][3]
            R_ref2foot = [
                cos(θ) -sin(θ) 0
                sin(θ) cos(θ) 0
                0 0 1
            ]
            d_ref2foot = hip_coord
            H = [R_ref2foot d_ref2foot; 0 0 0 1]
            H = inv(H)

            # Swing foot trajectories
            p_l_global = [stepL[:, time_idx] .+ [0, 0, d4]; 1]
            p_l_local = *(H, p_l_global)

            p_r_global = [stepR[:, time_idx] .+ [0, 0, d4]; 1]
            p_r_local = *(H, p_r_global)

            # Inverse Kinematics 
            qr = twoLinksInverseKinematics.(p_r_local[1], p_r_local[3])
            ql = twoLinksInverseKinematics.(p_l_local[1], p_l_local[3])
            qref = [ql[1]; qr[1]; ql[2]; qr[2]]

            if t >= br.previewTime
                # Update the next time
                time_idx = time_idx + 1
            end
            # desired_q = qref 
            desired_q = popfirst!(fifo_ref_states)
            if debugMode && t > 0
                println("+=========================DEBUG MODE===========================+")
                println("t = $(t)")
                println("tref[time_idx] = $(tref[time_idx])")
                println("stepNum = $(stepNum)")
                # println("fifo_ΔZMPx = $(fifo_ΔZMPx)")
                # println("fifo_states = $(fifo_states)")
                # println("fifo_ref_states = $(fifo_ref_states)")
                println("qactual = $(qactual)")
                println("desired_q = $(desired_q)")
                println("qref = $(qref)")
                println("ΔCoMx = $(ΔCoMx)")
                println("hip_coord = $(hip_coord)")
                println("Δhip_coord = $(Δhip_coord)")
                println("q_dk = $(q_dk)")
                println("xhip = $(xhip)")
                println("yhip = $(yhip)")
                println("p_l_local = $(p_l_local)")
                println("p_r_local = $(p_r_local)")

                println("qr = $(qr)")
                println("ql = $(ql)")

                println("τ = $(τ[end-3-ddl : end-ddl])")
                println("Press Enter to continue...")
                readline()  # Wait for user input
                println("Program resumed.")
            end
        end
        # Reset torque
        for joint in joints(mechanism)
            τ[velocity_range(state, joint)] .= 0 # no control
        end

        # Position control 
        v̇ = copy(velocity(state))
        actual_q = configuration(state)[(end - 3 - ddl):(end - ddl)]

        Δq = desired_q .- actual_q
        Δq̇ = vec(zeros(4, 1) .- velocity(state)[(end - 3 - ddl):(end - ddl)])
        v̇_new = pid_control!(pid, Δq, Δq̇, ∫edt, Δt)
        v̇[(end - 3 - ddl):(end - ddl)] .= v̇_new

        # Torque in the constrained space 
        newτ = inverse_dynamics(state, v̇) #- dynamics_bias(state)
        τ[(end - 3 - ddl):(end - ddl)] .= newτ[(end - 3 - ddl):(end - ddl)]

        # Get the torque and CoM measurements 
        if (t >= sim_index * Δt && t <= tref[end] + br.previewTime)
            sim_index = sim_index + 1
            com = center_of_mass(state).v
            push!(rs.CoM, com)
            push!(rs.torques, τ)
        end
        return nothing
    end
end

function measureZMP(
    rs::RobotSimulator,
    dynamics_results::DynamicsResult,
    state::MechanismState,
)
    RigidBodyDynamics.contact_dynamics!(dynamics_results, state)
    foot_link = findbody(rs.mechanism, "r_foot_link")
    sensor_right = RigidBodyDynamics.contact_wrench(dynamics_results, foot_link)
    foot_link = findbody(rs.mechanism, "l_foot_link")
    sensor_left = RigidBodyDynamics.contact_wrench(dynamics_results, foot_link)
    contact_torque_right = sensor_right.angular
    contact_force_right = sensor_right.linear
    contact_torque_left = sensor_left.angular
    contact_force_left = sensor_left.linear

    d = 0.0045 # half foot hight
    pr_x_right =
        (-contact_torque_right[2] - contact_force_right[1] * d) / contact_force_right[3]
    pr_y_right =
        (contact_torque_right[1] - contact_force_right[2] * d) / contact_force_right[3]

    pr_x_left =
        (-contact_torque_left[2] - contact_force_left[1] * d) / contact_force_left[3]
    pr_y_left = (contact_torque_left[1] - contact_force_left[2] * d) / contact_force_left[3]

    px =
        (pr_x_right * contact_force_right[3] + pr_x_left * contact_force_left[3]) /
        (contact_force_right[3] + contact_force_left[3])
    py =
        (pr_y_right * contact_force_right[3] + pr_y_left * contact_force_left[3]) /
        (contact_force_right[3] + contact_force_left[3])
    if (contact_force_left[3] <= 0 && contact_force_right[3] >= 0)
        p = [pr_x_right, pr_y_right]
    elseif (contact_force_left[3] >= 0 && contact_force_right[3] <= 0)
        p = [pr_x_left, pr_y_left]
    else
        p = [px, py]
    end
    push!(rs.ZMP, p)
    return p
end
