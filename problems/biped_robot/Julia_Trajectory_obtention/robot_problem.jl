module RobotProblem

using MathematicalSystems
using LinearAlgebra, StaticArrays
using RigidBodyDynamics
using Base.Threads

# include Dionysos
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const PB = DI.Problem
const ST = DI.System
const SY = DI.Symbolic

# include the tools for the simulator from src
include(joinpath(@__DIR__, "..", "src", "RS_tools.jl"))
import .RS_tools

function get_visualization_tool(;
    robot_urdf = joinpath(@__DIR__, "..", "deps/ZMP_2DBipedRobot_nodamping.urdf"),
)
    rs = RS_tools.RobotSimulator(;
        fileName = robot_urdf,
        symbolic = false,
        add_contact_points = true,
        add_gravity = true,
        add_flat_ground = true,
    )
    vis = RS_tools.set_visulalizer(; mechanism = rs.mechanism, fileName = robot_urdf)
    return rs, vis
end

function system(;
    tstep = 5e-1,
    robot_urdf = joinpath(@__DIR__, "..", "deps/ZMP_2DBipedRobot_nodamping.urdf"),
)
    rs = RS_tools.RobotSimulator(;
        fileName = robot_urdf,
        symbolic = false,
        add_contact_points = true,
        add_gravity = true,
        add_flat_ground = true,
    )

    mechanism = rs.mechanism
    tmp_state = MechanismState(mechanism) # just to inspect dimensions
    n_pos = num_positions(tmp_state)
    n_vel = num_velocities(tmp_state)
    Δt_simu = 1e-4           # Simulation step
    Δt_dionysos = tstep      # Time step used by Dionysos

    println("n_pos: ", n_pos)
    println("n_vel: ", n_vel)

     # --- One MechanismState per thread (thread-local mutable state) ---
    states_per_thread = [MechanismState(mechanism) for _ in 1:Threads.nthreads()]

    ## Motor parameters ##
    HGR = 353.5                   # Hip gear-ratio
    KGR = 212.6                   # Knee gear-ratio
    ktp  = 0.395 / HGR            # Torque constant with respect to the voltage [Nm/V] 
    Kvp  = 1.589 / (HGR * HGR)    # Viscous friction constant [Nm*s/rad] (linked to motor speed)
    τc_u  = 0.065 / HGR           # Dry friction torque [Nm]
    GR = SVector{4,Float64}(HGR, HGR, KGR, KGR)  # Gear ratios
    Kp = 900.0 / 128.0            # DXL controller gain
    ddl  = 2                      # constant used in indexing

    ## Robots geometry parameters ##
    Lthigh = 0.20125
    Lleg = 0.172
    Hip_offset = 0.04025
    Foot_height = 0.009
    Init_offset = -0.0006559432

    # --- Controller factory: all locals, no shared scratch ---
    function voltage_controller!(u::SVector{3,Float64}, q_ref::SVector{1,Float64})
        function controller!(τ, t, state)
            τ .= 0.0

            # indices: last 4 actuated joints plus offset ddl
            q  = configuration(state)
            qd = velocity(state)
            idx_lo = length(q) - 3 - ddl
            idx_hi = length(q) - ddl

            current_q  = @view q[idx_lo:idx_hi]
            current_qd  = @view qd[idx_lo:idx_hi]
            ω = current_qd .* GR

            # DXL controller on the right knee
            # (q_ref is 1x1 here, you might adapt if needed)
            PWM = (q_ref .- current_q[4]) .* (4095.0 / (2π) * Kp) # Only true because profile acceleration and profile velocity are null
            PWM_sat = clamp.(PWM, -885.0, 885.0) # Apply_saturation
            u_K = PWM_sat .* (12.0 / 885.0)

            # Total motor commands for 4 actuators
            U_tot = SVector{4,Float64}(u[1], u[2], u[3], u_K[1])

            τ_m = U_tot .* GR .* ktp .- ω .* GR .* Kvp .- sign.(ω) .* GR .* τc_u
            
            τ[idx_lo:idx_hi] .= τ_m
            return nothing
        end
    end

    # --- Fill state (pure) ---
    function fill_state!(x)
        # x = [LH RH LK RK LA RA] / [position, velocity] for 3 actuated joints
        q  = @SVector [0.0, 0.0, x[1], x[2], x[3], 0.0, 0.0, 0.0]
        qd = @SVector [0.0, 0.0, x[4], x[5], x[6], 0.0, 0.0, 0.0]

        # heights of the two legs (double pendulum)
        zl = Lthigh * cos(q[3]) + Lleg * cos(q[5] + q[3])
        zr = Lthigh * cos(q[4]) + Lleg * cos(q[6] + q[4])

        # boom z position: most extended leg is in contact
        q2 = max(zl, zr) - Lthigh - Lleg + Init_offset
        q  = SVector{8,Float64}(q[1], q2, q[3], q[4], q[5], q[6], q[7], q[8])

        # identify contact leg
        i1 = zl > zr ? 3 : 4
        i2 = zl > zr ? 5 : 6

        # speed equations of the double pendulum
        xboom = Lthigh * sin(q[i1]) + Lleg * sin(q[i2] + q[i1])
        ẋboom = Lthigh * qd[i1] * cos(q[i1]) + Lleg * (qd[i1] + qd[i2]) * cos(q[i1] + q[i2])
        żboom = -(Lthigh * qd[i1] * sin(q[i1]) + Lleg * (qd[i1] + qd[i2]) * sin(q[i1] + q[i2]))

        q1  = xboom
        qd1 = ẋboom
        qd2 = żboom

        # adjust the angular speed of the feet to remain mostly horizontal
        qd7 = -(qd[3] + qd[5])
        qd8 = -(qd[4] + qd[6])

        # The x position is set to 0, the feet are kept to the ground 
        q  = SVector{8,Float64}(q1, q2, q[3], q[4], q[5], q[6], -(q[3] + q[5]), -(q[4] + q[6]))
        qd = SVector{8,Float64}(qd1, qd2, qd[3], qd[4], qd[5], qd[6], qd7, qd8)

        return q, qd
    end

    # --- Thread-safe vector field: one MechanismState per thread ---
    function vectorFieldBipedRobot(x, u)
        tid = Threads.threadid()
        state = states_per_thread[tid]

        # Step 1: build full state
        q, q̇ = fill_state!(x)

        # Step 2: set mechanism state
        set_configuration!(state, q)
        set_velocity!(state, q̇)

        # Step 3: simulate
        q_ref = SVector{1}(0.0)
        controller! = voltage_controller!(u, q_ref)
        ts, qs, vs =
            RigidBodyDynamics.simulate(state, Δt_dionysos, controller!; Δt = Δt_simu)

        # Only final joint states are used
        q_end = qs[end]
        v_end = vs[end]

        x_next = SVector{6,Float64}(q_end[3], q_end[4], q_end[5],
                                    v_end[3], v_end[4], v_end[5])
        return x_next
    end

    # --- State and input spaces ---
    disc_steps = [fill(π/180, 3)..., fill(0.075, 3)...]
    state_lower_bounds = [-12*π/180, 0, 0, -0.6, -0.15, -0.15] .- disc_steps
    state_upper_bounds = [0, 12*π/180, 14*π/180, 0.15, 0.6, 0.6] .+ disc_steps
    state_space = UT.HyperRectangle(state_lower_bounds, state_upper_bounds)

    input_lower_bounds = [-3, -3, -3]
    input_upper_bounds = [3, 2, 3]
    input_space = UT.HyperRectangle(input_lower_bounds, input_upper_bounds)

    sys = MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem(
        vectorFieldBipedRobot,
        6, # state dimension (right knee not included)
        3, # input dimension (right knee not included)
        state_space,
        input_space,
    )

    return sys
end

function problem(;
    tstep = 1e-1,
    robot_urdf = joinpath(@__DIR__, "..", "deps/ZMP_2DBipedRobot_nodamping.urdf"),
)
    sys = system(; tstep = tstep, robot_urdf = robot_urdf)
    return PB.EmptyProblem(sys, nothing)
end

end # module
