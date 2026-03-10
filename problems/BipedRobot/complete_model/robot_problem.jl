module RobotProblem

using MathematicalSystems
using LinearAlgebra, StaticArrays
using RigidBodyDynamics
using Base.Threads

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const PR = DI.Problem

include(joinpath(@__DIR__, "..", "src", "RS_tools.jl"))
import .RS_tools

# ----------------------------------------------------------------------
# Visualization
# ----------------------------------------------------------------------
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

# ----------------------------------------------------------------------
# System
# ----------------------------------------------------------------------
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
    tmp_state = MechanismState(mechanism)
    println("n_pos: ", num_positions(tmp_state))
    println("n_vel: ", num_velocities(tmp_state))

    Δt_simu = 1e-4
    Δt_dionysos = tstep

    # --- One MechanismState per thread (thread-safe) ---
    states_per_thread = [MechanismState(mechanism) for _ in 1:Threads.nthreads()]

    # ---------------- Motor params ----------------
    HGR = 353.5
    KGR = 212.6
    ktp = 0.395 / HGR
    Kvp = 1.589 / (HGR * HGR)
    τc_u = 0.065 / HGR
    GR = SVector{4, Float64}(HGR, HGR, KGR, KGR)
    ddl = 2

    # Controller factory: no shared scratch
    function voltage_controller!(u::SVector{4, Float64})
        function controller!(τ, t, state)
            τ .= 0.0

            qd = velocity(state)
            idx_lo = length(qd) - 3 - ddl
            idx_hi = length(qd) - ddl

            current_qd = @view qd[idx_lo:idx_hi]
            ω = current_qd .* GR

            τ_m = u .* GR .* ktp .- ω .* GR .* Kvp .- sign.(ω) .* GR .* τc_u
            τ[idx_lo:idx_hi] .= τ_m
            return nothing
        end
    end

    # ---------------- Geometry params ----------------
    Lthigh = 0.20125
    Lleg = 0.172
    Init_offset = -0.0006559432

    # Pure mapping x(8) -> (q(8), qd(8))
    # x = [LH RH LK RK dLH dRH dLK dRK]
    function fill_state(x::SVector{8, Float64})
        q = @SVector [0.0, 0.0, x[1], x[2], x[3], x[4], 0.0, 0.0]
        qd = @SVector [0.0, 0.0, x[5], x[6], x[7], x[8], 0.0, 0.0]

        zl = Lthigh * cos(q[3]) + Lleg * cos(q[5] + q[3])
        zr = Lthigh * cos(q[4]) + Lleg * cos(q[6] + q[4])

        q2 = max(zl, zr) - Lthigh - Lleg + Init_offset

        # identify contact leg
        i1 = zl > zr ? 3 : 4
        i2 = zl > zr ? 5 : 6

        xboom = Lthigh * sin(q[i1]) + Lleg * sin(q[i2] + q[i1])
        ẋboom =
            Lthigh * qd[i1] * cos(q[i1]) + Lleg * (qd[i1] + qd[i2]) * cos(q[i1] + q[i2])
        żboom =
            -(Lthigh * qd[i1] * sin(q[i1]) + Lleg * (qd[i1] + qd[i2]) * sin(q[i1] + q[i2]))

        q1 = xboom
        qd1 = ẋboom
        qd2 = żboom

        # keep feet horizontal
        q7 = -(q[3] + q[5])
        q8 = -(q[4] + q[6])
        qd7 = -(qd[3] + qd[5])
        qd8 = -(qd[4] + qd[6])

        q = SVector{8, Float64}(q1, q2, q[3], q[4], q[5], q[6], q7, q8)
        qd = SVector{8, Float64}(qd1, qd2, qd[3], qd[4], qd[5], qd[6], qd7, qd8)

        return q, qd
    end

    # Thread-safe vector field
    function vectorFieldBipedRobot(x, u)
        tid = Threads.threadid()
        state = states_per_thread[tid]

        xS = SVector{8, Float64}(x)
        uS = SVector{4, Float64}(u)

        q, qd = fill_state(xS)
        set_configuration!(state, q)
        set_velocity!(state, qd)

        controller! = voltage_controller!(uS)
        _, qs, vs =
            RigidBodyDynamics.simulate(state, Δt_dionysos, controller!; Δt = Δt_simu)

        q_end = qs[end]
        v_end = vs[end]

        # return next reduced state (same as original intent)
        return SVector{8, Float64}(
            q_end[3],
            q_end[4],
            q_end[5],
            q_end[6],
            v_end[3],
            v_end[4],
            v_end[5],
            v_end[6],
        )
    end

    # State/input spaces (keep your existing bounds)
    state_lower_bounds = [-0.5, -0.5, -0.2, -0.2, -0.8, -0.8, -0.8, -0.8]
    state_upper_bounds = [0.5, 0.5, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8]
    state_space = UT.HyperRectangle(state_lower_bounds, state_upper_bounds)

    input_lower_bounds = [-3, -3, -3, -3]
    input_upper_bounds = [3, 3, 3, 3]
    input_space = UT.HyperRectangle(input_lower_bounds, input_upper_bounds)

    return MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem(
        vectorFieldBipedRobot,
        8,
        4,
        state_space,
        input_space,
    )
end

function problem(;
    tstep = 2e-2,
    robot_urdf = joinpath(@__DIR__, "..", "deps/ZMP_2DBipedRobot_nodamping.urdf"),
)
    sys = system(; tstep = tstep, robot_urdf = robot_urdf)
    return PR.EmptyProblem(sys, nothing)
end

end # module
