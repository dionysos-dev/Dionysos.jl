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
        ẋboom = Lthigh * qd[i1] * cos(q[i1]) + Lleg * (qd[i1] + qd[i2]) * cos(q[i1] + q[i2])
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
    state_space = UT.box(state_lower_bounds, state_upper_bounds)

    input_lower_bounds = [-3, -3, -3, -3]
    input_upper_bounds = [3, 3, 3, 3]
    input_space = UT.box(input_lower_bounds, input_upper_bounds)

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
    return PR.AlternatingSimulationProblem(sys, nothing)
end

#-----------------------------------------------------------
# - Additional pruning constraints for a biped gait tube   -
#-----------------------------------------------------------

const LTHIGH = 0.20125
const LLEG = 0.172

deg2rad(d) = d * π / 180

"""
    stance_foot(x; tol = 1e-9)

Determine which leg is the stance leg from a simple geometric estimate
of the foot height.

The vertical positions of the feet are approximated as

    zl = Lthigh*cos(LH) + Lleg*cos(LK + LH)
    zr = Lthigh*cos(RH) + Lleg*cos(RK + RH)

The leg whose foot is lower is assumed to be in stance.

Returns:
- `:left` if the left foot is clearly lower,
- `:right` if the right foot is clearly lower,
- `:ambiguous` if both feet have nearly equal height.
"""
@inline function stance_foot(x::SVector{8, Float64}; tol = 1e-9)
    LH, RH, LK, RK = x[1], x[2], x[3], x[4]

    zl = Lthigh * cos(LH) + Lleg * cos(LK + LH)
    zr = Lthigh * cos(RH) + Lleg * cos(RK + RH)

    if zl >= zr + tol
        return :left
    elseif zr >= zl + tol
        return :right
    else
        return :ambiguous
    end
end

"""
    in_gait_tube(x)

Return `true` if the state `x` belongs to a desired *walking tube* in the
state space of the robot.

This predicate is used as a pruning rule to discard states that do not
resemble a plausible walking configuration.

The following constraints are enforced:

1. **Hip symmetry**
   The sum of hip angles and hip angular velocities should remain small,
   encouraging an anti-symmetric gait (one leg forward while the other
   moves backward).

2. **Stance/swing configuration**
   One knee should be nearly straight (stance leg) while the other
   remains bent (swing leg), with a minimum angular gap between them.

3. **Velocity bounds**
   Knee angular velocities must remain within reasonable limits, with
   stricter bounds for the stance leg.

Returns `true` if all constraints are satisfied.
"""
function in_gait_tube(x::AbstractVector{<:Real})
    LH, RH, LK, RK, dLH, dRH, dLK, dRK = x

    hip_sum_tol = deg2rad(10)
    dhip_sum_tol = 0.4

    stance_knee_max = deg2rad(12)
    swing_knee_min = deg2rad(10)
    knee_gap_min = deg2rad(8)

    stance_dk_max = 0.3
    swing_dk_max = 0.5

    # Enforce approximate anti-symmetry of hips
    if abs(LH + RH) > hip_sum_tol
        ;
        return false;
    end
    if abs(dLH + dRH) > dhip_sum_tol
        ;
        return false;
    end

    # Detect stance/swing configuration
    left_stance =
        (LK <= stance_knee_max) && (RK >= swing_knee_min) && ((RK - LK) >= knee_gap_min)
    right_stance =
        (RK <= stance_knee_max) && (LK >= swing_knee_min) && ((LK - RK) >= knee_gap_min)
    if !(left_stance || right_stance)
        ;
        return false;
    end

    # Velocity limits on knees
    if left_stance
        if abs(dLK) > stance_dk_max || abs(dRK) > swing_dk_max
            ;
            return false;
        end
    else
        if abs(dRK) > stance_dk_max || abs(dLK) > swing_dk_max
            ;
            return false;
        end
    end

    return true
end

"""
    input_allowed(x, u)

Return `true` if control input `u` is admissible in state `x`.

The following heuristics are enforced:

1. **Avoid strong symmetric hip torques**
   Both hips should not push strongly in the same direction.

2. **Balance hip torques**
   The sum of hip torques should remain moderate.

3. **Swing/stance knee behavior**
   The swing leg should receive stronger torque than the stance leg.
   This encourages lifting the swing leg rather than pushing with the
   stance leg.

Returns `true` if the input is considered physically plausible for walking.
"""
function input_allowed(x::AbstractVector{<:Real}, u::AbstractVector{<:Real})
    LH, RH, LK, RK, dLH, dRH, dLK, dRK = x
    uLH, uRH, uLK, uRK = u

    # Avoid both hips pushing strongly in the same direction
    if (sign(uLH) == sign(uRH)) && (abs(uLH) >= 2) && (abs(uRH) >= 2)
        return false
    end

    # Keep hip torques roughly balanced
    if abs(uLH + uRH) > 2
        return false
    end

    # Encourage stronger torque on the swing knee than on the stance knee
    left_is_stance = LK < RK

    if left_is_stance
        if abs(uLK) > 1
            ;
            return false;
        end
        if abs(uRK) < 1
            ;
            return false;
        end
        if abs(uRK) < abs(uLK)
            ;
            return false;
        end
    else
        if abs(uRK) > 1
            ;
            return false;
        end
        if abs(uLK) < 1
            ;
            return false;
        end
        if abs(uLK) < abs(uRK)
            ;
            return false;
        end
    end

    return true
end

end # module
