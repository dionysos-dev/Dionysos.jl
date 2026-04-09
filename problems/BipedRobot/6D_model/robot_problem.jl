module RobotProblem

using MathematicalSystems
using LinearAlgebra, StaticArrays
using RigidBodyDynamics
using Base.Threads

# include Dionysos
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const SY = DI.Symbolic

# include the tools for the simulator from src
include(joinpath(@__DIR__, "..", "src", "RS_tools.jl"))
import .RS_tools

const _robot_cache = Ref{Any}(nothing)

function _get_robot_data(robot_urdf)
    if _robot_cache[] === nothing
        rs = RS_tools.RobotSimulator(;
            fileName = robot_urdf,
            symbolic = false,
            add_contact_points = true,
            add_gravity = true,
            add_flat_ground = true,
        )
        mechanism = rs.mechanism
        _robot_cache[] = (
            rs = rs,
            mechanism = mechanism,
            states_per_thread = [MechanismState(mechanism) for _ in 1:Threads.nthreads()],
        )
    end
    return _robot_cache[]
end

# The robot state is
#     x = [LH, RH, LK, RK, dLH, dRH, dLK, dRK]
# where:
# - `LH`, `RH` are the hip angles,
# - `LK`, `RK` are the knee angles,
# - `dLH`, `dRH`, `dLK`, `dRK` are the corresponding angular velocities.
# The control vector is
#     u = [uLH, uRH, uLK, uRK]
# where each component represents a torque applied at the corresponding joint.
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
    HGR = 353.5                  # Hip gear-ratio
    KGR = 212.6                  # Knee gear-ratio
    ktp = 0.395 / HGR            # Torque constant with respect to the voltage [Nm/V] 
    Kvp = 1.589 / (HGR * HGR)    # Viscous friction constant [Nm*s/rad] (linked to motor speed)
    τc_u = 0.065 / HGR           # Dry friction torque [Nm]
    GR = SVector{4, Float64}(HGR, HGR, KGR, KGR)  # Gear ratios
    Kp = 900.0 / 128.0           # DXL controller gain
    ddl = 2                      # constant used in indexing

    ## Robots geometry parameters ##
    Lthigh = 0.20125
    Lleg = 0.172
    Hip_offset = 0.04025
    Foot_height = 0.009
    Init_offset = -0.0006559432

    # --- Controller factory: all locals, no shared scratch ---
    function voltage_controller!(u::SVector{3, Float64}, q_ref::SVector{1, Float64})
        function controller!(τ, t, state)
            τ .= 0.0

            # indices: last 4 actuated joints plus offset ddl
            q = configuration(state)
            q̇ = velocity(state)
            idx_lo = length(q) - 3 - ddl
            idx_hi = length(q) - ddl

            current_q = @view q[idx_lo:idx_hi]
            current_q̇ = @view q̇[idx_lo:idx_hi]
            ω = current_q̇ .* GR

            # DXL controller on the right knee
            # (q_ref is 1x1 here, you might adapt if needed)
            PWM = (q_ref .- current_q[4]) .* (4095.0 / (2π) * Kp) # Only true because profile acceleration and profile velocity are null
            PWM_sat = clamp.(PWM, -885.0, 885.0) # Apply_saturation
            u_K = PWM_sat .* (12.0 / 885.0)

            # Total motor commands for 4 actuators
            U_tot = SVector{4, Float64}(u[1], u[2], u[3], u_K[1])

            τ_m = U_tot .* GR .* ktp .- ω .* GR .* Kvp .- sign.(ω) .* GR .* τc_u

            τ[idx_lo:idx_hi] .= τ_m
            return nothing
        end
    end

    """
    x ∈ R^8 = [LH, RH, LK, RK, dLH, dRH, dLK, dRK]
    Ordre des q (8) : [x z LH RH LK RK LA RA]
    """
    function fill_state!(x)
        q = vcat(zeros(2), x[1:3], zeros(3)) # q[3:6] = LH RH LK RK ; q[7:8]=LA RA
        q̇ = vcat(zeros(2), x[4:6], zeros(3)) # q̇[3:6] = dLH dRH dLK dRK

        # heights of the two legs (double pendulum)
        zl = Lthigh * cos(q[3]) + Lleg * cos(q[5] + q[3])
        zr = Lthigh * cos(q[4]) + Lleg * cos(q[6] + q[4])  # RK=0

        # boom z position: most extended leg is in contact
        q2 = max(zl, zr) - Lthigh - Lleg + Init_offset
        q = SVector{8, Float64}(q[1], q2, q[3], q[4], q[5], q[6], q[7], q[8])

        # identify contact leg
        i1 = zl > zr ? 3 : 4
        i2 = zl > zr ? 5 : 6

        # speed equations of the double pendulum
        xboom = Lthigh * sin(q[i1]) + Lleg * sin(q[i2] + q[i1])
        ẋboom =
            Lthigh * q̇[i1] * cos(q[i1]) + Lleg * (q̇[i1] + q̇[i2]) * cos(q[i1] + q[i2])
        żboom =
            -(Lthigh * q̇[i1] * sin(q[i1]) + Lleg * (q̇[i1] + q̇[i2]) * sin(q[i1] + q[i2]))

        q1 = xboom
        q̇1 = ẋboom
        q̇2 = żboom

        # adjust the angular speed of the feet to remain mostly horizontal
        q̇7 = -(q̇[3] + q̇[5])
        q̇8 = -(q̇[4] + q̇[6])

        # The x position is set to 0, the feet are kept to the ground 
        q = SVector{8, Float64}(
            q1,
            q2,
            q[3],
            q[4],
            q[5],
            q[6],
            -(q[3] + q[5]),
            -(q[4] + q[6]),
        )
        q̇ = SVector{8, Float64}(q̇1, q̇2, q̇[3], q̇[4], q̇[5], q̇[6], q̇7, q̇8)

        return q, q̇
    end

    function vectorFieldBipedRobot(x, u)
        cache = _get_robot_data(robot_urdf)
        state = cache.states_per_thread[Threads.threadid()]

        q, q̇ = fill_state!(x)
        set_configuration!(state, q)
        set_velocity!(state, q̇)

        q_ref = SVector{1}(0.0)
        controller! = voltage_controller!(u, q_ref)
        ts, qs, vs =
            RigidBodyDynamics.simulate(state, Δt_dionysos, controller!; Δt = Δt_simu)

        q_end = qs[end]
        v_end = vs[end]

        return SVector{6, Float64}(
            q_end[3],
            q_end[4],
            q_end[5],
            v_end[3],
            v_end[4],
            v_end[5],
        )
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
    return PR.AlternatingSimulationProblem(sys, nothing)
end

#-----------------------------------------------------------
# - Additional pruning constraints for the reduced 6D model -
#-----------------------------------------------------------

const LTHIGH = 0.20125
const LLEG = 0.172

deg2rad(d) = d * π / 180

"""
    stance_foot(x; tol = 1e-9)

Determine which leg is lower from a simple geometric estimate of foot height
for the reduced 6D model

    x = [LH, RH, LK, dLH, dRH, dLK]

This reduced model does not include the right knee explicitly, so the right knee
is assumed to remain at `RK = 0`.

The foot heights are approximated as

    zl = LTHIGH*cos(LH) + LLEG*cos(LK + LH)
    zr = LTHIGH*cos(RH) + LLEG*cos(RH)

Returns:
- `:left` if the left foot is clearly lower,
- `:right` if the right foot is clearly lower,
- `:ambiguous` otherwise.
"""
@inline function stance_foot(x::SVector{6, Float64}; tol = 1e-9)
    LH, RH, LK = x[1], x[2], x[3]

    zl = LTHIGH * cos(LH) + LLEG * cos(LK + LH)
    zr = LTHIGH * cos(RH) + LLEG * cos(RH)   # RK assumed to be 0

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

Return `true` if the reduced 6D state

    x = [LH, RH, LK, dLH, dRH, dLK]

belongs to a desired walking tube.

Since the reduced model includes only one knee (`LK`), this predicate no longer
tries to distinguish between two fully modeled stance/swing legs. Instead, it
enforces a simplified gait pattern:

1. **Hip anti-symmetry**
   The two hips should roughly mirror each other.

2. **Swing-knee configuration**
   The modeled knee should remain bent enough to represent a swing leg, but not
   excessively bent.

3. **Velocity bounds**
   Hip and knee angular velocities should remain within reasonable ranges.

This predicate is intended as a conservative pruning rule for abstraction.
"""
function in_gait_tube(x::AbstractVector{<:Real})
    LH, RH, LK, dLH, dRH, dLK = x

    hip_sum_tol = deg2rad(10)
    dhip_sum_tol = 0.4

    swing_knee_min = deg2rad(5)
    swing_knee_max = deg2rad(20)

    dhip_max = 0.8
    dknee_max = 0.5

    # Approximate anti-symmetry of hips
    if abs(LH + RH) > hip_sum_tol
        return false
    end
    if abs(dLH + dRH) > dhip_sum_tol
        return false
    end

    # The modeled knee should represent a plausible swing knee
    if LK < swing_knee_min || LK > swing_knee_max
        return false
    end

    # Velocity bounds
    if abs(dLH) > dhip_max || abs(dRH) > dhip_max || abs(dLK) > dknee_max
        return false
    end

    return true
end

"""
    input_allowed(x, u)

Return `true` if the reduced 3D control input

    u = [uLH, uRH, uLK]

is admissible at the reduced 6D state

    x = [LH, RH, LK, dLH, dRH, dLK].

This reduced model has no right-knee torque, so the admissibility conditions
are adapted accordingly.

The following heuristics are enforced:

1. **Avoid strong symmetric hip torques**
   Both hips should not push strongly in the same direction.

2. **Balance hip torques**
   The sum of hip torques should remain moderate.

3. **Swing-knee actuation**
   The modeled knee torque should remain large enough to support swing-leg
   motion, but not excessively large.
"""
function input_allowed(x::AbstractVector{<:Real}, u::AbstractVector{<:Real})
    LH, RH, LK, dLH, dRH, dLK = x
    uLH, uRH, uLK = u

    # Avoid both hips pushing strongly in the same direction
    if (sign(uLH) == sign(uRH)) && (abs(uLH) >= 2) && (abs(uRH) >= 2)
        return false
    end

    # Keep hip torques roughly balanced
    if abs(uLH + uRH) > 2
        return false
    end

    # Keep knee torque in a plausible swing range
    if abs(uLK) < 1 || abs(uLK) > 3
        return false
    end

    return true
end

end # module
