module RobotProblem

using StaticArrays
using MathematicalSystems
using RigidBodyDynamics
using Base.Threads

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem

include(joinpath(@__DIR__, "..", "src", "RSCore.jl"))
import .RSCore

const _robot_cache = Ref{Any}(nothing)

function _get_robot_data(robot_urdf::AbstractString)
    if _robot_cache[] === nothing || _robot_cache[].robot_urdf != robot_urdf
        rs = RSCore.RobotSimulator(;
            fileName = robot_urdf,
            symbolic = false,
            add_contact_points = true,
            add_gravity = true,
            add_flat_ground = true,
        )

        mechanism = rs.mechanism

        _robot_cache[] = (
            robot_urdf = robot_urdf,
            rs = rs,
            mechanism = mechanism,
            states_per_thread = [MechanismState(mechanism) for _ in 1:Threads.nthreads()],
        )
    end

    return _robot_cache[]
end

function clear_robot_cache!()
    _robot_cache[] = nothing
    return nothing
end

"""
    warmup_robot_problem!(; robot_urdf, tstep, Δt_simu)

Force robot-data initialization and one dynamics evaluation on the current process.
Useful before distributed benchmarks to separate setup/JIT from abstraction time.
"""
function warmup_robot_problem!(;
    robot_urdf = joinpath(@__DIR__, "..", "deps", "ZMP_2DBipedRobot_nodamping.urdf"),
    tstep = 1e-1,
    Δt_simu = 1e-4,
)
    sys = system(; tstep = tstep, robot_urdf = robot_urdf, Δt_simu = Δt_simu)

    x = SVector{6, Float64}(
        -0.10471975511965978,
        0.10471975511965978,
        0.10471975511965978,
        -0.45,
        0.0,
        0.0,
    )

    u = SVector{3, Float64}(0.0, 0.0, 1.0)

    sys.f(x, u)
    return nothing
end

struct VoltageController{U, Q, G}
    u::U
    q_ref::Q
    GR::G
    Kp::Float64
    ktp::Float64
    Kvp::Float64
    τc_u::Float64
    ddl::Int
end

function (vc::VoltageController)(τ, t, state)
    τ .= 0.0

    q = configuration(state)
    q̇ = velocity(state)

    idx_lo = length(q) - 3 - vc.ddl
    idx_hi = length(q) - vc.ddl

    current_q = @view q[idx_lo:idx_hi]
    current_q̇ = @view q̇[idx_lo:idx_hi]

    ω = current_q̇ .* vc.GR

    PWM = (vc.q_ref .- current_q[4]) .* (4095.0 / (2π) * vc.Kp)
    PWM_sat = clamp.(PWM, -885.0, 885.0)
    u_K = PWM_sat .* (12.0 / 885.0)

    U_tot = SVector{4, Float64}(vc.u[1], vc.u[2], vc.u[3], u_K[1])

    τ_m = U_tot .* vc.GR .* vc.ktp .- ω .* vc.GR .* vc.Kvp .- sign.(ω) .* vc.GR .* vc.τc_u

    τ[idx_lo:idx_hi] .= τ_m

    return nothing
end

function system(;
    tstep = 5e-1,
    robot_urdf = joinpath(@__DIR__, "..", "deps", "ZMP_2DBipedRobot_nodamping.urdf"),
    Δt_simu = 1e-4,
)
    Δt_dionysos = tstep

    HGR = 353.5
    KGR = 212.6
    ktp = 0.395 / HGR
    Kvp = 1.589 / (HGR * HGR)
    τc_u = 0.065 / HGR
    GR = SVector{4, Float64}(HGR, HGR, KGR, KGR)
    Kp = 900.0 / 128.0
    ddl = 2

    Lthigh = 0.20125
    Lleg = 0.172
    Init_offset = -0.0006559432

    q_ref = SVector{1, Float64}(0.0)

    function build_state(x::SVector{6, Float64})
        q3, q4, q5 = x[1], x[2], x[3]
        q̇3, q̇4, q̇5 = x[4], x[5], x[6]

        zl = Lthigh * cos(q3) + Lleg * cos(q5 + q3)
        zr = Lthigh * cos(q4) + Lleg * cos(q4)

        q2 = max(zl, zr) - Lthigh - Lleg + Init_offset

        if zl > zr
            xboom = Lthigh * sin(q3) + Lleg * sin(q5 + q3)

            ẋboom = Lthigh * q̇3 * cos(q3) + Lleg * (q̇3 + q̇5) * cos(q3 + q5)

            żboom = -(Lthigh * q̇3 * sin(q3) + Lleg * (q̇3 + q̇5) * sin(q3 + q5))
        else
            xboom = Lthigh * sin(q4) + Lleg * sin(q4)

            ẋboom = Lthigh * q̇4 * cos(q4) + Lleg * q̇4 * cos(q4)

            żboom = -(Lthigh * q̇4 * sin(q4) + Lleg * q̇4 * sin(q4))
        end

        q7 = -(q3 + q5)
        q8 = -q4

        q̇7 = -(q̇3 + q̇5)
        q̇8 = -q̇4

        q = SVector{8, Float64}(xboom, q2, q3, q4, q5, 0.0, q7, q8)

        q̇ = SVector{8, Float64}(ẋboom, żboom, q̇3, q̇4, q̇5, 0.0, q̇7, q̇8)

        return q, q̇
    end

    function vectorFieldBipedRobot(x::SVector{6, Float64}, u::SVector{3, Float64})
        cache = _get_robot_data(robot_urdf)
        state = cache.states_per_thread[Threads.threadid()]

        q, q̇ = build_state(x)

        set_configuration!(state, q)
        set_velocity!(state, q̇)

        controller! = VoltageController(u, q_ref, GR, Kp, ktp, Kvp, τc_u, ddl)

        _, qs, vs =
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

    disc_steps = [fill(π / 180, 3)..., fill(0.075, 3)...]

    state_lower_bounds = [-12 * π / 180, 0, 0, -0.6, -0.15, -0.15] .- disc_steps

    state_upper_bounds = [0, 12 * π / 180, 14 * π / 180, 0.15, 0.6, 0.6] .+ disc_steps

    state_space = UT.HyperRectangle(state_lower_bounds, state_upper_bounds)

    input_lower_bounds = [-3, -3, -3]
    input_upper_bounds = [3, 2, 3]
    input_space = UT.HyperRectangle(input_lower_bounds, input_upper_bounds)

    return MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem(
        vectorFieldBipedRobot,
        6,
        3,
        state_space,
        input_space,
    )
end

function problem(;
    tstep = 1e-1,
    robot_urdf = joinpath(@__DIR__, "..", "deps", "ZMP_2DBipedRobot_nodamping.urdf"),
    Δt_simu = 1e-4,
)
    sys = system(; tstep = tstep, robot_urdf = robot_urdf, Δt_simu = Δt_simu)
    return PR.AlternatingSimulationProblem(sys, nothing)
end

const LTHIGH = 0.20125
const LLEG = 0.172

deg2rad(d) = d * π / 180

@inline function stance_foot(x::SVector{6, Float64}; tol = 1e-9)
    LH, RH, LK = x[1], x[2], x[3]

    zl = LTHIGH * cos(LH) + LLEG * cos(LK + LH)
    zr = LTHIGH * cos(RH) + LLEG * cos(RH)

    if zl >= zr + tol
        return :left
    elseif zr >= zl + tol
        return :right
    else
        return :ambiguous
    end
end

function in_gait_tube(x::AbstractVector{<:Real})
    LH, RH, LK, dLH, dRH, dLK = x

    hip_sum_tol = deg2rad(10)
    dhip_sum_tol = 0.4
    swing_knee_min = deg2rad(5)
    swing_knee_max = deg2rad(20)
    dhip_max = 0.8
    dknee_max = 0.5

    abs(LH + RH) <= hip_sum_tol || return false
    abs(dLH + dRH) <= dhip_sum_tol || return false
    swing_knee_min <= LK <= swing_knee_max || return false
    abs(dLH) <= dhip_max || return false
    abs(dRH) <= dhip_max || return false
    abs(dLK) <= dknee_max || return false

    return true
end

function input_allowed(x::AbstractVector{<:Real}, u::AbstractVector{<:Real})
    uLH, uRH, uLK = u

    if (sign(uLH) == sign(uRH)) && (abs(uLH) >= 2) && (abs(uRH) >= 2)
        return false
    end

    abs(uLH + uRH) <= 2 || return false
    1 <= abs(uLK) <= 3 || return false

    return true
end

end # module
