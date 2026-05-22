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

include(joinpath(@__DIR__, "..", "src", "RBDCustomSimulation.jl"))
import .RBDCustomSimulation

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

Base.@kwdef struct RobotDomainConfig{X1, X2, U1, U2}
    x_lb::X1
    x_ub::X2
    u_lb::U1
    u_ub::U2
end

function default_robot_domain()
    disc_steps = [fill(π / 180, 3)..., fill(0.075, 3)...]

    return RobotDomainConfig(;
        x_lb = SVector(-12π / 180, 0.0, 0.0, -0.6, -0.15, -0.15) .- disc_steps,
        x_ub = SVector(0.0, 12π / 180, 14π / 180, 0.15, 0.6, 0.6) .+ disc_steps,
        u_lb = SVector(-3.0, -3.0, -3.0),
        u_ub = SVector(3.0, 2.0, 3.0),
    )
end

function warmup_robot_problem!(;
    robot_urdf = joinpath(@__DIR__, "..", "deps", "ZMP_2DBipedRobot_nodamping.urdf"),
    tstep = 1e-1,
    domain = default_robot_domain(),
    Δt_simu = 1e-4,
    simulator::Symbol = :history,
)
    sys = system(;
        robot_urdf = robot_urdf,
        tstep = tstep,
        domain = domain,
        Δt_simu = Δt_simu,
        simulator = simulator,
    )

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

mutable struct VoltageController{Q, G}
    u::SVector{3, Float64}
    q_ref::Q
    GR::G
    Kp::Float64
    ktp::Float64
    Kvp::Float64
    τc_u::Float64
    ddl::Int
end

@inline function (vc::VoltageController)(τ, t, state)
    τ .= 0.0

    q = configuration(state)
    q̇ = velocity(state)

    idx_lo = length(q) - 3 - vc.ddl

    q4 = q[idx_lo + 3]

    v1 = q̇[idx_lo]
    v2 = q̇[idx_lo + 1]
    v3 = q̇[idx_lo + 2]
    v4 = q̇[idx_lo + 3]

    GR1 = vc.GR[1]
    GR2 = vc.GR[2]
    GR3 = vc.GR[3]
    GR4 = vc.GR[4]

    ω1 = v1 * GR1
    ω2 = v2 * GR2
    ω3 = v3 * GR3
    ω4 = v4 * GR4

    pwm = (vc.q_ref[1] - q4) * (4095.0 / (2π) * vc.Kp)
    pwm_sat = clamp(pwm, -885.0, 885.0)
    uK = pwm_sat * (12.0 / 885.0)

    τ[idx_lo] = vc.u[1] * GR1 * vc.ktp - ω1 * GR1 * vc.Kvp - sign(ω1) * GR1 * vc.τc_u

    τ[idx_lo + 1] = vc.u[2] * GR2 * vc.ktp - ω2 * GR2 * vc.Kvp - sign(ω2) * GR2 * vc.τc_u

    τ[idx_lo + 2] = vc.u[3] * GR3 * vc.ktp - ω3 * GR3 * vc.Kvp - sign(ω3) * GR3 * vc.τc_u

    τ[idx_lo + 3] = uK * GR4 * vc.ktp - ω4 * GR4 * vc.Kvp - sign(ω4) * GR4 * vc.τc_u

    return nothing
end

function system(;
    robot_urdf = joinpath(@__DIR__, "..", "deps", "ZMP_2DBipedRobot_nodamping.urdf"),
    tstep = 1e-1,
    domain = default_robot_domain(),
    Δt_simu = 1e-4,
    simulator::Symbol = :history,
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

    controllers_per_thread = [
        VoltageController(
            SVector{3, Float64}(0.0, 0.0, 0.0),
            q_ref,
            GR,
            Kp,
            ktp,
            Kvp,
            τc_u,
            ddl,
        ) for _ in 1:Threads.nthreads()
    ]

    custom_sims_per_thread = Ref{Any}(nothing)

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

    function init_custom_sims!(cache)
        if custom_sims_per_thread[] === nothing
            custom_sims_per_thread[] = [
                RBDCustomSimulation.CachedFinalStateSimulator(
                    cache.states_per_thread[i],
                    controllers_per_thread[i],
                ) for i in 1:Threads.nthreads()
            ]
        end

        return custom_sims_per_thread[]
    end

    function vectorFieldBipedRobot(x::SVector{6, Float64}, u::SVector{3, Float64})
        tid = Threads.threadid()

        cache = _get_robot_data(robot_urdf)
        state = cache.states_per_thread[tid]
        controller! = controllers_per_thread[tid]

        controller!.u = u

        q, q̇ = build_state(x)

        set_configuration!(state, q)
        set_velocity!(state, q̇)
        set_additional_state!(state, zero(RigidBodyDynamics.additional_state(state)))

        if simulator == :custom
            sims = init_custom_sims!(cache)
            sim = sims[tid]

            q_end, v_end = RBDCustomSimulation.simulate_final_state_cached!(
                sim,
                Δt_dionysos,
                Δt_simu;
                check_finite = false,
            )
        else
            q_end, v_end = RBDCustomSimulation.simulate_final_state!(
                state,
                Δt_dionysos,
                controller!;
                Δt = Δt_simu,
                backend = :history,
            )
        end

        return SVector{6, Float64}(
            q_end[3],
            q_end[4],
            q_end[5],
            v_end[3],
            v_end[4],
            v_end[5],
        )
    end

    state_space = UT.HyperRectangle(domain.x_lb, domain.x_ub)
    input_space = UT.HyperRectangle(domain.u_lb, domain.u_ub)

    return MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem(
        vectorFieldBipedRobot,
        6,
        3,
        state_space,
        input_space,
    )
end

function problem(;
    robot_urdf = joinpath(@__DIR__, "..", "deps", "ZMP_2DBipedRobot_nodamping.urdf"),
    tstep = 1e-1,
    domain = default_robot_domain(),
    Δt_simu = 1e-4,
    simulator = :history,
)
    sys = system(;
        robot_urdf = robot_urdf,
        tstep = tstep,
        domain = domain,
        Δt_simu = Δt_simu,
        simulator = simulator,
    )

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
