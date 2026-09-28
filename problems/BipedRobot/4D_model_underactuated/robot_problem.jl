# Underactuated four-link biped, velocity-controlled design model.
#
# Unlike the fully actuated 4-D model (`../4D_model/robot_problem.jl`, `ẋ = u`), the
# grounded ankle `q₁` is *not* actuated: only the three motor joints (grounded knee
# `q₂`, hip `q₃`, swing knee `q₄`) track commanded velocities, and the ankle obeys
# the first, unactuated Euler–Lagrange row of the multibody model. The state is
# the velocity-controlled design state of `BipedFourLink`,
#
#     x = [q̇₁, q₁, q₂, q₃, q₄],   u = [v₂, v₃, v₄],
#
# with `ẋ = BipedFourLink.velocity_controlled_derivative(model, x, u)` (commands held
# constant over a sampling period, i.e. `input_accelerations = 0`). `q₁` falls like
# an inverted pendulum; the controller keeps it in the domain only through the
# configuration of the motor joints and the Coriolis coupling.
#
# The dynamics are nonlinear in `(q̇₁, q₁)`, so the exact-lattice argument of the
# fully actuated model holds for the three motor coordinates only: the abstraction
# built with `CENTER_SIMULATION` is **not** a certified over-approximation here.
# The closed-loop trajectory of the concrete system is what validates a controller.
module UnderactuatedRobotProblem

using LinearAlgebra
using StaticArrays
using MathematicalSystems
import LazySets
import Plots
import RigidBodyDynamics as RBD

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping

include(joinpath(@__DIR__, "BipedFourLink.jl"))
import .BipedFourLink as BFL

const URDF_FILE = joinpath(@__DIR__, "four_link_biped.urdf")

"""
    build_model(urdf_file = URDF_FILE)

The `BipedFourLink` locked-ankle design model of the four-link biped.
"""
build_model(urdf_file::AbstractString = URDF_FILE) = BFL.build_model(urdf_file)

# Index of each coordinate in the state `x = [q̇₁, q₁, q₂, q₃, q₄]`.
const IDX_Q̇1 = 1
const IDX_Q = SVector(2, 3, 4, 5)

joint_angles(x) = SVector(x[2], x[3], x[4], x[5])

# ============================================================
# Geometry & kinematics
# ============================================================

"""
    Geometry{T}

Segment lengths of the four-link biped, matching `four_link_biped.urdf`: tibias
(ankle → knee, knee → foot) of `tibia` and femurs (knee → hip) of `femur`.
[`check_against_module`](@ref) verifies the forward kinematics below against the
`RigidBodyDynamics` model.
"""
Base.@kwdef struct Geometry{T <: Real}
    tibia::T = 0.172
    femur::T = 0.202
end

default_geometry() = Geometry{Float64}()

"""
    cartesian_coordinates(q, geometry)

World `(x, z)` positions of the joints for `q = [q₁, q₂, q₃, q₄]`, grounded ankle
at the origin: `(grounded_foot, grounded_knee, hip, swing_knee, swing_foot)`.

Absolute segment angles from the upward vertical (positive towards `+x`): grounded
tibia `q₁`, grounded femur `q₁ + q₂`; the swing femur `φ₃ = q₁ + q₂ + q₃` and swing
tibia `φ₄ = φ₃ + q₄` hang downwards. `q = 0` is the robot standing straight with the
two legs side by side.
"""
function cartesian_coordinates(q::SVector{4}, geometry::Geometry)
    Lt, Lf = geometry.tibia, geometry.femur
    φ1 = q[1]
    φ2 = φ1 + q[2]
    φ3 = φ2 + q[3]
    φ4 = φ3 + q[4]
    grounded_foot = SVector(0.0, 0.0)
    grounded_knee = grounded_foot + Lt * SVector(sin(φ1), cos(φ1))
    hip = grounded_knee + Lf * SVector(sin(φ2), cos(φ2))
    swing_knee = hip - Lf * SVector(sin(φ3), cos(φ3))
    swing_foot = swing_knee - Lt * SVector(sin(φ4), cos(φ4))
    return (grounded_foot, grounded_knee, hip, swing_knee, swing_foot)
end

swing_foot_position(q::SVector{4}, geometry::Geometry) =
    cartesian_coordinates(q, geometry)[5]

# Per-coordinate gradient bounds of the swing-foot map on the state
# `[q̇₁, q₁, q₂, q₃, q₄]`: each joint rotates the rest of the chain about itself, so
# the foot moves by at most the pivot-to-foot distance per radian (ankle: both
# legs; grounded knee: femur + swing leg; hip: swing leg; swing knee: tibia). The
# foot does not depend on q̇₁.
function swing_foot_gradient_bound(geometry::Geometry)
    Lt, Lf = geometry.tibia, geometry.femur
    return SVector(0.0, 2 * (Lt + Lf), 2 * Lf + Lt, Lf + Lt, Lt)
end

"""
    swing_foot_deviation_bound(geometry, radius)

Upper bound on the displacement of the swing foot within a cell of half-widths
`radius` (in the 5-D state), from the swing foot of the cell center.
"""
swing_foot_deviation_bound(geometry::Geometry, radius::AbstractVector) =
    sum(swing_foot_gradient_bound(geometry) .* abs.(radius))

# ============================================================
# Velocity-controlled dynamics (from `BipedFourLink`)
# ============================================================

# Preallocated RigidBodyDynamics buffers for one thread.
struct _DynamicsBuffers{S, R, J}
    state::S
    result::R
    joints::J # (q₁, q₂, q₃, q₄) joints of the design mechanism
    i1::Int   # velocity index of q₁
end

function _buffers(model::BFL.BipedModel)
    mechanism = model.mechanism
    joints =
        (model.q1_joint, (RBD.findjoint(mechanism, name) for name in model.joint_names)...)
    state = RBD.MechanismState(mechanism)
    i1 = first(RBD.velocity_range(state, joints[1]))
    return _DynamicsBuffers(state, RBD.DynamicsResult(mechanism), joints, i1)
end

"""
    VelocityControlledDynamics(model)

Callable `(x, u) -> ẋ` for the velocity-controlled design state
`x = [q̇₁, q₁, q₂, q₃, q₄]` and commanded motor velocities `u = [v₂, v₃, v₄]`:

    q̈₁ = -c₁(q, [q̇₁, u]) / M₁₁(q),   q̇₁ = x₁,   (q̇₂, q̇₃, q̇₄) = u,

which is `BipedFourLink.velocity_controlled_derivative(model, x, u)` with its default
`input_accelerations = 0`. That function allocates a fresh mechanism state and
dynamics result on every call (≈ 170 µs); this one reuses per-thread buffers on the
same mechanism and returns an `SVector`, which is what makes an abstraction with
millions of RK4 roll-outs affordable. [`check_against_module`](@ref) compares the two.
"""
struct VelocityControlledDynamics{M, B}
    model::M
    buffers::Vector{B}
end

function VelocityControlledDynamics(model::BFL.BipedModel)
    # `maxthreadid` also covers interactive threads; the call never yields, so a
    # task cannot migrate between reading its thread id and releasing the buffer.
    return VelocityControlledDynamics(
        model,
        [_buffers(model) for _ in 1:Threads.maxthreadid()],
    )
end

function (f::VelocityControlledDynamics)(x, u)
    buf = f.buffers[Threads.threadid()]
    state, result, joints = buf.state, buf.result, buf.joints
    RBD.set_configuration!(state, joints[1], x[2])
    RBD.set_velocity!(state, joints[1], x[1])
    for i in 1:3
        RBD.set_configuration!(state, joints[i + 1], x[i + 2])
        RBD.set_velocity!(state, joints[i + 1], u[i])
    end
    RBD.mass_matrix!(result.massmatrix, state)
    RBD.dynamics_bias!(
        result.dynamicsbias,
        result.accelerations,
        result.jointwrenches,
        state,
    )
    i1 = buf.i1
    q̈1 = -parent(result.dynamicsbias)[i1] / result.massmatrix[i1, i1]
    return SVector(q̈1, x[1], u[1], u[2], u[3])
end

"""
    check_against_module(model; geometry, n = 50, atol = 1e-9)

Compare, at `n` random states and inputs, [`VelocityControlledDynamics`](@ref) with
`BipedFourLink.velocity_controlled_derivative` and [`cartesian_coordinates`](@ref)
with the `RigidBodyDynamics` joint positions. Errors on any mismatch.
"""
function check_against_module(
    model::BFL.BipedModel;
    geometry::Geometry = default_geometry(),
    n::Int = 50,
    atol = 1e-9,
)
    f = VelocityControlledDynamics(model)
    for k in 1:n
        # Deterministic, well-spread samples without a dependency on `Random`.
        s(j) = sin(12.9898 * k + 78.233 * j)
        x = SVector(ntuple(j -> s(j), 5))
        u = SVector(ntuple(j -> 2 * s(j + 5), 3))
        ẋ_ref = BFL.velocity_controlled_derivative(model, collect(x), collect(u))
        isapprox(f(x, u), ẋ_ref; atol = atol) ||
            error("fast dynamics disagree with BipedFourLink at x = $x, u = $u")

        q = joint_angles(x)
        rbd_state = BFL.state_from_coordinates(model, collect(q))
        _, knee, hip, swing_knee, swing_foot = cartesian_coordinates(q, geometry)
        for (p, ref) in (
            (knee, BFL.world_position(model, rbd_state, "grounded_knee")),
            (hip, BFL.world_position(model, rbd_state, "hip")),
            (swing_knee, BFL.world_position(model, rbd_state, "swing_knee")),
            (swing_foot, BFL.foot_position(model, rbd_state; foot = :swing)),
        )
            isapprox(p, SVector(ref[1], ref[3]); atol = atol) ||
                error("forward kinematics disagree with the URDF at q = $q")
        end
    end
    return true
end

# ============================================================
# Domain carving: Cartesian obstacle + ground → state-space cells
# ============================================================

_obstacle_list(::Nothing) = LazySets.LazySet[]
_obstacle_list(obstacle::LazySets.LazySet) = [obstacle]
_obstacle_list(obstacles::AbstractVector) = obstacles

"""
    infeasible_cells(X, state_grid, obstacle, geometry)

The `MP.CellUnion` of the cells of `X` (a 5-D hyperrectangle) to exclude, as in the
fully actuated model: cells where the swing foot *might* touch the Cartesian
`obstacle` (sound, `MP.image_blocked_cells`), and cells where the swing foot is
surely below the ground (tolerant, so the footholds on the ground stay).
"""
function infeasible_cells(
    X::LazySets.AbstractHyperrectangle,
    state_grid::MP.GridFree,
    obstacle,
    geometry::Geometry,
)
    grad = swing_foot_gradient_bound(geometry)
    dev = swing_foot_deviation_bound(geometry, MP.get_h(state_grid) ./ 2)

    # The swing foot does not depend on q̇₁: decide on the 4-D angle cells once.
    q_grid = MP.GridFree(MP.get_origin(state_grid)[IDX_Q], MP.get_h(state_grid)[IDX_Q])
    Q = LazySets.Hyperrectangle(;
        low = LazySets.low(X)[IDX_Q],
        high = LazySets.high(X)[IDX_Q],
    )
    q_foot(q) = swing_foot_position(SVector{4}(q), geometry)
    below_ground = MP.cells_where(q_grid, Q) do pos
        return q_foot(MP.get_coord_by_pos(q_grid, pos))[2] < -dev
    end
    obstacles = _obstacle_list(obstacle)
    blocked_q =
        isempty(obstacles) ? Set{NTuple{4, Int}}() :
        MP.image_blocked_cells(q_foot, grad[IDX_Q], q_grid, Q, obstacles).positions
    removed_q = union(below_ground.positions, blocked_q)

    # Lift to 5-D along the q̇₁ axis.
    h = MP.get_h(state_grid)
    k_lo = round(
        Int,
        (LazySets.low(X)[IDX_Q̇1] - MP.get_origin(state_grid)[IDX_Q̇1]) / h[IDX_Q̇1],
    )
    k_hi = round(
        Int,
        (LazySets.high(X)[IDX_Q̇1] - MP.get_origin(state_grid)[IDX_Q̇1]) / h[IDX_Q̇1],
    )
    removed = Set{NTuple{5, Int}}()
    sizehint!(removed, length(removed_q) * (k_hi - k_lo + 1))
    for pos in removed_q, k in k_lo:k_hi
        push!(removed, (k, pos...))
    end
    return MP.CellUnion(state_grid, removed)
end

# ============================================================
# Target set: configurations placing the swing foot at a foothold
# ============================================================

# Two-link inverse kinematics of the swing leg: absolute angles (φ₃, φ₄) placing
# the swing foot at `hip + d`, with the knee bent backwards (`q₄ = φ₄ − φ₃ ≥ 0`, the
# human-like knee for walking towards `+x`), or `nothing` if unreachable.
function _swing_leg_ik(d::SVector{2}, geometry::Geometry)
    L1, L2 = geometry.femur, geometry.tibia
    c = (d[1]^2 + d[2]^2 - L1^2 - L2^2) / (2 * L1 * L2)
    -1.0 <= c <= 1.0 || return nothing
    # Same planar two-link problem as the fully actuated model with ψ = −φ.
    ψ4_rel = -acos(c)
    A = L1 + L2 * cos(ψ4_rel)
    B = L2 * sin(ψ4_rel)
    ψ3 = atan(d[1], -d[2]) - atan(B, A)
    φ3 = -ψ3
    φ4 = -(ψ3 + ψ4_rel)
    return φ3, φ4
end

"""
    target_cells(state_grid, foothold, geometry, X; nsub = 3)

5-D grid cells containing a configuration with the grounded foot at the origin
and the swing foot at the Cartesian point `foothold`. `(q₁, q₂)` are swept over
their cells (`nsub` samples per axis per cell), `(q₃, q₄)` follow by inverse
kinematics, and every q̇₁ cell of `X` is accepted.
"""
function target_cells(
    state_grid::MP.GridFree,
    foothold::SVector{2},
    geometry::Geometry,
    X::LazySets.AbstractHyperrectangle;
    nsub::Int = 3,
)
    h = MP.get_h(state_grid)
    lo = SVector{5}(LazySets.low(X))
    hi = SVector{5}(LazySets.high(X))
    q̇1_values = range(
        MP.get_coord_by_pos(state_grid, MP.get_pos_by_coord(state_grid, lo))[1],
        hi[1];
        step = h[1],
    )

    found = Set{NTuple{5, Int}}()
    for q1 in range(lo[2], hi[2]; step = h[2] / nsub)
        for q2 in range(lo[3], hi[3]; step = h[3] / nsub)
            hip = cartesian_coordinates(SVector(q1, q2, 0.0, 0.0), geometry)[3]
            ik = _swing_leg_ik(foothold - hip, geometry)
            ik === nothing && continue
            φ3, φ4 = ik
            q3 = φ3 - (q1 + q2)
            q4 = φ4 - φ3
            for q̇1 in q̇1_values
                x = SVector(q̇1, q1, q2, q3, q4)
                all(lo .<= x .<= hi) || continue
                push!(found, MP.get_pos_by_coord(state_grid, x))
            end
        end
    end
    return collect(found)
end

"""
    target_set(state_grid, foothold, geometry, X; kwargs...)

The `MP.CellUnion` of [`target_cells`](@ref).
"""
function target_set(
    state_grid::MP.GridFree,
    foothold::SVector{2},
    geometry::Geometry,
    X::LazySets.AbstractHyperrectangle;
    kwargs...,
)
    cells = target_cells(state_grid, foothold, geometry, X; kwargs...)
    isempty(cells) && error("no reachable target configuration for foothold $foothold")
    return MP.CellUnion(state_grid, cells)
end

# ============================================================
# Discretization
# ============================================================

"""
    default_discretization(; dx = 0.1, dq1 = 0.05, dq̇1 = 0.1, tstep = 0.1)

5-D state grid (steps `(dq̇1, dq1, dx, dx, dx)`), 3-D input grid of step
`du = dx / tstep` and the sampling time. `tstep · du = dx` keeps the three motor
coordinates on the lattice (each input moves them by exactly one cell); the ankle
coordinates `(q̇₁, q₁)` follow the nonlinear dynamics and are rounded to the grid.

Returns `(; state_grid, input_grid, tstep, u_max, du)` with `u_max = du`.
"""
function default_discretization(; dx = 0.1, dq1 = 0.05, dq̇1 = 0.1, tstep = 0.1)
    du = dx / tstep
    state_grid =
        MP.GridFree(SVector(0.0, 0.0, 0.0, 0.0, 0.0), SVector(dq̇1, dq1, dx, dx, dx))
    input_grid = MP.GridFree(SVector(0.0, 0.0, 0.0), SVector(du, du, du))
    return (; state_grid, input_grid, tstep, u_max = du, du)
end

# ============================================================
# System / problem factories
# ============================================================

Base.@kwdef struct RobotDomainConfig{X1, X2, U1, U2}
    x_lb::X1
    x_ub::X2
    u_lb::U1
    u_ub::U2
end

"""
    system(; model, domain, removed_cells = nothing)

Continuous-time velocity-controlled system `ẋ = f(x, u)` with
`f = VelocityControlledDynamics(model)`, state domain the box of `domain` minus the
`removed_cells` (a precomputed [`infeasible_cells`](@ref)). The abstraction
discretizes it with its own `time_step` (RK4).
"""
function system(;
    model::BFL.BipedModel = build_model(),
    domain::RobotDomainConfig,
    removed_cells = nothing,
)
    X = LazySets.Hyperrectangle(; low = domain.x_lb, high = domain.x_ub)
    if removed_cells !== nothing && !isempty(removed_cells)
        X = UT.set_minus(X, removed_cells)
    end
    U = LazySets.Hyperrectangle(; low = domain.u_lb, high = domain.u_ub)
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        VelocityControlledDynamics(model),
        5,
        3,
        X,
        U,
    )
end

"""
    balanced_start(f; q2 = 0.0, φ3 = 0.2, q4 = 0.0, q1_range = (-0.5, 0.5))

Resting configuration `x = [0, q₁, q₂, q₃, q₄]` whose ankle is in static equilibrium
(`q̈₁ = 0` at zero velocities) for the dynamics `f` (a
[`VelocityControlledDynamics`](@ref)): the grounded knee at `q₂`, the swing femur at
absolute angle `φ₃` (positive = behind the hip) and the swing knee at `q₄`; `q₁` is
found by bisection on `q1_range` and `q₃ = φ₃ − q₁ − q₂`.

The equilibrium is unstable (the ankle is an inverted pendulum), but it is the one
start from which the ankle does not already fall at the first sample.
"""
function balanced_start(
    f::VelocityControlledDynamics;
    q2 = 0.0,
    φ3 = 0.2,
    q4 = 0.0,
    q1_range = (-0.5, 0.5),
)
    config(q1) = SVector(0.0, q1, q2, φ3 - q1 - q2, q4)
    q̈1(q1) = f(config(q1), SVector(0.0, 0.0, 0.0))[1]
    a, b = q1_range
    sign(q̈1(a)) != sign(q̈1(b)) ||
        error("no ankle equilibrium in q1 ∈ $q1_range for this posture")
    for _ in 1:60
        m = (a + b) / 2
        sign(q̈1(m)) == sign(q̈1(a)) ? (a = m) : (b = m)
    end
    return config((a + b) / 2)
end

"""
    step_problem(concrete_system, state_grid; x0, foothold, geometry)

One footstep as a reach specification: from `x0` (a small box strictly inside its
grid cell), reach a configuration placing the swing foot at `foothold`.
"""
function step_problem(
    concrete_system,
    state_grid;
    x0::SVector{5},
    foothold::SVector{2} = SVector(0.2, 0.0),
    geometry::Geometry = default_geometry(),
)
    h = MP.get_h(state_grid)
    _I_ = LazySets.Hyperrectangle(x0, 0.25 .* h)
    X_box = UT._outer_box(concrete_system.X)
    _T_ = target_set(state_grid, foothold, geometry, X_box)
    return PR.OptimalControlProblem(
        concrete_system,
        _I_,
        _T_,
        nothing,
        nothing,
        PR.Infinity(),
    )
end

# ============================================================
# Visualization
# ============================================================

"""
    draw_robot!(fig, x; geometry, alpha)

Draw the grounded (blue) and swing (red) legs of the state `x = [q̇₁, q₁, q₂, q₃, q₄]`.
"""
function draw_robot!(fig, x; geometry::Geometry = default_geometry(), alpha::Real = 1.0)
    gf, gk, hip, sk, sf = cartesian_coordinates(joint_angles(x), geometry)
    Plots.plot!(
        fig,
        [gf[1], gk[1], hip[1]],
        [gf[2], gk[2], hip[2]];
        lw = 4,
        marker = :circle,
        color = :steelblue,
        alpha = alpha,
        label = "",
    )
    Plots.plot!(
        fig,
        [hip[1], sk[1], sf[1]],
        [hip[2], sk[2], sf[2]];
        lw = 4,
        marker = :circle,
        color = :indianred,
        alpha = alpha,
        label = "",
    )
    return fig
end

"""
    system_plot!(; geometry, obstacle, foothold, xlims, ylims)

Closure `(fig, x, u) -> fig` drawing the robot, the ground, the obstacle and the
target foothold, for the trajectory dashboards.
"""
function system_plot!(;
    geometry::Geometry = default_geometry(),
    obstacle = nothing,
    foothold = nothing,
    xlims = (-0.5, 0.5),
    ylims = (-0.05, 0.45),
)
    return function (fig, x, u)
        _ = u
        Plots.plot!(
            fig,
            [xlims[1], xlims[2]],
            [0.0, 0.0];
            color = :black,
            lw = 1,
            label = "",
        )
        for ob in _obstacle_list(obstacle)
            Plots.plot!(fig, ob; color = :black, alpha = 0.8, label = "")
        end
        foothold !== nothing && Plots.scatter!(
            fig,
            [foothold[1]],
            [foothold[2]];
            marker = :xcross,
            color = :green,
            label = "",
        )
        draw_robot!(fig, x; geometry = geometry)
        Plots.plot!(fig; xlims = xlims, ylims = ylims, aspect_ratio = :equal)
        return fig
    end
end

end # module
