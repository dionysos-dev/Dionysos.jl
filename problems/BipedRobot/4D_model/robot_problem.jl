# 4-D velocity-controlled biped model.
#
# Instead of controlling the full multibody robot in torque (6D/8D models, one
# RigidBodyDynamics roll-out per transition), the robot is controlled at the
# *velocity* level: states are the 4 joint angles (stance/swing hip and knee),
# inputs are the joint angular velocities, and the dynamics reduce to the pure
# integrator `ẋ = u` — the motors' low-level controllers are trusted to track
# the commanded velocities. This halves the abstraction dimension and makes
# `f(x, u)` trivial.
#
# The discretization is chosen *commensurable*: with `tstep · du = hx`, every
# input translates the state grid exactly onto itself, so the abstraction is
# deterministic and exact (a bisimulation — zero spurious transitions), which
# is what makes the `CENTER_SIMULATION` approximation mode sound here (checked
# by `MP.is_lattice_exact`; see `default_discretization`).
#
# Physical constraints (a Cartesian obstacle the swing foot must clear, the
# ground) are pulled back into joint space by removing grid cells whose image
# might touch them (`infeasible_cells`); the removal is a sound OUTER carving
# via a Lipschitz bound of the forward kinematics, closing the "obstacle jump"
# gap together with the one-cell-per-step displacement bound
# (`MP.intersample_safe_time_step`).
module RobotProblem

using StaticArrays
using MathematicalSystems
import LazySets
import Plots

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const PR = DI.Problem
const MP = DI.Mapping

# ============================================================
# Geometry & kinematics
# ============================================================

"""
    Geometry{T}

Planar two-leg geometry: torso/hip mounting angle `theta_hip` and the two
segment lengths. Defaults follow the canonical robot constants used by the
6D/8D models (`Lthigh = 0.20125`, `Lleg = 0.172`).
"""
Base.@kwdef struct Geometry{T <: Real}
    theta_hip::T = 0.0
    hip_to_knee::T = 0.20125
    knee_to_foot::T = 0.172
end

default_geometry() = Geometry{Float64}()

"""
    cartesian_coordinates(θ, geometry, grounded_left_foot)

Joint positions from the 4 joint angles `θ = (left hip, left knee, right hip,
right knee)`, with the grounded foot at the origin. Returns the tuple
`(left_hip, left_knee, left_foot, right_hip, right_knee, right_foot)`.
"""
function cartesian_coordinates(θ::SVector{4}, geometry::Geometry, grounded_left_foot::Bool)
    (θ1, θ2, θ3, θ4) =
        grounded_left_foot ? (θ[1], θ[2], θ[3], θ[4]) : (θ[3], θ[4], θ[1], θ[2])
    α = geometry.theta_hip

    grounded_foot = SVector(0.0, 0.0)
    grounded_knee =
        grounded_foot + geometry.knee_to_foot * SVector(-sin(α + θ1 + θ2), cos(α + θ1 + θ2))
    grounded_hip = grounded_knee + geometry.hip_to_knee * SVector(-sin(α + θ1), cos(α + θ1))
    swing_hip = grounded_hip
    swing_knee = swing_hip + geometry.hip_to_knee * SVector(sin(α + θ3), -cos(α + θ3))
    swing_foot =
        swing_knee + geometry.knee_to_foot * SVector(sin(α + θ3 + θ4), -cos(α + θ3 + θ4))

    return grounded_left_foot ?
           (grounded_hip, grounded_knee, grounded_foot, swing_hip, swing_knee, swing_foot) :
           (swing_hip, swing_knee, swing_foot, grounded_hip, grounded_knee, grounded_foot)
end

"""
    angular_coordinates(cartesian, geometry)

Inverse of [`cartesian_coordinates`](@ref) on the joint positions tuple
`(left_hip, left_knee, left_foot, right_hip, right_knee, right_foot)`. Returns
`(θ, grounded_left_foot)`. Valid on the working domain where every segment
angle stays within `(-π/2, π/2)` of the vertical.
"""
function angular_coordinates(cartesian::Tuple, geometry::Geometry)
    lh, lk, lf, rh, rk, _ = cartesian
    grounded_left_foot = lf == SVector(0.0, 0.0)
    α = geometry.theta_hip

    θ1 = asin((lk[1] - lh[1]) / geometry.hip_to_knee) - α
    θ2 = asin((lf[1] - lk[1]) / geometry.knee_to_foot) - α - θ1
    θ3 = asin((rk[1] - rh[1]) / geometry.hip_to_knee) - α
    θ4 = asin((cartesian[6][1] - rk[1]) / geometry.knee_to_foot) - α - θ3

    return SVector(θ1, θ2, θ3, θ4), grounded_left_foot
end

swing_foot_position(θ::SVector{4}, geometry::Geometry, grounded_left_foot::Bool) =
    cartesian_coordinates(θ, geometry, grounded_left_foot)[grounded_left_foot ? 6 : 3]

# Per-joint gradient bounds of the swing-foot map: the foot position moves by at
# most (L1 + L2) per radian of a hip angle and L2 per radian of a knee angle, so
# |foot(x) − foot(c)| ≤ Σᵢ gᵢ |xᵢ − cᵢ| on any set.
function swing_foot_gradient_bound(geometry::Geometry)
    g_hip = geometry.hip_to_knee + geometry.knee_to_foot
    g_knee = geometry.knee_to_foot
    return SVector(g_hip, g_knee, g_hip, g_knee)
end

"""
    swing_foot_deviation_bound(geometry, radius)

Upper bound on how far the swing foot of any configuration inside a cell of
half-widths `radius` can be from the swing foot of the cell center (Lipschitz
bound of the forward kinematics). This is what makes the center-based cell
tests of [`infeasible_cells`](@ref) sound.
"""
swing_foot_deviation_bound(geometry::Geometry, radius::AbstractVector) =
    sum(swing_foot_gradient_bound(geometry) .* abs.(radius))

# ============================================================
# Domain carving: Cartesian obstacle + ground → joint-space cells
# ============================================================

# Boxes passed to `UT.set_minus` are shrunk by a relative ε around their own
# center: a full cell box shares its faces with the neighboring cells, and the
# (closed-boundary) OUTER hole discretization would then also swallow every
# face-adjacent neighbor — a 3ⁿ over-removal. The shrink only excludes the
# intended cell. (This replaces an origin-anchored `0.95 * bounds` shrink in
# the original code, which shifted the removed boxes away from their cells.)
function _shrunk_cell(grid, pos; ε = 1e-8)
    cell = MP.get_elem_by_pos(grid, pos)
    return LazySets.Hyperrectangle(LazySets.center(cell), (1 - ε) .* cell.radius)
end

function _inflated_cell(grid, pos; ε = 1e-8)
    cell = MP.get_elem_by_pos(grid, pos)
    return LazySets.Hyperrectangle(LazySets.center(cell), (1 + ε) .* cell.radius)
end

# `nothing`, a single set, or a vector of sets → vector of sets.
_obstacle_list(::Nothing) = LazySets.LazySet[]
_obstacle_list(obstacle::LazySets.LazySet) = [obstacle]
_obstacle_list(obstacles::AbstractVector) = obstacles

"""
    infeasible_cells(X, state_grid, obstacle, geometry; grounded_left_foot = true)

Positions of the grid cells of `X` (a hyperrectangle) that a sound abstraction
must exclude:

- **obstacle** (sound, OUTER): a cell is removed as soon as *some* configuration
  in it *might* place the swing foot inside the Cartesian `obstacle` — tested
  as `foot(center)` within the Lipschitz deviation bound of the obstacle, so a
  kept cell is certified obstacle-free. Accepts a single `LazySet`, a vector of
  them (e.g. a step plus a ceiling), or `nothing` to skip.
- **ground** (tolerant): a cell is removed only when *every* configuration in
  it surely puts the swing foot strictly below the ground. The ground `y ≥ 0`
  is a closed contact constraint — the target foothold lies *on* it — so sound
  OUTER carving would remove the target cells themselves; boundary cells are
  deliberately kept, accurate to the grid resolution.

Returns a `Vector` of grid positions (integer tuples).
"""
function infeasible_cells(
    X::LazySets.AbstractHyperrectangle,
    state_grid::MP.GridFree,
    obstacle,
    geometry::Geometry;
    grounded_left_foot::Bool = true,
)
    dev = swing_foot_deviation_bound(geometry, MP.get_h(state_grid) ./ 2)
    obstacles = _obstacle_list(obstacle)
    # Cheap prefilter: the swing foot must be within `dev` of an obstacle's
    # bounding box for the exact disjointness test to be worth running.
    obstacle_boxes = [UT._outer_box(ob) for ob in obstacles]

    removed = NTuple{4, Int}[]
    for pos in MP.get_pos_from_set(state_grid, X, MP.INNER)
        center = MP.get_coord_by_pos(state_grid, pos)
        foot = swing_foot_position(SVector{4}(center), geometry, grounded_left_foot)

        if foot[2] < -dev # surely below ground everywhere in the cell
            push!(removed, pos)
            continue
        end

        for (ob, box) in zip(obstacles, obstacle_boxes)
            near_box = all(
                abs(foot[i] - LazySets.center(box)[i]) <= dev + box.radius[i] for i in 1:2
            )
            # Box over-approximation of the deviation ball: conservative, and
            # box vs set disjointness is robust for any LazySet obstacle.
            if near_box &&
               !UT.is_disjoint(LazySets.Hyperrectangle(foot, SVector(dev, dev)), ob)
                push!(removed, pos)
                break
            end
        end
    end
    return removed
end

"""
    remove_cells(X, state_grid, positions)

State domain with the cells at `positions` removed: `X ∖ ⋃ cells`. The removed
boxes are slightly shrunk so that only the intended cells drop out of the grid
mapping.
"""
function remove_cells(
    X::LazySets.AbstractHyperrectangle,
    state_grid::MP.GridFree,
    positions,
)
    isempty(positions) && return X
    holes = [_shrunk_cell(state_grid, pos) for pos in positions]
    return UT.set_minus(X, UT.set_union(holes))
end

"""
    carve_domain(X, state_grid, obstacle, geometry; grounded_left_foot = true)

State domain with the infeasible cells removed:
`X ∖ ⋃ cells(infeasible_cells(...))` (see [`remove_cells`](@ref)).
"""
function carve_domain(
    X::LazySets.AbstractHyperrectangle,
    state_grid::MP.GridFree,
    obstacle,
    geometry::Geometry;
    grounded_left_foot::Bool = true,
)
    removed = infeasible_cells(
        X,
        state_grid,
        obstacle,
        geometry;
        grounded_left_foot = grounded_left_foot,
    )
    return remove_cells(X, state_grid, removed)
end

# ============================================================
# Target set: configurations placing the swing foot at a foothold
# ============================================================

# Two-link inverse kinematics of the swing leg: joint angles placing the swing
# foot at `hip + d`, elbow-down (knee angle ≤ 0), or `nothing` if unreachable.
function _swing_leg_ik(d::SVector{2}, geometry::Geometry)
    L1 = geometry.hip_to_knee
    L2 = geometry.knee_to_foot
    c4 = (d[1]^2 + d[2]^2 - L1^2 - L2^2) / (2 * L1 * L2)
    -1.0 <= c4 <= 1.0 || return nothing
    θ4 = -acos(c4)
    A = L1 + L2 * cos(θ4)
    B = L2 * sin(θ4)
    θ3 = atan(d[1], -d[2]) - atan(B, A) - geometry.theta_hip
    return θ3, θ4
end

"""
    target_cells(state_grid, foothold, geometry, X; grounded_left_foot = true, nsub = 3)

Grid cells containing a configuration that keeps the grounded foot at the
origin and places the swing foot at the Cartesian point `foothold` — the target
manifold of one step. Stance angles `(θ1, θ2)` are swept over their grid cells
(`nsub` samples per axis per cell), the swing angles follow by two-link inverse
kinematics, and the visited cells are deduplicated.

Returns a `Vector` of grid positions.
"""
function target_cells(
    state_grid::MP.GridFree,
    foothold::SVector{2},
    geometry::Geometry,
    X::LazySets.AbstractHyperrectangle;
    grounded_left_foot::Bool = true,
    nsub::Int = 3,
)
    h = MP.get_h(state_grid)
    lo = SVector{4}(LazySets.low(X))
    hi = SVector{4}(LazySets.high(X))

    stance_lo = grounded_left_foot ? (lo[1], lo[2]) : (lo[3], lo[4])
    stance_hi = grounded_left_foot ? (hi[1], hi[2]) : (hi[3], hi[4])
    h12 = grounded_left_foot ? (h[1], h[2]) : (h[3], h[4])

    found = Set{NTuple{4, Int}}()
    for θ1 in range(stance_lo[1], stance_hi[1]; step = h12[1] / nsub)
        for θ2 in range(stance_lo[2], stance_hi[2]; step = h12[2] / nsub)
            stance =
                grounded_left_foot ? SVector(θ1, θ2, 0.0, 0.0) : SVector(0.0, 0.0, θ1, θ2)
            hip =
                cartesian_coordinates(stance, geometry, grounded_left_foot)[grounded_left_foot ?
                                                                            1 : 4]
            ik = _swing_leg_ik(foothold - hip, geometry)
            ik === nothing && continue
            θ3, θ4 = ik
            θ = grounded_left_foot ? SVector(θ1, θ2, θ3, θ4) : SVector(θ3, θ4, θ1, θ2)
            all(lo .<= θ .<= hi) || continue
            push!(found, MP.get_pos_by_coord(state_grid, θ))
        end
    end
    return collect(found)
end

"""
    target_set(state_grid, foothold, geometry, X; kwargs...)

Union of the (slightly inflated, so that INNER discretization robustly recovers
them) cells of [`target_cells`](@ref).
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
    return UT.set_union([_inflated_cell(state_grid, pos) for pos in cells])
end

# ============================================================
# Discretization: commensurable grids (the exact-lattice design)
# ============================================================

"""
    default_discretization(; dx = 0.1, tstep = 0.1, speed_levels = 1, swept_transitions = false)

State grid, input grid and sampling time such that `tstep * du = dx`: every
admissible input translates the state grid exactly onto itself, so the
abstraction is deterministic and exact (checked with `MP.is_lattice_exact`).

With `speed_levels = 1` each step moves the state by at most one cell per
axis, which closes the inter-sample soundness gap over carved domains by
itself (`MP.intersample_safe_time_step`). Higher `speed_levels` (a richer
velocity alphabet, up to `speed_levels · du` per joint) move up to
`speed_levels` cells per axis per step and therefore **require** swept-cell
transition validation: pass `swept_transitions = true` to acknowledge it and
wire [`swept_state_input_filter`](@ref) into the abstraction.

Returns `(; state_grid, input_grid, tstep, u_max)`.
"""
function default_discretization(;
    dx = 0.1,
    tstep = 0.1,
    speed_levels::Int = 1,
    swept_transitions::Bool = false,
)
    du = dx / tstep
    state_grid = MP.GridFree(SVector(0.0, 0.0, 0.0, 0.0), SVector(dx, dx, dx, dx))
    input_grid = MP.GridFree(SVector(0.0, 0.0, 0.0, 0.0), SVector(du, du, du, du))
    u_max = speed_levels * du

    MP.is_lattice_exact(state_grid, input_grid, tstep) || error(
        "state grid, input grid and tstep are not commensurable: " *
        "CENTER_SIMULATION would be unsound. Use tstep * du = dx.",
    )
    if !swept_transitions
        tstep <= MP.intersample_safe_time_step(state_grid, SVector(fill(u_max, 4)...)) ||
            error(
                "tstep exceeds one cell of displacement per axis per step: " *
                "inter-sample obstacle crossings become possible. Either use " *
                "speed_levels = 1, or validate multi-cell steps with " *
                "`swept_state_input_filter` and pass `swept_transitions = true`.",
            )
    end

    return (; state_grid, input_grid, tstep, u_max, du)
end

"""
    swept_state_input_filter(state_grid, tstep, removed)

`(x, u) -> Bool` filter for the abstraction (`state_input_filter` attribute):
keeps a transition only when **every** grid cell crossed by the inter-sample
segment from `x` to `x + tstep · u` avoids the `removed` cells (an iterable of
grid positions, e.g. from [`infeasible_cells`](@ref)).

For the integrator dynamics `ẋ = u` the inter-sample trajectory *is* that
segment, so this check restores inter-sample soundness for multi-cell steps
(`speed_levels ≥ 2`), lifting the one-cell-per-step speed cap. Set the
`intersample_checked` attribute of the abstraction solver to `true` alongside,
so the multi-cell jump warning knows the sweep is covered.
"""
function swept_state_input_filter(state_grid::MP.GridFree, tstep::Real, removed)
    removed_set = removed isa Set ? removed : Set(removed)
    return function (x, u)
        for pos in MP.cells_on_segment(state_grid, x, x .+ tstep .* u)
            pos in removed_set && return false
        end
        return true
    end
end

# ============================================================
# System / problem factories (6D-compatible surface)
# ============================================================

Base.@kwdef struct RobotDomainConfig{X1, X2, U1, U2}
    x_lb::X1
    x_ub::X2
    u_lb::U1
    u_ub::U2
end

"""
    default_robot_domain(; u_max = 1.0)

Joint angles within ±π/3, joint velocities within ±`u_max` (which must match
the input grid: `u_max = speed_levels * du`).
"""
function default_robot_domain(; u_max = 1.0)
    x_bar = π / 3
    return RobotDomainConfig(;
        x_lb = SVector(-x_bar, -x_bar, -x_bar, -x_bar),
        x_ub = SVector(x_bar, x_bar, x_bar, x_bar),
        u_lb = SVector(-u_max, -u_max, -u_max, -u_max),
        u_ub = SVector(u_max, u_max, u_max, u_max),
    )
end

"""
    system(; tstep, domain, obstacle, state_grid, geometry, grounded_left_foot, kwargs...)

Discrete-time velocity-controlled system `x⁺ = x + tstep · u` (exact for
`ẋ = u`). When `obstacle` and `state_grid` are given, the state domain is
carved (see [`carve_domain`](@ref)); the extra keyword arguments (`robot_urdf`,
`Δt_simu`, `simulator`) are accepted for drop-in compatibility with the 6D/8D
execution drivers and ignored — the model is purely kinematic.
"""
function system(;
    robot_urdf = nothing,
    tstep = 0.1,
    domain = default_robot_domain(),
    Δt_simu = nothing,
    simulator::Symbol = :none,
    obstacle = nothing,
    state_grid = nothing,
    geometry::Geometry = default_geometry(),
    grounded_left_foot::Bool = true,
    removed_cells = nothing, # precomputed `infeasible_cells` result, to carve once
)
    _ = (robot_urdf, Δt_simu, simulator)

    X = LazySets.Hyperrectangle(; low = domain.x_lb, high = domain.x_ub)
    if removed_cells !== nothing
        state_grid === nothing && error("`removed_cells` requires `state_grid`")
        X = remove_cells(X, state_grid, removed_cells)
    elseif obstacle !== nothing || state_grid !== nothing
        state_grid === nothing &&
            error("carving the obstacle out of the domain requires `state_grid`")
        X = carve_domain(
            X,
            state_grid,
            obstacle,
            geometry;
            grounded_left_foot = grounded_left_foot,
        )
    end
    U = LazySets.Hyperrectangle(; low = domain.u_lb, high = domain.u_ub)

    step(x, u) = x + tstep * u
    return MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem(step, 4, 4, X, U)
end

"""
    problem(; kwargs...)

Abstraction-time specification ([`Dionysos.Problem.AlternatingSimulationProblem`](@ref))
over the velocity-controlled system; same keyword surface as [`system`](@ref).
"""
function problem(; kwargs...)
    sys = system(; kwargs...)
    return PR.AlternatingSimulationProblem(sys, nothing)
end

"""
    warmup_robot_problem!(; kwargs...)

One dynamics call to force compilation (kept for drop-in compatibility with the
execution drivers; the kinematic model has nothing expensive to warm up).
"""
function warmup_robot_problem!(; kwargs...)
    sys = system(; kwargs...)
    sys.f(SVector(0.0, 0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0, 0.0))
    return nothing
end

"""
    step_problem(concrete_system, state_grid; x0, foothold, geometry, grounded_left_foot)

One footstep as a reach specification: from the configuration `x0` (a small box
strictly inside its grid cell, so the abstract initial set is that single
cell), reach a configuration placing the swing foot at `foothold` (see
[`target_set`](@ref)).
"""
function step_problem(
    concrete_system,
    state_grid;
    x0::SVector{4} = SVector(0.2, 0.0, -0.2, 0.0),
    foothold::SVector{2} = SVector(0.2, 0.0),
    geometry::Geometry = default_geometry(),
    grounded_left_foot::Bool = true,
)
    h = MP.get_h(state_grid)
    _I_ = LazySets.Hyperrectangle(x0, 0.25 .* h)
    X_box = UT._outer_box(concrete_system.X)
    _T_ = target_set(
        state_grid,
        foothold,
        geometry,
        X_box;
        grounded_left_foot = grounded_left_foot,
    )
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

function robot_segments(θ::SVector{4}, geometry::Geometry, grounded_left_foot::Bool)
    lh, lk, lf, rh, rk, rf = cartesian_coordinates(θ, geometry, grounded_left_foot)
    return ([lf[1], lk[1], lh[1]], [lf[2], lk[2], lh[2]]),
    ([rf[1], rk[1], rh[1]], [rf[2], rk[2], rh[2]])
end

"""
    draw_robot!(fig, θ; geometry, grounded_left_foot)

Draw the two legs (and the hip joint) of configuration `θ` on `fig`.
"""
function draw_robot!(
    fig,
    θ::SVector{4};
    geometry::Geometry = default_geometry(),
    grounded_left_foot::Bool = true,
)
    (lx, ly), (rx, ry) = robot_segments(θ, geometry, grounded_left_foot)
    Plots.plot!(fig, lx, ly; lw = 4, marker = :circle, color = :steelblue, label = "")
    Plots.plot!(fig, rx, ry; lw = 4, marker = :circle, color = :indianred, label = "")
    Plots.plot!(
        fig,
        [lx[end], rx[end]],
        [ly[end], ry[end]];
        lw = 4,
        color = :gray,
        label = "",
    )
    return fig
end

"""
    system_plot!(; geometry, grounded_left_foot, obstacle, foothold, xlims, ylims)

Closure `(fig, x, u) -> fig` drawing the robot in the Cartesian plane, with the
ground, the obstacle (black, house convention) and the target foothold; for the
trajectory dashboards.
"""
function system_plot!(;
    geometry::Geometry = default_geometry(),
    grounded_left_foot::Bool = true,
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
        draw_robot!(
            fig,
            SVector{4}(x);
            geometry = geometry,
            grounded_left_foot = grounded_left_foot,
        )
        Plots.plot!(fig; xlims = xlims, ylims = ylims, aspect_ratio = :equal)
        return fig
    end
end

end # module
