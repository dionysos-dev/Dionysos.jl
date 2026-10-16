# Underactuated four-link biped: mechanism, normal-form design model and footstep problem.
#
# The robot is defined **once**, symbolically, with `Multibody2D` (`symbolic_model`): everything
# else in this file is derived from that definition (the forward kinematics, the normal-form design
# model) or is the footstep scenario built on top of it (obstacle carving, target set, grids,
# start state, closed loop, plots). It is the same robot as `four_link_biped.urdf`, which
# `BipedFourLink.jl` simulates with RigidBodyDynamics (checked in
# `tests/test_biped_four_link.jl`); no RigidBodyDynamics is needed here.
#
# Unlike the fully actuated 4-D model (`../4D_model/robot_problem.jl`, `ẋ = u`), the grounded ankle
# `q₁` is *not* actuated: only the three motor joints (grounded knee `q₂`, hip `q₃`, swing knee
# `q₄`) track commanded velocities `u = [v₂, v₃, v₄]`, held constant over a sampling period (the
# actuated accelerations are zero), and the ankle obeys the first, unactuated Euler–Lagrange row of
# the multibody model: `q₁` falls like an inverted pendulum, and the controller keeps it in the
# domain only through the configuration of the motor joints and the Coriolis coupling. The design
# model (`reduced_dynamics`) describes the ankle by its generalized momentum `p₁ = ∂L/∂q̇₁`, in the
# first slot of the 5-D state:
#
#     x = [p₁, q₁, q₂, q₃, q₄],     u = [v₂, v₃, v₄]
#
# Angle convention: all joint angles, velocities and the ankle coordinate are **positive
# counter-clockwise** in the `(x, z)` plane (`x` pointing forwards, `z` up), the convention of
# `Multibody2D.jl`. A positive `q₁` leans the robot backwards, a positive `q₃` swings the swing
# foot forwards, and the knee bends backwards for a negative `q₄`.
#
# The dynamics are nonlinear in the ankle coordinates, so the abstraction has to be a certified
# over-approximation: the local growth bound of `LocalGrowthAbstraction.jl` is. The closed-loop
# trajectory of the concrete system is the end-to-end check of a controller (`closed_loop`).
module UnderactuatedRobotProblem

using StaticArrays
import LazySets
import MathOptInterface as MOI
import Plots
import Symbolics
using Symbolics: Num

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping

# Only what is unambiguous is imported (`Multibody2D` exports `position`, `mass_matrix`, … which
# clash with other packages): the rest is qualified (`Multibody2D.position`, …).
include(joinpath(@__DIR__, "Multibody2D.jl"))
using .Multibody2D: MultibodySystem, add_frame!, add_point!, normal_form

# Index of each coordinate in the state `x = [ankle, q₁, q₂, q₃, q₄]`.
const IDX_ANKLE = 1
const IDX_Q = SVector(2, 3, 4, 5)

joint_angles(x) = SVector(x[2], x[3], x[4], x[5])

# ============================================================
# Mechanism: the symbolic model and its numbers
# ============================================================

"""
    symbolic_model() -> NamedTuple

The four-link biped as a [`MultibodySystem`](@ref Multibody2D.MultibodySystem) with symbolic
parameters, and its normal form.

Frames (all translations in the parent frame, all angles relative to the parent, counter-clockwise):

* `RF1`: rotated by `q₁` about the world, no translation (grounded tibia, ankle at the origin);
* `RF2`: `[0, L2]` from `RF1`, rotated by `q₂` (grounded knee);
* `RF3`: `[0, L1]` from `RF2`, rotated by `q₃` (hip);
* `RF4`: `[0, -L1]` from `RF3`, rotated by `q₄` (swing knee).

Bodies: a tibia (`m2`, length `L2`) and a femur (`m1`, length `L1`) on each leg, modelled as thin
rods `I = m L² / 12` about their centres; motor point masses `mM` at the grounded and swing knees and
`2 mM` at the hip (in the frames `RF2`, `RF4`, `RF3`), with no rotor inertia; the swing foot is a
massless point at the end of the swing tibia. This is exactly the model of `four_link_biped.urdf`.

The input matrix is `B = [0 0 0; I₃]`: the last three coordinates are actuated, `q₁` is not.

Returns `(; sys, nf, B, parameters)` with `parameters` the `NamedTuple` of symbolic variables
`(L1, L2, m1, m2, mM, g)`. Building it takes a few seconds (the normal form is the costly part):
build it once and pass it to the functions below.
"""
function symbolic_model()
    Symbolics.@variables q[1:4] qd[1:4] L1 L2 m1 m2 mM g
    sys = MultibodySystem(collect(q), collect(qd); gravity = g)

    rf1 = add_frame!(sys, :RF1, sys.world, [0, 0], q[1])
    rf2 = add_frame!(sys, :RF2, rf1, [0, L2], q[2])
    rf3 = add_frame!(sys, :RF3, rf2, [0, L1], q[3])
    rf4 = add_frame!(sys, :RF4, rf3, [0, -L1], q[4])

    # Links: position = center of mass, inertia of a thin rod about its center.
    add_point!(sys, :link1, rf1, [0, L2 / 2]; mass = m2, inertia = m2 * L2^2 / 12)
    add_point!(sys, :link2, rf2, [0, L1 / 2]; mass = m1, inertia = m1 * L1^2 / 12)
    add_point!(sys, :link3, rf3, [0, -L1 / 2]; mass = m1, inertia = m1 * L1^2 / 12)
    add_point!(sys, :link4, rf4, [0, -L2 / 2]; mass = m2, inertia = m2 * L2^2 / 12)
    # Motors: point masses at the joints.
    add_point!(sys, :grounded_knee, rf2, [0, 0]; mass = mM, inertia = 0)
    add_point!(sys, :hip, rf3, [0, 0]; mass = 2 * mM, inertia = 0)
    add_point!(sys, :swing_knee, rf4, [0, 0]; mass = mM, inertia = 0)
    add_point!(sys, :swing_foot, rf4, [0, -L2]; mass = 0, inertia = 0)

    B = Num[
        0 0 0
        1 0 0
        0 1 0
        0 0 1
    ]
    nf = normal_form(sys, B)
    return (; sys, nf, B, parameters = (; L1, L2, m1, m2, mM, g))
end

"""
    default_parameters() -> NamedTuple

Numerical values of the parameters of [`symbolic_model`](@ref), those of `four_link_biped.urdf`:
femur `L1 = 0.202` m and tibia `L2 = 0.172` m (different but similar), femur/tibia masses
`m1 = 0.0202`, `m2 = 0.0172` kg and motor mass `mM = 0.4` kg, which is by far the dominant mass (the
hip carries `2 mM`), `g = 9.81`.
"""
default_parameters() =
    (; L1 = 0.202, L2 = 0.172, m1 = 0.0202, m2 = 0.0172, mM = 0.4, g = 9.81)

"""
    substitute_parameters(model, expressions, values = default_parameters())

Replace the symbolic parameters of `model` (`model.parameters`) by the numbers `values` (a
`NamedTuple` with the same keys) in `expressions`, a symbolic expression or an array of them. What
is left is symbolic in the coordinates (and the inputs) only: this is the step "build the model
symbolically, put the numbers in later".
"""
function substitute_parameters(model, expressions, values = default_parameters())
    assignments = Dict(
        getfield(model.parameters, k) => float(getfield(values, k)) for
        k in keys(model.parameters)
    )
    expressions isa AbstractArray &&
        return Symbolics.substitute.(expressions, Ref(assignments))
    return Symbolics.substitute(expressions, assignments)
end

# Compile a numerical symbolic expression (no free parameter) into a Julia function of `args`.
function _compile(expression, args...)
    f = Symbolics.build_function(expression, args...; expression = Val(false), cse = true)
    return f isa Tuple ? first(f) : f
end

# ============================================================
# Geometry & kinematics (derived from the symbolic model)
# ============================================================

"""
    Geometry(model = symbolic_model(), values = default_parameters())

Segment lengths of the four-link biped (`tibia = values.L2`, `femur = values.L1`) and the forward
kinematics of its joints, **compiled from the symbolic positions** of `model`: there is no second
copy of the kinematics to keep consistent with the dynamics. `default_geometry()` is the one of
`default_parameters()`.
"""
struct Geometry{T <: Real, F}
    tibia::T
    femur::T
    joint_positions::F
end

function Geometry(model = symbolic_model(), values = default_parameters())
    sys = model.sys
    names = (:grounded_knee, :hip, :swing_knee, :swing_foot)
    positions = reduce(
        vcat,
        [substitute_parameters(model, Multibody2D.position(sys[n]), values) for n in names],
    )
    return Geometry(float(values.L2), float(values.L1), _compile(positions, sys.q))
end

default_geometry() = Geometry()

"""
    cartesian_coordinates(q, geometry)

World `(x, z)` positions of the joints for `q = [q₁, q₂, q₃, q₄]`, grounded ankle
at the origin: `(grounded_foot, grounded_knee, hip, swing_knee, swing_foot)`.

Absolute segment angles from the upward vertical (positive counter-clockwise, i.e.
towards `-x`): grounded tibia `q₁`, grounded femur `q₁ + q₂`; the swing femur
`φ₃ = q₁ + q₂ + q₃` and swing tibia `φ₄ = φ₃ + q₄` hang downwards. `q = 0` is the robot
standing straight with the two legs side by side.
"""
function cartesian_coordinates(q::SVector{4}, geometry::Geometry)
    p = geometry.joint_positions(q)
    return (
        SVector(0.0, 0.0),
        SVector(p[1], p[2]),
        SVector(p[3], p[4]),
        SVector(p[5], p[6]),
        SVector(p[7], p[8]),
    )
end

swing_foot_position(q::SVector{4}, geometry::Geometry) =
    cartesian_coordinates(q, geometry)[5]

# Per-coordinate gradient bounds of the swing-foot map on the state
# `[ankle, q₁, q₂, q₃, q₄]`: each joint rotates the rest of the chain about itself, so
# the foot moves by at most the pivot-to-foot distance per radian (ankle: both
# legs; grounded knee: femur + swing leg; hip: swing leg; swing knee: tibia). The
# foot does not depend on the ankle velocity or momentum.
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
# The design model, in normal form
# ============================================================

"""
    reduced_dynamics(model) -> NamedTuple

The 5-dimensional velocity-controlled model extracted from the normal form `x = (q, v_a, p_u)` of
`model`, with the **ankle momentum first**: the state is `x = (p₁, q₁, q₂, q₃, q₄)` (`p₁ = ∂L/∂q̇₁`, the
momentum of the unactuated ankle), the input `u = v_a = (q̇₂, q̇₃, q̇₄)` and

    ẋ = (ṗ₁, q̇₁, v₂, v₃, v₄),     q̇₁ = D₁₁⁻¹ (p₁ - D₁ₐ v_a),     ṗ₁ = ∂L/∂q₁,

i.e. the `p_u` and `q̇` rows of the normal form. The rows of the actuated velocities are dropped: the
commanded velocities are held constant over each sampling period. `ṗ₁` is independent of the
torques, so no torque appears in the result.

Returns `(; f, x, u)`, all symbolic, with the parameters of `model` still free: put the numbers in
with [`substitute_parameters`](@ref).
"""
function reduced_dynamics(model)
    nf = model.nf
    nf.actuated == [2, 3, 4] && nf.unactuated == [1] ||
        error("expected the ankle q₁ to be the only unactuated coordinate")
    # nf.state = [q₁ q₂ q₃ q₄ | v₂ v₃ v₄ | p₁]; `drift` is ẋ at τ = 0 (its rows do not depend on τ).
    order = [8, 1, 2, 3, 4]
    x = nf.state[order]
    u = nf.state[5:7]
    f = nf.drift[order]
    for e in f
        isempty(intersect(Symbolics.get_variables(e), nf.input)) ||
            error("the reduced dynamics must not depend on the torques")
    end
    return (; f, x, u)
end

# ============================================================
# Domain carving: Cartesian obstacle + ground → state-space cells
# ============================================================

_obstacle_list(::Nothing) = LazySets.LazySet[]
_obstacle_list(obstacle::LazySets.LazySet) = [obstacle]
_obstacle_list(obstacles::AbstractVector) = obstacles

"The name under which the ground constraint (swing foot above `z = 0`) appears in a [`Carving`](@ref)."
const GROUND_NAME = "ground"

"""
    Carving

The result of [`infeasible_cells`](@ref): `cells`, the `MP.CellUnion` of the 5-D cells to exclude from
the state domain, and `causes`, the same cells grouped by the constraint that removed them, as
`name => positions` pairs (`"ground"` first, then the obstacles by their `obstacle_name`). A cell can
be removed by several constraints. The positions are those of the 4-D angle cells
`(q₁, q₂, q₃, q₄)`: the swing foot, hence the carving, does not depend on the ankle coordinate.
"""
struct Carving
    cells::MP.CellUnion
    causes::Vector{Pair{String, Set{NTuple{4, Int}}}}
end

"""
    carved_by(carving, pos) -> Vector{String}

The names of the constraints of `carving` that removed the 5-D cell `pos`.
"""
function carved_by(carving::Carving, pos::NTuple{5, Int})
    return [name for (name, cells) in carving.causes if pos[2:5] in cells]
end

# The obstacles grouped by name, in order of appearance: `names` is one name for all of them, or
# one name per obstacle.
function _obstacle_groups(obstacles, names::AbstractString)
    return isempty(obstacles) ? Pair{String, typeof(obstacles)}[] :
           [String(names) => obstacles]
end

function _obstacle_groups(obstacles, names::AbstractVector)
    length(names) == length(obstacles) || error(
        "`obstacle_name` must be one name or one name per obstacle: got " *
        "$(length(names)) names for $(length(obstacles)) obstacles",
    )
    groups = Pair{String, typeof(obstacles)}[]
    for (obstacle, name) in zip(obstacles, names)
        k = findfirst(g -> first(g) == name, groups)
        if k === nothing
            push!(groups, String(name) => eltype(obstacles)[obstacle])
        else
            push!(last(groups[k]), obstacle)
        end
    end
    return groups
end

"""
    infeasible_cells(X, state_grid, obstacle, geometry; obstacle_name = "obstacle") -> Carving

The cells of `X` (a 5-D hyperrectangle) to exclude, as in the fully actuated model: cells where the
swing foot *might* touch the Cartesian `obstacle` (a `LazySet` or a vector of them; sound,
`MP.image_blocked_cells`), and cells where the swing foot is surely below the ground (tolerant, so
the footholds on the ground stay). Returns a [`Carving`](@ref): the `cells` to remove (what
[`system`](@ref) and `UT.set_minus` take), and which constraint removed which, so that
[`start_state`](@ref) can say what carved the start cell.

The constraints are named: the ground constraint `"ground"` ([`GROUND_NAME`](@ref)), and the
obstacles by `obstacle_name`, a string for all of them (the default `"obstacle"`) or one string per
obstacle, e.g. `"lift constraint"` for the obstacle at the origin that forces the swing foot to lift.
"""
function infeasible_cells(
    X::LazySets.AbstractHyperrectangle,
    state_grid::MP.GridFree,
    obstacle,
    geometry::Geometry;
    obstacle_name = "obstacle",
)
    grad = swing_foot_gradient_bound(geometry)
    dev = swing_foot_deviation_bound(geometry, MP.get_h(state_grid) ./ 2)

    # The swing foot does not depend on the ankle velocity/momentum: decide on the 4-D angle cells once.
    q_grid = MP.GridFree(MP.get_origin(state_grid)[IDX_Q], MP.get_h(state_grid)[IDX_Q])
    Q = LazySets.Hyperrectangle(;
        low = LazySets.low(X)[IDX_Q],
        high = LazySets.high(X)[IDX_Q],
    )
    q_foot(q) = swing_foot_position(SVector{4}(q), geometry)
    below_ground = MP.cells_where(q_grid, Q) do pos
        return q_foot(MP.get_coord_by_pos(q_grid, pos))[2] < -dev
    end
    causes = Pair{String, Set{NTuple{4, Int}}}[GROUND_NAME => below_ground.positions]
    for (name, group) in _obstacle_groups(_obstacle_list(obstacle), obstacle_name)
        blocked = MP.image_blocked_cells(q_foot, grad[IDX_Q], q_grid, Q, group)
        push!(causes, name => blocked.positions)
    end
    removed_q = reduce(union, last.(causes))

    # Lift to 5-D along the ankle axis.
    h = MP.get_h(state_grid)
    k_lo = round(
        Int,
        (LazySets.low(X)[IDX_ANKLE] - MP.get_origin(state_grid)[IDX_ANKLE]) / h[IDX_ANKLE],
    )
    k_hi = round(
        Int,
        (LazySets.high(X)[IDX_ANKLE] - MP.get_origin(state_grid)[IDX_ANKLE]) / h[IDX_ANKLE],
    )
    removed = Set{NTuple{5, Int}}()
    sizehint!(removed, length(removed_q) * (k_hi - k_lo + 1))
    for pos in removed_q, k in k_lo:k_hi
        push!(removed, (k, pos...))
    end
    return Carving(MP.CellUnion(state_grid, removed), causes)
end

# ============================================================
# Target set: configurations placing the swing foot at a foothold
# ============================================================

# Two-link inverse kinematics of the swing leg: absolute angles (φ₃, φ₄), counter-clockwise,
# placing the swing foot at `hip + d`, with the knee bent backwards (`q₄ = φ₄ − φ₃ ≤ 0`, the
# human-like knee for walking towards `+x`), or `nothing` if unreachable.
function _swing_leg_ik(d::SVector{2}, geometry::Geometry)
    L1, L2 = geometry.femur, geometry.tibia
    c = (d[1]^2 + d[2]^2 - L1^2 - L2^2) / (2 * L1 * L2)
    -1.0 <= c <= 1.0 || return nothing
    # Same planar two-link problem as the fully actuated model, whose angles ψ are the
    # counter-clockwise ones used here.
    ψ4_rel = -acos(c)
    A = L1 + L2 * cos(ψ4_rel)
    B = L2 * sin(ψ4_rel)
    ψ3 = atan(d[1], -d[2]) - atan(B, A)
    return ψ3, ψ3 + ψ4_rel
end

"""
    target_cells(state_grid, foothold, geometry, X; nsub = 3, velocity = nothing, tolerance = 0.0)

5-D grid cells containing a configuration with the grounded foot at the origin
and the swing foot at the Cartesian point `foothold`. `(q₁, q₂)` are swept over
their cells (`nsub` samples per axis per cell), `(q₃, q₄)` follow by inverse
kinematics. The first state coordinate (the ankle momentum `p₁`) is free: every cell of `X` along it is
accepted, or, with `velocity = (vlo, vhi)`, every cell lying entirely inside that range.

With `tolerance > 0` the foot only has to land *near* the foothold: within `tolerance` metres of
it along `x`, and between the foothold height and `tolerance` above it (never below the ground
level of the foothold), sampled every centimetre (so that the targets of increasing tolerances,
multiples of 1 cm, are nested). Without a tolerance the target is a thin
sliver of the angle space: every step the over-approximated ankle spreads over several cells, and a
controller that is certified against the worst case cannot land inside a sliver (see the experiment
document).
"""
function target_cells(
    state_grid::MP.GridFree,
    foothold::SVector{2},
    geometry::Geometry,
    X::LazySets.AbstractHyperrectangle;
    nsub::Int = 3,
    velocity = nothing,
    tolerance::Real = 0.0,
)
    tolerance >= 0 || error("`tolerance` must be non-negative, got $tolerance")
    h = MP.get_h(state_grid)
    lo = SVector{5}(LazySets.low(X))
    hi = SVector{5}(LazySets.high(X))
    ankle_values = range(
        MP.get_coord_by_pos(state_grid, MP.get_pos_by_coord(state_grid, lo))[1],
        hi[1];
        step = h[1],
    )
    if velocity !== nothing
        vlo, vhi = velocity
        tol = 1e-9 * h[1]
        ankle_values = [
            v for
            v in ankle_values if vlo - tol <= v - h[1] / 2 && v + h[1] / 2 <= vhi + tol
        ]
        isempty(ankle_values) && error(
            "the velocity range $velocity contains no whole cell of width $(h[1]); " *
            "the target must be made of cells entirely inside it",
        )
    end

    # Offsets of the foothold, every centimetre up to the tolerance: the foot within `tolerance` metres
    # along `x`, and from the foothold height up to `tolerance` above it.
    n = ceil(Int, tolerance / 0.01 - 1e-9)
    footholds =
        tolerance == 0 ? [foothold] :
        [
            foothold + SVector(dx, dz) for
            dx in range(-tolerance, tolerance; length = 2n + 1) for
            dz in range(0, tolerance; length = n + 1)
        ]
    found = Set{NTuple{5, Int}}()
    for fh in footholds, q1 in range(lo[2], hi[2]; step = h[2] / nsub)
        for q2 in range(lo[3], hi[3]; step = h[3] / nsub)
            hip = cartesian_coordinates(SVector(q1, q2, 0.0, 0.0), geometry)[3]
            ik = _swing_leg_ik(fh - hip, geometry)
            ik === nothing && continue
            φ3, φ4 = ik
            q3 = φ3 - (q1 + q2)
            q4 = φ4 - φ3
            for a in ankle_values
                x = SVector(a, q1, q2, q3, q4)
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
coordinates follow the nonlinear dynamics and are rounded to the grid. `dq̇1` is the step of the
first state coordinate, the ankle momentum `p₁`.

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

# The range of the first state coordinate (the ankle momentum) of `velocity`, clipped to the box `X_box`;
# `nothing` means the whole interval of definition.
function _velocity_range(velocity, X_box, what)
    lo, hi = LazySets.low(X_box)[IDX_ANKLE], LazySets.high(X_box)[IDX_ANKLE]
    velocity === nothing && return (lo, hi)
    vlo, vhi = velocity
    vlo <= vhi ||
        error("`$what` must be an interval (vlo, vhi) with vlo ≤ vhi, got $velocity")
    (vhi >= lo && vlo <= hi) ||
        error("`$what` = $velocity does not intersect the domain interval ($lo, $hi)")
    return (max(vlo, lo), min(vhi, hi))
end

"""
    step_problem(concrete_system, state_grid; q0, foothold, geometry,
                 initial_velocity = nothing, target_velocity = nothing, foothold_tolerance = 0.0)

One footstep as a reach specification: from the configuration `q0 = [q₁, q₂, q₃, q₄]`, with
the ankle momentum `p₁` (the first state coordinate) anywhere in `initial_velocity`, reach a
configuration placing the swing foot at `foothold`, with the ankle momentum anywhere in
`target_velocity`.

Nothing is imposed on the ankle at the start or at the end: the swing foot has no
reaction force from the ground, so requiring rest is needlessly conservative. Both
`initial_velocity` and `target_velocity` are `(vlo, vhi)` ranges on the first state
coordinate and default to `nothing`, the whole interval of definition of the state domain.
`foothold_tolerance` (metres) lets the foot land near the foothold instead of exactly on it, see
[`target_cells`](@ref).
The initial set is the box around the cell of `q0` (snapped to the grid, the box lies
strictly inside the cell) times the initial range: the controller must work from
**every** cell of the first coordinate the range touches. The target is the set of cells lying
entirely inside the target range.

Errors if the cell of `q0` is not in the state domain (outside its bounds or carved out): the
initial set would be empty, which the synthesis can only report as a failure, without saying why.
"""
function step_problem(
    concrete_system,
    state_grid;
    q0::SVector{4},
    foothold::SVector{2} = SVector(0.2, 0.0),
    geometry::Geometry = default_geometry(),
    initial_velocity = nothing,
    target_velocity = nothing,
    foothold_tolerance::Real = 0.0,
)
    h = MP.get_h(state_grid)
    X_box = UT._outer_box(concrete_system.X)
    v_init = _velocity_range(initial_velocity, X_box, "initial_velocity")
    v_target = _velocity_range(target_velocity, X_box, "target_velocity")

    # Snap the configuration to the center of its cell so that the initial box is inside it.
    snapped = MP.get_coord_by_pos(
        state_grid,
        MP.get_pos_by_coord(state_grid, SVector(zero(eltype(q0)), q0...)),
    )
    qc = snapped[IDX_Q]
    r = 0.25 .* h[IDX_Q]
    _I_ = LazySets.Hyperrectangle(;
        low = SVector(v_init[1], (qc .- r)...),
        high = SVector(v_init[2], (qc .+ r)...),
    )
    SVector((v_init[1] + v_init[2]) / 2, qc...) in concrete_system.X || error(
        "the cell of the initial configuration q0 = $q0 is not in the state domain (outside " *
        "its bounds, or carved out: the swing foot is within the obstacle/ground margin); " *
        "choose another configuration or refine the grid",
    )

    _T_ = target_set(
        state_grid,
        foothold,
        geometry,
        X_box;
        velocity = target_velocity === nothing ? nothing : v_target,
        tolerance = foothold_tolerance,
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
# Size of the abstraction
# ============================================================

"CPU time of the local growth bound per `(cell, input)` pair, in seconds (measured on the biped: ≈ 0.2 ms)."
const CPU_SECONDS_PER_PAIR = 2.0e-4

"Average number of transitions per `(cell, input)` pair of the biped abstraction (measured: ≈ 4)."
const TRANSITIONS_PER_PAIR = 4

"""
    abstraction_size(X, state_grid, removed = nothing; ninputs = 27,
                     threads = Threads.nthreads()) -> NamedTuple

How big the abstraction of the box `X` will be, before computing it: `cells_in_box` is the number of
cells of `state_grid` lying entirely inside `X` (the cells the abstraction keeps, `MP.INNER`), `cells`
the number left after removing the carved ones (`removed`: a [`Carving`](@ref) or its `cells`),
`pairs = cells · ninputs` the number of `(cell, input)` pairs whose image must be computed, and
`transitions` and `minutes` rough estimates of the number of transitions (`TRANSITIONS_PER_PAIR` per
pair) and of the wall time on `threads` threads (`CPU_SECONDS_PER_PAIR` per pair).

The cost is linear in `pairs`, and `cells` is the *product* of the cell counts along the five axes:
widening every bound of the state box by a factor 2.5 multiplies the cells by up to 100.
"""
function abstraction_size(
    X::LazySets.AbstractHyperrectangle,
    state_grid,
    removed = nothing;
    ninputs::Int = 27,
    threads::Int = Threads.nthreads(),
)
    o, h = MP.get_origin(state_grid), MP.get_h(state_grid)
    lo, hi = LazySets.low(X), LazySets.high(X)
    # The cells whose box `c ± h/2` lies inside `[lo, hi]`.
    kmin = ntuple(i -> ceil(Int, (lo[i] - o[i]) / h[i] + 0.5 - 1e-9), 5)
    kmax = ntuple(i -> floor(Int, (hi[i] - o[i]) / h[i] - 0.5 + 1e-9), 5)
    cells_in_box = prod(max.(kmax .- kmin .+ 1, 0))
    cells = cells_in_box
    if removed !== nothing
        positions = removed isa Carving ? removed.cells.positions : removed.positions
        cells -= count(pos -> all(i -> kmin[i] <= pos[i] <= kmax[i], 1:5), positions)
    end
    pairs = cells * ninputs
    return (;
        cells_in_box,
        cells,
        pairs,
        transitions = TRANSITIONS_PER_PAIR * pairs,
        minutes = pairs * CPU_SECONDS_PER_PAIR / threads / 60,
    )
end

# ============================================================
# Start state and closed loop
# ============================================================

"""
    start_state(state_grid, q0, X; velocity = nothing, removed = nothing,
                geometry = default_geometry()) -> SVector{5}

Concrete start of a closed-loop run: the configuration `q0 = [q₁, q₂, q₃, q₄]` with the first
state coordinate (the ankle momentum) in the middle of `velocity` (a `(vlo, vhi)`
range clipped to the box `X`; by default the whole interval of `X`), snapped to the centre of its
grid cell so that the run starts where the abstraction does.

With `removed` (the [`Carving`](@ref) of [`infeasible_cells`](@ref), or just its `cells`), it errors,
before any expensive abstraction, if the start cell was carved out. With a `Carving` the error names
the constraints that did it (`"ground"`, or the `obstacle_name`) and tells how the swing foot at the
start violates each, so that the culprit is visible in the output; a finer grid or another `q0` is
needed.
"""
function start_state(
    state_grid,
    q0::SVector{4},
    X::LazySets.AbstractHyperrectangle;
    velocity = nothing,
    removed = nothing,
    geometry::Geometry = default_geometry(),
)
    v = _velocity_range(velocity, X, "initial_velocity")
    x0 = MP.get_coord_by_pos(
        state_grid,
        MP.get_pos_by_coord(state_grid, SVector((v[1] + v[2]) / 2, q0...)),
    )
    pos = MP.get_pos_by_coord(state_grid, x0)
    cells = removed isa Carving ? removed.cells : removed
    if cells !== nothing && pos in cells
        margin = swing_foot_deviation_bound(geometry, MP.get_h(state_grid) ./ 2)
        foot = swing_foot_position(SVector{4}(x0[2:5]), geometry)
        error(_carved_start_message(removed, pos, foot, margin))
    end
    return x0
end

function _carved_start_message(removed, pos, foot, margin)
    cm = round(100 * margin; digits = 1)
    where_ =
        "the swing foot at the start is at (x, z) = ($(round(foot[1]; digits = 3)), " *
        "$(round(foot[2]; digits = 3))) m"
    advice = "choose another q0 or refine the grid (the margin scales with the cell sizes)"
    removed isa Carving ||
        return "the start cell is carved out (the swing foot at x0 is within " *
               "the $cm cm margin of the obstacle or below the ground): " *
               "refine the grid"
    reasons = map(carved_by(removed, pos)) do name
        return name == GROUND_NAME ?
               "\"$name\" (the swing foot is more than $cm cm below the ground)" :
               "\"$name\" (the swing foot is within $cm cm of it)"
    end
    return "the start cell is carved out by " *
           join(reasons, " and ") *
           ": $where_; $advice"
end

"""
    closed_loop(optimizer, problem, x0; horizon = 400, ankle_name = "ankle") -> Trajectory

Simulate the concrete (RK4) system from the state `x0` under the controller that `optimizer`
synthesized for `problem` (the `PR.OptimalControlProblem` of [`step_problem`](@ref)), for at most
`horizon` steps. The run stops when the target set of `problem` is reached, or when the controller
has no action at the current state, which happens when the concrete ankle has drifted away from the
rounded abstract one. It prints the number of steps, whether the target was reached, whether the
run stopped outside the controller's domain, and the ranges of the ankle angle `q₁` and of the first
state coordinate (named `ankle_name`, e.g. `"p₁"`) along the run.

This is the end-to-end check of a controller and of the abstraction under it. The returned
trajectory goes to [`plot_footstep`](@ref).
"""
function closed_loop(optimizer, problem, x0; horizon::Int = 400, ankle_name = "ankle")
    controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
    system = MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system"))
    reached(x) = x ∈ problem.target_set
    traj =
        ST.get_closed_loop_trajectory(system, controller, x0, horizon; stopping = reached)
    xs = collect(ST.states(traj))
    println("steps: ", length(xs) - 1, ", reached: ", reached(xs[end]))
    reached(xs[end]) ||
        length(xs) > horizon ||
        println("stopped at x = ", xs[end], ": outside the controller's domain")
    println(
        "ankle along the run: q₁ ∈ ",
        extrema(x[2] for x in xs),
        ", ",
        ankle_name,
        " ∈ ",
        extrema(x[1] for x in xs),
    )
    return traj
end

# ============================================================
# Visualization
# ============================================================

"""
    draw_robot!(fig, x; geometry, alpha)

Draw the grounded (blue) and swing (red) legs of the state `x = [ankle, q₁, q₂, q₃, q₄]`.
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

# The slice, at the current `(q₁, q₂)`, of the carved cells in the `(q₃, q₄)` plane, as a background
# `(panel, x) -> panel` for the state panel of the dashboard: the carving does not depend on the
# ankle coordinate.
function _blocked_slice(state_grid, removed, x_lb, x_ub)
    h = MP.get_h(state_grid)
    k3_range = round(Int, x_lb[4] / h[4]):round(Int, x_ub[4] / h[4])
    k4_range = round(Int, x_lb[5] / h[5]):round(Int, x_ub[5] / h[5])
    return function (panel, x)
        pos = MP.get_pos_by_coord(state_grid, x)
        shapes = Plots.Shape[]
        for k3 in k3_range, k4 in k4_range
            (pos[1], pos[2], pos[3], k3, k4) in removed || continue
            c = MP.get_coord_by_pos(state_grid, (pos[1], pos[2], pos[3], k3, k4))
            push!(
                shapes,
                Plots.Shape(
                    [c[4] - h[4] / 2, c[4] + h[4] / 2, c[4] + h[4] / 2, c[4] - h[4] / 2],
                    [c[5] - h[5] / 2, c[5] - h[5] / 2, c[5] + h[5] / 2, c[5] + h[5] / 2],
                ),
            )
        end
        isempty(shapes) ||
            Plots.plot!(panel, shapes; color = :black, alpha = 0.35, lw = 0, label = "")
        return panel
    end
end

"""
    plot_footstep(traj, state_grid, removed; x_lb, x_ub, tstep, file_prefix,
                  geometry = default_geometry(), obstacle = nothing, foothold = nothing,
                  ankle_label = "ankle", title, xlims = (-0.6, 0.6), ylims = (-0.05, 0.55),
                  snapshot_every = 5, dynamics = nothing, subframes = 1, fps = nothing)

Save the three figures of a closed-loop footstep `traj` (a `Dionysos.System.Trajectory` with
inputs, in the state `[ankle, q₁, q₂, q₃, q₄]`) and return their paths as
`(; snapshot, ankle, animation)`, or `nothing` if the trajectory has no input (the controller had no
action at the start: there is nothing to draw):

* `file_prefix * ".png"`: the robot every `snapshot_every` steps, with the ground, the `obstacle` and
  the `foothold`;
* `file_prefix * "_ankle.png"`: the ankle angle `q₁` and the ankle coordinate (labelled
  `ankle_label`, e.g. `"p₁ (N m s)"`) against time, `tstep` being the sampling
  period;
* `file_prefix * ".gif"`: the animated dashboard (`Dionysos.animate_trajectory_dashboard`) with the
  robot, the state panel `(q₃, q₄)` over the box `[x_lb, x_ub]` with the slice of the carved cells
  `removed` ([`infeasible_cells`](@ref)) at the current `(q₁, q₂)` as a grey background, and the
  input panel `(v₃, v₄)`; `title` is its title.

A closed loop has one state per sampling period, which makes a short run (a few steps) a flash. With
`dynamics`, the continuous dynamics `(x, u) -> ẋ` of the design model (e.g. `lgb.dynamics.f`), and
`subframes > 1`, each step is re-integrated (RK4, the input held constant) at `subframes` intermediate
instants for the ankle plot and the animation, which then has `subframes` frames per step. `fps`
(frames per second of the gif) defaults to the value that makes it last about four seconds, between 5
and 20.
"""
function plot_footstep(
    traj,
    state_grid,
    removed;
    x_lb,
    x_ub,
    tstep,
    file_prefix,
    geometry::Geometry = default_geometry(),
    obstacle = nothing,
    foothold = nothing,
    ankle_label = "ankle",
    title = "Underactuated four-link biped — footstep",
    xlims = (-0.6, 0.6),
    ylims = (-0.05, 0.55),
    snapshot_every::Int = 5,
    dynamics = nothing,
    subframes::Int = 1,
    fps = nothing,
)
    ST.inputs(traj) === nothing && return nothing
    xs = collect(ST.states(traj))
    us = collect(ST.inputs(traj))
    isempty(us) && return nothing
    plot_robot = system_plot!(; geometry, obstacle, foothold, xlims, ylims)
    snapshot = file_prefix * ".png"
    ankle = file_prefix * "_ankle.png"
    animation = file_prefix * ".gif"

    fig = Plots.plot(; aspect_ratio = :equal, legend = false)
    rest = zero(us[1])
    for (k, x) in enumerate(xs)
        (k - 1) % snapshot_every == 0 || k == length(xs) || continue
        plot_robot(fig, x, k <= length(us) ? us[k] : rest)
    end
    Plots.savefig(fig, snapshot)
    println("saved ", snapshot)

    # With the dynamics, the intermediate instants of every step, for a smoother ankle plot and animation.
    if dynamics !== nothing && subframes > 1
        frames = typeof(xs[1])[]
        frame_inputs = eltype(us)[]
        for k in eachindex(us), j in 0:(subframes - 1)
            push!(
                frames,
                ST.runge_kutta4(dynamics, xs[k], us[k], tstep * j / subframes, max(2, j)),
            )
            push!(frame_inputs, us[k])
        end
        push!(frames, xs[end])
        traj = ST.Trajectory(frames; inputs = frame_inputs)
        xs, us, tstep = frames, frame_inputs, tstep / subframes
    end

    # The unactuated ankle: angle and ankle coordinate against time.
    ts = tstep .* (0:(length(xs) - 1))
    Plots.savefig(
        Plots.plot(
            ts,
            [[x[2] for x in xs], [x[1] for x in xs]];
            label = permutedims(["q₁ (rad)", ankle_label]),
            xlabel = "time (s)",
            title = "Unactuated ankle",
        ),
        ankle,
    )
    println("saved ", ankle)

    DI.animate_trajectory_dashboard(
        plot_robot,
        traj;
        Δt = tstep,
        xdims = (4, 5),
        udims = (2, 3),
        xlabel_state = "q₃ (hip)",
        ylabel_state = "q₄ (swing knee)",
        xlabel_input = "v₃ (hip)",
        ylabel_input = "v₄ (swing knee)",
        xlims_state = (x_lb[4], x_ub[4]),
        ylims_state = (x_lb[5], x_ub[5]),
        state_background! = _blocked_slice(state_grid, removed, x_lb, x_ub),
        title,
        fps = fps === nothing ? clamp(round(Int, length(xs) / 4), 5, 20) : fps,
        filename = animation,
    )
    println("saved ", animation)
    return (; snapshot, ankle, animation)
end

end # module
