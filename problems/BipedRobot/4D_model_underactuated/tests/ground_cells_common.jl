# The classification of the cells of a 4-D box by the ground constraint (swing foot above `z = 0`),
# shared by the scripts of the ground-cells study (`ground_cells_classify.jl`,
# `check_ground_cells_extremes.jl`).
#
# Not a script: it loads `robot_problem.jl` (as `RP`) and `foot_preimage_sweep.jl`, and defines
# functions. The range `[z_min, z_max]` of the foot height over a cell is the one of
# `RP.infeasible_cells`, by one of its bounds `RP.FOOT_BOUNDS` (`:linear`, `:second_order_global`,
# `:second_order_cell`, `:corners`; their maths is in the docstring of `RP.FOOT_BOUNDS`). The former
# rule of `infeasible_cells` is kept here as a reference, `:margin`: the range `z_c ± dev`, with
# `dev = Σᵢ Gᵢ rᵢ` the Lipschitz margin, one worst-case constant for every cell, so that its `:below`
# cells are exactly the cells the former rule removed. A cell is then
#
# * `:below` (code 1, `z_max < 0`): removed by `infeasible_cells`;
# * `:across` (code 2, `z_min < 0 ≤ z_max`): ambiguous, removed only with `ground_ambiguous = :remove`;
# * `:above` (code 3, `z_min ≥ 0`): kept.
#
# Only the class of each cell is stored (one byte, in an array indexed by position), so that the
# classification scales to tens of millions of cells.

using LinearAlgebra
using Random
using StaticArrays
import Dionysos
import LazySets

const MP = Dionysos.Mapping

include(joinpath(@__DIR__, "..", "robot_problem.jl"))
import .UnderactuatedRobotProblem as RP

# The sweep of the pre-image of a point, and the cells it crosses.
include(joinpath(@__DIR__, "foot_preimage_sweep.jl"))

# The ranges of the study: the former rule and the bounds of `RP.infeasible_cells`.
const STUDY_BOUNDS = (:margin, RP.FOOT_BOUNDS...)

function check_study_bound(bound)
    bound in STUDY_BOUNDS || error(
        "the ground bound must be one of $(join(repr.(STUDY_BOUNDS), ", ")), got $(repr(bound))",
    )
    return nothing
end

# A label of the range, for the printed reports and the plots.
function bound_label(bound)
    return bound == :margin ? "worst-case margin (former rule)" :
           bound == :linear ? "linear range" :
           bound == :second_order_global ? "linear range + global second-order remainder" :
           bound == :second_order_cell ? "per-cell second order" : "corner values"
end

# The range of the foot height over the cell `q ± r`: the one of `RP.infeasible_cells`, or the former
# rule `z_c ± dev`. `R2, R3` are the remainder bounds of `RP._foot_remainders`.
@inline function cell_height_range(q, r, geometry, bound, R2, R3, dev)
    if bound == :margin
        z = RP._joints(q, geometry)[5][2]
        return z - dev, z + dev
    end
    return RP._height_range(q, r, geometry, bound, R2, R3)
end

# Classify every cell of the box `[q_lb, q_ub]` (`INNER`, like `RP.infeasible_cells`) as `:below`,
# `:across` or `:above` by the range `foot_bound` (one of `STUDY_BOUNDS`). Returns the class array
# (`kmin:kmax` along each axis), the number of cells of each class, the mean and largest half-width of
# the range, the grid, the half-widths `r`, the remainder bounds `R2` (`remainder`) and `R3` (`third`),
# and the margin `dev` of the former rule.
function ground_cells_by_bound(
    q_lb,
    q_ub,
    cell_size,
    geometry;
    foot_bound = :second_order_cell,
)
    check_study_bound(foot_bound)
    q_grid = MP.GridFree(zero(cell_size), cell_size)
    Q = LazySets.Hyperrectangle(; low = q_lb, high = q_ub)
    ranges = MP.get_pos_from_set(q_grid, Q, MP.INNER).iterators
    kmin = first.(ranges)
    kmax = last.(ranges)
    r = SVector{4, Float64}(cell_size ./ 2)
    R2, R3 = RP._foot_remainders(geometry, r)
    # The ankle has no weight in the margin (its gradient bound is 0), hence the 0.0 in first place.
    dev = RP.swing_foot_deviation_bound(geometry, SVector(0.0, r...))

    class = Array{UInt8, 4}(undef, length.(ranges)...)
    counts = zeros(Int, 3)
    sum_half = max_half = 0.0
    for I in CartesianIndices(class)
        pos = Tuple(I) .+ kmin .- 1
        q = SVector{4, Float64}(MP.get_coord_by_pos(q_grid, pos))
        lo, hi = cell_height_range(q, r, geometry, foot_bound, R2, R3, dev)
        code = hi < 0 ? 0x01 : lo >= 0 ? 0x03 : 0x02
        class[I] = code
        counts[code] += 1
        half = (hi - lo) / 2
        sum_half += half
        max_half = max(max_half, half)
    end
    return (;
        class,
        kmin,
        kmax,
        counts,
        mean_half = sum_half / length(class),
        max_half,
        q_grid,
        r,
        remainder = R2,
        third = R3,
        dev,
        foot_bound,
    )
end

# A second pass over the cells of `cells` (the result of `ground_cells_by_bound`): how many cells the
# former rule removes (the foot of the centre below `z = -dev`) and how many of them are not `:below`,
# and the window of the plane that holds the image of the cells (the first-order range of the foot
# around the foot of every centre, widened by the second-order remainder).
function former_rule_and_window(cells, geometry)
    r, dev, R2 = cells.r, cells.dev, cells.remainder
    old_removed = old_not_below = 0
    xlo = zlo = Inf
    xhi = zhi = -Inf
    for I in CartesianIndices(cells.class)
        pos = Tuple(I) .+ cells.kmin .- 1
        q = SVector{4, Float64}(MP.get_coord_by_pos(cells.q_grid, pos))
        joints = RP._joints(q, geometry)
        foot = joints[5]
        if foot[2] < -dev
            old_removed += 1
            cells.class[I] != 0x01 && (old_not_below += 1)
        end
        # A rotation `dθ` about the joint `a` moves the foot `p` by `dθ ⋅ (-(p_z - a_z), p_x - a_x)`.
        half_x = sum(i -> abs(joints[i][2] - foot[2]) * r[i], 1:4)
        half_z = sum(i -> abs(foot[1] - joints[i][1]) * r[i], 1:4)
        xlo, xhi = min(xlo, foot[1] - half_x - R2), max(xhi, foot[1] + half_x + R2)
        zlo, zhi = min(zlo, foot[2] - half_z - R2), max(zhi, foot[2] + half_z + R2)
    end
    return (; old_removed, old_not_below, window = (xlo, xhi, zlo, zhi))
end

# The largest amount by which the foot height leaves the range of its cell, on the 16 corners and
# `npoints` random points of `ncells` random cells (zero, up to rounding, for the certified bounds).
function check_range(cells, geometry; ncells = 300, npoints = 200)
    rng = MersenneTwister(1)
    r = cells.r
    corners = [r .* SVector(s) for s in Iterators.product(ntuple(_ -> (-1.0, 1.0), 4)...)]
    worst_exit = 0.0
    for _ in 1:ncells
        pos = ntuple(a -> cells.kmin[a] + rand(rng, 0:(size(cells.class, a) - 1)), 4)
        q = SVector{4, Float64}(MP.get_coord_by_pos(cells.q_grid, pos))
        lo, hi = cell_height_range(
            q,
            r,
            geometry,
            cells.foot_bound,
            cells.remainder,
            cells.third,
            cells.dev,
        )
        points = [
            vec(corners);
            [r .* SVector(ntuple(_ -> 2rand(rng) - 1, 4)) for _ in 1:npoints]
        ]
        for δ in points
            z = RP.swing_foot_position(q + δ, geometry)[2]
            worst_exit = max(worst_exit, z - hi, lo - z)
        end
    end
    return (; worst_exit)
end
