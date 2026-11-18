# Where does the ground constraint of `RP.infeasible_cells` remove cells, seen from the swing foot?
#
# The carving decides on the 4-D angle cells `(q₁, q₂, q₃, q₄)` (the swing foot does not depend on
# the ankle momentum). This script runs `RP.infeasible_cells` on a 4-D box (with no obstacle) and maps
# every cell of its `"ground"` constraint to a shape in the Cartesian plane of the swing foot (`x`
# forwards, `z` up, the stance foot at the origin), filled in red with some transparency. The shapes
# of different cells overlap, so the opacity of a region tells how many removed cells map to it.
#
# The shape of a cell encloses the foot of every configuration of the cell: along `z`, the range
# `[z_min, z_max]` of the foot height over the cell that `infeasible_cells` uses
# (`RP.swing_foot_height_range`, with the same bound); along `x`, the foot of the cell centre plus or
# minus the Lipschitz margin `dev` (`RP.swing_foot_deviation_bound`). With a certified bound and the
# ambiguous cells kept, a cell is removed when `z_max < 0`, so every shape lies entirely below the
# ground, `z = 0`. Another shape (e.g. the exact image of the cell) is a change of `foot_shape` only.
#
# Set `CELL_SIZE`, `STATE_BOUNDS = (q_lb, q_ub)`, `FOOT_BOUND` (one of `RP.FOOT_BOUNDS`) and
# `GROUND_AMBIGUOUS` (`:keep` or `:remove`) before including this file to change them (or edit the
# defaults).
#
#   julia --project=problems/BipedRobot tests/plot_ground_cells.jl

using StaticArrays
import Dionysos
import LazySets
import Plots

const MP = Dionysos.Mapping

include(joinpath(@__DIR__, "..", "robot_problem.jl"))
import .UnderactuatedRobotProblem as RP

# ------------------------------------------------------------
# Parameters (4-D: q₁, q₂, q₃, q₄, in radians)
# ------------------------------------------------------------

# Size of the cells of the grid, centred on the origin: cell `pos` is centred on `cell_size .* pos`.
cell_size = @isdefined(CELL_SIZE) ? CELL_SIZE : SVector(0.1, 0.1, 0.1, 0.1)

# Lower and upper bounds of the state box. Like the abstraction, only the cells lying entirely inside
# it are examined. The default is the box of `examples/BipedRobot/biped_4d_normal_form.jl`, widened so
# that the ground constraint acts on a grid of 0.1 rad (its margin is then 9.35 cm). The number of
# overlaps grows quickly with the box: about 2000 for `[-1, 1]⁴`.
q_lb, q_ub =
    @isdefined(STATE_BOUNDS) ? STATE_BOUNDS :
    (SVector(-0.5, 0.0, -0.8, -1.0), SVector(0.5, 1.0, 0.8, 0.0))

# The options of `RP.infeasible_cells`: how the range of the foot height over a cell is bounded, and
# what to do with the cells across the ground.
foot_bound = @isdefined(FOOT_BOUND) ? FOOT_BOUND : :second_order_cell
ground_ambiguous = @isdefined(GROUND_AMBIGUOUS) ? GROUND_AMBIGUOUS : :keep

# Optional: the opacity of one shape (`nothing`: chosen so that the most overlapped region is nearly
# opaque, but not below `MIN_ALPHA`) and the file of the plot (next to this script by default).
fill_alpha = @isdefined(FILL_ALPHA) ? FILL_ALPHA : nothing
# Below this opacity the 8-bit blending of the image can no longer tell overlaps apart.
const MIN_ALPHA = 0.02
output_file =
    @isdefined(OUTPUT_FILE) ? OUTPUT_FILE :
    joinpath(@__DIR__, "ground_cells_foot_plane.png")

# ------------------------------------------------------------
# The ground constraint of `RP.infeasible_cells`
# ------------------------------------------------------------

# The shape drawn for a removed cell, from the swing foot `foot = (x, z)` of its centre, the range
# `(z_min, z_max)` of the foot height over it and the margin `dev`, as the vertex coordinates
# `(xs, zs)`.
function foot_shape(foot, range, dev)
    return (
        [foot[1] - dev, foot[1] + dev, foot[1] + dev, foot[1] - dev],
        [range[1], range[1], range[2], range[2]],
    )
end

# All the shapes as one pair of vectors, `NaN` between two shapes: Plots draws them as one series
# (one legend entry), each polygon filled on its own, so that overlaps accumulate opacity.
function nan_separated(shapes)
    xs, zs = Float64[], Float64[]
    for (x, z) in shapes
        append!(xs, x)
        push!(xs, NaN)
        append!(zs, z)
        push!(zs, NaN)
    end
    return xs, zs
end

# The `"ground"` cells of `RP.infeasible_cells` on the box `[q_lb, q_ub]` (on a grid with one cell
# along the ankle axis, which the ground constraint ignores). Returns the number of cells examined,
# the positions of the removed cells, the swing foot of their centres, the range of the foot height
# over them, and the margin `dev`.
function ground_cells(q_lb, q_ub, cell_size, geometry; foot_bound, ground_ambiguous)
    r = cell_size ./ 2
    # The ankle has no weight in the margin (its gradient bound is 0), hence the 0.0 in first place.
    dev = RP.swing_foot_deviation_bound(geometry, SVector(0.0, r...))
    q_grid = MP.GridFree(zero(cell_size), cell_size)
    Q = LazySets.Hyperrectangle(; low = q_lb, high = q_ub)
    ncells = length(MP.get_pos_from_set(q_grid, Q, MP.INNER))
    h = SVector(cell_size[1], cell_size...)
    X = LazySets.Hyperrectangle(;
        low = SVector(-h[1], q_lb...),
        high = SVector(h[1], q_ub...),
    )
    carving = RP.infeasible_cells(
        X,
        MP.GridFree(zero(h), h),
        nothing,
        geometry;
        foot_bound,
        ground_ambiguous,
    )
    positions = sort!(collect(carving.causes[1].second))
    centres = [SVector{4}(MP.get_coord_by_pos(q_grid, pos)) for pos in positions]
    feet = [RP.swing_foot_position(q, geometry) for q in centres]
    ranges = [RP.swing_foot_height_range(q, r, geometry; foot_bound) for q in centres]
    return (; ncells, positions, feet, ranges, dev)
end

# Largest number of rectangles `[xlo, xhi] × [zlo, zhi]` that cover a common point, measured on a
# raster of step `step` (a difference array, so linear in the number of rectangles).
function max_overlap(rectangles, step)
    xmin = minimum(b -> b[1], rectangles)
    zmin = minimum(b -> b[3], rectangles)
    index(v, vmin) = floor(Int, (v - vmin) / step) + 1
    nx = index(maximum(b -> b[2], rectangles), xmin) + 1
    nz = index(maximum(b -> b[4], rectangles), zmin) + 1
    diff = zeros(Int, nx, nz)
    for (xlo, xhi, zlo, zhi) in rectangles
        i0, i1 = index(xlo, xmin), index(xhi, xmin) + 1
        j0, j1 = index(zlo, zmin), index(zhi, zmin) + 1
        diff[i0, j0] += 1
        diff[i1, j0] -= 1
        diff[i0, j1] -= 1
        diff[i1, j1] += 1
    end
    return maximum(cumsum(cumsum(diff; dims = 1); dims = 2))
end

# ------------------------------------------------------------
# Run, plot
# ------------------------------------------------------------

println("Building the symbolic model and its kinematics…")
geometry = RP.default_geometry()

removed = ground_cells(q_lb, q_ub, cell_size, geometry; foot_bound, ground_ambiguous)
println(
    "Ground constraint (foot_bound = $(repr(foot_bound)), ground_ambiguous = ",
    "$(repr(ground_ambiguous))): ",
    length(removed.positions),
    " of ",
    removed.ncells,
    " cells removed",
)

if isempty(removed.positions)
    println("No cell is removed: nothing to plot (a finer grid or a wider box?).")
else
    shapes = [foot_shape(f, z, removed.dev) for (f, z) in zip(removed.feet, removed.ranges)]
    n_max = max_overlap(
        [
            (f[1] - removed.dev, f[1] + removed.dev, z[1], z[2]) for
            (f, z) in zip(removed.feet, removed.ranges)
        ],
        removed.dev / 10,
    )
    # `n` overlapping shapes of opacity α give 1 - (1 - α)ⁿ: choose α so that `n_max` give 90%.
    α = fill_alpha === nothing ? max(1 - 0.1^(1 / n_max), MIN_ALPHA) : fill_alpha
    # The number of overlaps from which the region looks 90% opaque.
    n_sat = round(Int, log(0.1) / log(1 - α))
    println(
        "Shapes overlap up to ≈ ",
        n_max,
        " times; opacity of one shape: ",
        round(α; sigdigits = 3),
        " (90% opaque from ≈ ",
        n_sat,
        " overlaps)",
    )

    fig = Plots.plot(;
        aspect_ratio = :equal,
        size = (900, 760),
        legend = :outerbottom,
        xlabel = "swing foot x [m]",
        ylabel = "swing foot z [m]",
        title = "Cells removed by the ground constraint, in the swing-foot plane\n" *
                "$(length(removed.positions)) of $(removed.ncells) cells, cell size " *
                "$(round.(cell_size; sigdigits = 3)) rad, up to ≈ $n_max overlaps",
        titlefontsize = 10,
    )
    xs, zs = nan_separated(shapes)
    Plots.plot!(
        fig,
        xs,
        zs;
        seriestype = :shape,
        color = :red,
        fillalpha = α,
        linealpha = 0,
        linewidth = 0,
        label = "removed cell: x of the centre ± $(round(100 * removed.dev; digits = 1)) cm, " *
                "z over its range (opacity $(round(α; sigdigits = 2)) each)",
    )
    Plots.hline!(fig, [0.0]; color = :black, label = "ground, z = 0")
    Plots.scatter!(
        fig,
        [0.0],
        [0.0];
        color = :black,
        markershape = :diamond,
        label = "stance foot",
    )

    Plots.savefig(fig, output_file)
    println("Plot saved to ", output_file)
    isinteractive() && display(fig)
end
