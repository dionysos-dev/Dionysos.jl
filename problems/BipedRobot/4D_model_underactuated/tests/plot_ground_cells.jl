# Where does the ground constraint of `RP.infeasible_cells` remove cells, seen from the swing foot?
#
# The carving decides on the 4-D angle cells `(q₁, q₂, q₃, q₄)` (the swing foot does not depend on
# the ankle momentum). This script runs the ground test of `infeasible_cells` on every cell of a 4-D
# box, but removes nothing: a cell that would be removed is mapped to a shape in the Cartesian plane of
# the swing foot (`x` forwards, `z` up, the stance foot at the origin) and drawn there, filled in red
# with some transparency. The shapes of different cells overlap, so the opacity of a region tells how
# many removed cells map to it.
#
# The shape of a cell is the square the carving itself reasons with: the swing foot of the cell
# centre, plus or minus the Lipschitz margin `dev` (`RP.swing_foot_deviation_bound`) along `x` and
# `z`, which encloses the foot of every configuration of the cell. A cell is removed when the foot of
# its centre is below `z = -dev`, so every shape lies entirely below the ground, `z = 0`. Another
# shape (e.g. the exact image of the cell) is a change of `foot_shape` only.
#
# The two parameters are the cell size and the state bounds, in 4-D. Set `CELL_SIZE` and
# `STATE_BOUNDS = (q_lb, q_ub)` before including this file to change them (or edit the defaults).
#
# Run it with an environment that has Dionysos, Symbolics, Plots, LazySets, StaticArrays and
# MathOptInterface:
#   julia --project=<env> tests/plot_ground_cells.jl

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

# Optional: the opacity of one shape (`nothing`: chosen so that the most overlapped region is nearly
# opaque, but not below `MIN_ALPHA`) and the file of the plot (next to this script by default).
fill_alpha = @isdefined(FILL_ALPHA) ? FILL_ALPHA : nothing
# Below this opacity the 8-bit blending of the image can no longer tell overlaps apart.
const MIN_ALPHA = 0.02
output_file =
    @isdefined(OUTPUT_FILE) ? OUTPUT_FILE :
    joinpath(@__DIR__, "ground_cells_foot_plane.png")

# ------------------------------------------------------------
# The ground test of `RP.infeasible_cells`, without the removal
# ------------------------------------------------------------

# The shape drawn for a removed cell whose centre has the swing foot at `foot = (x, z)`: the square
# of half-width `dev` around it, as the vertex coordinates `(xs, zs)`.
function foot_shape(foot, dev)
    return (
        [foot[1] - dev, foot[1] + dev, foot[1] + dev, foot[1] - dev],
        [foot[2] - dev, foot[2] - dev, foot[2] + dev, foot[2] + dev],
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

# Run the ground test of `RP.infeasible_cells` (same margin, same cells, same test) on every cell of
# the box `[q_lb, q_ub]`. Returns the number of cells examined, the positions of the cells that
# `infeasible_cells` would remove, the swing foot of their centres, and the margin `dev`.
function ground_cells(q_lb, q_ub, cell_size, geometry)
    # `infeasible_cells`: `dev = swing_foot_deviation_bound(geometry, h ./ 2)`, with `h` the 5-D cell
    # size; the ankle has no weight in it (its gradient bound is 0), hence the 0.0 in first place.
    dev = RP.swing_foot_deviation_bound(geometry, SVector(0.0, (cell_size ./ 2)...))
    # `infeasible_cells`: the 4-D grid of the angles and its box `Q`.
    q_grid = MP.GridFree(zero(cell_size), cell_size)
    Q = LazySets.Hyperrectangle(; low = q_lb, high = q_ub)
    q_foot(q) = RP.swing_foot_position(SVector{4}(q), geometry)

    ncells = 0
    positions = NTuple{4, Int}[]
    feet = SVector{2, Float64}[]
    # `MP.cells_where` evaluates its predicate on these cells (`INNER` by default).
    for pos in MP.get_pos_from_set(q_grid, Q, MP.INNER)
        ncells += 1
        foot = q_foot(MP.get_coord_by_pos(q_grid, pos))
        # The ground test of `infeasible_cells`: the foot of the centre is surely below the ground.
        foot[2] < -dev || continue
        # `infeasible_cells` would remove this cell: here it is only recorded.
        push!(positions, NTuple{4, Int}(pos))
        push!(feet, foot)
    end
    return (; ncells, positions, feet, dev)
end

# Largest number of squares of half-width `dev`, centred on `feet`, that cover a common point,
# measured on a raster of step `dev / 10` (a difference array, so linear in the number of squares).
function max_overlap(feet, dev; step = dev / 10)
    side = round(Int, 2 * dev / step)
    xmin, zmin = minimum(first, feet) - dev, minimum(last, feet) - dev
    nx = floor(Int, (maximum(first, feet) - minimum(first, feet)) / step) + side + 2
    nz = floor(Int, (maximum(last, feet) - minimum(last, feet)) / step) + side + 2
    diff = zeros(Int, nx, nz)
    for foot in feet
        i = floor(Int, (foot[1] - dev - xmin) / step) + 1
        j = floor(Int, (foot[2] - dev - zmin) / step) + 1
        diff[i, j] += 1
        diff[i + side, j] -= 1
        diff[i, j + side] -= 1
        diff[i + side, j + side] += 1
    end
    return maximum(cumsum(cumsum(diff; dims = 1); dims = 2))
end

# ------------------------------------------------------------
# Run, check against `infeasible_cells`, plot
# ------------------------------------------------------------

println("Building the symbolic model and its kinematics…")
geometry = RP.default_geometry()

removed = ground_cells(q_lb, q_ub, cell_size, geometry)
println(
    "Ground constraint: ",
    length(removed.positions),
    " of ",
    removed.ncells,
    " cells removed (margin dev = ",
    round(100 * removed.dev; digits = 2),
    " cm)",
)

# The copy above must agree with the real thing: compare with the `"ground"` cells of
# `RP.infeasible_cells` on the same grid (its ankle axis is one cell, which the ground test ignores).
let h = SVector(cell_size[1], cell_size...)
    state_grid = MP.GridFree(zero(h), h)
    X = LazySets.Hyperrectangle(;
        low = SVector(-h[1], q_lb...),
        high = SVector(h[1], q_ub...),
    )
    carving = RP.infeasible_cells(X, state_grid, nothing, geometry)
    same = carving.causes[1].second == Set(removed.positions)
    same ? println("Same cells as RP.infeasible_cells: true") :
    @warn "the removed cells differ from the \"ground\" cells of RP.infeasible_cells"
end

if isempty(removed.positions)
    println("No cell is removed: nothing to plot (a finer grid or a wider box?).")
else
    shapes = [foot_shape(foot, removed.dev) for foot in removed.feet]
    n_max = max_overlap(removed.feet, removed.dev)
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
        label = "removed cell: foot of the centre ± $(round(100 * removed.dev; digits = 1)) cm " *
                "(opacity $(round(α; sigdigits = 2)) each)",
    )
    Plots.hline!(fig, [0.0]; color = :black, label = "ground, z = 0")
    Plots.hline!(
        fig,
        [-removed.dev];
        color = :gray,
        linestyle = :dash,
        label = "z = -dev: removal threshold of the centre",
    )
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
