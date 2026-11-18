# For every point of the swing-foot plane: how many 4-D cells does its pre-image cross?
#
# The swing foot is a map from the four joint angles to the plane, so the pre-image of a point `p` is
# a 2-D surface in 4-D (not a finite set of points, and no pseudo-inverse gives it). It is known in
# closed form, though: the hip depends on `(q₁, q₂)` only, and the swing leg is a two-link inverse-
# kinematics problem (`RP._swing_leg_ik`), so the pre-image is parametrized by `(q₁, q₂)` and the
# sign of the knee angle `q₄`:
#
#     hip = hip(q₁, q₂),   d = p - hip,   cos q₄ = (|d|² - L₁² - L₂²) / (2 L₁ L₂)   (needs |cos q₄| ≤ 1),
#     q₄ = ±acos(cos q₄),   q₃ = atan(dₓ, -d_z) - atan(L₂ sin q₄, L₁ + L₂ cos q₄) - (q₁ + q₂).
#
# For each pixel `p` of the plane, this script sweeps `(q₁, q₂)` over the state domain, computes
# `(q₃, q₄)` for both knee branches, and collects the 4-D cells that the points `(q₁, q₂, q₃, q₄)`
# fall in (the sweep is in `foot_preimage_sweep.jl`). The cells are split by the `"ground"`
# constraint of `RP.infeasible_cells` (with the bound `FOOT_BOUND` and the choice `GROUND_AMBIGUOUS`
# for the cells across the ground, as in `plot_ground_cells.jl`):
#
# * left: the number of **admitted** cells crossed by the pre-image;
# * right: the number of cells that the ground constraint **removes** that it crosses, i.e. the exact
#   counterpart of the overlap plotted by `plot_ground_cells.jl`, which draws an enclosing rectangle
#   for each removed cell instead of its image.
#
# The count is exact up to the sweep: it is a lower bound that converges as `sweep_step` decreases
# (a cell that the pre-image only grazes needs a fine sweep), and the script re-counts a few pixels at
# half the step to show how far it has converged. The kinematics are the ones of `RP`: the
# pre-image points are checked to map back to the pixel with `RP.swing_foot_position`.
#
# The parameters are the cell size and the state bounds, in 4-D, and the options of the ground
# constraint, as in `plot_ground_cells.jl`. Set `CELL_SIZE`, `STATE_BOUNDS = (q_lb, q_ub)`,
# `FOOT_BOUND` and `GROUND_AMBIGUOUS` before including this file to change them (or edit the
# defaults). Angles must lie in `[-π, π]`.
#
# The cost is the number of pixels times the number of `(q₁, q₂)` samples: use threads (the defaults
# take about two minutes on 10 threads; a larger `pixel_size` or `sweep_step` makes it faster and
# coarser).
#   julia -t auto --project=problems/BipedRobot tests/plot_ground_preimage_counts.jl

using LinearAlgebra
using Random
using StaticArrays
import Dionysos
import LazySets
import Plots

const MP = Dionysos.Mapping

include(joinpath(@__DIR__, "..", "robot_problem.jl"))
import .UnderactuatedRobotProblem as RP

# The sweep of the pre-image of a point, and the cells it crosses.
include(joinpath(@__DIR__, "foot_preimage_sweep.jl"))

# ------------------------------------------------------------
# Parameters (4-D: q₁, q₂, q₃, q₄, in radians)
# ------------------------------------------------------------

# Size of the cells of the grid, centred on the origin: cell `pos` is centred on `cell_size .* pos`.
cell_size = @isdefined(CELL_SIZE) ? CELL_SIZE : SVector(0.1, 0.1, 0.1, 0.1)

# Lower and upper bounds of the state box. Like the abstraction, only the cells lying entirely inside
# it are examined (the domain is the union of those cells). The default is the box of
# `plot_ground_cells.jl`.
q_lb, q_ub =
    @isdefined(STATE_BOUNDS) ? STATE_BOUNDS :
    (SVector(-0.5, 0.0, -0.8, -1.0), SVector(0.5, 1.0, 0.8, 0.0))

# The options of `RP.infeasible_cells`: how the range of the foot height over a cell is bounded, and
# what to do with the cells across the ground.
foot_bound = @isdefined(FOOT_BOUND) ? FOOT_BOUND : :second_order_cell
ground_ambiguous = @isdefined(GROUND_AMBIGUOUS) ? GROUND_AMBIGUOUS : :keep

# Optional: the size of a pixel of the plane (m), the step of the `(q₁, q₂)` sweep (rad: the smaller,
# the closer the count to the exact one, the longer the run) and the file of the plot (next to this
# script by default).
pixel_size = @isdefined(PIXEL_SIZE) ? PIXEL_SIZE : 0.005
sweep_step = @isdefined(SWEEP_STEP) ? SWEEP_STEP : 0.001
output_file =
    @isdefined(OUTPUT_FILE) ? OUTPUT_FILE : joinpath(@__DIR__, "ground_preimage_counts.png")

all(abs.(q_lb) .<= π) && all(abs.(q_ub) .<= π) ||
    error("the angle bounds must lie in [-π, π]: got $q_lb, $q_ub")

# ------------------------------------------------------------
# The ground constraint of `RP.infeasible_cells`, on every cell of the box
# ------------------------------------------------------------

# The cells of the box (`INNER`, like `RP.infeasible_cells`) and the subset `removed` that the
# `"ground"` constraint of `RP.infeasible_cells` removes, on a grid with one cell along the ankle axis
# (which the ground constraint ignores).
function classify_cells(q_lb, q_ub, cell_size, geometry; foot_bound, ground_ambiguous)
    q_grid = MP.GridFree(zero(cell_size), cell_size)
    Q = LazySets.Hyperrectangle(; low = q_lb, high = q_ub)
    positions = [NTuple{4, Int}(pos) for pos in MP.get_pos_from_set(q_grid, Q, MP.INNER)]
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
    return (; positions, removed = carving.causes[1].second, q_grid)
end

# ------------------------------------------------------------
# Run, check, plot
# ------------------------------------------------------------

println("Building the symbolic model and its kinematics…")
geometry = RP.default_geometry()

cells = classify_cells(q_lb, q_ub, cell_size, geometry; foot_bound, ground_ambiguous)
# Class codes: 1 admitted, 2 removed.
class, kmin, kmax = cell_class_array(
    cells.positions,
    [pos in cells.removed ? 0x02 : 0x01 for pos in cells.positions],
)
println(
    "Domain: ",
    length(cells.positions),
    " cells, ",
    length(cells.removed),
    " removed by the ground constraint (foot_bound = $(repr(foot_bound)), ",
    "ground_ambiguous = $(repr(ground_ambiguous)))",
)

h = SVector(cell_size)
sweep = build_sweep(geometry, class, kmin, kmax, h, sweep_step)
println(
    "Sweep: ",
    size(sweep.hx, 1),
    " × ",
    size(sweep.hx, 2),
    " samples of (q₁, q₂) per knee branch (",
    length(sweep.branches),
    " branch(es)), step ≈ ",
    sweep_step,
    " rad",
)

# The pre-image points map back to the pixel under the kinematics of `RP` (random hips and points).
let rng = MersenneTwister(1), err = 0.0, n = 0
    lo = SVector((kmin[1] - 0.5) * h[1], (kmin[2] - 0.5) * h[2])
    hi = SVector((kmax[1] + 0.5) * h[1], (kmax[2] + 0.5) * h[2])
    for _ in 1:20_000
        q1, q2 = lo .+ rand(rng, 2) .* (hi .- lo)
        p = SVector(0.6 * (2rand(rng) - 1), 0.6 * (2rand(rng) - 1))
        for sgn in (1.0, -1.0)
            q = preimage_point(geometry, q1, q2, p, sgn)
            q === nothing && continue
            err = max(err, norm(RP.swing_foot_position(q, geometry) - p))
            n += 1
        end
    end
    println(
        "Pre-image points map back to the point: max error ",
        err,
        " over ",
        n,
        " points",
    )
end

xmin, xmax, zmin, zmax =
    foot_window(cells.positions, cells.q_grid, geometry, 3 * pixel_size)
xs = range(xmin + pixel_size / 2, xmax; step = pixel_size)
zs = range(zmin + pixel_size / 2, zmax; step = pixel_size)
println(
    "Plane: x ∈ [",
    round(xmin; digits = 3),
    ", ",
    round(xmax; digits = 3),
    "], z ∈ [",
    round(zmin; digits = 3),
    ", ",
    round(zmax; digits = 3),
    "], ",
    length(xs),
    " × ",
    length(zs),
    " pixels of ",
    pixel_size,
    " m, ",
    Threads.nthreads(),
    " thread(s)",
)
counts = @time count_plane(sweep, xs, zs, 2)
admitted, removed = counts[:, :, 1], counts[:, :, 2]

reached = count(>(0), admitted .+ removed)
println(
    "Pixels reached by the image of the domain: ",
    reached,
    " of ",
    length(admitted),
    "; admitted cells crossed: up to ",
    maximum(admitted),
    "; removed cells crossed: up to ",
    maximum(removed),
)
border =
    [admitted[1, :]; admitted[end, :]; admitted[:, 1]; admitted[:, end]] .+
    [removed[1, :]; removed[end, :]; removed[:, 1]; removed[:, end]]
any(>(0), border) &&
    @warn "the image of the domain touches the border of the window: it is cut"

# How far the sweep has converged: re-count a few pixels at half the step.
let rng = MersenneTwister(2)
    fine = build_sweep(geometry, class, kmin, kmax, h, sweep_step / 2)
    nonempty = findall(>(0), admitted .+ removed)
    seen = fill(false, length(class))
    touched = Int[]
    finer = zeros(Int, 2)
    println(
        "Convergence (step $sweep_step → $(sweep_step / 2)), admitted / removed cells crossed:",
    )
    for I in rand(rng, nonempty, 3)
        i, j = Tuple(I)
        preimage_counts!(finer, seen, touched, fine, xs[i], zs[j])
        println(
            "  p = (",
            round(xs[i]; digits = 3),
            ", ",
            round(zs[j]; digits = 3),
            "): ",
            (admitted[i, j], removed[i, j]),
            " → ",
            Tuple(finer),
        )
    end
end

if reached == 0
    println("The pre-images cross no cell: nothing to plot.")
else
    function count_panel(counts, label, color, title)
        data = [n == 0 ? NaN : Float64(n) for n in counts]
        panel = Plots.heatmap(
            xs,
            zs,
            permutedims(data);
            aspect_ratio = :equal,
            color,
            clims = (1, max(maximum(counts), 2)),
            xlims = (xmin, xmax),
            ylims = (zmin, zmax),
            xlabel = "swing foot x [m]",
            ylabel = "swing foot z [m]",
            colorbar_title = label,
            title,
            titlefontsize = 10,
        )
        Plots.hline!(panel, [0.0]; color = :black, label = "")
        Plots.scatter!(
            panel,
            [0.0],
            [0.0];
            color = :black,
            markershape = :diamond,
            label = "",
        )
        return panel
    end

    setup = "cell size $(round.(cell_size; sigdigits = 3)) rad, $(length(cells.positions)) cells"
    fig = Plots.plot(
        count_panel(
            admitted,
            "admitted cells crossed",
            :viridis,
            "Admitted cells crossed by the pre-image of each point\n$setup",
        ),
        count_panel(
            removed,
            "removed cells crossed",
            :Reds,
            "Removed cells crossed by the pre-image of each point\n$(length(cells.removed)) removed cells (ground constraint)",
        );
        layout = (1, 2),
        size = (1600, 620),
        left_margin = 12 * Plots.mm,
        bottom_margin = 6 * Plots.mm,
    )
    Plots.savefig(fig, output_file)
    println("Plot saved to ", output_file)
    isinteractive() && display(fig)
end
