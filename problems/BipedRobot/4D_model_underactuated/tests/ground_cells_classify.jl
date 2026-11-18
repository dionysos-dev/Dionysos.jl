# The ground constraint of `RP.infeasible_cells` on every cell of a 4-D box: how many cells are below,
# across and above the ground, and where they map in the plane of the swing foot.
#
# The cells are classified by `ground_cells_common.jl`, with the range of the foot height over a cell
# of `RP.infeasible_cells` (one of `RP.FOOT_BOUNDS`) or the former rule of `infeasible_cells` as a
# reference (`:margin`):
#
# * `:below`  (`z_max < 0`): the whole cell is below the ground: removed;
# * `:across` (`z_min < 0 ≤ z_max`): the cell straddles the ground, the "ambiguous" cells, kept by
#   `infeasible_cells` (`ground_ambiguous = :keep`, the default) or removed (`:remove`);
# * `:above`  (`z_min ≥ 0`): the whole cell is above the ground: kept.
#
# The script prints how many cells are in each class, compares with the former rule and with
# `RP.infeasible_cells` itself, checks the range by sampling points in random cells, and times every
# stage.
#
# The plot, in the style of `plot_ground_preimage_counts.jl`: for every pixel `p` of the swing-foot
# plane, the number of cells of each class that the pre-image of `p` crosses (the sweep is in
# `foot_preimage_sweep.jl`), in three panels: `:below`, `:across` and `:above`. Pixels crossed by
# exactly one cell are solid grey and the outline of each region is drawn. A pixel that is dark in the
# first panel is a foot position that many removed cells can produce.
#
# Set before including this file (or edit the defaults): `CELL_SIZE`, `STATE_BOUNDS = (q_lb, q_ub)`
# (angles in `[-π, π]`), `FOOT_BOUND` (one of `RP.FOOT_BOUNDS`, or `:margin`), `GROUND_AMBIGUOUS`
# (`:keep` or `:remove`), `PIXEL_SIZE` and `SWEEP_STEP` (the resolution of the plot: a smaller
# `SWEEP_STEP` makes the counts more exact and the run longer), `OUTPUT_FILE` (the plot) and
# `RESULTS_FILE` (optional: the counts and the timings, saved with `Serialization` for
# `compare_ground_cells.jl`). `run_ground_cells_matrix.sh` runs the configurations of the presentation.
#
#   julia -t auto --project=problems/BipedRobot tests/ground_cells_classify.jl

using Serialization
import Plots

include(joinpath(@__DIR__, "ground_cells_common.jl"))

# ------------------------------------------------------------
# Parameters (4-D: q₁, q₂, q₃, q₄, in radians)
# ------------------------------------------------------------

# Size of the cells of the grid, centred on the origin: cell `pos` is centred on `cell_size .* pos`.
cell_size = @isdefined(CELL_SIZE) ? CELL_SIZE : SVector(0.1, 0.1, 0.1, 0.1)

# Lower and upper bounds of the state box. Like the abstraction, only the cells lying entirely inside
# it are examined. The default is the box of `plot_ground_cells.jl` (the "step box").
q_lb, q_ub =
    @isdefined(STATE_BOUNDS) ? STATE_BOUNDS :
    (SVector(-0.5, 0.0, -0.8, -1.0), SVector(0.5, 1.0, 0.8, 0.0))

foot_bound = @isdefined(FOOT_BOUND) ? FOOT_BOUND : :second_order_cell
ground_ambiguous = @isdefined(GROUND_AMBIGUOUS) ? GROUND_AMBIGUOUS : :keep
check_study_bound(foot_bound)
ground_ambiguous in (:keep, :remove) ||
    error("`GROUND_AMBIGUOUS` must be `:keep` or `:remove`, got $(repr(ground_ambiguous))")
all(abs.(q_lb) .<= π) && all(abs.(q_ub) .<= π) ||
    error("the angle bounds must lie in [-π, π]: got $q_lb, $q_ub")
method_label = bound_label(foot_bound)

pixel_size = @isdefined(PIXEL_SIZE) ? PIXEL_SIZE : 0.005
sweep_step = @isdefined(SWEEP_STEP) ? SWEEP_STEP : 0.002
output_file =
    @isdefined(OUTPUT_FILE) ? OUTPUT_FILE :
    joinpath(@__DIR__, "ground_cells_$(foot_bound).png")
results_file = @isdefined(RESULTS_FILE) ? RESULTS_FILE : nothing

# ------------------------------------------------------------
# Classification, checks
# ------------------------------------------------------------

# The time of every stage of the run: `timed(f, label)` runs `f()` and records how long it takes (the
# first call of a function includes its compilation, which is recorded apart). With `repeat = n`, `f`
# runs `n` times and the mean is recorded. `compare_ground_cells.jl` reads the stages by label.
stage_times = NamedTuple{(:label, :time, :compile), Tuple{String, Float64, Float64}}[]
function timed(f, label; repeat = 1)
    stats = @timed begin
        result = nothing
        for _ in 1:repeat
            result = f()
        end
        result
    end
    push!(
        stage_times,
        (; label, time = stats.time / repeat, compile = stats.compile_time / repeat),
    )
    return stats.value
end

run_start = time()
println("Building the symbolic model and its kinematics…")
geometry = timed("symbolic model and kinematics") do
    return RP.default_geometry()
end

classify() = ground_cells_by_bound(q_lb, q_ub, cell_size, geometry; foot_bound)
cells = timed("classification of the cells") do
    return classify()
end
total = length(cells.class)
# Again, now that everything is compiled: the cost of the range of a cell, without the compilation
# (mean of several runs when it takes milliseconds, and not at all for very large grids).
nrepeat = total <= 50_000 ? 20 : total <= 500_000 ? 3 : 1
repeated = total <= 3_000_000
if repeated
    timed("classification again (compiled, mean of $nrepeat)"; repeat = nrepeat) do
        return classify()
    end
end
let first_call =
        stage_times[findfirst(s -> s.label == "classification of the cells", stage_times)]
    println(
        "Classification of $total cells, range: $method_label: ",
        round(first_call.time; sigdigits = 3),
        " s the first time, of which ",
        round(first_call.compile; sigdigits = 3),
        " s of compilation",
    )
    repeated && println(
        "  once compiled: ",
        round(1e3 * stage_times[end].time; sigdigits = 3),
        " ms, ",
        round(1e6 * stage_times[end].time / total; sigdigits = 3),
        " µs per cell",
    )
end
former = timed("former rule and plot window") do
    return former_rule_and_window(cells, geometry)
end
cm(x) = round(100 * x; digits = 2)

check = timed("check by sampling") do
    return check_range(cells, geometry)
end
println(
    "Range of the foot height over a cell, $method_label (random points of random cells): ",
    "z outside [z_min, z_max] by up to ",
    cm(check.worst_exit),
    " cm",
)
foot_bound in (:second_order_global, :second_order_cell) &&
    check.worst_exit > 1e-12 &&
    @warn "a certified range is left by the foot height"

n_below, n_across, n_above = cells.counts
# The classification must agree with the real thing: the `"ground"` cells of `RP.infeasible_cells` on
# the same grid (its ankle axis is one cell, which the ground test ignores), for both choices of the
# ambiguous cells (only for grids small enough for `infeasible_cells` to hold the cells in sets).
if foot_bound != :margin && total <= 3_000_000
    agree = timed("RP.infeasible_cells (same bound)") do
        h = SVector(cell_size[1], cell_size...)
        state_grid = MP.GridFree(zero(h), h)
        X = LazySets.Hyperrectangle(;
            low = SVector(-h[1], q_lb...),
            high = SVector(h[1], q_ub...),
        )
        removed(a) = length(
            RP.infeasible_cells(
                X,
                state_grid,
                nothing,
                geometry;
                foot_bound,
                ground_ambiguous = a,
            ).causes[1].second,
        )
        return removed(:keep) == n_below && removed(:remove) == n_below + n_across
    end
    agree ? println("Same removed cells as RP.infeasible_cells (:keep and :remove): true") :
    @warn "the classes differ from the \"ground\" cells of RP.infeasible_cells"
end

function report(label, n)
    println(
        "  ",
        rpad(label, 56),
        lpad(n, 9),
        "  (",
        lpad(round(100 * n / total; digits = 1), 5),
        " %)",
    )
    return nothing
end

println("\nGround constraint, range: ", method_label, ", ", total, " cells:")
report("below the ground (z_max < 0)", n_below)
report("across the ground (z_min < 0 ≤ z_max), ambiguous", n_across)
report("above the ground (z_min ≥ 0)", n_above)
println("\nCells removed:")
report("ground_ambiguous = :keep   (the :below cells)", n_below)
report("ground_ambiguous = :remove (the :below and :across cells)", n_below + n_across)
report("former rule of infeasible_cells (z_centre < -dev)", former.old_removed)
println(
    "\nMean half-width of the range of z over a cell: ",
    cm(cells.mean_half),
    " cm (largest ",
    cm(cells.max_half),
    " cm), against the margin dev = ",
    cm(cells.dev),
    " cm of the former rule",
)
println(
    "Cells removed by the former rule that are not :below: ",
    former.old_not_below,
    "; cells that are :below and not removed by the former rule: ",
    n_below - former.old_removed + former.old_not_below,
)

# ------------------------------------------------------------
# The plot: the cells of each class that the pre-image of every pixel crosses
# ------------------------------------------------------------

class, kmin, kmax = cells.class, cells.kmin, cells.kmax
h = SVector(cell_size)
sweep = timed("pre-image sweep set-up") do
    return build_sweep(geometry, class, kmin, kmax, h, sweep_step)
end
xmin, xmax, zmin, zmax = former.window .+ (-3, 3, -3, 3) .* pixel_size
xs = range(xmin + pixel_size / 2, xmax; step = pixel_size)
zs = range(zmin + pixel_size / 2, zmax; step = pixel_size)
println(
    "\nPlane: ",
    length(xs),
    " × ",
    length(zs),
    " pixels of ",
    pixel_size,
    " m, sweep of ",
    size(sweep.hx, 1),
    " × ",
    size(sweep.hx, 2),
    " samples of (q₁, q₂) per knee branch (",
    length(sweep.branches),
    " branch(es)), ",
    Threads.nthreads(),
    " thread(s)",
)
counts = timed("pixel counts (pre-image sweep)") do
    return count_plane(sweep, xs, zs, 3)
end
below, across, above = counts[:, :, 1], counts[:, :, 2], counts[:, :, 3]
total_counts = below .+ across .+ above
println(
    "Pixels reached by the image of the domain: ",
    count(>(0), total_counts),
    " of ",
    length(below),
    "; cells crossed up to ",
    maximum(below),
    " below, ",
    maximum(across),
    " across, ",
    maximum(above),
    " above",
)
any(
    >(0),
    [total_counts[1, :]; total_counts[end, :]; total_counts[:, 1]; total_counts[:, end]],
) && @warn "the image of the domain touches the border of the window: it is cut"

# How far the sweep has converged: re-count a few pixels at half the step.
convergence = timed("convergence check (half step)") do
    rng = MersenneTwister(2)
    fine = build_sweep(geometry, class, kmin, kmax, h, sweep_step / 2)
    seen = fill(false, length(class))
    touched = Int[]
    finer = zeros(Int, 3)
    println(
        "Convergence (step $sweep_step → $(sweep_step / 2)), below / across / above cells crossed:",
    )
    rows = NamedTuple[]
    for I in rand(rng, findall(>(0), total_counts), 3)
        i, j = Tuple(I)
        preimage_counts!(finer, seen, touched, fine, xs[i], zs[j])
        coarse = (below[i, j], across[i, j], above[i, j])
        println(
            "  p = (",
            round(xs[i]; digits = 3),
            ", ",
            round(zs[j]; digits = 3),
            "): ",
            coarse,
            " → ",
            Tuple(finer),
        )
        push!(rows, (; coarse, fine = Tuple(finer)))
    end
    return rows
end

# The outline of the pixels of `mask` (`mask[i, j]` for the pixel centred on `(xs[i], zs[j])`): the
# edges between a pixel of the mask and one outside it, as `NaN`-separated segments `(x, z)`.
function outline_segments(mask, xs, zs)
    dx, dz = step(xs) / 2, step(zs) / 2
    nx, nz = size(mask)
    inside(i, j) = 1 <= i <= nx && 1 <= j <= nz && mask[i, j]
    X, Z = Float64[], Float64[]
    function edge!(x0, x1, z0, z1)
        append!(X, (x0, x1, NaN))
        append!(Z, (z0, z1, NaN))
        return nothing
    end
    for j in 1:nz, i in 1:nx
        mask[i, j] || continue
        x, z = xs[i], zs[j]
        inside(i - 1, j) || edge!(x - dx, x - dx, z - dz, z + dz)
        inside(i + 1, j) || edge!(x + dx, x + dx, z - dz, z + dz)
        inside(i, j - 1) || edge!(x - dx, x + dx, z - dz, z - dz)
        inside(i, j + 1) || edge!(x - dx, x + dx, z + dz, z + dz)
    end
    return X, Z
end

# A colormap that starts at a visible tint instead of white.
function tinted(name)
    base = Plots.cgrad(name)
    return Plots.cgrad([base[t] for t in range(0.15, 1.0; length = 32)])
end

# One panel: the number of cells crossed at every pixel. Pixels crossed by exactly one cell are solid
# grey (they are the rim of the region), the others use the colormap, and the outline of the region
# (at least one cell) is drawn in black.
function count_panel(counts, colormap, title)
    more = [n >= 2 ? Float64(n) : NaN for n in counts]
    single = [n == 1 ? 1.0 : NaN for n in counts]
    panel = Plots.heatmap(
        xs,
        zs,
        permutedims(more);
        aspect_ratio = :equal,
        color = tinted(colormap),
        clims = (2, max(maximum(counts), 3)),
        xlims = (xmin, xmax),
        ylims = (zmin, zmax),
        xlabel = "swing foot x [m]",
        ylabel = "swing foot z [m]",
        colorbar_title = "cells crossed (grey: 1)",
        title,
        titlefontsize = 10,
    )
    Plots.heatmap!(
        panel,
        xs,
        zs,
        permutedims(single);
        color = Plots.cgrad([:gray70, :gray70]),
        clims = (0, 1),
        colorbar_entry = false,
    )
    X, Z = outline_segments(counts .>= 1, xs, zs)
    isempty(X) || Plots.plot!(panel, X, Z; color = :black, linewidth = 1, label = "")
    Plots.hline!(panel, [0.0]; color = :black, label = "")
    Plots.scatter!(panel, [0.0], [0.0]; color = :black, markershape = :diamond, label = "")
    return panel
end

setup = "cell size $(round.(cell_size; sigdigits = 3)) rad, $total cells"
fate = ground_ambiguous == :remove ? "removed" : "kept"
fig = timed("plot and save") do
    fig = Plots.plot(
        count_panel(below, :Reds, "Cells below the ground (removed): $n_below\n$setup"),
        count_panel(
            across,
            :Oranges,
            "Cells across the ground (ambiguous, $fate): $n_across\nrange: $method_label",
        ),
        count_panel(
            above,
            :Blues,
            "Cells above the ground (kept): $n_above\npre-image of each point",
        );
        layout = (1, 3),
        size = (2100, 480),
        left_margin = 12 * Plots.mm,
        bottom_margin = 8 * Plots.mm,
    )
    Plots.savefig(fig, output_file)
    return fig
end
println("Plot saved to ", output_file)
isinteractive() && display(fig)

println("\nExecution time (without loading the packages):")
for stage in stage_times
    println(
        "  ",
        rpad(stage.label, 46),
        lpad(round(stage.time; sigdigits = 3), 9),
        " s",
        stage.compile > 0.05 ?
        "  ($(round(stage.compile; sigdigits = 3)) s of compilation)" : "",
    )
end
total_time = time() - run_start
println("  ", rpad("total", 46), lpad(round(total_time; sigdigits = 3), 9), " s")

if results_file !== nothing
    Serialization.serialize(
        results_file,
        (;
            config = (;
                cell_size = collect(cell_size),
                q_lb = collect(q_lb),
                q_ub = collect(q_ub),
                foot_bound,
                ground_ambiguous,
                pixel_size,
                sweep_step,
                method_label,
            ),
            ncells = total,
            class_counts = (n_below, n_across, n_above),
            old_removed = former.old_removed,
            old_not_below = former.old_not_below,
            mean_half = cells.mean_half,
            max_half = cells.max_half,
            dev = cells.dev,
            remainder = cells.remainder,
            check,
            xs = collect(xs),
            zs = collect(zs),
            counts,
            convergence,
            stage_times = [(; s.label, s.time, s.compile) for s in stage_times],
            total_time,
            threads = Threads.nthreads(),
        ),
    )
    println("Results saved to ", results_file)
end
