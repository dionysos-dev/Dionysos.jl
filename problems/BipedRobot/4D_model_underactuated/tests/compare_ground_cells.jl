# Comparison of the runs of `run_ground_cells_matrix.sh` (`ground_cells_linearized.jl` on several cell
# sizes and state boxes, with five ranges of the foot height).
#
# Reads the results files of `ground_cells_plots/` and writes, in the same folder (or `OUTPUT_DIR`):
#
# * `<configuration>_grid.png`: one figure per configuration, three rows (the cells below the ground,
#   across it, above it) by five columns (the current rule, the first-order range, the first order plus
#   the worst-case second-order remainder, the per-cell second order, the corner values): for every point
#   of the swing-foot plane, the number of cells of the class whose pre-image contains it. Pixels with
#   exactly one cell are grey, the outline of each region is black, and the colour scale is shared along
#   a row, so the methods can be compared by eye;
# * `old_rule.png`: the current rule alone, removed against admitted cells, for the first configuration;
# * bar charts of the numbers: `bars_classes.png` (the share of removed / ambiguous / kept cells of
#   every method, for every configuration), `bars_removed.png` and `bars_ambiguous.png` (the removed
#   and the ambiguous shares), each for the first three boxes and, with the suffix `_knee`, for the knee
#   boxes, `bars_time.png` (the cost of the classification) and, if the results of
#   `check_ground_cells_extremes.jl` are in the folder, `bars_reliability.png` and `bars_truth.png`
#   (the classes against the true range of the foot height);
# * `summary_table.tex`, `timing_table.tex` and `overview.txt`: the shares of cells in each class and the
#   times (LaTeX tables for the presentation, and a text overview).
#
# Run it with an environment that has Plots (and a run of the matrix before):
#   julia --project=<env> tests/compare_ground_cells.jl

using Printf
using Serialization
import Plots
import Plots: grid  # used by `Plots.@layout`

dir = @isdefined(RESULTS_DIR) ? RESULTS_DIR : joinpath(@__DIR__, "ground_cells_plots")
out_dir = @isdefined(OUTPUT_DIR) ? OUTPUT_DIR : dir

# The state boxes of the runs (the suffix of the configuration name), with a name (in text and in
# LaTeX) and the bounds. The step box is the one of the forward-step example `biped_4d_normal_form.jl`, widened
# so that the ground constraint acts on a grid of 0.1 rad; then two cubes, and three knee boxes, in which the
# knees bend backwards (`q₂ ≥ 0`, `q₄ ≤ 0`).
BOXES = [
    (
        key = "box",
        name = "step box",
        tex_name = "step box",
        bounds = "q₁∈[-0.5,0.5], q₂∈[0,1], q₃∈[-0.8,0.8], q₄∈[-1,0]",
    ),
    (
        key = "b15",
        name = "cube [-1.5, 1.5]⁴",
        tex_name = "cube \\([-1.5,1.5]^4\\)",
        bounds = "q₁, q₂, q₃, q₄ ∈ [-1.5,1.5]",
    ),
    (
        key = "b20",
        name = "cube [-2, 2]⁴",
        tex_name = "cube \\([-2,2]^4\\)",
        bounds = "q₁, q₂, q₃, q₄ ∈ [-2,2]",
    ),
    (
        key = "k10",
        name = "knee box 1",
        tex_name = "knee box 1",
        bounds = "q₁, q₃ ∈ [-1,1], q₂ ∈ [0,1], q₄ ∈ [-1,0]",
    ),
    (
        key = "k15",
        name = "knee box 1.5",
        tex_name = "knee box 1.5",
        bounds = "q₁, q₃ ∈ [-1.5,1.5], q₂ ∈ [0,1.5], q₄ ∈ [-1.5,0]",
    ),
    (
        key = "k20",
        name = "knee box 2",
        tex_name = "knee box 2",
        bounds = "q₁, q₃ ∈ [-2,2], q₂ ∈ [0,2], q₄ ∈ [-2,0]",
    ),
]
# The boxes of the bar charts, in groups of three panels, and the suffix of the file of each group.
BOX_SETS = [("", BOXES[1:3]), ("_knee", BOXES[4:6])]
# The cell sizes (the prefix of the configuration name), in rad.
CELLS = [("c010", 0.1), ("c005", 0.05)]
# Configuration name and its description (the box first, as in the slides).
CONFIGS = [("$(c)_$(b.key)", "cell $h rad, $(b.name)") for b in BOXES for (c, h) in CELLS]
box_of(name) = BOXES[findfirst(b -> endswith(name, b.key), BOXES)]
cell_of(name) = CELLS[findfirst(c -> startswith(name, c[1]), CELLS)][2]
# Method name, its label, a label on two lines (for the plots), and its colour in the bar charts.
METHODS = [
    ("margin", "current rule", "current\nrule", Plots.RGB(0.55, 0.55, 0.55)),
    ("first_order", "first order", "first\norder", Plots.RGB(0.72, 0.64, 0.85)),
    (
        "second_order",
        "second order, global",
        "2nd order\nglobal",
        Plots.RGB(0.50, 0.75, 0.48),
    ),
    (
        "second_order_cell",
        "second order, per cell",
        "2nd order\nper cell",
        Plots.RGB(0.11, 0.47, 0.22),
    ),
    ("corners", "corners", "corner\nvalues", Plots.RGB(0.46, 0.16, 0.51)),
]
# Class name, what happens to it, and its colormap (the rows of the grids).
CLASSES = [
    ("below", "removed", :Reds),
    ("across", "ambiguous", :Oranges),
    ("above", "kept", :Blues),
]

# ------------------------------------------------------------
# Panels
# ------------------------------------------------------------

# A colormap that starts at a visible tint instead of white.
function tinted(name)
    base = Plots.cgrad(name)
    return Plots.cgrad([base[t] for t in range(0.15, 1.0; length = 32)])
end

# The outline of the pixels of `mask` (`mask[i, j]` for the pixel centred on `(xs[i], zs[j])`): the
# edges between a pixel of the mask and one outside it, as `NaN`-separated segments `(x, z)`.
function outline_segments(mask, xs, zs)
    dx, dz = (xs[2] - xs[1]) / 2, (zs[2] - zs[1]) / 2
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

# One panel of counts: grey where exactly one cell is crossed, the colormap (limits `clims`) beyond, the
# outline of the region in black, the ground and the stance foot. `compact` drops the colorbar.
function panel_of_counts(xs, zs, counts, colormap, clims; compact = true, kwargs...)
    more = [n >= 2 ? Float64(n) : NaN for n in counts]
    single = [n == 1 ? 1.0 : NaN for n in counts]
    xlims = (xs[1] - (xs[2] - xs[1]) / 2, xs[end] + (xs[2] - xs[1]) / 2)
    zlims = (zs[1] - (zs[2] - zs[1]) / 2, zs[end] + (zs[2] - zs[1]) / 2)
    panel = Plots.heatmap(
        xs,
        zs,
        permutedims(more);
        aspect_ratio = :equal,
        color = tinted(colormap),
        clims,
        xlims,
        ylims = zlims,
        colorbar = !compact,
        kwargs...,
    )
    Plots.heatmap!(
        panel,
        xs,
        zs,
        permutedims(single);
        color = Plots.cgrad([:gray70, :gray70]),
        clims,
        colorbar_entry = false,
    )
    outline_x, outline_z = outline_segments(counts .>= 1, xs, zs)
    isempty(outline_x) || Plots.plot!(
        panel,
        outline_x,
        outline_z;
        color = :black,
        linewidth = 1,
        label = "",
        colorbar_entry = false,
    )
    Plots.hline!(panel, [0.0]; color = :black, linewidth = 0.8, label = "")
    Plots.scatter!(
        panel,
        [0.0],
        [0.0];
        color = :black,
        markershape = :diamond,
        markersize = 3,
        label = "",
    )
    return panel
end

# Round tick positions every `step` between `lo` and `hi` (a step of 0.5 m for windows over 1 m wide).
ticks(lo, hi, step) = collect(range(ceil(lo / step) * step, floor(hi / step) * step; step))
tick_step(lo, hi) = hi - lo > 1.0 ? 0.5 : 0.2

# A number of cells for a panel title: exact up to 99 999, then in thousands or millions.
function compact_count(n)
    n < 10^5 && return string(n)
    n < 10^6 && return "$(round(Int, n / 1e3)) k"
    return "$(round(n / 1e6; digits = 2)) M"
end

# The grid of a configuration: rows = classes, columns = methods, one colour scale per row.
function grid_figure(results)
    xs, zs = results["margin"].xs, results["margin"].zs
    xticks = ticks(xs[1], xs[end], tick_step(xs[1], xs[end]))
    yticks = ticks(zs[1], zs[end], tick_step(zs[1], zs[end]))
    row_max = [
        maximum(maximum(results[m].counts[:, :, row]) for (m, _) in METHODS) for row in 1:3
    ]
    panels = []
    for (row, (class, role, colormap)) in enumerate(CLASSES),
        (col, (method, _, plot_label, _)) in enumerate(METHODS)

        r = results[method]
        n = r.class_counts[row]
        title =
            (row == 1 ? "$plot_label\n" : "") *
            "$(compact_count(n)) ($(round(100 * n / r.ncells; digits = 1))%)"
        push!(
            panels,
            panel_of_counts(
                xs,
                zs,
                r.counts[:, :, row],
                colormap,
                (2, max(row_max[row], 3));
                title,
                titlefontsize = 22,
                tickfontsize = 18,
                guidefontsize = 20,
                ylabel = col == 1 ? "$class\n($role)" : "",
                xlabel = row == 3 ? "x [m]" : "",
                xticks,
                yticks,
                xformatter = row == 3 ? :auto : (_ -> ""),
                yformatter = col == 1 ? :auto : (_ -> ""),
            ),
        )
    end
    return Plots.plot(
        panels...;
        layout = (3, length(METHODS)),
        size = (2000, 960),
        left_margin = 13 * Plots.mm,
        bottom_margin = 8 * Plots.mm,
    )
end

# The current rule alone: the cells it removes against the ones it admits (the `:across` and `:above`
# cells together), with a colorbar each.
function old_rule_figure(results)
    r = results["margin"]
    below = r.counts[:, :, 1]
    admitted = r.counts[:, :, 2] .+ r.counts[:, :, 3]
    n_removed = r.class_counts[1]
    n_admitted = r.ncells - n_removed
    common = (;
        compact = false,
        titlefontsize = 15,
        tickfontsize = 12,
        guidefontsize = 14,
        xlabel = "swing foot x [m]",
        ylabel = "swing foot z [m]",
        colorbar_title = "cells crossed (grey: 1)",
    )
    return Plots.plot(
        panel_of_counts(
            r.xs,
            r.zs,
            below,
            :Reds,
            (2, max(maximum(below), 3));
            title = "Removed by the current rule: $n_removed cells ($(round(100 * n_removed / r.ncells; digits = 1))%)",
            common...,
        ),
        panel_of_counts(
            r.xs,
            r.zs,
            admitted,
            :Blues,
            (2, max(maximum(admitted), 3));
            title = "Admitted: $n_admitted cells ($(round(100 * n_admitted / r.ncells; digits = 1))%)",
            common...,
        );
        layout = (1, 2),
        size = (1500, 560),
        left_margin = 8 * Plots.mm,
        bottom_margin = 6 * Plots.mm,
    )
end

# ------------------------------------------------------------
# Bar charts
# ------------------------------------------------------------

# The title of a panel of a state box: its name, and the bounds on a second line.
box_title(box) = "$(box.name)\n$(box.bounds)"

const CLASS_COLORS = [
    Plots.RGB(190 / 255, 35 / 255, 35 / 255),
    Plots.RGB(215 / 255, 105 / 255, 10 / 255),
    Plots.RGB(30 / 255, 100 / 255, 170 / 255),
]
const TRUTH_COLOR = Plots.RGB(0.15, 0.15, 0.15)
share(r, k) = 100 * r.class_counts[k] / r.ncells
fmt1(v) = string(round(v; digits = 1))

# A strip with the legend of the methods (and of anything else, `extra`), to put above the panels.
function legend_strip(;
    extra = [],
    columns = length(METHODS) + length(extra),
    fontsize = 15,
)
    strip = Plots.plot(;
        framestyle = :none,
        legend = :top,
        legend_columns = columns,
        legendfontsize = fontsize,
        grid = false,
        xlims = (0, 1),
        ylims = (0, 1),
    )
    for (_, label, _, color) in METHODS
        Plots.bar!(strip, Float64[], Float64[]; color, label, linecolor = :black)
    end
    for (label, color) in extra
        Plots.bar!(strip, Float64[], Float64[]; color, label, linecolor = :black)
    end
    return strip
end

# Bars grouped by `groups` (the x positions `1:length(groups)`) with one bar per series of `values` (a
# vector of vectors, one per series), in the colours `colors`; the value is written above every bar.
function grouped_bars(
    values,
    colors,
    groups;
    ylims,
    title = "",
    ylabel = "",
    digits = 1,
    labelsize = 13,
    width = 0.8,
    rotate = false,
    kwargs...,
)
    nseries = length(values)
    w = width / nseries
    panel = Plots.plot(;
        xticks = (1:length(groups), groups),
        xlims = (0.4, length(groups) + 0.6),
        ylims,
        title,
        ylabel,
        legend = false,
        titlefontsize = 15,
        tickfontsize = 15,
        guidefontsize = 16,
        kwargs...,
    )
    labels = Tuple{Float64, Float64, Plots.PlotText}[]
    for (s, (vals, color)) in enumerate(zip(values, colors))
        xs = (1:length(groups)) .+ (s - (nseries + 1) / 2) * w
        Plots.bar!(
            panel,
            xs,
            vals;
            color,
            bar_width = 0.92 * w,
            linecolor = :black,
            label = "",
        )
        for (x, v) in zip(xs, vals)
            str = digits == 0 ? string(round(Int, v)) : string(round(v; digits))
            offset = 0.012 * (ylims[2] - ylims[1])
            text =
                rotate ? Plots.text(str, labelsize, :black, :left, :center, 90.0) :
                Plots.text(str, labelsize, :black, v >= 0 ? :bottom : :top)
            push!(labels, (x, v >= 0 ? v + offset : v - offset, text))
        end
    end
    Plots.annotate!(panel, labels)
    return panel
end

# The share of the removed / ambiguous / kept cells of every method (stacked), for every configuration
# (rows: cell sizes, columns: state boxes), with the shares written on the bars.
function classes_figure(all_results, boxes)
    panels = []
    for (row, (c, h)) in enumerate(CELLS), (col, box) in enumerate(boxes)
        res = get(all_results, "$(c)_$(box.key)", nothing)
        panel = Plots.plot(;
            xticks = (1:length(METHODS), [m[3] for m in METHODS]),
            xlims = (0.4, length(METHODS) + 0.6),
            ylims = (0, 100),
            ylabel = col == 1 ? "cell $h rad\nshare of the cells [%]" : "",
            title = row == 1 ? box_title(box) : "",
            legend = false,
            titlefontsize = 15,
            tickfontsize = 14,
            guidefontsize = 16,
            top_margin = row == 1 ? 9 * Plots.mm : 2 * Plots.mm,
            yformatter = col == 1 ? :auto : (_ -> ""),
            xformatter = :auto,
        )
        res === nothing && (push!(panels, panel); continue)
        xs = 1:length(METHODS)
        below = [share(res[m[1]], 1) for m in METHODS]
        across = [share(res[m[1]], 2) for m in METHODS]
        above = [share(res[m[1]], 3) for m in METHODS]
        Plots.bar!(
            panel,
            xs,
            above .+ across .+ below;
            color = CLASS_COLORS[3],
            bar_width = 0.78,
            linecolor = :white,
            label = "",
        )
        Plots.bar!(
            panel,
            xs,
            across .+ below;
            color = CLASS_COLORS[2],
            bar_width = 0.78,
            linecolor = :white,
            label = "",
        )
        Plots.bar!(
            panel,
            xs,
            below;
            color = CLASS_COLORS[1],
            bar_width = 0.78,
            linecolor = :white,
            label = "",
        )
        labels = Tuple{Float64, Float64, Plots.PlotText}[]
        for i in xs
            below[i] >= 3.5 &&
                push!(labels, (i, below[i] / 2, Plots.text(fmt1(below[i]), 13, :white)))
            below[i] < 3.5 && push!(
                labels,
                (i, below[i] + 4, Plots.text(fmt1(below[i]), 13, CLASS_COLORS[1])),
            )
            across[i] >= 3.5 && push!(
                labels,
                (i, below[i] + across[i] / 2, Plots.text(fmt1(across[i]), 13, :white)),
            )
            above[i] >= 3.5 && push!(
                labels,
                (
                    i,
                    below[i] + across[i] + above[i] / 2,
                    Plots.text(fmt1(above[i]), 13, :white),
                ),
            )
        end
        Plots.annotate!(panel, labels)
        push!(panels, panel)
    end
    strip = Plots.plot(;
        framestyle = :none,
        legend = :top,
        legend_columns = 3,
        legendfontsize = 17,
        grid = false,
        xlims = (0, 1),
        ylims = (0, 1),
    )
    for (k, (name, role, _)) in enumerate(CLASSES)
        Plots.bar!(
            strip,
            Float64[],
            Float64[];
            color = CLASS_COLORS[k],
            label = "$name ($role)",
            linecolor = :black,
        )
    end
    return Plots.plot(
        strip,
        panels...;
        layout = Plots.@layout([a{0.07h}; grid(2, 3)]),
        size = (2000, 1000),
        left_margin = 10 * Plots.mm,
        bottom_margin = 6 * Plots.mm,
    )
end

# One value per method, as grouped bars: groups = cell sizes, one panel per state box, with the numbers
# on the bars. `value(result)` is the plotted value, `ylabel` its label.
function methods_figure(all_results, boxes, value, ylabel; ymax)
    panels = []
    for (col, box) in enumerate(boxes)
        groups = ["cell $h rad" for (_, h) in CELLS]
        values = [
            [
                haskey(all_results, "$(c)_$(box.key)") ?
                value(all_results["$(c)_$(box.key)"][m[1]]) : 0.0 for (c, _) in CELLS
            ] for m in METHODS
        ]
        push!(
            panels,
            grouped_bars(
                values,
                [m[4] for m in METHODS],
                groups;
                ylims = (0, ymax),
                title = box_title(box),
                labelsize = 12,
                ylabel = col == 1 ? ylabel : "",
                yformatter = col == 1 ? :auto : (_ -> ""),
            ),
        )
    end
    return Plots.plot(
        legend_strip(),
        panels...;
        layout = Plots.@layout([a{0.08h}; grid(1, 3)]),
        size = (2000, 850),
        left_margin = 10 * Plots.mm,
        bottom_margin = 6 * Plots.mm,
    )
end

# ------------------------------------------------------------
# Tables
# ------------------------------------------------------------

# An integer with thin spaces between the thousands, for LaTeX.
function latex_int(n)
    digits = string(n)
    parts = [digits[max(1, i - 2):i] for i in reverse(length(digits):-3:1)]
    return join(parts, "\\,")
end
percent(n, total) = @sprintf("%.1f", 100 * n / total)

# The time of the classification of a run (the first call, without its compilation), in seconds.
function classification_time(r)
    stage = r.stage_times[findfirst(
        s -> s.label == "classification of the cells",
        r.stage_times,
    )]
    return stage.time - stage.compile
end
stage_time(r, label) = r.stage_times[findfirst(s -> s.label == label, r.stage_times)].time

# The time of the classification per cell, in µs: the compiled repetition of the run when there is one
# (grids up to 3 million cells), else the first call without its compilation.
function time_per_cell(r)
    again = findfirst(s -> startswith(s.label, "classification again"), r.stage_times)
    t = again === nothing ? classification_time(r) : r.stage_times[again].time
    return 1e6 * t / r.ncells
end

function format_time(t)
    t < 0.01 && return @sprintf("%.1f ms", 1000 * t)
    t < 10 && return @sprintf("%.2f s", t)
    return @sprintf("%.0f s", t)
end

# The cost of the classification: per cell (µs, logarithmic; the bar is the median over the
# configurations, the dots are the configurations) against the cost of the abstraction of a cell
# (about 0.2 ms per (cell, input) pair, 27 inputs), and the time of the largest grid.
function time_figure(all_results)
    nm = length(METHODS)
    per_cell = [[time_per_cell(res[m[1]]) for res in values(all_results)] for m in METHODS]
    median_of(v) = sort(v)[(length(v) + 1) ÷ 2]
    ylo, yhi = 1e-2, 3e4
    left = Plots.plot(;
        xticks = (1:nm, [m[3] for m in METHODS]),
        xlims = (0.4, nm + 0.6),
        ylims = (ylo, yhi),
        yscale = :log10,
        yticks = 10.0 .^ (-2:4),
        ylabel = "classification time per cell [µs]",
        title = "classification (one thread)",
        legend = :topleft,
        titlefontsize = 17,
        tickfontsize = 15,
        guidefontsize = 16,
        legendfontsize = 14,
    )
    for (i, (_, _, _, color)) in enumerate(METHODS)
        Plots.bar!(
            left,
            [i],
            [median_of(per_cell[i])];
            fillto = ylo,
            color,
            bar_width = 0.7,
            linecolor = :black,
            label = "",
        )
        Plots.scatter!(
            left,
            fill(i, length(per_cell[i])) .+
            range(-0.18, 0.18; length = length(per_cell[i])),
            per_cell[i];
            color = :white,
            markerstrokecolor = :black,
            markersize = 5,
            label = i == 1 ? "one configuration" : "",
        )
        Plots.annotate!(
            left,
            [(
                i,
                3 * median_of(per_cell[i]),
                Plots.text(
                    @sprintf("%.2g µs", median_of(per_cell[i])),
                    14,
                    :black,
                    :bottom,
                ),
            )],
        )
    end
    Plots.hline!(
        left,
        [5000.0];
        color = CLASS_COLORS[1],
        linestyle = :dash,
        linewidth = 2,
        label = "abstraction of a cell (27 inputs): about 5 ms",
    )
    largest = argmax(name -> all_results[name]["margin"].ncells, collect(keys(all_results)))
    n = all_results[largest]["margin"].ncells
    seconds = [classification_time(all_results[largest][m[1]]) for m in METHODS]
    right = Plots.plot(;
        xticks = (1:nm, [m[3] for m in METHODS]),
        xlims = (0.4, nm + 0.6),
        ylims = (0, 1.25 * maximum(seconds)),
        ylabel = "classification time [s]",
        title = "the largest grid: $(latex_int(n) |> x -> replace(x, "\\," => " ")) cells",
        legend = false,
        titlefontsize = 17,
        tickfontsize = 15,
        guidefontsize = 16,
    )
    for (i, (_, _, _, color)) in enumerate(METHODS)
        Plots.bar!(
            right,
            [i],
            [seconds[i]];
            color,
            bar_width = 0.7,
            linecolor = :black,
            label = "",
        )
        Plots.annotate!(
            right,
            [(
                i,
                seconds[i] + 0.02 * maximum(seconds),
                Plots.text(format_time(seconds[i]), 14, :black, :bottom),
            )],
        )
    end
    return Plots.plot(
        left,
        right;
        layout = (1, 2),
        size = (2000, 850),
        left_margin = 12 * Plots.mm,
        bottom_margin = 8 * Plots.mm,
    )
end

# Reliability on the first configuration with the results of `check_ground_cells_extremes.jl`: the
# cells whose class is wrong (removed cells that reach the ground side `z ≥ 0`, kept cells that go below
# it) and the highest point of the removed cells.
function reliability_figure(extremes, description)
    res = extremes.results
    nm = length(METHODS)
    wrong_below = [res[m[1]].below_contradicted for m in METHODS]
    wrong_above = [res[m[1]].above_contradicted for m in METHODS]
    highest = [100 * res[m[1]].highest_below for m in METHODS]
    ymax = max(maximum(wrong_below), maximum(wrong_above), 1) * 1.2
    left = grouped_bars(
        [wrong_below, wrong_above],
        [CLASS_COLORS[1], CLASS_COLORS[3]],
        [m[3] for m in METHODS];
        ylims = (0, ymax),
        title = "cells with a wrong class ($description, $(latex_int(extremes.ncells) |> x -> replace(x, "\\," => " ")) cells)",
        ylabel = "number of cells",
        digits = 0,
        labelsize = 15,
        width = 0.7,
    )
    Plots.bar!(
        left,
        Float64[],
        Float64[];
        color = CLASS_COLORS[1],
        label = "removed cell that reaches z ≥ 0",
        linecolor = :black,
    )
    Plots.bar!(
        left,
        Float64[],
        Float64[];
        color = CLASS_COLORS[3],
        label = "kept cell that goes below z = 0",
        linecolor = :black,
    )
    Plots.plot!(left; legend = :topright, legendfontsize = 14)
    lo, hi = min(minimum(highest), 0.0) * 1.15, max(maximum(highest), 0.0) + 1.5
    right = Plots.plot(;
        xticks = (1:nm, [m[3] for m in METHODS]),
        xlims = (0.4, nm + 0.6),
        ylims = (lo, hi),
        ylabel = "highest point of a removed cell [cm]",
        title = "how far above the ground a removed cell reaches",
        legend = false,
        titlefontsize = 17,
        tickfontsize = 15,
        guidefontsize = 16,
    )
    for i in 1:nm
        Plots.bar!(
            right,
            [i],
            [highest[i]];
            color = highest[i] >= 0 ? CLASS_COLORS[1] : METHODS[i][4],
            bar_width = 0.7,
            linecolor = :black,
            label = "",
        )
        Plots.annotate!(
            right,
            [(
                i,
                highest[i] + (highest[i] >= 0 ? 0.25 : -0.25),
                Plots.text(
                    @sprintf("%.2f", highest[i]),
                    14,
                    :black,
                    highest[i] >= 0 ? :bottom : :top,
                ),
            )],
        )
    end
    Plots.hline!(right, [0.0]; color = :black, linewidth = 2, label = "")
    return Plots.plot(
        left,
        right;
        layout = (1, 2),
        size = (2000, 850),
        left_margin = 12 * Plots.mm,
        bottom_margin = 8 * Plots.mm,
    )
end

# How far the cells removed and the cells kept by every method are from the true range of the foot height
# (`check_ground_cells_extremes.jl`), in % of the cells, for the configurations that have it: the number
# of cells the method classifies as `:below` (`:above`) minus the number of cells that really are. A
# positive bar is a class that is too large (the first order removes cells that touch the feasible
# side), a negative bar are cells that could be removed (kept) and are left as ambiguous. The current
# rule (about 17 points below the true range) is left out.
function truth_figure(extremes_by_config)
    names = [n for (n, _) in CONFIGS if haskey(extremes_by_config, n)]
    groups = ["$(cell_of(n)) rad\n$(box_of(n).name)" for n in names]
    methods = METHODS[2:end]
    panels = []
    for (field, title) in (
        (:below, "removed cells: method minus true range"),
        (:above, "kept cells: method minus true range"),
    )
        values = [
            [
                100 * (
                    getfield(extremes_by_config[n].results[m[1]], field) -
                    getfield(extremes_by_config[n].true_range, field)
                ) / extremes_by_config[n].ncells for n in names
            ] for m in methods
        ]
        push!(
            panels,
            grouped_bars(
                values,
                [m[4] for m in methods],
                groups;
                ylims = field == :below ? (-3.3, 2.7) : (-10.5, 1.5),
                title,
                ylabel = "% of the cells",
                labelsize = 15,
                digits = 2,
                width = 0.84,
            ),
        )
    end
    strip = Plots.plot(;
        framestyle = :none,
        legend = :top,
        legend_columns = length(methods),
        legendfontsize = 15,
        grid = false,
        xlims = (0, 1),
        ylims = (0, 1),
    )
    for (_, label, _, color) in methods
        Plots.bar!(strip, Float64[], Float64[]; color, label, linecolor = :black)
    end
    for panel in panels
        Plots.hline!(panel, [0.0]; color = :black, linewidth = 1.5, label = "")
    end
    return Plots.plot(
        strip,
        panels...;
        layout = Plots.@layout([a{0.08h}; grid(1, 2)]),
        size = (2000, 850),
        left_margin = 12 * Plots.mm,
        bottom_margin = 8 * Plots.mm,
    )
end

# ------------------------------------------------------------
# Run
# ------------------------------------------------------------

all_results = Dict{String, Dict{String, Any}}()
for (name, _) in CONFIGS
    files = [joinpath(dir, "$(name)_$(m[1]).jls") for m in METHODS]
    if !all(isfile, files)
        @warn "results of $name missing: skipped" missing = filter(!isfile, files)
        continue
    end
    all_results[name] = Dict(m[1] => deserialize(f) for (m, f) in zip(METHODS, files))
end
isempty(all_results) &&
    error("no complete configuration in $dir: run run_ground_cells_matrix.sh")
extremes_by_config = Dict{String, Any}()
for (name, _) in CONFIGS
    file = joinpath(dir, "extremes_$name.jls")
    isfile(file) && (extremes_by_config[name] = deserialize(file))
end

mkpath(out_dir)
for (name, description) in CONFIGS
    haskey(all_results, name) || continue
    file = joinpath(out_dir, "$(name)_grid.png")
    Plots.savefig(grid_figure(all_results[name]), file)
    println("Grid of $name ($description) saved to $file")
end
first_name = first(name for (name, _) in CONFIGS if haskey(all_results, name))
Plots.savefig(old_rule_figure(all_results[first_name]), joinpath(out_dir, "old_rule.png"))
println("Old rule figure saved (", first_name, ")")

# The bar charts, for every group of boxes that has results.
for (suffix, boxes) in BOX_SETS
    any(haskey(all_results, "$(c)_$(b.key)") for b in boxes, (c, _) in CELLS) || continue
    # The scale of the charts of one value per method: a bit above the largest value.
    top(value) =
        5 * ceil(
            1.12 * maximum(
                value(all_results["$(c)_$(b.key)"][m[1]]) for
                b in boxes, (c, _) in CELLS, m in METHODS if
                haskey(all_results, "$(c)_$(b.key)")
            ) / 5,
        )
    Plots.savefig(
        classes_figure(all_results, boxes),
        joinpath(out_dir, "bars_classes$suffix.png"),
    )
    Plots.savefig(
        methods_figure(
            all_results,
            boxes,
            r -> share(r, 1),
            "removed cells [% of the cells]";
            ymax = top(r -> share(r, 1)),
        ),
        joinpath(out_dir, "bars_removed$suffix.png"),
    )
    Plots.savefig(
        methods_figure(
            all_results,
            boxes,
            r -> share(r, 2),
            "ambiguous cells [% of the cells]";
            ymax = top(r -> share(r, 2)),
        ),
        joinpath(out_dir, "bars_ambiguous$suffix.png"),
    )
end
Plots.savefig(time_figure(all_results), joinpath(out_dir, "bars_time.png"))
if haskey(extremes_by_config, first_name)
    Plots.savefig(
        reliability_figure(extremes_by_config[first_name], CONFIGS[1][2]),
        joinpath(out_dir, "bars_reliability.png"),
    )
end
length(extremes_by_config) > 0 &&
    Plots.savefig(truth_figure(extremes_by_config), joinpath(out_dir, "bars_truth.png"))
println("Bar charts saved")

# The shares of cells of each class, in % (removed and ambiguous, for every method).
open(joinpath(out_dir, "summary_table.tex"), "w") do io
    println(io, "\\begin{tabular}{lr|", "r"^length(METHODS), "|", "r"^length(METHODS), "}")
    println(io, "\\toprule")
    println(
        io,
        " & & \\multicolumn{$(length(METHODS))}{c|}{removed (\\Below)} & \\multicolumn{$(length(METHODS))}{c}{ambiguous (\\Across)} \\\\",
    )
    heads = join([replace(m[3], "\n" => " ") for m in METHODS], " & ")
    println(io, "configuration & cells & ", heads, " & ", heads, " \\\\")
    println(io, "\\midrule")
    for (name, description) in CONFIGS
        haskey(all_results, name) || continue
        res = all_results[name]
        n = res["margin"].ncells
        row = [
            "$(cell_of(name)) rad, $(box_of(name).tex_name)",
            latex_int(n),
            [percent(res[m[1]].class_counts[1], n) for m in METHODS]...,
            [percent(res[m[1]].class_counts[2], n) for m in METHODS]...,
        ]
        println(io, join(row, " & "), " \\\\")
    end
    println(io, "\\bottomrule")
    return println(io, "\\end{tabular}")
end

# The times of the classification of the cells.
open(joinpath(out_dir, "timing_table.tex"), "w") do io
    println(io, "\\begin{tabular}{lr", "r"^length(METHODS), "}")
    println(io, "\\toprule")
    println(
        io,
        "configuration & cells & ",
        join([replace(m[3], "\n" => " ") for m in METHODS], " & "),
        " \\\\",
    )
    println(io, "\\midrule")
    for (name, description) in CONFIGS
        haskey(all_results, name) || continue
        res = all_results[name]
        row = [
            "$(cell_of(name)) rad, $(box_of(name).tex_name)",
            latex_int(res["margin"].ncells),
            [format_time(classification_time(res[m[1]])) for m in METHODS]...,
        ]
        println(io, join(row, " & "), " \\\\")
    end
    println(io, "\\bottomrule")
    return println(io, "\\end{tabular}")
end

# A plain-text overview of everything, for the notes.
open(joinpath(out_dir, "overview.txt"), "w") do io
    for (name, description) in CONFIGS
        haskey(all_results, name) || continue
        res = all_results[name]
        println(
            io,
            "== $name: $description ($(box_of(name).bounds)), $(res["margin"].ncells) cells, dev = ",
            round(100 * res["margin"].dev; digits = 2),
            " cm, R = ",
            round(100 * res["margin"].remainder; digits = 2),
            " cm",
        )
        for (m, label, _, _) in METHODS
            r = res[m]
            c = r.class_counts
            println(
                io,
                rpad("  $label", 26),
                "below $(c[1]) ($(percent(c[1], r.ncells))%)  across $(c[2]) ($(percent(c[2], r.ncells))%)  above $(c[3]) ($(percent(c[3], r.ncells))%)",
                "  | mean half-width ",
                round(100 * r.mean_half; digits = 2),
                " cm, max ",
                round(100 * r.max_half; digits = 2),
                " cm | z leaves the range by up to ",
                round(100 * r.check.worst_exit; digits = 2),
                " cm | classification ",
                format_time(classification_time(r)),
                " (",
                @sprintf("%.3g", time_per_cell(r)),
                " µs/cell), pixel counts ",
                format_time(stage_time(r, "pixel counts (pre-image sweep)")),
                ", total ",
                format_time(r.total_time),
                " | max counts ",
                join([maximum(r.counts[:, :, k]) for k in 1:3], "/"),
            )
        end
        r = res["first_order"]
        println(
            io,
            "  current rule removed: $(r.old_removed); not :below under first order: $(r.old_not_below)",
        )
        if haskey(extremes_by_config, name)
            e = extremes_by_config[name]
            println(
                io,
                "  true range (multi-start): below $(e.true_range.below), above $(e.true_range.above), across $(e.true_range.across); contradicted classes: ",
                join(
                    [
                        "$label $(e.results[m].below_contradicted)/$(e.results[m].above_contradicted)"
                        for (m, label, _, _) in METHODS
                    ],
                    "; ",
                ),
            )
        end
    end
end
println(read(joinpath(out_dir, "overview.txt"), String))
