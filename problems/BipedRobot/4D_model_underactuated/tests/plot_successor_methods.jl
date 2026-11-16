# Figures and tables of the presentation "Four ways to compute the successors" (the two decks
# `successor_methods_detailed.tex` and `successor_methods_visual.tex` in
# `Documentation/transitions_comparison/`), from the results of `compare_successor_methods.jl`.
#
# The study has eight variants per system: the global and the cell-wise growth bound, and the random and
# the uniform sample points, each for three values of N (`study.Ns`). In every figure the growth bounds
# are red and orange, the random points are blue and the uniform points green, in three shades each:
# the lightest for the fewest points, the darkest for the most.
#
# It writes, in `OUTPUT_DIR` (by default the `figures/` folder of the presentations; every file starts
# with `sm_`):
#
# * `sm_bars_<system>_{transitions,time}.png`: the four statistics (total and average number of
#   transitions, total and average computation time) of the eight variants for every sampling time and
#   cell size, as bars; `sm_ref_<system>.png`: the same for the reference configuration, one big bar per
#   variant; `sm_sens_<system>.png`: the effect of the sampling time and of the cell size;
# * `sm_hist_<system>.png`: the histogram of the number of successors of a (cell, input) pair, for the
#   reference configuration;
# * `sm_table_<system>_{transitions,pair,time,upair,ref}.tex`: the same numbers as tables;
# * `sm_overlay_simple.png`, `sm_overlay_double_{small,large}.png`: the successors of one cell and one input
#   with the four methods (at the reference N) on top of each other (the cells that each method lists, the
#   growth-bound boxes, the images of the sample points, the true image), for the cells where the
#   cell-wise bound helps the least and the most; `sm_overlayN_simple.png`, `sm_overlayN_double.png`: the
#   cells listed by the random and the uniform points for each N, against those of the cell-wise bound;
#   `sm_gain_simple.png` and `sm_gain_double.png`: where the cell-wise bound helps, over the state space
#   (simple pendulum) and as a distribution (double pendulum);
# * `sm_coverage.png` and `sm_coverage.tex`: how much of what a fresh point of a cell does each variant
#   covers, and how many cells it lists for each cell that is actually reached, against N.
#
# Run it after `compare_successor_methods.jl`, with the same environment:
#
#   julia --project=<env> tests/plot_successor_methods.jl

using Printf
using Random
using Serialization
using Statistics
import Plots
import Plots: grid  # used by `Plots.@layout`

# The definitions of the study (not its run).
include(joinpath(@__DIR__, "compare_successor_methods.jl"))

const OUT =
    @isdefined(OUTPUT_DIR) ? OUTPUT_DIR :
    joinpath(@__DIR__, "..", "Documentation", "transitions_comparison", "figures")
const RESULTS = joinpath(@__DIR__, "successor_methods_results", "results.jls")

# The configuration that is looked at in detail, for every system.
const REFERENCE = (
    simple = (tau = 0.1, size = [0.1, 0.2]),
    double = (tau = 0.01, size = [0.3, 0.3, 0.6, 0.6]),
)

const STUDIES = studies()
study_of(system) = only(s for s in STUDIES if s.name == system)
# A variant is `(method, N)`, with `N = 0` for the growth bounds; see `variants` of the study script.
variants_of(system) = variants(study_of(system))
# The four methods at the reference number of points.
function reference_of(system)
    N = study_of(system).N
    return [(:global, 0), (:cellwise, 0), (:random, N), (:uniform, N)]
end

# ------------------------------------------------------------------------------------------
# Colors and names
# ------------------------------------------------------------------------------------------

rgb(r, g, b) = Plots.RGB(r / 255, g / 255, b / 255)

# The colors of the deck (`preamble.tex`): red, orange, blue, green; the blue and the green of the deck
# are the middle shades. Three shades for the three N: the lightest for the fewest points.
const GROWTH_RGB = Dict(:global => (190, 35, 35), :cellwise => (215, 105, 10))
const SHADE_RGB = Dict(
    :random => [(125, 175, 220), (30, 100, 170), (10, 45, 100)],
    :uniform => [(120, 195, 140), (27, 120, 55), (8, 62, 28)],
)
const RANKS = ("low", "mid", "high")

rank_of(system, N) = findfirst(==(N), study_of(system).Ns)
shade(m, k) = rgb(SHADE_RGB[m][k]...)

function color_of(system, v)
    m, N = v
    haskey(GROWTH_RGB, m) && return rgb(GROWTH_RGB[m]...)
    return shade(m, rank_of(system, N))
end

const FAMILY = Dict(
    :global => "global growth bound",
    :cellwise => "cell-wise growth bound",
    :random => "random points",
    :uniform => "uniform points",
)

function short_of(v)
    m, N = v
    m === :global && return "global"
    m === :cellwise && return "cell-wise"
    return "$m $N"
end

function tick_of(v)
    m, N = v
    m === :global && return "global"
    m === :cellwise && return "cell-\nwise"
    return "N=$N"
end

function tex_color(system, v)
    m, N = v
    m === :global && return "cbelow"
    m === :cellwise && return "cacross"
    return (m === :random ? "crand" : "cunif") * RANKS[rank_of(system, N)]
end

function tex_label(system, v)
    m, N = v
    text = m in SAMPLED ? "$(FAMILY[m]), \$N=$N\$" : FAMILY[m]
    return "\\textcolor{$(tex_color(system, v))}{$text}"
end

const SYSTEM_TITLES = Dict(:simple => "simple pendulum", :double => "double pendulum")

Plots.default(;
    titlefontsize = 18,
    guidefontsize = 16,
    tickfontsize = 14,
    legendfontsize = 15,
    framestyle = :box,
    grid = true,
    gridalpha = 0.25,
    dpi = 100,
)

mkpath(OUT)

# ------------------------------------------------------------------------------------------
# The results
# ------------------------------------------------------------------------------------------

const DATA = deserialize(RESULTS)
const RECORDS = DATA.records
const INFO = DATA.info

records_of(system) = [r for r in RECORDS if r.system == system]

# The configurations of a system, in the order of the study.
function configurations(system)
    seen = Tuple{Float64, Vector{Float64}}[]
    for r in records_of(system)
        (r.tau, r.size) in seen || push!(seen, (r.tau, r.size))
    end
    return seen
end

record(system, v, tau, size) = only(
    r for r in RECORDS if r.system == system &&
        r.method == v[1] &&
        r.N == v[2] &&
        r.tau == tau &&
        r.size == size
)

fmt_g(x) = @sprintf("%g", x)

# The cell size of a configuration, in words and in a short label: the pendulum has (angle, velocity),
# the double pendulum (angle, angle, velocity, velocity) with equal angles and equal velocities.
size_label(system, size) =
    system == :simple ? "$(fmt_g(size[1]))/$(fmt_g(size[2]))" :
    "$(fmt_g(size[1]))/$(fmt_g(size[3]))"
config_label(system, tau, size) = "τ=$(fmt_g(tau))\n$(size_label(system, size))"

# ------------------------------------------------------------------------------------------
# Bars
# ------------------------------------------------------------------------------------------

# Ticks of a logarithmic axis that read as numbers (1, 2, 5, 10, ...), between `lo` and `hi`.
function log_ticks(lo, hi)
    ticks = Float64[]
    for e in floor(Int, log10(lo)):ceil(Int, log10(hi)), m in (1, 2, 5)
        t = m * 10.0^e
        lo <= t <= hi && push!(ticks, t)
    end
    return ticks
end

# The label of a tick: a count (1k, 10k, 1M), a time (1 ms, 0.1 s) or a plain number.
function tick_label(v, kind)
    short(x) = x == round(x) ? string(Int(round(x))) : @sprintf("%.3g", x)
    if kind == :count
        return v >= 1e6 ? short(v / 1e6) * "M" : v >= 1e3 ? short(v / 1e3) * "k" : short(v)
    elseif kind == :seconds
        return v >= 1 ? short(v) * " s" :
               v >= 1e-3 ? short(v * 1e3) * " ms" : short(v * 1e6) * " µs"
    end
    return short(v)
end

# Bars: one group per configuration (`configs`, labelled by `labelf`), one bar per variant. `values` maps
# every variant to its numbers, in the order of `configs`. `taus`: also write the sampling time above
# every run of equal τ.
function bar_panel(
    system,
    values,
    ylabel;
    configs = configurations(system),
    labelf = (t, s) -> size_label(system, s),
    logscale = true,
    kind = :plain,
    taus = false,
    title = "",
    annotate = nothing,
    ymax = nothing,
    xtickfontsize = nothing,
)
    vs = variants_of(system)
    nv = length(vs)
    n = length(configs)
    positive = [x for v in vs for x in values[v] if x > 0]
    base = logscale ? minimum(positive) / 2.2 : 0.0
    top =
        ymax !== nothing ? ymax : maximum(positive) * (logscale ? (taus ? 4.0 : 1.8) : 1.15)
    ticks = logscale ? log_ticks(base, top) : nothing
    p = Plots.plot(;
        title,
        ylabel,
        yscale = logscale ? :log10 : :identity,
        ylims = (base, top),
        yticks = logscale ? (ticks, [tick_label(t, kind) for t in ticks]) : :auto,
        xlims = (0.4, n + 0.6),
        xticks = (1:n, [labelf(t, s) for (t, s) in configs]),
        legend = false,
        xtickfontsize = xtickfontsize !== nothing ? xtickfontsize : (n > 5 ? 11 : 14),
        left_margin = 9Plots.mm,
        bottom_margin = 3Plots.mm,
    )
    width = 0.84 / nv
    for (m, v) in enumerate(vs), g in 1:n
        y = values[v][g]
        x0 = g + (m - (nv + 1) / 2) * width
        if y > 0
            Plots.plot!(
                p,
                Plots.Shape(
                    [x0 - width / 2, x0 + width / 2, x0 + width / 2, x0 - width / 2],
                    [base, base, y, y],
                );
                fillcolor = color_of(system, v),
                linecolor = :black,
                linewidth = 0.4,
                label = "",
            )
            annotate === nothing || Plots.annotate!(
                p,
                x0,
                y * (logscale ? 1.25 : 1.0) + (logscale ? 0.0 : 0.05 * top),
                Plots.text(annotate(y), 11, :black),
            )
        else
            Plots.annotate!(
                p,
                x0,
                base * (logscale ? 1.6 : 1.0),
                Plots.text("0", 10, color_of(system, v)),
            )
        end
    end
    if taus
        # The configurations come in runs of equal τ: separators, and the value above each run.
        ts = [t for (t, _) in configs]
        start = 1
        for g in 1:n
            if g == n || ts[g + 1] != ts[g]
                Plots.annotate!(
                    p,
                    (start + g) / 2,
                    top * (logscale ? 0.5 : 0.96),
                    Plots.text("τ = $(fmt_g(ts[g]))", 14, :black),
                )
                g < n && Plots.vline!(
                    p,
                    [g + 0.5];
                    linecolor = :gray,
                    linestyle = :dash,
                    linewidth = 1,
                    label = "",
                )
                start = g + 1
            end
        end
    end
    return p
end

# The variants in the order of the legend, which has four columns filled row by row: the growth bound
# and the three N of the random points on the first row, the cell-wise bound and the uniform points on
# the second.
function legend_order(system)
    vs = variants_of(system)
    return vcat(
        [(:global, 0)],
        [v for v in vs if v[1] === :random],
        [(:cellwise, 0)],
        [v for v in vs if v[1] === :uniform],
    )
end

# A strip with the legend of the variants of a system (a plot without axes), on two rows.
function legend_strip(system; fontsize = 14)
    p = Plots.plot(;
        framestyle = :none,
        legend = :top,
        legend_columns = 4,
        legendfontsize = fontsize,
        foreground_color_legend = nothing,
    )
    for v in legend_order(system)
        Plots.bar!(
            p,
            [NaN],
            [NaN];
            fillcolor = color_of(system, v),
            linecolor = :black,
            label = short_of(v),
        )
    end
    return p
end

function values_of(system, field; configs = configurations(system))
    return Dict(
        v => Float64[field(record(system, v, t, s)) for (t, s) in configs] for
        v in variants_of(system)
    )
end

# The four statistics, for the nine configurations of a system: one figure for the transitions
# (total, per pair) and one for the time, the two panels one above the other.
function bars_figure(system)
    vs = variants_of(system)
    transitions = values_of(system, r -> r.transitions)
    average = values_of(system, r -> r.transitions / r.pairs)
    total_time = values_of(system, r -> r.time)
    average_time = values_of(system, r -> 1e6 * r.time / r.pairs)
    spread =
        maximum(average[(:global, 0)]) /
        max(minimum(x for v in vs for x in average[v] if x > 0), 1e-9)
    for (name, panels) in (
        (
            "transitions",
            [
                bar_panel(
                    system,
                    transitions,
                    "total transitions";
                    kind = :count,
                    taus = true,
                ),
                bar_panel(
                    system,
                    average,
                    "transitions per pair";
                    logscale = spread > 12,
                    taus = true,
                ),
            ],
        ),
        (
            "time",
            [
                bar_panel(system, total_time, "total time"; kind = :seconds, taus = true),
                bar_panel(system, average_time, "time per pair (µs)"; taus = true),
            ],
        ),
    )
        fig = Plots.plot(
            legend_strip(system),
            panels...;
            layout = Plots.@layout([l{0.1h}; grid(2, 1)]),
            size = (1500, 960),
        )
        Plots.savefig(fig, joinpath(OUT, "sm_bars_$(system)_$name.png"))
    end
    return nothing
end

# The four statistics of the reference configuration, one big bar per variant.
function ref_figure(system)
    ref = REFERENCE[system]
    vs = variants_of(system)
    rs = Dict(v => record(system, v, ref.tau, ref.size) for v in vs)
    panels = []
    for (name, field, kind) in (
        ("total transitions", r -> r.transitions, :count),
        ("transitions per pair", r -> r.transitions / r.pairs, :plain),
        ("total time", r -> r.time, :seconds),
        ("time per pair (µs)", r -> 1e6 * r.time / r.pairs, :plain),
    )
        vals = [Float64(field(rs[v])) for v in vs]
        p = Plots.bar(
            1:length(vs),
            vals;
            fillcolor = [color_of(system, v) for v in vs],
            linecolor = :black,
            linewidth = 0.6,
            bar_width = 0.8,
            label = "",
            xticks = (1:length(vs), [tick_of(v) for v in vs]),
            ylims = (0, maximum(vals) * 1.22),
            yformatter = y -> tick_label(y, kind),
            title = name,
            left_margin = 4Plots.mm,
            bottom_margin = 7Plots.mm,
            xtickfontsize = 13,
        )
        for (i, y) in enumerate(vals)
            Plots.annotate!(
                p,
                i,
                y + 0.07 * maximum(vals),
                Plots.text(tick_label(y, kind), 12, :black),
            )
        end
        push!(panels, p)
    end
    fig = Plots.plot(
        legend_strip(system),
        panels...;
        layout = Plots.@layout([l{0.1h}; grid(2, 2)]),
        size = (1400, 820),
    )
    Plots.savefig(fig, joinpath(OUT, "sm_ref_$system.png"))
    return nothing
end

# The effect of the sampling time (at the reference cell size) and of the cell size (at the reference
# sampling time) on the two statistics per pair.
function sensitivity_figure(system)
    ref = REFERENCE[system]
    vs = variants_of(system)
    all = configurations(system)
    by_tau = [(t, s) for (t, s) in all if s == ref.size]
    by_size = [(t, s) for (t, s) in all if t == ref.tau]
    average(c) = values_of(system, r -> r.transitions / r.pairs; configs = c)
    micro(c) = values_of(system, r -> 1e6 * r.time / r.pairs; configs = c)
    spread(v) =
        maximum(maximum(v[m]) for m in vs) /
        max(minimum(x for m in vs for x in v[m] if x > 0), 1e-9)
    tau_label = (t, s) -> "τ=$(fmt_g(t))"
    size_label_ = (t, s) -> size_label(system, s)
    panels = [
        bar_panel(
            system,
            average(by_tau),
            "transitions per pair";
            configs = by_tau,
            labelf = tau_label,
            logscale = spread(average(by_tau)) > 12,
            title = "effect of τ",
        ),
        bar_panel(
            system,
            micro(by_tau),
            "time per pair (µs)";
            configs = by_tau,
            labelf = tau_label,
            title = "effect of τ",
        ),
        bar_panel(
            system,
            average(by_size),
            "transitions per pair";
            configs = by_size,
            labelf = size_label_,
            logscale = spread(average(by_size)) > 12,
            title = "effect of the cell size",
        ),
        bar_panel(
            system,
            micro(by_size),
            "time per pair (µs)";
            configs = by_size,
            labelf = size_label_,
            title = "effect of the cell size",
        ),
    ]
    fig = Plots.plot(
        legend_strip(system; fontsize = 13),
        panels...;
        layout = Plots.@layout([l{0.1h}; grid(2, 2)]),
        size = (1300, 800),
    )
    Plots.savefig(fig, joinpath(OUT, "sm_sens_$system.png"))
    return nothing
end

# ------------------------------------------------------------------------------------------
# Histograms of the number of successors of a pair
# ------------------------------------------------------------------------------------------

# Index `k + 1` of `hist`: the number of pairs with `k` successors. Bins: every `k` up to `exact`, then
# powers of two (a pair of the double pendulum can have hundreds).
function bins_of(hist; exact)
    edges = Tuple{Int, Int}[(k, k) for k in 0:min(exact, length(hist) - 1)]
    lo = exact + 1
    while lo <= length(hist) - 1
        hi = 2 * (lo - 1)                 # powers of two: 5–8, 9–16, 17–32, ...
        push!(edges, (lo, hi))
        lo = hi + 1
    end
    return edges
end

# One panel per method; the sampled ones have one bar per N in a bin, in the three shades.
function histogram_figure(system)
    ref = REFERENCE[system]
    exact = system == :simple ? 14 : 4
    study = study_of(system)
    vs = variants_of(system)
    all_bins = Dict(
        v => begin
            hist = record(system, v, ref.tau, ref.size).hist
            edges = bins_of(hist; exact)
            [sum(hist[(a + 1):min(b + 1, length(hist))]) for (a, b) in edges]
        end for v in vs
    )
    longest = maximum(length(all_bins[v]) for v in vs)
    edges = bins_of(zeros(Int, 2^(longest) + 4); exact)[1:longest]
    labels = [a == b ? string(a) : "$(a)–$(b)" for (a, b) in edges]
    pairs = record(system, (:global, 0), ref.tau, ref.size).pairs
    ymax = maximum(maximum(all_bins[v]) for v in vs) / pairs * 100 * 1.12
    families = [
        (:global, [(:global, 0)]),
        (:cellwise, [(:cellwise, 0)]),
        (:random, [(:random, N) for N in study.Ns]),
        (:uniform, [(:uniform, N) for N in study.Ns]),
    ]
    panels = []
    for (i, (family, members)) in enumerate(families)
        p = Plots.plot(;
            title = FAMILY[family],
            xticks = (1:longest, labels),
            xrotation = 45,
            xtickfontsize = 13,
            xlims = (0.3, longest + 0.7),
            ylims = (0, ymax),
            ylabel = i in (1, 3) ? "share of the pairs (%)" : "",
            xlabel = i in (3, 4) ? "number of successors" : "",
            legend = length(members) > 1 ? :topright : false,
            legendfontsize = 14,
            left_margin = 11Plots.mm,
            bottom_margin = 8Plots.mm,
        )
        k = length(members)
        width = 0.9 / k
        for (j, v) in enumerate(members)
            counts = vcat(all_bins[v], zeros(Int, longest - length(all_bins[v])))
            Plots.bar!(
                p,
                (1:longest) .+ (j - (k + 1) / 2) * width,
                counts ./ pairs .* 100;
                fillcolor = color_of(system, v),
                linecolor = :black,
                linewidth = 0.3,
                bar_width = width,
                label = k > 1 ? "N = $(v[2])" : "",
            )
        end
        push!(panels, p)
    end
    fig = Plots.plot(
        panels...;
        layout = (2, 2),
        size = (1240, 740),
        plot_title = "τ = $(fmt_g(ref.tau)), cells $(join(map(fmt_g, ref.size), "×"))",
        plot_titlefontsize = 18,
    )
    Plots.savefig(fig, joinpath(OUT, "sm_hist_$system.png"))
    return nothing
end

# ------------------------------------------------------------------------------------------
# Tables
# ------------------------------------------------------------------------------------------

function fmt_int(n)
    s = string(n)
    out = ""
    while length(s) > 3
        out = "\\,$(s[(end - 2):end])" * out
        s = s[1:(end - 3)]
    end
    return s * out
end

function fmt_time(t)
    t < 1e-3 && return @sprintf("%.2f ms", 1e3 * t)
    t < 0.1 && return @sprintf("%.1f ms", 1e3 * t)
    t < 10 && return @sprintf("%.2f s", t)
    return @sprintf("%.1f s", t)
end

fmt_avg(x) =
    x >= 100 ? @sprintf("%.0f", x) : x >= 10 ? @sprintf("%.1f", x) : @sprintf("%.2f", x)
fmt_us(x) =
    x >= 100 ? @sprintf("%.0f", x) : x >= 10 ? @sprintf("%.1f", x) : @sprintf("%.2f", x)

cell_text(system, size) =
    system == :simple ? "\$$(fmt_g(size[1]))\\times$(fmt_g(size[2]))\$" :
    "\$$(fmt_g(size[1]))^2\\times$(fmt_g(size[3]))^2\$"

# One statistic of the eight variants, one row per configuration: τ, cell size, number of cells, and a
# column per variant under its group (growth bounds, random points N, uniform points N).
function stat_table(system, f)
    study = study_of(system)
    vs = variants_of(system)
    nN = length(study.Ns)
    spec = "rrr" * "r"^length(vs)
    io = IOBuffer()
    println(io, "\\begin{tabular}{$spec}")
    println(io, "\\toprule")
    println(
        io,
        " & & & \\multicolumn{2}{c}{growth bounds} & \\multicolumn{$nN}{c}{\\textcolor{crandmid}{random points}, \$N\$} & \\multicolumn{$nN}{c}{\\textcolor{cunifmid}{uniform points}, \$N\$} \\\\",
    )
    println(
        io,
        "\\cmidrule(lr){4-5} \\cmidrule(lr){6-$(5 + nN)} \\cmidrule(lr){$(6 + nN)-$(5 + 2nN)}",
    )
    heads = join(
        (
            "\\textcolor{$(tex_color(system, v))}{$(v[1] in SAMPLED ? v[2] : v[1] === :global ? "global" : "cell-wise")}"
            for v in vs
        ),
        " & ",
    )
    println(io, "\$\\tau\$ & cell & cells & $heads \\\\")
    println(io, "\\midrule")
    last_tau = nothing
    for (tau, size) in configurations(system)
        last_tau !== nothing && tau != last_tau && println(io, "\\addlinespace")
        last_tau = tau
        rs = [record(system, v, tau, size) for v in vs]
        println(
            io,
            "$(fmt_g(tau)) & $(cell_text(system, size)) & $(fmt_int(rs[1].n_cells)) & $(join((f(r) for r in rs), " & ")) \\\\",
        )
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
    return String(take!(io))
end

function write_tables(system)
    for (name, f) in (
        ("transitions", r -> fmt_int(r.transitions)),
        (
            "pair",
            r -> @sprintf(
                "%s / %.0f\\%%",
                fmt_avg(r.transitions / r.pairs),
                100 * r.enabled_pairs / r.pairs
            ),
        ),
        ("time", r -> fmt_time(r.time)),
        ("upair", r -> fmt_us(1e6 * r.time / r.pairs)),
    )
        write(joinpath(OUT, "sm_table_$(system)_$name.tex"), stat_table(system, f))
    end
    # The reference configuration: one row per variant.
    ref = REFERENCE[system]
    io = IOBuffer()
    println(io, "\\begin{tabular}{lrrrrrr}")
    println(io, "\\toprule")
    println(
        io,
        "method & transitions & per pair & per enabled pair & enabled pairs & time & \$\\mu\$s per pair \\\\",
    )
    println(io, "\\midrule")
    for v in variants_of(system)
        r = record(system, v, ref.tau, ref.size)
        per_enabled = r.enabled_pairs == 0 ? "--" : fmt_avg(r.transitions / r.enabled_pairs)
        println(
            io,
            "$(tex_label(system, v)) & $(fmt_int(r.transitions)) & $(fmt_avg(r.transitions / r.pairs)) & $per_enabled & $(@sprintf("%.1f\\%%", 100 * r.enabled_pairs / r.pairs)) & $(fmt_time(r.time)) & $(fmt_us(1e6 * r.time / r.pairs)) \\\\",
        )
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
    write(joinpath(OUT, "sm_table_$(system)_ref.tex"), String(take!(io)))
    return nothing
end

# ------------------------------------------------------------------------------------------
# The abstractions of the reference configuration, to look at them
# ------------------------------------------------------------------------------------------

function reference_abstractions(system; domain = nothing)
    study = study_of(system)
    # A domain of its own, for looking at one cell: the pendulum of the study ends at θ = π, which
    # disables every input of the cells next to it.
    domain === nothing ||
        (study = merge(study, (; system = SimplePendulum.system(; _X_ = domain))))
    ref = REFERENCE[system]
    dynamics = SymbolicDynamics(study.system)
    grids = grids_of(study, ref.size)
    opts = Dict(
        (m, N) => build(m, study, dynamics, grids, ref.tau, N == 0 ? study.N : N) for
        (m, N) in variants(study)
    )
    return (; study, dynamics, grids, tau = ref.tau, opts)
end

# The set of successors of every pair that has one.
function successor_sets(opt)
    sets = Dict{Tuple{Int, Int}, Set{Int}}()
    for (target, source, input) in SY.enum_transitions(abstraction_of(opt))
        push!(get!(sets, (source, input), Set{Int}()), target)
    end
    return sets
end

# The box of a growth bound (the centre's image with its radius), before the cells are looked up.
function growth_box(study, jacobian_bound, tau, cell, u)
    approx = ST.discretize(ST.ContinuousTimeGrowthBound(study.system; jacobian_bound), tau;)
    return ST.get_over_approximation_map(approx)(cell, u)
end

# ------------------------------------------------------------------------------------------
# Drawing the successor sets
# ------------------------------------------------------------------------------------------

box_shape(lo, hi) =
    Plots.Shape([lo[1], hi[1], hi[1], lo[1], lo[1]], [lo[2], lo[2], hi[2], hi[2], lo[2]])

# One plane of a (possibly higher-dimensional) system: the coordinates `plane` of everything.
struct Scene
    plane::Tuple{Int, Int}
    source::Any                       # the cell
    cells::Dict{Tuple{Symbol, Int}, Vector{Any}}   # variant => the cells that it lists
    boxes::Dict{Symbol, Any}          # growth bounds => the box
    points::Dict{Tuple{Symbol, Int}, Vector{Any}}  # sampled variants => the images of the sample points
    truth::Vector{Any}                # images of many points of the cell
    outline::Union{Nothing, Vector{Any}}   # the image of the boundary of the cell, in order (2-D only)
end

with_plane(s::Scene, plane) =
    Scene(plane, s.source, s.cells, s.boxes, s.points, s.truth, nothing)

function scene_of(ctx, source_cell, input_index, q, u; n_truth = 3000, plane = (1, 2))
    sym = abstraction_of(ctx.opts[(:global, 0)])
    study, tau = ctx.study, ctx.tau
    system = study.system
    approx = SampledSimulation(system, tau; nsamples = 1)
    F = ST.get_system_map(approx)
    cells = Dict{Tuple{Symbol, Int}, Vector{Any}}()
    for v in variants(study)
        s = abstraction_of(ctx.opts[v])
        cells[v] = [SY.get_concrete_elem(s, t) for t in SY.post(s, q, input_index)]
    end
    boxes = Dict(
        :global => growth_box(study, study.jacobian_bound, tau, source_cell, u),
        :cellwise =>
            growth_box(study, cellwise_bound(ctx.dynamics), tau, source_cell, u),
    )
    points = Dict{Tuple{Symbol, Int}, Vector{Any}}()
    for (m, N) in variants(study)
        m in SAMPLED || continue
        a = SampledSimulation(system, tau; nsamples = N, sampling = m, seed = 1)
        image = ST.get_over_approximation_map(a)(source_cell, u)
        points[(m, N)] = collect(image.points)
    end
    c, ρ = LazySets.center(source_cell), LazySets.radius_hyperrectangle(source_cell)
    n = length(c)
    rng = MersenneTwister(7)
    truth = [F(c .+ ρ .* (2 .* rand(rng, n) .- 1), u) for _ in 1:n_truth]
    outline = nothing
    if n == 2
        edge = range(-1, 1; length = 120)
        boundary = vcat(
            [SVector(t, -1.0) for t in edge],
            [SVector(1.0, t) for t in edge],
            [SVector(t, 1.0) for t in reverse(edge)],
            [SVector(-1.0, t) for t in reverse(edge)],
        )
        outline = [F(c .+ ρ .* s, u) for s in boundary]
    end
    return Scene(plane, source_cell, cells, boxes, points, truth, outline)
end

# The projected squares of the cells of a method, each once.
function projected_cells(cells, plane)
    seen = Set{Tuple{Float64, Float64}}()
    out = []
    for cell in cells
        c, r = LazySets.center(cell), LazySets.radius_hyperrectangle(cell)
        key = (round(c[plane[1]]; digits = 9), round(c[plane[2]]; digits = 9))
        key in seen && continue
        push!(seen, key)
        push!(
            out,
            (
                SVector(c[plane[1]] - r[plane[1]], c[plane[2]] - r[plane[2]]),
                SVector(c[plane[1]] + r[plane[1]], c[plane[2]] + r[plane[2]]),
            ),
        )
    end
    return out
end

function scene_axes(scene::Scene; title, xlabel, ylabel, window, titlefontsize = 20)
    return Plots.plot(;
        title,
        xlabel,
        ylabel,
        xlims = (window[1], window[2]),
        ylims = (window[3], window[4]),
        legend = false,
        titlefontsize,
        left_margin = 11Plots.mm,
        bottom_margin = 6Plots.mm,
        top_margin = 4Plots.mm,
    )
end

# The true image of the cell: the image of its boundary (2-D), or of many of its points.
function draw_truth!(p, scene::Scene)
    i, j = scene.plane
    if scene.outline !== nothing
        Plots.plot!(
            p,
            Plots.Shape([q[1] for q in scene.outline], [q[2] for q in scene.outline]);
            fillcolor = :gray,
            fillalpha = 0.35,
            linecolor = :black,
            linewidth = 2.0,
            label = "",
        )
    else
        Plots.scatter!(
            p,
            [q[i] for q in scene.truth],
            [q[j] for q in scene.truth];
            markersize = 2.0,
            markerstrokewidth = 0,
            markercolor = :gray45,
            markeralpha = 0.5,
            label = "",
        )
    end
    return p
end

# The cell itself.
function draw_source!(p, scene::Scene)
    i, j = scene.plane
    c, r = LazySets.center(scene.source), LazySets.radius_hyperrectangle(scene.source)
    Plots.plot!(
        p,
        box_shape(SVector(c[i] - r[i], c[j] - r[j]), SVector(c[i] + r[i], c[j] + r[j]));
        fillalpha = 0,
        linecolor = :black,
        linewidth = 3.2,
        label = "",
    )
    return p
end

# The four methods at the reference N (`vs`, in the order global, cell-wise, random, uniform).
function draw_scene(system, scene::Scene, vs; title, xlabel, ylabel, window)
    i, j = scene.plane
    p = scene_axes(scene; title, xlabel, ylabel, window)
    draw_truth!(p, scene)
    # The cells that each method lists, the largest set first and with the thickest outline, so that
    # the nested sets show as stripes where their edges coincide.
    for (v, width) in zip(vs, (5.5, 4, 2.8, 1.6))
        for (lo, hi) in projected_cells(scene.cells[v], scene.plane)
            Plots.plot!(
                p,
                box_shape(lo, hi);
                fillcolor = color_of(system, v),
                fillalpha = 0.07,
                linecolor = color_of(system, v),
                linewidth = width,
                label = "",
            )
        end
    end
    # The growth-bound boxes (before the cells are looked up).
    for v in vs
        v[1] in (:global, :cellwise) || continue
        box = scene.boxes[v[1]]
        c, r = LazySets.center(box), LazySets.radius_hyperrectangle(box)
        Plots.plot!(
            p,
            box_shape(SVector(c[i] - r[i], c[j] - r[j]), SVector(c[i] + r[i], c[j] + r[j]));
            fillalpha = 0,
            linecolor = color_of(system, v),
            linestyle = :dash,
            linewidth = 2.5,
            label = "",
        )
    end
    # The images of the sample points.
    for (v, shape, size, stroke) in ((vs[3], :circle, 8, 1.5), (vs[4], :diamond, 9, 2))
        Plots.scatter!(
            p,
            [q[i] for q in scene.points[v]],
            [q[j] for q in scene.points[v]];
            markershape = shape,
            markersize = size,
            markercolor = color_of(system, v),
            markerstrokecolor = :white,
            markerstrokewidth = stroke,
            label = "",
        )
    end
    draw_source!(p, scene)
    return p
end

# One sampled variant against the cell-wise bound: the cells that it lists (filled, in its shade), the
# cells of the cell-wise bound (dashed orange outline) as the reference, and the images of its points.
function draw_scene_N(system, scene::Scene, v; title, xlabel, ylabel, window)
    i, j = scene.plane
    p = scene_axes(scene; title, xlabel, ylabel, window, titlefontsize = 17)
    draw_truth!(p, scene)
    for (lo, hi) in projected_cells(scene.cells[(:cellwise, 0)], scene.plane)
        Plots.plot!(
            p,
            box_shape(lo, hi);
            fillalpha = 0,
            linecolor = color_of(system, (:cellwise, 0)),
            linestyle = :dash,
            linewidth = 2.5,
            label = "",
        )
    end
    color = color_of(system, v)
    for (lo, hi) in projected_cells(scene.cells[v], scene.plane)
        Plots.plot!(
            p,
            box_shape(lo, hi);
            fillcolor = color,
            fillalpha = 0.3,
            linecolor = color,
            linewidth = 3,
            label = "",
        )
    end
    Plots.scatter!(
        p,
        [q[i] for q in scene.points[v]],
        [q[j] for q in scene.points[v]];
        markershape = v[1] === :random ? :circle : :diamond,
        markersize = 8,
        markercolor = color,
        markerstrokecolor = :white,
        markerstrokewidth = 1.5,
        label = "",
    )
    draw_source!(p, scene)
    return p
end

# The key to the drawings of the successor sets: a plot without axes.
function overlay_legend(system)
    vs = reference_of(system)
    p = Plots.plot(;
        framestyle = :none,
        legend = :top,
        legend_columns = 3,
        legendfontsize = 15,
        foreground_color_legend = nothing,
    )
    for v in vs
        name = v[1] in SAMPLED ? "$(FAMILY[v[1]]) (N=$(v[2]))" : FAMILY[v[1]]
        Plots.plot!(
            p,
            [NaN],
            [NaN];
            color = color_of(system, v),
            linewidth = 5,
            label = name,
        )
    end
    Plots.plot!(
        p,
        [NaN],
        [NaN];
        color = :black,
        linestyle = :dash,
        linewidth = 4,
        label = "growth-bound box",
    )
    Plots.scatter!(
        p,
        [NaN],
        [NaN];
        markershape = :circle,
        markersize = 9,
        markercolor = color_of(system, vs[3]),
        label = "image: random pt",
    )
    Plots.scatter!(
        p,
        [NaN],
        [NaN];
        markershape = :diamond,
        markersize = 9,
        markercolor = color_of(system, vs[4]),
        label = "image: uniform pt",
    )
    Plots.plot!(p, [NaN], [NaN]; color = :gray, linewidth = 6, label = "true image")
    Plots.plot!(p, [NaN], [NaN]; color = :black, linewidth = 5, label = "the cell")
    return p
end

# The key to the drawings of the cells for each N.
function overlay_legend_N(system)
    study = study_of(system)
    p = Plots.plot(;
        framestyle = :none,
        legend = :top,
        legend_columns = 3,
        legendfontsize = 14,
        foreground_color_legend = nothing,
    )
    Plots.plot!(
        p,
        [NaN],
        [NaN];
        color = color_of(system, (:cellwise, 0)),
        linestyle = :dash,
        linewidth = 4,
        label = "cells of the cell-wise bound",
    )
    Plots.plot!(
        p,
        [NaN],
        [NaN];
        color = shade(:random, 2),
        linewidth = 6,
        label = "cells of the random points",
    )
    Plots.plot!(
        p,
        [NaN],
        [NaN];
        color = shade(:uniform, 2),
        linewidth = 6,
        label = "cells of the uniform points",
    )
    Plots.plot!(p, [NaN], [NaN]; color = :gray, linewidth = 6, label = "true image")
    Plots.scatter!(
        p,
        [NaN],
        [NaN];
        markershape = :circle,
        markersize = 9,
        markercolor = shade(:random, 2),
        label = "image of a point",
    )
    Plots.plot!(p, [NaN], [NaN]; color = :black, linewidth = 5, label = "the cell")
    return p
end

# A window that contains everything that is drawn for the variants `vs`, with a margin.
function window_of(scene::Scene, vs)
    i, j = scene.plane
    xs = Float64[]
    ys = Float64[]
    for v in vs, cell in scene.cells[v]
        c, r = LazySets.center(cell), LazySets.radius_hyperrectangle(cell)
        push!(xs, c[i] - r[i], c[i] + r[i])
        push!(ys, c[j] - r[j], c[j] + r[j])
    end
    for q in scene.truth
        push!(xs, q[i])
        push!(ys, q[j])
    end
    for v in vs
        v[1] in (:global, :cellwise) || continue
        box = scene.boxes[v[1]]
        bc, br = LazySets.center(box), LazySets.radius_hyperrectangle(box)
        push!(xs, bc[i] - br[i], bc[i] + br[i])
        push!(ys, bc[j] - br[j], bc[j] + br[j])
    end
    mx, my = 0.05 * (maximum(xs) - minimum(xs)), 0.05 * (maximum(ys) - minimum(ys))
    return (minimum(xs) - mx, maximum(xs) + mx, minimum(ys) - my, maximum(ys) + my)
end

# The cells listed for each N of the two sampled methods, against the cell-wise bound: random points in
# the top row, uniform points in the bottom row, N growing from left to right.
function overlay_N_figure(system, scene::Scene, file; xlabel, ylabel)
    study = study_of(system)
    shown = vcat([(:cellwise, 0)], [(m, N) for m in SAMPLED for N in study.Ns])
    window = window_of(scene, shown)
    panels = []
    for m in SAMPLED, (k, N) in enumerate(study.Ns)
        v = (m, N)
        push!(
            panels,
            draw_scene_N(
                system,
                scene,
                v;
                title = "$m, N = $N: $(length(scene.cells[v])) cells",
                xlabel = m === :uniform ? xlabel : "",
                ylabel = k == 1 ? ylabel : "",
                window,
            ),
        )
    end
    fig = Plots.plot(
        panels...,
        overlay_legend_N(system);
        layout = Plots.@layout([grid(2, 3); l{0.14h}]),
        size = (1400, 780),
    )
    Plots.savefig(fig, joinpath(OUT, file))
    return nothing
end

# ------------------------------------------------------------------------------------------
# The simple pendulum: the cells at θ ≈ π and θ ≈ π/2, and where the cell-wise bound helps
# ------------------------------------------------------------------------------------------

function overlay_simple()
    domain =
        LazySets.Hyperrectangle(; low = SVector(-3.95, -5.0), high = SVector(3.95, 5.0))
    ctx = reference_abstractions(:simple; domain)
    vs = reference_of(:simple)
    sym = abstraction_of(ctx.opts[(:global, 0)])
    sets = Dict(v => successor_sets(ctx.opts[v]) for v in variants_of(:simple))
    u = SVector(0.0)
    input = SY.get_abstract_input(sym, u)
    # The cells just below θ = π (the domain ends there) and around θ = π/2, at zero velocity.
    chosen = [("θ ≈ π", SVector(3.1, 0.0)), ("θ ≈ π/2", SVector(1.6, 0.0))]
    panels = []
    notes = String[]
    pi2 = nothing
    for (k, (name, x)) in enumerate(chosen)
        q = SY.get_abstract_state(sym, x)
        cell = SY.get_concrete_elem(sym, q)
        scene = scene_of(ctx, cell, input, q, u)
        k == 2 && (pi2 = scene)
        counts = join(
            ("$(short_of(v)): $(length(scene.cells[v]))" for v in variants_of(:simple)),
            ", ",
        )
        push!(notes, "$name: $counts")
        push!(
            panels,
            draw_scene(
                :simple,
                scene,
                vs;
                title = "$name, u = 0: cells " *
                        join(("$(length(scene.cells[v]))" for v in vs), " / "),
                xlabel = "θ (rad)",
                ylabel = k == 1 ? "ω (rad/s)" : "",
                window = window_of(scene, vs),
            ),
        )
    end
    fig = Plots.plot(
        panels...,
        overlay_legend(:simple);
        layout = Plots.@layout([grid(1, 2); l{0.2h}]),
        size = (1400, 680),
    )
    Plots.savefig(fig, joinpath(OUT, "sm_overlay_simple.png"))
    println("simple pendulum overlay: ", join(notes, "\n  "))
    overlay_N_figure(
        :simple,
        pi2,
        "sm_overlayN_simple.png";
        xlabel = "θ (rad)",
        ylabel = "ω (rad/s)",
    )

    # Where does the cell-wise bound help? The number of cells that the global bound lists more, for
    # every cell and u = 0, and its average over the velocities as a function of θ.
    thetas = Float64[]
    omegas = Float64[]
    for q in 1:SY.get_n_state(sym)
        x = SY.get_concrete_state(sym, q)
        push!(thetas, x[1])
        push!(omegas, x[2])
    end
    θs, ωs =
        sort(unique(round.(thetas; digits = 9))), sort(unique(round.(omegas; digits = 9)))
    gain = fill(NaN, length(ωs), length(θs))
    for q in 1:SY.get_n_state(sym)
        x = SY.get_concrete_state(sym, q)
        g = get(sets[(:global, 0)], (q, input), nothing)
        c = get(sets[(:cellwise, 0)], (q, input), nothing)
        (g === nothing || c === nothing) && continue
        gain[
            searchsortedfirst(ωs, round(x[2]; digits = 9)),
            searchsortedfirst(θs, round(x[1]; digits = 9)),
        ] = length(g) - length(c)
    end
    heat = Plots.heatmap(
        θs,
        ωs,
        gain;
        color = Plots.cgrad(:YlOrRd),
        xlabel = "θ (rad)",
        ylabel = "ω (rad/s)",
        title = "extra cells of the global bound, u = 0",
        colorbar_title = "extra cells",
        clims = (0, maximum(filter(!isnan, gain))),
        left_margin = 11Plots.mm,
        bottom_margin = 8Plots.mm,
        titlefontsize = 20,
    )
    for (_, x) in chosen
        Plots.scatter!(
            heat,
            [x[1]],
            [x[2]];
            markershape = :star5,
            markersize = 13,
            markercolor = :white,
            markerstrokecolor = :black,
            markerstrokewidth = 2,
            label = "",
        )
    end
    mean_gain = [mean(filter(!isnan, gain[:, k])) for k in eachindex(θs)]
    profile = Plots.plot(
        θs,
        mean_gain;
        color = color_of(:simple, (:global, 0)),
        linewidth = 4,
        label = "average over ω",
        xlabel = "θ (rad)",
        ylabel = "extra cells of the global bound",
        title = "average over the velocities",
        legend = :top,
        ylims = (0, 1.35 * maximum(mean_gain)),
        left_margin = 11Plots.mm,
        bottom_margin = 8Plots.mm,
        titlefontsize = 20,
    )
    # |cos θ| is what the bound of the Jacobian depends on: the gain is largest where it is smallest.
    Plots.plot!(
        profile,
        θs,
        maximum(mean_gain) .* (1 .- abs.(cos.(θs)));
        color = :black,
        linestyle = :dash,
        linewidth = 3,
        label = "shape of 1 − |cos θ|",
    )
    for (_, x) in chosen
        Plots.vline!(profile, [x[1]]; color = :gray, linewidth = 2, label = "")
    end
    fig = Plots.plot(heat, profile; layout = (1, 2), size = (1240, 500))
    Plots.savefig(fig, joinpath(OUT, "sm_gain_simple.png"))
    println(
        "simple pendulum gain: θ ≈ π: ",
        mean_gain[argmin(abs.(θs .- 3.0))],
        ", θ ≈ π/2: ",
        mean_gain[argmin(abs.(θs .- 1.6))],
        ", maximum at θ = ",
        θs[argmax(mean_gain)],
        ", minimum at θ = ",
        θs[argmin(mean_gain)],
    )
    return nothing
end

# ------------------------------------------------------------------------------------------
# The double pendulum: the pair where the cell-wise bound helps the least and the most
# ------------------------------------------------------------------------------------------

function overlay_double()
    ctx = reference_abstractions(:double)
    vs = reference_of(:double)
    sym = abstraction_of(ctx.opts[(:global, 0)])
    sets = Dict(v => successor_sets(ctx.opts[v]) for v in variants_of(:double))
    G, C = sets[(:global, 0)], sets[(:cellwise, 0)]
    both = [p for p in keys(G) if haskey(C, p)]
    ratio(p) = length(G[p]) / length(C[p])
    println("double pendulum: $(length(both)) pairs are enabled by both growth bounds")
    # Least: no difference, among the pairs with a good number of successors; most: the largest ratio.
    equal = [p for p in both if length(G[p]) == length(C[p])]
    least =
        isempty(equal) ? both[argmin(ratio.(both))] :
        equal[argmax([length(G[p]) for p in equal])]
    most = both[argmax(ratio.(both))]
    panels = []
    largest = nothing
    for (name, (q, input)) in (("smallest gain", least), ("largest gain", most))
        u = SY.get_concrete_input(sym, input)
        cell = SY.get_concrete_elem(sym, q)
        counts = join(("$(length(get(sets[v], (q, input), [])))" for v in vs), " / ")
        println(
            "  $name: cell centre $(SY.get_concrete_state(sym, q)), u = $(u[1]); cells listed: ",
            join(
                (
                    "$(short_of(v)): $(length(get(sets[v], (q, input), [])))" for
                    v in variants_of(:double)
                ),
                ", ",
            ),
        )
        scene13 = scene_of(ctx, cell, input, q, u; plane = (1, 2))
        scene24 = with_plane(scene13, (3, 4))
        name == "largest gain" && (largest = scene24)
        push!(
            panels,
            draw_scene(
                :double,
                scene13,
                vs;
                title = "$name: angles (cells $counts)",
                xlabel = "θ1 (rad)",
                ylabel = "θ2 (rad)",
                window = window_of(scene13, vs),
            ),
        )
        push!(
            panels,
            draw_scene(
                :double,
                scene24,
                vs;
                title = "$name: velocities",
                xlabel = "ω1 (rad/s)",
                ylabel = "ω2 (rad/s)",
                window = window_of(scene24, vs),
            ),
        )
    end
    # One figure per pair: the angle plane and the velocity plane side by side.
    for (tag, two) in (("small", panels[1:2]), ("large", panels[3:4]))
        fig = Plots.plot(
            two...,
            overlay_legend(:double);
            layout = Plots.@layout([grid(1, 2); l{0.2h}]),
            size = (1400, 640),
        )
        Plots.savefig(fig, joinpath(OUT, "sm_overlay_double_$tag.png"))
    end
    # The same pair, in the velocity plane, for every N.
    overlay_N_figure(
        :double,
        largest,
        "sm_overlayN_double.png";
        xlabel = "ω1 (rad/s)",
        ylabel = "ω2 (rad/s)",
    )

    # The distribution of the ratio global / cell-wise over the pairs that both enable.
    r = ratio.(both)
    edges = [1.0, 1.0001, 1.5, 2, 3, 4, 6, 10, 20, 50, maximum(r) + 1]
    counts = [count(x -> a <= x < b, r) for (a, b) in zip(edges[1:(end - 1)], edges[2:end])]
    labels = ["1", "1–1.5", "1.5–2", "2–3", "3–4", "4–6", "6–10", "10–20", "20–50", "> 50"]
    keep = 1:length(counts)
    p = Plots.bar(
        keep,
        counts[keep] ./ length(r) .* 100;
        fillcolor = color_of(:double, (:cellwise, 0)),
        linecolor = :black,
        linewidth = 0.4,
        bar_width = 0.9,
        label = "",
        xticks = (keep, labels[keep]),
        xlabel = "cells listed by the global bound / by the cell-wise bound",
        ylabel = "share of the pairs (%)",
        title = "$(length(r)) pairs enabled by both growth bounds",
        size = (800, 470),
        ylims = (0, maximum(counts) / length(r) * 100 * 1.15),
        left_margin = 11Plots.mm,
        top_margin = 4Plots.mm,
        bottom_margin = 10Plots.mm,
    )
    Plots.savefig(p, joinpath(OUT, "sm_gain_double.png"))
    println(
        "  ratio global / cell-wise: median $(median(r)), mean $(mean(r)), maximum $(maximum(r))",
    )
    return nothing
end

# ------------------------------------------------------------------------------------------
# Coverage and tightness, for the reference configuration
# ------------------------------------------------------------------------------------------

# For pairs that all the variants enable: a fresh, independent sample of the cell is simulated.
# `covered`: the share of those points whose image goes to a cell that the variant lists; `listed`:
# the cells that the variant lists, for each cell that the sample actually reaches.
function coverage(system; pairs_used = 1500, fresh = 300)
    ctx = reference_abstractions(system)
    vs = variants_of(system)
    sym = abstraction_of(ctx.opts[(:global, 0)])
    sets = Dict(v => successor_sets(ctx.opts[v]) for v in vs)
    common = collect(intersect((keys(sets[v]) for v in vs)...))
    sort!(common)
    shuffle!(MersenneTwister(3), common)
    n_common = length(common)
    common = common[1:min(pairs_used, length(common))]
    approx = SampledSimulation(ctx.study.system, ctx.tau; nsamples = 1)
    F = ST.get_system_map(approx)
    covered = Dict(v => Float64[] for v in vs)
    listed = Dict(v => Float64[] for v in vs)
    for (q, input) in common
        cell = SY.get_concrete_elem(sym, q)
        u = SY.get_concrete_input(sym, input)
        key = SampledAbstraction._cell_key(
            999,
            SVector{LazySets.dim(cell), Float64}(LazySets.center(cell)),
        )
        reached = Int[]
        for x in RandomSampler()(cell, fresh, key)
            t = SY.get_abstract_state(sym, F(x, u))
            t === nothing || push!(reached, t)
        end
        isempty(reached) && continue
        distinct = length(unique(reached))
        for v in vs
            s = sets[v][(q, input)]
            push!(covered[v], count(t -> t in s, reached) / length(reached))
            push!(listed[v], length(s) / distinct)
        end
    end
    return (;
        used = length(covered[(:global, 0)]),
        covered,
        listed,
        common = length(common),
        n_common,
    )
end

# The figure: the two statistics against N, one panel per system and statistic; the two growth bounds
# are the horizontal dashed lines.
function coverage_figure(results)
    panels = []
    for (row, (field, ylabel)) in
        enumerate(((:covered, "points covered (%)"), (:listed, "listed / reached")))
        for (col, system) in enumerate((:simple, :double))
            study = study_of(system)
            Ns = collect(study.Ns)
            scale = field === :covered ? 100 : 1
            value(v) = scale * mean(getfield(results[system], field)[v])
            ys = [value((m, N)) for m in SAMPLED for N in Ns]
            lo, hi = extrema(vcat(ys, [value((g, 0)) for g in (:global, :cellwise)]))
            p = Plots.plot(;
                title = row == 1 ? SYSTEM_TITLES[system] : "",
                xlabel = row == 2 ? "points per cell N" : "",
                ylabel = col == 1 ? ylabel : "",
                xscale = :log10,
                xticks = (Ns, string.(Ns)),
                xlims = (Ns[1] / 1.5, Ns[end] * 1.5),
                yscale = field === :covered ? :identity : :log10,
                ylims = field === :covered ? (max(0, 5 * floor(lo / 5) - 5), 103) :
                        (lo / 2.2, hi * 1.8),
                yticks = field === :covered ? :auto :
                         (
                    log_ticks(lo / 2.2, hi * 1.8),
                    [tick_label(t, :plain) for t in log_ticks(lo / 2.2, hi * 1.8)],
                ),
                legend = false,
                left_margin = 11Plots.mm,
                bottom_margin = 7Plots.mm,
            )
            for g in (:global, :cellwise)
                Plots.hline!(
                    p,
                    [value((g, 0))];
                    color = color_of(system, (g, 0)),
                    linestyle = :dash,
                    linewidth = 3.5,
                    label = "",
                )
            end
            for m in SAMPLED
                y = [value((m, N)) for N in Ns]
                Plots.plot!(p, Ns, y; color = shade(m, 2), linewidth = 3.5, label = "")
                Plots.scatter!(
                    p,
                    Ns,
                    y;
                    markercolor = [shade(m, k) for k in eachindex(Ns)],
                    markerstrokecolor = :black,
                    markersize = 12,
                    label = "",
                )
            end
            # The values: at each N, the higher of the two series above its marker, the lower below.
            for N in Ns
                r, u = value((:random, N)), value((:uniform, N))
                for (m, y, other) in ((:random, r, u), (:uniform, u, r))
                    above = y > other || (y == other && m === :uniform)
                    label =
                        field === :covered ? (y > 99.95 ? "100" : @sprintf("%.1f", y)) :
                        (y < 10 ? @sprintf("%.2f", y) : @sprintf("%.1f", y))
                    Plots.annotate!(
                        p,
                        N,
                        field === :covered ? y + (above ? 3.4 : -3.4) :
                        y * (above ? 1.5 : 1 / 1.5),
                        Plots.text(label, 12, shade(m, 3)),
                    )
                end
            end
            push!(panels, p)
        end
    end
    legend = Plots.plot(;
        framestyle = :none,
        legend = :top,
        legend_columns = 2,
        legendfontsize = 14,
        foreground_color_legend = nothing,
    )
    Plots.plot!(
        legend,
        [NaN],
        [NaN];
        color = shade(:random, 2),
        linewidth = 5,
        label = "random points",
    )
    Plots.plot!(
        legend,
        [NaN],
        [NaN];
        color = shade(:uniform, 2),
        linewidth = 5,
        label = "uniform points",
    )
    Plots.plot!(
        legend,
        [NaN],
        [NaN];
        color = rgb(GROWTH_RGB[:global]...),
        linestyle = :dash,
        linewidth = 5,
        label = "global growth bound",
    )
    Plots.plot!(
        legend,
        [NaN],
        [NaN];
        color = rgb(GROWTH_RGB[:cellwise]...),
        linestyle = :dash,
        linewidth = 5,
        label = "cell-wise growth bound",
    )
    fig = Plots.plot(
        legend,
        panels...;
        layout = Plots.@layout([l{0.12h}; grid(2, 2)]),
        size = (1250, 700),
    )
    Plots.savefig(fig, joinpath(OUT, "sm_coverage.png"))
    return nothing
end

function coverage_outputs()
    results = Dict(s => coverage(s) for s in (:simple, :double))
    io = IOBuffer()
    println(io, "\\begin{tabular}{llrr}")
    println(io, "\\toprule")
    println(
        io,
        "system & method & fresh points covered & cells listed per cell reached \\\\",
    )
    println(io, "\\midrule")
    for s in (:simple, :double)
        for (k, v) in enumerate(variants_of(s))
            name = k == 1 ? SYSTEM_TITLES[s] : ""
            println(
                io,
                "$name & $(tex_label(s, v)) & $(@sprintf("%.1f\\%%", 100 * mean(results[s].covered[v]))) & $(@sprintf("%.2f", mean(results[s].listed[v]))) \\\\",
            )
        end
        s == :simple && println(io, "\\midrule")
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
    write(joinpath(OUT, "sm_coverage.tex"), String(take!(io)))
    coverage_figure(results)
    for s in (:simple, :double)
        println(
            "coverage $s ($(results[s].used) pairs used, $(results[s].n_common) enabled by all variants): ",
        )
        for v in variants_of(s)
            println(
                "    $(rpad(short_of(v), 12)) $(round(100 * mean(results[s].covered[v]); digits = 1))% / $(round(mean(results[s].listed[v]); digits = 2))",
            )
        end
    end
    return nothing
end

# ------------------------------------------------------------------------------------------

# The colors and the numbers of points, for the LaTeX of the decks.
function write_setup()
    open(joinpath(OUT, "sm_setup.tex"), "w") do io
        println(io, "% generated by plot_successor_methods.jl")
        for m in SAMPLED, (k, rank) in enumerate(RANKS)
            c = SHADE_RGB[m][k]
            println(
                io,
                "\\definecolor{$(m === :random ? "crand" : "cunif")$rank}{RGB}{$(c[1]),$(c[2]),$(c[3])}",
            )
        end
        for system in (:simple, :double), (k, rank) in enumerate(RANKS)
            println(io, "\\newcommand{\\smN$system$rank}{$(study_of(system).Ns[k])}")
        end
        println(
            io,
            "\\newcommand{\\smsetupsimple}{$(@sprintf("%.2f", INFO.setups[:simple]))}",
        )
        println(
            io,
            "\\newcommand{\\smsetupdouble}{$(@sprintf("%.2f", INFO.setups[:double]))}",
        )
        println(io, "\\newcommand{\\smsetupcold}{$(@sprintf("%.1f", INFO.cold_symbolic))}")
        return println(
            io,
            "\\newcommand{\\smmachine}{$(replace(INFO.cpu, "_" => "\\_")), Julia $(INFO.julia)}",
        )
    end
    return nothing
end

const PARTS = (
    "bars",
    "ref",
    "sens",
    "hist",
    "tables",
    "overlay_simple",
    "overlay_double",
    "coverage",
)

# All the figures and tables, or only the parts named on the command line (`PARTS`), e.g.
# `julia --project=<env> tests/plot_successor_methods.jl coverage ref`.
function main_plots(parts = ARGS)
    unknown = setdiff(parts, PARTS)
    isempty(unknown) || error("unknown parts $unknown; the parts are $(join(PARTS, ", "))")
    run(part) = isempty(parts) || part in parts
    for system in (:simple, :double)
        run("bars") && bars_figure(system)
        run("ref") && ref_figure(system)
        run("sens") && sensitivity_figure(system)
        run("hist") && histogram_figure(system)
        run("tables") && write_tables(system)
    end
    run("overlay_simple") && overlay_simple()
    run("overlay_double") && overlay_double()
    run("coverage") && coverage_outputs()
    write_setup()
    println("written to $(abspath(OUT))")
    return nothing
end

main_plots()
