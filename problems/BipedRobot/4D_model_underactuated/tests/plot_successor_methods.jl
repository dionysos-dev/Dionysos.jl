# Figures and tables of the presentation "Four ways to compute the successors" (the two decks
# `successor_methods_detailed.tex` and `successor_methods_visual.tex` in
# `Documentation/ground_cells_presentation/`), from the results of `compare_successor_methods.jl`.
#
# It writes, in `OUTPUT_DIR` (by default the `figures/` folder of the presentations; every file starts
# with `sm_`):
#
# * `sm_bars_<system>.png`: the four statistics (total and average number of transitions, total and
#   average computation time) of the four methods for every sampling time and cell size, as bars;
# * `sm_hist_<system>.png`: the histogram of the number of successors of a (cell, input) pair, for the
#   reference configuration;
# * `sm_table_<system>_transitions.tex`, `sm_table_<system>_time.tex`, `sm_table_<system>_ref.tex`: the
#   same numbers as tables;
# * `sm_overlay_simple.png`, `sm_overlay_double.png`: the successors of one cell and one input with the four
#   methods on top of each other (the cells that each method lists, the growth-bound boxes, the images of
#   the sample points, the true image), for the cells where the cell-wise bound helps the least and the
#   most; `sm_gain_simple.png` and `sm_gain_double.png`: where the cell-wise bound helps, over the state
#   space (simple pendulum) and as a distribution (double pendulum);
# * `sm_coverage.png` and `sm_coverage.tex`: how much of what a fresh point of a cell does each method
#   covers, and how many cells it lists for each cell that is actually reached.
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
    joinpath(@__DIR__, "..", "Documentation", "ground_cells_presentation", "figures")
const RESULTS = joinpath(@__DIR__, "successor_methods_results", "results.jls")

# The configuration that is looked at in detail, for every system.
const REFERENCE = (
    simple = (tau = 0.1, size = [0.1, 0.2]),
    double = (tau = 0.01, size = [0.3, 0.3, 0.6, 0.6]),
)

# The colors of the deck (`preamble.tex`): red, orange, blue, green.
const COLORS = Dict(
    :global => Plots.RGB(190 / 255, 35 / 255, 35 / 255),
    :cellwise => Plots.RGB(215 / 255, 105 / 255, 10 / 255),
    :random => Plots.RGB(30 / 255, 100 / 255, 170 / 255),
    :uniform => Plots.RGB(27 / 255, 120 / 255, 55 / 255),
)
const LABELS = Dict(
    :global => "global growth bound",
    :cellwise => "cell-wise growth bound",
    :random => "random points",
    :uniform => "uniform points",
)
const SHORT = Dict(
    :global => "global bound",
    :cellwise => "cell-wise bound",
    :random => "random points",
    :uniform => "uniform points",
)
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

record(system, method, tau, size) = only(
    r for r in RECORDS if
    r.system == system && r.method == method && r.tau == tau && r.size == size
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

# Ticks of a logarithmic axis that read as numbers (1, 2, 5, 10, ...), between `lo` and `hi`.
function log_ticks(lo, hi)
    ticks = Float64[]
    for e in floor(Int, log10(lo)):ceil(Int, log10(hi)), m in (1, 2, 5)
        t = m * 10.0^e
        lo <= t <= hi && push!(ticks, t)
    end
    return ticks
end

# Bars: one group per configuration (`configs`, labelled by `labelf`), one bar per method. `values` maps
# every method to its numbers, in the order of `configs`. `taus`: also write the sampling time above
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
)
    n = length(configs)
    positive = [v for m in METHODS for v in values[m] if v > 0]
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
        xtickfontsize = n > 5 ? 11 : 14,
        left_margin = 9Plots.mm,
        bottom_margin = 3Plots.mm,
    )
    width = 0.2
    for (m, method) in enumerate(METHODS), g in 1:n
        v = values[method][g]
        x0 = g + (m - 2.5) * width
        if v > 0
            Plots.plot!(
                p,
                Plots.Shape(
                    [x0 - width / 2, x0 + width / 2, x0 + width / 2, x0 - width / 2],
                    [base, base, v, v],
                );
                fillcolor = COLORS[method],
                linecolor = :black,
                linewidth = 0.4,
                label = "",
            )
            annotate === nothing || Plots.annotate!(
                p,
                x0,
                v * (logscale ? 1.25 : 1.0) + (logscale ? 0.0 : 0.05 * top),
                Plots.text(annotate(v), 11, :black),
            )
        else
            Plots.annotate!(
                p,
                x0,
                base * (logscale ? 1.6 : 1.0),
                Plots.text("0", 11, COLORS[method]),
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

# A strip with the legend of the four methods (a plot without axes).
function legend_strip(; fontsize = 15)
    p = Plots.plot(;
        framestyle = :none,
        legend = :top,
        legend_columns = 4,
        legendfontsize = fontsize,
        foreground_color_legend = nothing,
    )
    for method in METHODS
        Plots.bar!(
            p,
            [NaN],
            [NaN];
            fillcolor = COLORS[method],
            linecolor = :black,
            label = SHORT[method],
        )
    end
    return p
end

function values_of(system, field; configs = configurations(system))
    return Dict(
        m => Float64[field(record(system, m, t, s)) for (t, s) in configs] for m in METHODS
    )
end

# The four statistics, for the nine configurations of a system: one figure for the transitions
# (total, per pair) and one for the time.
function bars_figure(system)
    transitions = values_of(system, r -> r.transitions)
    average = values_of(system, r -> r.transitions / r.pairs)
    total_time = values_of(system, r -> r.time)
    average_time = values_of(system, r -> 1e6 * r.time / r.pairs)
    spread =
        maximum(average[:global]) /
        max(minimum(v for m in METHODS for v in average[m] if v > 0), 1e-9)
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
            legend_strip(),
            panels...;
            layout = Plots.@layout([l{0.11h}; grid(1, 2)]),
            size = (1500, 700),
        )
        Plots.savefig(fig, joinpath(OUT, "sm_bars_$(system)_$name.png"))
    end
    return nothing
end

# The four statistics of the reference configuration, one big bar per method.
function ref_figure(system)
    ref = REFERENCE[system]
    rs = Dict(m => record(system, m, ref.tau, ref.size) for m in METHODS)
    panels = []
    for (name, field, kind, logscale) in (
        ("total transitions", r -> r.transitions, :count, false),
        ("transitions per pair", r -> r.transitions / r.pairs, :plain, false),
        ("total time", r -> r.time, :seconds, false),
        ("time per pair (µs)", r -> 1e6 * r.time / r.pairs, :plain, false),
    )
        vals = [Float64(field(rs[m])) for m in METHODS]
        p = Plots.bar(
            1:4,
            vals;
            fillcolor = [COLORS[m] for m in METHODS],
            linecolor = :black,
            linewidth = 0.6,
            bar_width = 0.72,
            label = "",
            xticks = (1:4, ["global", "cell-wise", "random", "uniform"]),
            ylims = (0, maximum(vals) * 1.22),
            yformatter = v -> tick_label(v, kind),
            title = name,
            left_margin = 4Plots.mm,
            xtickfontsize = 15,
        )
        for (i, v) in enumerate(vals)
            Plots.annotate!(
                p,
                i,
                v + 0.07 * maximum(vals),
                Plots.text(tick_label(v, kind == :plain ? :plain : kind), 15, :black),
            )
        end
        push!(panels, p)
    end
    fig = Plots.plot(
        panels...;
        layout = (2, 2),
        size = (1300, 640),
        plot_title = "$(SYSTEM_TITLES[system]): τ = $(fmt_g(ref.tau)), cells $(join(map(fmt_g, ref.size), "×"))",
        plot_titlefontsize = 18,
    )
    Plots.savefig(fig, joinpath(OUT, "sm_ref_$system.png"))
    return nothing
end

# The effect of the sampling time (at the reference cell size) and of the cell size (at the reference
# sampling time) on the two statistics per pair.
function sensitivity_figure(system)
    ref = REFERENCE[system]
    all = configurations(system)
    by_tau = [(t, s) for (t, s) in all if s == ref.size]
    by_size = [(t, s) for (t, s) in all if t == ref.tau]
    average(c) = values_of(system, r -> r.transitions / r.pairs; configs = c)
    micro(c) = values_of(system, r -> 1e6 * r.time / r.pairs; configs = c)
    spread(v) =
        maximum(maximum(v[m]) for m in METHODS) /
        max(minimum(x for m in METHODS for x in v[m] if x > 0), 1e-9)
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
        legend_strip(),
        panels...;
        layout = Plots.@layout([l{0.08h}; grid(2, 2)]),
        size = (1200, 720),
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

function histogram_figure(system)
    ref = REFERENCE[system]
    exact = system == :simple ? 14 : 4
    all_bins = Dict(
        m => begin
            hist = record(system, m, ref.tau, ref.size).hist
            edges = bins_of(hist; exact)
            [sum(hist[(a + 1):min(b + 1, length(hist))]) for (a, b) in edges]
        end for m in METHODS
    )
    longest = maximum(length(all_bins[m]) for m in METHODS)
    edges = bins_of(zeros(Int, 2^(longest) + 4); exact)[1:longest]
    labels = [a == b ? string(a) : "$(a)–$(b)" for (a, b) in edges]
    panels = []
    pairs = record(system, :global, ref.tau, ref.size).pairs
    ymax = maximum(maximum(all_bins[m]) for m in METHODS) / pairs * 100 * 1.12
    for (i, method) in enumerate(METHODS)
        counts = vcat(all_bins[method], zeros(Int, longest - length(all_bins[method])))
        share = counts ./ pairs .* 100
        p = Plots.bar(
            1:longest,
            share;
            fillcolor = COLORS[method],
            linecolor = :black,
            linewidth = 0.4,
            bar_width = 0.9,
            label = "",
            title = LABELS[method],
            xticks = (1:longest, labels),
            xrotation = 45,
            xtickfontsize = 13,
            ylims = (0, ymax),
            ylabel = i in (1, 3) ? "share of the pairs (%)" : "",
            xlabel = i in (3, 4) ? "number of successors" : "",
            left_margin = 11Plots.mm,
            bottom_margin = 8Plots.mm,
        )
        mean_k = sum((a + b) / 2 * c for ((a, b), c) in zip(edges, counts)) / pairs
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

function header_rows(columns)
    cols = join(("c" for _ in 1:(2 + 1 + length(columns) * length(METHODS)));)
    return cols
end

function table_config(system, body_columns)
    # body_columns: a vector of (header, formatter(record) -> String); one group per method.
    k = length(body_columns)
    spec = "rrr" * join(("r" ^ k for _ in METHODS))
    io = IOBuffer()
    println(io, "\\begin{tabular}{$spec}")
    println(io, "\\toprule")
    heads = join(
        (
            "\\multicolumn{$k}{c}{\\textcolor{$(m == :global ? "cbelow" : m == :cellwise ? "cacross" : m == :random ? "cabove" : "cgood")}{$(LABELS[m])}}"
            for m in METHODS
        ),
        " & ",
    )
    println(io, " & & & $heads \\\\")
    start = 4
    rules = join(
        (
            "\\cmidrule(lr){$(start + (i - 1) * k)-$(start + i * k - 1)}" for
            i in 1:length(METHODS)
        ),
        " ",
    )
    println(io, rules)
    sub = join((join((h for (h, _) in body_columns), " & ") for _ in METHODS), " & ")
    println(io, "\$\\tau\$ & cell & cells & $sub \\\\")
    println(io, "\\midrule")
    configs = configurations(system)
    last_tau = nothing
    for (tau, size) in configs
        last_tau !== nothing && tau != last_tau && println(io, "\\addlinespace")
        last_tau = tau
        rs = [record(system, m, tau, size) for m in METHODS]
        cells = fmt_int(rs[1].n_cells)
        vals = join((join((f(r) for (_, f) in body_columns), " & ") for r in rs), " & ")
        println(io, "$(fmt_g(tau)) & $(cell_text(system, size)) & $cells & $vals \\\\")
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
    return String(take!(io))
end

function write_tables(system)
    write(
        joinpath(OUT, "sm_table_$(system)_transitions.tex"),
        table_config(
            system,
            [
                ("total", r -> fmt_int(r.transitions)),
                ("per pair", r -> fmt_avg(r.transitions / r.pairs)),
                ("enabled", r -> @sprintf("%.0f\\%%", 100 * r.enabled_pairs / r.pairs)),
            ],
        ),
    )
    write(
        joinpath(OUT, "sm_table_$(system)_time.tex"),
        table_config(
            system,
            [
                ("total", r -> fmt_time(r.time)),
                ("\$\\mu\$s per pair", r -> fmt_us(1e6 * r.time / r.pairs)),
            ],
        ),
    )
    # The reference configuration: one row per method.
    ref = REFERENCE[system]
    io = IOBuffer()
    println(io, "\\begin{tabular}{lrrrrrr}")
    println(io, "\\toprule")
    println(
        io,
        "method & transitions & per pair & per enabled pair & enabled pairs & time & \$\\mu\$s per pair \\\\",
    )
    println(io, "\\midrule")
    for m in METHODS
        r = record(system, m, ref.tau, ref.size)
        per_enabled = r.enabled_pairs == 0 ? "--" : fmt_avg(r.transitions / r.enabled_pairs)
        println(
            io,
            "\\textcolor{$(m == :global ? "cbelow" : m == :cellwise ? "cacross" : m == :random ? "cabove" : "cgood")}{$(LABELS[m])} & $(fmt_int(r.transitions)) & $(fmt_avg(r.transitions / r.pairs)) & $per_enabled & $(@sprintf("%.1f\\%%", 100 * r.enabled_pairs / r.pairs)) & $(fmt_time(r.time)) & $(fmt_us(1e6 * r.time / r.pairs)) \\\\",
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
    study = only(s for s in studies() if s.name == system)
    # A domain of its own, for looking at one cell: the pendulum of the study ends at θ = π, which
    # disables every input of the cells next to it.
    domain === nothing ||
        (study = merge(study, (; system = SimplePendulum.system(; _X_ = domain))))
    ref = REFERENCE[system]
    dynamics = SymbolicDynamics(study.system)
    grids = grids_of(study, ref.size)
    opts = Dict(m => build(m, study, dynamics, grids, ref.tau) for m in METHODS)
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
    cells::Dict{Symbol, Vector{Any}}  # method => the cells that it lists (hyperrectangles)
    boxes::Dict{Symbol, Any}          # growth bounds => the box
    points::Dict{Symbol, Vector{Any}} # sampled methods => the images of the sample points
    truth::Vector{Any}                # images of many points of the cell
    outline::Union{Nothing, Vector{Any}}   # the image of the boundary of the cell, in order (2-D only)
end

function scene_of(ctx, source_cell, input_index, q, u; n_truth = 3000, plane = (1, 2))
    sym = abstraction_of(ctx.opts[:global])
    study, tau = ctx.study, ctx.tau
    system = study.system
    approx = SampledSimulation(system, tau; nsamples = 1)
    F = ST.get_system_map(approx)
    cells = Dict{Symbol, Vector{Any}}()
    for m in METHODS
        s = abstraction_of(ctx.opts[m])
        cells[m] = [SY.get_concrete_elem(s, t) for t in SY.post(s, q, input_index)]
    end
    boxes = Dict(
        :global => growth_box(study, study.jacobian_bound, tau, source_cell, u),
        :cellwise =>
            growth_box(study, cellwise_bound(ctx.dynamics), tau, source_cell, u),
    )
    points = Dict{Symbol, Vector{Any}}()
    for (m, sampling) in ((:random, :random), (:uniform, :uniform))
        a = SampledSimulation(system, tau; nsamples = study.N, sampling, seed = 1)
        image = ST.get_over_approximation_map(a)(source_cell, u)
        points[m] = collect(image.points)
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

function draw_scene(scene::Scene; title, xlabel, ylabel, window)
    i, j = scene.plane
    p = Plots.plot(;
        title,
        xlabel,
        ylabel,
        xlims = (window[1], window[2]),
        ylims = (window[3], window[4]),
        legend = false,
        titlefontsize = 20,
        left_margin = 11Plots.mm,
        bottom_margin = 6Plots.mm,
    )
    # The true image of the cell: the image of its boundary (2-D), or of many of its points.
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
    # The cells that each method lists, the largest set first and with the thickest outline, so that
    # the nested sets show as stripes where their edges coincide.
    for (m, width) in zip(METHODS, (5.5, 4, 2.8, 1.6))
        for (lo, hi) in projected_cells(scene.cells[m], scene.plane)
            Plots.plot!(
                p,
                box_shape(lo, hi);
                fillcolor = COLORS[m],
                fillalpha = 0.07,
                linecolor = COLORS[m],
                linewidth = width,
                label = "",
            )
        end
    end
    # The growth-bound boxes (before the cells are looked up).
    for m in (:global, :cellwise)
        box = scene.boxes[m]
        box === nothing && continue
        c, r = LazySets.center(box), LazySets.radius_hyperrectangle(box)
        Plots.plot!(
            p,
            box_shape(SVector(c[i] - r[i], c[j] - r[j]), SVector(c[i] + r[i], c[j] + r[j]));
            fillalpha = 0,
            linecolor = COLORS[m],
            linestyle = :dash,
            linewidth = 2.5,
            label = "",
        )
    end
    # The images of the sample points.
    Plots.scatter!(
        p,
        [q[i] for q in scene.points[:random]],
        [q[j] for q in scene.points[:random]];
        markershape = :circle,
        markersize = 8,
        markercolor = COLORS[:random],
        markerstrokecolor = :white,
        markerstrokewidth = 1.5,
        label = "",
    )
    Plots.scatter!(
        p,
        [q[i] for q in scene.points[:uniform]],
        [q[j] for q in scene.points[:uniform]];
        markershape = :diamond,
        markersize = 9,
        markercolor = COLORS[:uniform],
        markerstrokecolor = :white,
        markerstrokewidth = 2,
        label = "",
    )
    # The cell itself.
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

# The key to the drawings of the successor sets: a plot without axes.
function overlay_legend()
    p = Plots.plot(;
        framestyle = :none,
        legend = :top,
        legend_columns = 3,
        legendfontsize = 15,
        foreground_color_legend = nothing,
    )
    for m in METHODS
        Plots.plot!(p, [NaN], [NaN]; color = COLORS[m], linewidth = 5, label = SHORT[m])
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
        markercolor = COLORS[:random],
        label = "image: random pt",
    )
    Plots.scatter!(
        p,
        [NaN],
        [NaN];
        markershape = :diamond,
        markersize = 9,
        markercolor = COLORS[:uniform],
        label = "image: uniform pt",
    )
    Plots.plot!(p, [NaN], [NaN]; color = :gray, linewidth = 6, label = "true image")
    Plots.plot!(p, [NaN], [NaN]; color = :black, linewidth = 5, label = "the cell")
    return p
end

# A window that contains everything that is drawn, with a margin.
function window_of(scene::Scene)
    i, j = scene.plane
    xs = Float64[]
    ys = Float64[]
    for m in METHODS, cell in scene.cells[m]
        c, r = LazySets.center(cell), LazySets.radius_hyperrectangle(cell)
        push!(xs, c[i] - r[i], c[i] + r[i])
        push!(ys, c[j] - r[j], c[j] + r[j])
    end
    for q in scene.truth
        push!(xs, q[i])
        push!(ys, q[j])
    end
    for m in (:global, :cellwise)
        box = scene.boxes[m]
        box === nothing && continue
        bc, br = LazySets.center(box), LazySets.radius_hyperrectangle(box)
        push!(xs, bc[i] - br[i], bc[i] + br[i])
        push!(ys, bc[j] - br[j], bc[j] + br[j])
    end
    mx, my = 0.05 * (maximum(xs) - minimum(xs)), 0.05 * (maximum(ys) - minimum(ys))
    return (minimum(xs) - mx, maximum(xs) + mx, minimum(ys) - my, maximum(ys) + my)
end

# ------------------------------------------------------------------------------------------
# The simple pendulum: the cells at θ ≈ π and θ ≈ π/2, and where the cell-wise bound helps
# ------------------------------------------------------------------------------------------

function overlay_simple()
    domain =
        LazySets.Hyperrectangle(; low = SVector(-3.95, -5.0), high = SVector(3.95, 5.0))
    ctx = reference_abstractions(:simple; domain)
    sym = abstraction_of(ctx.opts[:global])
    sets = Dict(m => successor_sets(ctx.opts[m]) for m in METHODS)
    u = SVector(0.0)
    input = SY.get_abstract_input(sym, u)
    # The cells just below θ = π (the domain ends there) and around θ = π/2, at zero velocity.
    chosen = [("θ ≈ π", SVector(3.1, 0.0)), ("θ ≈ π/2", SVector(1.6, 0.0))]
    panels = []
    notes = String[]
    for (k, (name, x)) in enumerate(chosen)
        q = SY.get_abstract_state(sym, x)
        cell = SY.get_concrete_elem(sym, q)
        scene = scene_of(ctx, cell, input, q, u)
        counts = join(("$(m): $(length(scene.cells[m]))" for m in METHODS), ", ")
        push!(notes, "$name: $counts")
        c = LazySets.center(cell)
        push!(
            panels,
            draw_scene(
                scene;
                title = "$name, u = 0: cells " *
                        join(("$(length(scene.cells[m]))" for m in METHODS), " / "),
                xlabel = "θ (rad)",
                ylabel = k == 1 ? "ω (rad/s)" : "",
                window = window_of(scene),
            ),
        )
    end
    fig = Plots.plot(
        panels...,
        overlay_legend();
        layout = Plots.@layout([grid(1, 2); l{0.2h}]),
        size = (1400, 680),
    )
    Plots.savefig(fig, joinpath(OUT, "sm_overlay_simple.png"))
    println("simple pendulum overlay: ", join(notes, "  |  "))

    # Where does the cell-wise bound help? The number of cells that the global bound lists more, for
    # every cell and u = 0, and its average over the velocities as a function of θ.
    mapping = SY.get_state_mapping(sym)
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
        g = get(sets[:global], (q, input), nothing)
        c = get(sets[:cellwise], (q, input), nothing)
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
        color = COLORS[:global],
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
    sym = abstraction_of(ctx.opts[:global])
    sets = Dict(m => successor_sets(ctx.opts[m]) for m in METHODS)
    both = [p for p in keys(sets[:global]) if haskey(sets[:cellwise], p)]
    ratio(p) = length(sets[:global][p]) / length(sets[:cellwise][p])
    println("double pendulum: $(length(both)) pairs are enabled by both growth bounds")
    # Least: no difference, among the pairs with a good number of successors; most: the largest ratio.
    equal = [p for p in both if length(sets[:global][p]) == length(sets[:cellwise][p])]
    least =
        isempty(equal) ? both[argmin(ratio.(both))] :
        equal[argmax([length(sets[:global][p]) for p in equal])]
    most = both[argmax(ratio.(both))]
    panels = []
    for (name, (q, input)) in (("smallest gain", least), ("largest gain", most))
        u = SY.get_concrete_input(sym, input)
        cell = SY.get_concrete_elem(sym, q)
        counts = join(("$(length(get(sets[m], (q, input), [])))" for m in METHODS), " / ")
        println(
            "  $name: cell centre $(SY.get_concrete_state(sym, q)), u = $(u[1]); cells listed (global, cell-wise, random, uniform): ",
            [length(get(sets[m], (q, input), [])) for m in METHODS],
        )
        scene13 = scene_of(ctx, cell, input, q, u; plane = (1, 2))
        scene24 = Scene(
            (3, 4),
            scene13.source,
            scene13.cells,
            scene13.boxes,
            scene13.points,
            scene13.truth,
            nothing,
        )
        push!(
            panels,
            draw_scene(
                scene13;
                title = "$name: angles (cells $counts)",
                xlabel = "θ1 (rad)",
                ylabel = "θ2 (rad)",
                window = window_of(scene13),
            ),
        )
        push!(
            panels,
            draw_scene(
                scene24;
                title = "$name: velocities",
                xlabel = "ω1 (rad/s)",
                ylabel = "ω2 (rad/s)",
                window = window_of(scene24),
            ),
        )
    end
    # One figure per pair: the angle plane and the velocity plane side by side.
    for (tag, two) in (("small", panels[1:2]), ("large", panels[3:4]))
        fig = Plots.plot(
            two...,
            overlay_legend();
            layout = Plots.@layout([grid(1, 2); l{0.2h}]),
            size = (1400, 640),
        )
        Plots.savefig(fig, joinpath(OUT, "sm_overlay_double_$tag.png"))
    end

    # The distribution of the ratio global / cell-wise over the pairs that both enable.
    r = ratio.(both)
    edges = [1.0, 1.0001, 1.5, 2, 3, 4, 6, 10, 20, 50, maximum(r) + 1]
    counts = [count(x -> a <= x < b, r) for (a, b) in zip(edges[1:(end - 1)], edges[2:end])]
    labels = ["1", "1–1.5", "1.5–2", "2–3", "3–4", "4–6", "6–10", "10–20", "20–50", "> 50"]
    keep = 1:length(counts)
    p = Plots.bar(
        keep,
        counts[keep] ./ length(r) .* 100;
        fillcolor = COLORS[:cellwise],
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

# For pairs that all four methods enable: a fresh, independent sample of the cell is simulated.
# `covered`: the share of those points whose image goes to a cell that the method lists; `listed`:
# the cells that the method lists, for each cell that the sample actually reaches.
function coverage(system; pairs_used = 1500, fresh = 300)
    ctx = reference_abstractions(system)
    sym = abstraction_of(ctx.opts[:global])
    sets = Dict(m => successor_sets(ctx.opts[m]) for m in METHODS)
    common = collect(intersect((keys(sets[m]) for m in METHODS)...))
    sort!(common)
    shuffle!(MersenneTwister(3), common)
    common = common[1:min(pairs_used, length(common))]
    approx = SampledSimulation(ctx.study.system, ctx.tau; nsamples = 1)
    F = ST.get_system_map(approx)
    covered = Dict(m => Float64[] for m in METHODS)
    listed = Dict(m => Float64[] for m in METHODS)
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
        for m in METHODS
            s = sets[m][(q, input)]
            push!(covered[m], count(t -> t in s, reached) / length(reached))
            push!(listed[m], length(s) / distinct)
        end
    end
    return (; used = length(covered[:global]), covered, listed, common = length(common))
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
        for (k, m) in enumerate(METHODS)
            name = k == 1 ? SYSTEM_TITLES[s] : ""
            println(
                io,
                "$name & \\textcolor{$(m == :global ? "cbelow" : m == :cellwise ? "cacross" : m == :random ? "cabove" : "cgood")}{$(LABELS[m])} & $(@sprintf("%.1f\\%%", 100 * mean(results[s].covered[m]))) & $(@sprintf("%.2f", mean(results[s].listed[m]))) \\\\",
            )
        end
        s == :simple && println(io, "\\midrule")
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
    write(joinpath(OUT, "sm_coverage.tex"), String(take!(io)))
    systems = (:simple, :double)
    configs = [(1, nothing), (2, nothing)]
    group_label = (k, _) -> k == 1 ? "simple" : "double"
    covered =
        Dict(m => [100 * mean(results[s].covered[m]) for s in systems] for m in METHODS)
    listed = Dict(m => [mean(results[s].listed[m]) for s in systems] for m in METHODS)
    panels = [
        bar_panel(
            :simple,
            covered,
            "fresh points covered (%)";
            configs,
            labelf = group_label,
            logscale = false,
            annotate = v -> v > 99.95 ? "100" : @sprintf("%.1f", v),
            ymax = 118,
        ),
        bar_panel(
            :simple,
            listed,
            "cells listed / cells reached";
            configs,
            labelf = group_label,
            logscale = true,
            annotate = v -> @sprintf("%.2g", v),
            ymax = 60,
        ),
    ]
    fig = Plots.plot(
        legend_strip(; fontsize = 12),
        panels...;
        layout = Plots.@layout([l{0.13h}; grid(1, 2)]),
        size = (1000, 480),
    )
    Plots.savefig(fig, joinpath(OUT, "sm_coverage.png"))
    for s in (:simple, :double)
        println(
            "coverage $s ($(results[s].used) pairs): ",
            join(
                (
                    "$(m): $(round(100 * mean(results[s].covered[m]); digits = 1))% / $(round(mean(results[s].listed[m]); digits = 2))"
                    for m in METHODS
                ),
                ", ",
            ),
        )
    end
    return nothing
end

# ------------------------------------------------------------------------------------------

function main_plots()
    for system in (:simple, :double)
        bars_figure(system)
        ref_figure(system)
        sensitivity_figure(system)
        histogram_figure(system)
        write_tables(system)
    end
    overlay_simple()
    overlay_double()
    coverage_outputs()
    open(joinpath(OUT, "sm_setup.tex"), "w") do io
        println(io, "% generated by plot_successor_methods.jl")
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
    println("written to $(abspath(OUT))")
    return nothing
end

main_plots()
