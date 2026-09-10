# Schematic of the two lifts, drawn exactly as the code builds them:
#   ClockLift: every base transition (source → target) is replicated as
#              (source, p) → (target, p + 1), for p = 1 … ntime-1.
#   ModeLift:  one copy of the base automaton per mode, plus the switch
#              transitions carried by the automaton's edges.
# Output: presentation/gallery/lift_principle.png
#
# Run: julia --project=test presentation/asset_lift_principle.jl

ENV["GKSwstype"] = "100"
using Plots
using Plots.PlotUtils.Colors: RGBA

const NODE = :steelblue
const EDGE = :gray35
const LIFT = :darkorange

node!(p, x, y, label; c = NODE) = begin
    scatter!(
        p,
        [x],
        [y];
        markersize = 17,
        color = c,
        markerstrokecolor = :black,
        markerstrokewidth = 1,
        label = false,
    )
    annotate!(p, x, y, text(label, 6, :white))
end

arrow!(p, x1, y1, x2, y2; c = EDGE, w = 1.6, shrink = 0.13) = begin
    dx, dy = x2 - x1, y2 - y1
    L = sqrt(dx^2 + dy^2)
    ux, uy = dx / L, dy / L
    plot!(
        p,
        [x1 + shrink * ux, x2 - shrink * ux],
        [y1 + shrink * uy, y2 - shrink * uy];
        color = c,
        linewidth = w,
        arrow = arrow(:closed, 0.6, 0.6),
        label = false,
    )
end

## ---- panel A: the base abstraction -------------------------------------
pa = plot(;
    title = "base abstraction",
    titlefontsize = 9,
    axis = false,
    grid = false,
    legend = false,
    xlims = (-0.5, 2.5),
    ylims = (-1.2, 2.2),
)
for (i, q) in enumerate(["q₁", "q₂", "q₃"])
    node!(pa, i - 1, 0.5, q)
end
arrow!(pa, 0, 0.5, 1, 0.5)
arrow!(pa, 1, 0.5, 2, 0.5)

## ---- panel B: the clock lift -------------------------------------------
pb = plot(;
    title = "ClockLift: × one slice per step",
    titlefontsize = 9,
    axis = false,
    grid = false,
    legend = false,
    xlims = (-0.6, 2.6),
    ylims = (-0.6, 2.6),
)
for p in 0:2, i in 0:2
    node!(pb, i, p, "")
end
for p in 0:2
    annotate!(pb, -0.45, p, text("t$(p+1)", 7, :gray30))
end
## every base transition is replicated one slice up
for p in 0:1
    arrow!(pb, 0, p, 1, p + 1; c = LIFT)
    arrow!(pb, 1, p, 2, p + 1; c = LIFT)
end
annotate!(pb, 1.0, -0.45, text("last slice: no successor", 6, :gray30))

## ---- panel C: the mode lift --------------------------------------------
pc = plot(;
    title = "ModeLift: × one copy per mode",
    titlefontsize = 9,
    axis = false,
    grid = false,
    legend = false,
    xlims = (-0.6, 2.6),
    ylims = (-0.6, 2.6),
)
for k in (0, 1.6), i in 0:2
    node!(pc, i, k, ""; c = k == 0 ? NODE : :seagreen)
end
annotate!(pc, -0.45, 0, text("OFF", 7, :gray30))
annotate!(pc, -0.45, 1.6, text("ON", 7, :gray30))
for k in (0, 1.6)
    arrow!(pc, 0, k, 1, k)
    arrow!(pc, 1, k, 2, k)
end
## switches: one per direction
arrow!(pc, 1, 0.0, 1, 1.6; c = :firebrick, w = 1.3)
arrow!(pc, 2, 1.6, 2, 0.0; c = :firebrick, w = 1.3)
annotate!(pc, 1.55, 0.8, text("switch", 6, :firebrick))

## ---- panel D: both lifts, the product of three factors ------------------
pd = plot(;
    title = "both: state × slice × mode",
    titlefontsize = 9,
    legend = false,
    camera = (38, 24),
    grid = true,
    xlabel = "state",
    ylabel = "time",
    zlabel = "mode",
    xguidefontsize = 7,
    yguidefontsize = 7,
    zguidefontsize = 7,
    xticks = (1:3, ["q₁", "q₂", "q₃"]),
    yticks = (1:3, ["t₁", "t₂", "t₃"]),
    zticks = (1:2, ["OFF", "ON"]),
    tickfontsize = 6,
)
## GR draws no arrowhead in 3-D, so the direction is carried by a short,
## thicker stub drawn over the last third of each segment.
function seg3!(p, a, b; c, w = 2.0)
    plot!(
        p,
        [a[1], b[1]],
        [a[2], b[2]],
        [a[3], b[3]];
        color = c,
        linewidth = w,
        label = false,
    )
    m = ntuple(i -> a[i] + 0.78 * (b[i] - a[i]), 3)
    return plot!(
        p,
        [m[1], b[1]],
        [m[2], b[2]],
        [m[3], b[3]];
        color = c,
        linewidth = 2.1w,
        label = false,
    )
end
for k in 1:2, pslice in 1:3, i in 1:3
    scatter3d!(
        pd,
        [i],
        [pslice],
        [k];
        markersize = 7,
        color = k == 1 ? NODE : :seagreen,
        markerstrokecolor = :black,
        markerstrokewidth = 0.6,
        label = false,
    )
end
for k in 1:2, pslice in 1:2, i in 1:2      ## time advances with every step
    seg3!(pd, (i, pslice, k), (i + 1, pslice + 1, k); c = LIFT)
end
## switches: same two as the ModeLift panel (up at q₂, down at q₃), and they
## are available at every slice
for pslice in 1:3
    seg3!(pd, (2, pslice, 1), (2, pslice, 2); c = :firebrick, w = 1.4)
    seg3!(pd, (3, pslice, 2), (3, pslice, 1); c = :firebrick, w = 1.4)
end

fig =
    plot(pa, pb, pc, pd; layout = (1, 4), size = (1500, 330), dpi = 200, margin = 3Plots.mm)
out = joinpath(@__DIR__, "gallery", "lift_principle.png")
savefig(fig, out)
println("saved ", out)
