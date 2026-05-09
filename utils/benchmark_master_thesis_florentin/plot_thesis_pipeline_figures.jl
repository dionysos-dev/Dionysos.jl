# Generates conceptual thesis Figures 4.2--4.4 for the trajectory-guided
# certified abstraction pipeline. Outputs are written as PDF and PNG files to:
#   utils/benchmark_master_thesis_florentin/outputs/plots/thesis_pipeline_figures
#
# Run from the repository root with:
#   julia --project=. utils/benchmark_master_thesis_florentin/plot_thesis_pipeline_figures.jl

using LaTeXStrings
using LinearAlgebra
using Plots
using Random

const OUTDIR = joinpath(@__DIR__, "outputs", "plots", "thesis_pipeline_figures")

const XLIMS = (0.0, 10.0)
const YLIMS = (0.0, 5.0)
const ETA = 0.25

const COLORS = (
    background = RGB(1.0, 1.0, 1.0),
    frame = RGB(0.55, 0.55, 0.55),
    grid = RGB(0.82, 0.82, 0.82),
    obstacle = RGB(0.07, 0.07, 0.07),
    target = RGB(0.62, 0.12, 0.10),
    initial = RGB(0.10, 0.48, 0.27),
    abstract = RGB(0.20, 0.23, 0.26),
    simulated = RGB(0.50, 0.50, 0.50),
    quantization = RGB(0.25, 0.25, 0.25),
    abstract_cell = RGB(0.55, 0.66, 0.72),
    mppi = RGB(0.91, 0.36, 0.08),
    sample_light = RGB(0.98, 0.82, 0.48),
    sample_mid = RGB(0.98, 0.56, 0.16),
    ellipsoid = RGB(0.00, 0.48, 0.56),
    cert_highlight = RGB(0.00, 0.36, 0.42),
    annotation = RGB(0.10, 0.10, 0.10),
    legend_border = RGB(0.82, 0.82, 0.82),
)

const OBSTACLE_RECTS = [
    (0.0, 2.65, 6.20, 0.55),
    (5.65, 1.25, 0.55, 1.95),
]
const TARGET_RECT = (0.52, 3.73, 1.10, 0.78)
const INITIAL_CENTER = (1.00, 0.70)
const INITIAL_RADIUS = 0.34

avg(values) = sum(values) / length(values)

function rectshape(x, y, w, h)
    return Shape([x, x + w, x + w, x], [y, y, y + h, y + h])
end

function circleshape(cx, cy, r; n = 96)
    theta = range(0, 2pi; length = n)
    return Shape(cx .+ r .* cos.(theta), cy .+ r .* sin.(theta))
end

function ellipseshape(cx, cy, a, b, phi; n = 120)
    theta = range(0, 2pi; length = n)
    c, s = cos(phi), sin(phi)
    x = cx .+ a .* cos.(theta) .* c .- b .* sin.(theta) .* s
    y = cy .+ a .* cos.(theta) .* s .+ b .* sin.(theta) .* c
    return Shape(x, y)
end

function in_obstacle(pt)
    x, y = pt
    return any(xr <= x <= xr + w && yr <= y <= yr + h for (xr, yr, w, h) in OBSTACLE_RECTS)
end

function base_plot(; size = (980, 520), show_grid = false)
    p = plot(;
        size = size,
        dpi = 300,
        xlims = XLIMS,
        ylims = YLIMS,
        aspect_ratio = :equal,
        legend = false,
        ticks = false,
        grid = false,
        framestyle = :box,
        foreground_color_border = COLORS.frame,
        background_color = COLORS.background,
        background_color_inside = COLORS.background,
        margin = 4Plots.mm,
        fontfamily = "Computer Modern",
    )
    show_grid && draw_uniform_grid!(p)
    return p
end

function add_direct_label!(p, x, y, label; color = COLORS.annotation, size = 10, halign = :center)
    annotate!(p, x, y, text(label, size, color, halign))
    return p
end

function draw_environment!(p)
    plot!(p, rectshape(XLIMS[1], YLIMS[1], XLIMS[2] - XLIMS[1], YLIMS[2] - YLIMS[1]);
        fillcolor = :white, fillalpha = 0.0, linecolor = COLORS.frame, linewidth = 0.85, label = "")

    for (x, y, w, h) in OBSTACLE_RECTS
        plot!(p, rectshape(x, y, w, h);
            fillcolor = COLORS.obstacle, fillalpha = 1.0, linecolor = COLORS.obstacle,
            linewidth = 0.0, label = "")
    end
    annotate!(p, 3.02, 2.93, text(L"\mathcal{X}_O", 14, :white, :center))

    plot!(p, rectshape(TARGET_RECT...);
        fillcolor = COLORS.target, fillalpha = 0.94, linecolor = COLORS.target,
        linewidth = 1.0, label = "")
    annotate!(p, TARGET_RECT[1] + TARGET_RECT[3] / 2, TARGET_RECT[2] + TARGET_RECT[4] / 2,
        text(L"\mathcal{X}_T", 14, :white, :center))

    plot!(p, circleshape(INITIAL_CENTER..., INITIAL_RADIUS);
        fillcolor = COLORS.initial, fillalpha = 0.90, linecolor = COLORS.initial,
        linewidth = 1.0, label = "")
    annotate!(p, INITIAL_CENTER[1], INITIAL_CENTER[2],
        text(L"\mathcal{X}_I", 13, :white, :center))
    return p
end

function draw_uniform_grid!(p; eta = ETA)
    for x in XLIMS[1]:eta:XLIMS[2]
        plot!(p, [x, x], [YLIMS[1], YLIMS[2]];
            color = COLORS.grid, alpha = 0.22, linewidth = 0.36, label = "")
    end
    for y in YLIMS[1]:eta:YLIMS[2]
        plot!(p, [XLIMS[1], XLIMS[2]], [y, y];
            color = COLORS.grid, alpha = 0.22, linewidth = 0.36, label = "")
    end
    return p
end

function abstract_centers()
    return [
        (1.125, 0.625), (1.625, 0.875), (2.375, 1.125), (3.125, 1.375),
        (4.125, 1.375), (5.125, 1.375), (6.625, 1.625), (7.125, 2.625),
        (6.625, 3.625), (5.125, 4.125), (3.375, 4.125), (1.125, 4.125),
    ]
end

function draw_visited_cells!(p, centers; eta = ETA)
    for (x, y) in centers
        plot!(p, rectshape(x - eta / 2, y - eta / 2, eta, eta);
            fillcolor = COLORS.abstract_cell, fillalpha = 0.16, linecolor = COLORS.abstract_cell,
            linealpha = 0.30, linewidth = 0.55, label = "")
    end
    return p
end

function draw_discrete_arrows!(p, pts; color = COLORS.abstract, linewidth = 1.5, alpha = 1.0, arrowsize = 0.12)
    for i in 1:(length(pts) - 1)
        x0, y0 = pts[i]
        x1, y1 = pts[i + 1]
        dx, dy = x1 - x0, y1 - y0
        quiver!(p, [x0 + 0.10dx], [y0 + 0.10dy];
            quiver = ([0.80dx], [0.80dy]), color = color, linewidth = linewidth,
            alpha = alpha, arrow = arrow(:closed, arrowsize), label = "")
    end
    return p
end

function draw_centers!(p, centers; markersize = 3.6, alpha = 1.0)
    scatter!(p, first.(centers), last.(centers);
        markercolor = :white, markerstrokecolor = COLORS.abstract, markersize = markersize,
        markerstrokewidth = 1.0, alpha = alpha, label = "")
    return p
end

function simulated_images_from_centers(centers; eta = ETA)
    xsim = Tuple{Float64,Float64}[]
    offsets = [0.075, -0.065, 0.080, -0.070, 0.060, -0.075, 0.080, -0.085, 0.070, -0.060, 0.070]
    backsteps = [0.055, 0.045, 0.060, 0.050, 0.055, 0.060, 0.050, 0.060, 0.050, 0.045, 0.050]
    for i in 1:(length(centers) - 1)
        x0, y0 = centers[i]
        x1, y1 = centers[i + 1]
        dx, dy = x1 - x0, y1 - y0
        nrm = max(sqrt(dx^2 + dy^2), eps())
        px, py = -dy / nrm, dx / nrm
        candidate = (x1 - backsteps[i] * dx / nrm + offsets[i] * px,
                     y1 - backsteps[i] * dy / nrm + offsets[i] * py)
        candidate = (clamp(candidate[1], x1 - eta / 2 + 0.015, x1 + eta / 2 - 0.015),
                     clamp(candidate[2], y1 - eta / 2 + 0.015, y1 + eta / 2 - 0.015))
        push!(xsim, candidate)
    end
    return xsim
end

function draw_simulated_transition_arrows!(p, centers, xsim = simulated_images_from_centers(centers);
    alpha = 0.62)
    for i in 1:length(xsim)
        x0, y0 = centers[i]
        xs, ys = xsim[i]
        quiver!(p, [x0], [y0];
            quiver = ([xs - x0], [ys - y0]), color = COLORS.simulated,
            alpha = alpha, linewidth = 1.05, arrow = arrow(:closed, 0.085), label = "")
    end
    return p
end

function draw_quantization_mismatch!(p, centers, xsim; highlight_indices = [3, 7, 9])
    highlight = Set(highlight_indices)
    for i in 1:length(xsim)
        xq, yq = centers[i + 1]
        xs, ys = xsim[i]
        is_highlight = i in highlight
        plot!(p, [xs, xq], [ys, yq];
            color = COLORS.quantization, alpha = is_highlight ? 0.84 : 0.42,
            linewidth = is_highlight ? 1.25 : 0.85, linestyle = :dot, label = "")
        scatter!(p, [xs], [ys];
            markershape = :diamond, markercolor = COLORS.simulated,
            markerstrokecolor = COLORS.simulated, markersize = is_highlight ? 4.6 : 3.5,
            alpha = is_highlight ? 0.88 : 0.62, label = "")
    end
    return p
end

function draw_abstraction_quantization_inset!(p)
    x, y, w, h = 3.18, 0.15, 3.24, 1.00
    draw_legend_box!(p, x, y + h, w, h; alpha = 0.92)

    cell = 0.42
    src = (x + 0.62, y + 0.58)
    dst = (x + 2.42, y + 0.58)
    sim = (dst[1] - 0.18, dst[2] + 0.10)

    for c in (src, dst)
        plot!(p, rectshape(c[1] - cell / 2, c[2] - cell / 2, cell, cell);
            fillcolor = COLORS.abstract_cell, fillalpha = 0.16,
            linecolor = COLORS.abstract_cell, linealpha = 0.42, linewidth = 0.55, label = "")
    end
    quiver!(p, [src[1] + 0.07], [src[2] + 0.01];
        quiver = ([sim[1] - src[1] - 0.11], [sim[2] - src[2] - 0.02]),
        color = COLORS.simulated, linewidth = 1.15, alpha = 0.82,
        arrow = arrow(:closed, 0.095), label = "")
    plot!(p, [sim[1], dst[1]], [sim[2], dst[2]];
        color = COLORS.quantization, linewidth = 1.15, linestyle = :dot, alpha = 0.92, label = "")

    scatter!(p, [src[1], dst[1]], [src[2], dst[2]];
        markercolor = :white, markerstrokecolor = COLORS.abstract,
        markerstrokewidth = 0.90, markersize = 3.8, label = "")
    scatter!(p, [sim[1]], [sim[2]];
        markershape = :diamond, markercolor = COLORS.simulated,
        markerstrokecolor = COLORS.simulated, markersize = 4.2, alpha = 0.90, label = "")

    add_direct_label!(p, src[1], src[2] + 0.31, L"c(q_k)"; size = 7)
    add_direct_label!(p, dst[1] + 0.13, dst[2] - 0.31, L"c(q_{k+1})"; size = 7)
    add_direct_label!(p, x + w / 2, y + 0.15,
        L"q_{k+1}=Q_h\!\left(f(c(q_k),\bar u_k^{\mathrm{abs}},0)\right)";
        size = 7)
    return p
end

function draw_quantization_reference!(p, centers; idx = 4)
    x0, y0 = centers[idx]
    x1, y1 = centers[idx + 1]
    cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
    plot!(p, ellipseshape(cx, cy, 0.58, 0.20, atan(y1 - y0, x1 - x0));
        fillcolor = :white, fillalpha = 0.0, linecolor = COLORS.simulated,
        linealpha = 0.32, linewidth = 0.65, label = "")
    plot!(p, [cx, 4.72], [cy - 0.20, 1.14];
        color = COLORS.simulated, alpha = 0.25, linewidth = 0.50, label = "")
    return p
end

function bezier_path(control; n_per_segment = 28)
    pts = Tuple{Float64,Float64}[]
    for i in 1:(length(control) - 1)
        p0 = control[max(i - 1, 1)]
        p1 = control[i]
        p2 = control[i + 1]
        p3 = control[min(i + 2, length(control))]
        for t in range(0, 1; length = n_per_segment)
            t2, t3 = t^2, t^3
            x = 0.5 * ((2p1[1]) + (-p0[1] + p2[1]) * t +
                (2p0[1] - 5p1[1] + 4p2[1] - p3[1]) * t2 +
                (-p0[1] + 3p1[1] - 3p2[1] + p3[1]) * t3)
            y = 0.5 * ((2p1[2]) + (-p0[2] + p2[2]) * t +
                (2p0[2] - 5p1[2] + 4p2[2] - p3[2]) * t2 +
                (-p0[2] + 3p1[2] - 3p2[2] + p3[2]) * t3)
            push!(pts, (x, y))
        end
    end
    return pts
end

function abstract_control_rollout_path()
    return bezier_path([
        (1.00, 0.70), (2.30, 0.98), (3.80, 1.12), (5.38, 1.18),
        (6.18, 1.48), (6.05, 2.15), (5.82, 2.82), (5.10, 3.25),
        (3.70, 3.55), (2.05, 3.58),
    ]; n_per_segment = 18)
end

function final_mppi_path()
    return bezier_path([
        (1.00, 0.72), (2.15, 0.90), (3.75, 0.82), (5.35, 0.92),
        (6.72, 1.26), (7.22, 2.18), (6.72, 3.45), (5.24, 4.08),
        (3.10, 4.22), (1.20, 4.20),
    ]; n_per_segment = 24)
end

function smooth_noise(rng, s, phase)
    return 0.72sin(2pi * s + phase) + 0.28sin(4pi * s + 0.6phase) + 0.10randn(rng)
end

function noisy_sample(nominal, seed, quality)
    rng = MersenneTwister(seed)
    m = length(nominal)
    phase_x = 2pi * rand(rng)
    phase_y = 2pi * rand(rng)
    amp = 0.14 + 0.78 * (1 - quality)
    pts = Tuple{Float64,Float64}[]
    for (i, (x, y)) in enumerate(nominal)
        s = (i - 1) / max(m - 1, 1)
        taper = sin(pi * s)^0.62
        nx = amp * taper * smooth_noise(rng, s, phase_x)
        ny = 1.28amp * taper * smooth_noise(rng, s, phase_y)
        push!(pts, (clamp(x + nx, XLIMS[1] + 0.08, XLIMS[2] - 0.08),
                    clamp(y + ny, YLIMS[1] + 0.08, YLIMS[2] - 0.08)))
    end
    return pts
end

function draw_path!(p, pts; color, linewidth, alpha = 1.0, linestyle = :solid)
    plot!(p, first.(pts), last.(pts);
        color = color, linewidth = linewidth, alpha = alpha, linestyle = linestyle, label = "")
    return p
end

function draw_sparse_markers_on_path!(p, pts; every = 24, color = :white,
    strokecolor = COLORS.simulated, markersize = 2.6)
    idxs = 1:every:length(pts)
    scatter!(p, first.(pts[idxs]), last.(pts[idxs]);
        markercolor = color, markerstrokecolor = strokecolor, markerstrokewidth = 0.9,
        markersize = markersize, label = "")
    return p
end

function draw_mppi_samples!(p, nominal; nsamples = 150)
    qualities = range(0.03, 0.98; length = nsamples)
    palette = cgrad([RGB(1.00, 0.90, 0.60), COLORS.sample_light, COLORS.sample_mid, COLORS.mppi])
    for (j, q) in enumerate(qualities)
        pts = noisy_sample(nominal, 900 + j, q)
        draw_path!(p, pts; color = palette[q], linewidth = 0.32 + 0.50q, alpha = 0.035 + 0.20q)
    end
    return p
end

function sampled_ellipsoids()
    path = final_mppi_path()
    idxs = round.(Int, range(5, length(path) - 5; length = 14))
    ells = []
    for (j, idx) in enumerate(idxs)
        x, y = path[idx]
        xnext, ynext = path[min(idx + 4, length(path))]
        phi = atan(ynext - y, xnext - x)
        open_space = min(abs(x - 5.90) + abs(y - 2.30), 2.0) / 2.0
        a = 0.34 + 0.15open_space
        b = 0.16 + 0.07open_space
        if 4.7 < x < 7.1 && 0.9 < y < 3.55
            a *= 0.70
            b *= 0.70
        end
        if j == length(idxs)
            x, y = (1.13, 4.13)
            a, b = 0.30, 0.15
            phi = 0.00
        end
        push!(ells, (x, y, a, b, phi, j))
    end
    return ells
end

function draw_ellipsoid_tuple!(p, ell; color = COLORS.ellipsoid, fillalpha = 0.12,
    linealpha = 0.90, linewidth = 1.3)
    x, y, a, b, phi, _ = ell
    plot!(p, ellipseshape(x, y, a, b, phi);
        fillcolor = color, fillalpha = fillalpha, linecolor = color,
        linealpha = linealpha, linewidth = linewidth, label = "")
    return p
end

function draw_ellipsoids!(p; full = true)
    shown = full ? sampled_ellipsoids() : sampled_ellipsoids()[1:2:end]
    for ell in shown
        draw_ellipsoid_tuple!(p, ell)
    end
    return p
end

function draw_ellipsoid_centers!(p, ells)
    scatter!(p, [ell[1] for ell in ells], [ell[2] for ell in ells];
        markercolor = COLORS.cert_highlight, markerstrokecolor = :white,
        markerstrokewidth = 0.65, markersize = 3.2, alpha = 0.95, label = "")
    return p
end

function draw_certified_initial_set!(p, ell)
    shape = certified_initial_intersection_shape(ell)
    overlap_fill = RGB(0.30, 0.82, 0.76)
    plot!(p, shape;
        fillcolor = overlap_fill, fillalpha = 0.58,
        linecolor = COLORS.cert_highlight, linealpha = 0.98,
        linewidth = 1.75, label = "")
    cx = avg(shape.x)
    cy = avg(shape.y)
    plot!(p, [cx + 0.05, cx + 0.46], [cy - 0.02, cy - 0.20];
        color = COLORS.cert_highlight, alpha = 0.62, linewidth = 0.70, label = "")
    add_direct_label!(p, cx + 0.54, cy - 0.25, L"\mathcal{X}_I^{\mathrm{cert}}";
        color = COLORS.cert_highlight, size = 11, halign = :left)
    return p
end

function ellipse_contains(pt, ell)
    x, y = pt
    cx, cy, a, b, phi, _ = ell
    c, s = cos(phi), sin(phi)
    dx, dy = x - cx, y - cy
    xr = c * dx + s * dy
    yr = -s * dx + c * dy
    return (xr / a)^2 + (yr / b)^2 <= 1.0 + 1e-8
end

function circle_contains(pt, center, radius)
    x, y = pt
    cx, cy = center
    return (x - cx)^2 + (y - cy)^2 <= radius^2 + 1e-8
end

function certified_initial_intersection_shape(ell; n = 360)
    pts = Tuple{Float64,Float64}[]
    for theta in range(0, 2pi; length = n)
        p = (INITIAL_CENTER[1] + INITIAL_RADIUS * cos(theta),
             INITIAL_CENTER[2] + INITIAL_RADIUS * sin(theta))
        ellipse_contains(p, ell) && push!(pts, p)
    end
    cx, cy, a, b, phi, _ = ell
    c, s = cos(phi), sin(phi)
    for theta in range(0, 2pi; length = n)
        p = (cx + a * cos(theta) * c - b * sin(theta) * s,
             cy + a * cos(theta) * s + b * sin(theta) * c)
        circle_contains(p, INITIAL_CENTER, INITIAL_RADIUS) && push!(pts, p)
    end
    if length(pts) < 4
        return ellipseshape(INITIAL_CENTER[1] + 0.10, INITIAL_CENTER[2] - 0.02, 0.18, 0.10, 0.05)
    end
    mx = avg(first.(pts))
    my = avg(last.(pts))
    sort!(pts; by = p -> atan(p[2] - my, p[1] - mx))
    return Shape(first.(pts), last.(pts))
end

function draw_legend_box!(p, x, y, w, h; alpha = 0.90)
    plot!(p, rectshape(x, y - h, w, h);
        fillcolor = :white, fillalpha = alpha, linecolor = COLORS.legend_border,
        linewidth = 0.60, label = "")
    return p
end

function draw_legend_line_item!(p, x, y, label; color, linewidth,
    linestyle = :solid, alpha = 1.0, textcolor = COLORS.annotation, size = 9)
    plot!(p, [x, x + 0.56], [y, y];
        color = color, linewidth = linewidth, linestyle = linestyle, alpha = alpha, label = "")
    add_direct_label!(p, x + 0.72, y, label; color = textcolor, size = size, halign = :left)
    return p
end

function draw_legend_marker_item!(p, x, y, label; marker = :circle, color = COLORS.simulated,
    strokecolor = color)
    scatter!(p, [x + 0.24], [y]; markershape = marker, markercolor = color,
        markerstrokecolor = strokecolor, markersize = 3.8, markerstrokewidth = 0.8, label = "")
    add_direct_label!(p, x + 0.72, y, label; size = 9, halign = :left)
    return p
end

function draw_legend_ellipse_item!(p, x, y, label; color = COLORS.ellipsoid, linewidth = 1.3)
    plot!(p, ellipseshape(x + 0.28, y, 0.25, 0.09, 0.0);
        fillcolor = color, fillalpha = 0.11, linecolor = color,
        linewidth = linewidth, label = "")
    add_direct_label!(p, x + 0.72, y, label; color = color, size = 9, halign = :left)
    return p
end

function draw_manual_legend!(p, kind)
    if kind == :abstraction
        x, y = 6.42, 4.78
        draw_legend_box!(p, x, y, 3.15, 0.86)
        quiver!(p, [x + 0.16], [y - 0.18]; quiver = ([0.42], [0.0]),
            color = COLORS.abstract, linewidth = 1.5, arrow = arrow(:closed, 0.10), label = "")
        add_direct_label!(p, x + 0.82, y - 0.18, "abstract transition"; size = 9, halign = :left)
        draw_legend_line_item!(p, x + 0.16, y - 0.45, "simulated image";
            color = COLORS.simulated, linewidth = 1.1, size = 9)
        plot!(p, [x + 0.16, x + 0.72], [y - 0.72, y - 0.72];
            color = COLORS.quantization, linewidth = 0.9, linestyle = :dot, alpha = 0.80, label = "")
        add_direct_label!(p, x + 0.82, y - 0.72, "quantization to center"; size = 9, halign = :left)
    elseif kind == :mppi
        x, y = 6.20, 4.78
        draw_legend_box!(p, x, y, 3.48, 0.88)
        draw_legend_line_item!(p, x + 0.16, y - 0.18, "rollout induced by abstract controls";
            color = COLORS.simulated, linewidth = 1.7, size = 9)
        draw_legend_line_item!(p, x + 0.16, y - 0.45, "MPPI sampled rollouts";
            color = COLORS.sample_mid, linewidth = 0.75, alpha = 0.55, size = 9)
        draw_legend_line_item!(p, x + 0.16, y - 0.72, "refined nominal trajectory";
            color = COLORS.mppi, linewidth = 2.8, textcolor = COLORS.mppi, size = 9)
    elseif kind == :certification
        x, y = 6.38, 4.78
        draw_legend_box!(p, x, y, 3.00, 0.62)
        draw_legend_line_item!(p, x + 0.16, y - 0.18, "nominal trajectory";
            color = COLORS.mppi, linewidth = 2.8, textcolor = COLORS.mppi, size = 9)
        draw_legend_ellipse_item!(p, x + 0.16, y - 0.48, "certified ellipsoids";
            color = COLORS.ellipsoid, linewidth = 1.3)
    end
    return p
end

draw_top_left_legend!(p, kind) = draw_manual_legend!(p, kind)

function draw_grid_scale!(p)
    x0, y0 = 3.05, 0.38
    plot!(p, [x0, x0 + ETA], [y0, y0]; color = COLORS.annotation, linewidth = 1.0, label = "")
    add_direct_label!(p, x0 + ETA / 2, y0 + 0.15, L"\eta"; size = 10)
    return p
end

function draw_local_certificate_transition!(p, ells; k = 9)
    Ek = ells[k]
    Enext = ells[k + 1]
    draw_ellipsoid_tuple!(p, Ek; color = COLORS.cert_highlight, fillalpha = 0.20,
        linealpha = 0.98, linewidth = 2.10)
    draw_ellipsoid_tuple!(p, Enext; color = COLORS.cert_highlight, fillalpha = 0.20,
        linealpha = 0.98, linewidth = 2.10)

    x0, y0 = Ek[1], Ek[2]
    x1, y1 = Enext[1], Enext[2]
    add_direct_label!(p, x0 + 0.52, y0 - 0.34, L"\mathcal{E}_k";
        color = COLORS.cert_highlight, size = 11)
    add_direct_label!(p, x1 + 0.52, y1 + 0.36, L"\mathcal{E}_{k+1}";
        color = COLORS.cert_highlight, size = 11)
    return p
end

function draw_local_certificate_inset!(p)
    x, y, w, h = 6.80, 0.12, 3.05, 1.34
    draw_legend_box!(p, x, y + h, w, h; alpha = 0.93)
    add_direct_label!(p, x + 0.16, y + h - 0.14, "one-step certificate"; size = 8, halign = :left)

    source = (x + 0.48, y + 0.78, 0.34, 0.15, 0.03, 1)
    succ = (x + 2.23, y + 0.78, 0.62, 0.28, 0.04, 2)
    image = (x + 2.05, y + 0.78, 0.25, 0.10, 0.08, 3)

    draw_ellipsoid_tuple!(p, source; color = COLORS.cert_highlight, fillalpha = 0.09,
        linealpha = 0.95, linewidth = 1.15)
    scatter!(p, [source[1]], [source[2]];
        markercolor = COLORS.cert_highlight, markerstrokecolor = :white,
        markerstrokewidth = 0.55, markersize = 2.5, label = "")
    draw_ellipsoid_tuple!(p, succ; color = COLORS.cert_highlight, fillalpha = 0.07,
        linealpha = 0.96, linewidth = 1.20)
    draw_ellipsoid_tuple!(p, image; color = COLORS.mppi, fillalpha = 0.18,
        linealpha = 0.95, linewidth = 1.00)

    quiver!(p, [source[1] + 0.36], [source[2]];
        quiver = ([image[1] - source[1] - 0.64], [image[2] - source[2]]),
        color = COLORS.annotation, linewidth = 1.15, alpha = 0.90,
        arrow = arrow(:closed, 0.10), label = "")

    add_direct_label!(p, source[1], y + 1.02, L"\xi";
        color = COLORS.cert_highlight, size = 7)
    add_direct_label!(p, image[1], y + 0.55, L"q(\xi)";
        color = COLORS.mppi, size = 7)
    add_direct_label!(p, succ[1] + 0.35, y + 1.04, L"\xi^+";
        color = COLORS.cert_highlight, size = 7)
    add_direct_label!(p, x + 1.24, y + 0.93, L"\kappa";
        color = COLORS.annotation, size = 7)
    add_direct_label!(p, x + w / 2, y + 0.18, L"q(\xi)\subseteq\xi^+";
        color = COLORS.annotation, size = 8)
    return p
end

function save_both(p, basename)
    mkpath(OUTDIR)
    savefig(p, joinpath(OUTDIR, basename * ".pdf"))
    savefig(p, joinpath(OUTDIR, basename * ".png"))
end

function figure_4_2()
    p = base_plot(; show_grid = true)
    centers = abstract_centers()
    xsim = simulated_images_from_centers(centers)
    draw_visited_cells!(p, centers)
    draw_environment!(p)
    draw_simulated_transition_arrows!(p, centers, xsim)
    draw_quantization_mismatch!(p, centers, xsim; highlight_indices = [3, 7, 9])
    draw_discrete_arrows!(p, centers; linewidth = 1.5, alpha = 0.98)
    draw_centers!(p, centers; markersize = 4.3)
    draw_quantization_reference!(p, centers; idx = 4)
    draw_abstraction_quantization_inset!(p)
    draw_manual_legend!(p, :abstraction)
    return p
end

function figure_4_3()
    p = base_plot()
    draw_environment!(p)
    rollout = abstract_control_rollout_path()
    final = final_mppi_path()
    draw_mppi_samples!(p, final; nsamples = 175)
    draw_path!(p, rollout; color = COLORS.simulated, linewidth = 1.7, alpha = 0.85)
    draw_sparse_markers_on_path!(p, rollout; every = 24, markersize = 2.4)
    draw_path!(p, final; color = COLORS.mppi, linewidth = 3.0, alpha = 0.98)
    scatter!(p, [first(final)[1], last(final)[1]], [first(final)[2], last(final)[2]];
        markercolor = :white, markerstrokecolor = COLORS.mppi, markersize = 4.0,
        markerstrokewidth = 1.0, label = "")
    add_direct_label!(p, 4.78, 4.07, L"\tau_{\mathrm{nom}}";
        color = COLORS.mppi, size = 9)
    draw_manual_legend!(p, :mppi)
    return p
end

function figure_4_4()
    p = base_plot()
    draw_environment!(p)
    final = final_mppi_path()
    ells = sampled_ellipsoids()
    draw_path!(p, final; color = COLORS.mppi, linewidth = 2.9, alpha = 0.98)
    draw_ellipsoids!(p; full = true)
    draw_ellipsoid_centers!(p, ells)
    draw_certified_initial_set!(p, ells[1])
    draw_local_certificate_transition!(p, ells; k = 9)
    add_direct_label!(p, 1.30, 4.62, L"\mathcal{E}_N\subseteq\mathcal{X}_T";
        color = COLORS.ellipsoid, size = 9)
    draw_local_certificate_inset!(p)
    draw_manual_legend!(p, :certification)
    return p
end

function main()
    default(; linewidth = 1.4, guidefont = font(11), tickfont = font(9), titlefont = font(12))
    figures = [
        ("fig_4_2_centered_simulation_abstraction", figure_4_2()),
        ("fig_4_3_mppi_refinement", figure_4_3()),
        ("fig_4_4_backward_certification", figure_4_4()),
    ]
    for (name, fig) in figures
        save_both(fig, name)
    end
    println("Saved thesis pipeline figures to: ", OUTDIR)
end

main()
