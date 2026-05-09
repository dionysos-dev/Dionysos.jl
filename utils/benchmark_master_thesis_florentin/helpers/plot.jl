using Plots
using LaTeXStrings
import LinearAlgebra as LA
import StaticArrays: SVector
import Dionysos

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System

const BENCH_DOMAIN_COLOR = :slategray
const BENCH_OBSTACLE_COLOR = :indianred3
const BENCH_INITIAL_COLOR = :seagreen3
const BENCH_TARGET_COLOR = :steelblue3
const BENCH_TRAJ_COLOR = :navy
const BENCH_ELLIPSOID_COLOR = :darkgoldenrod2

const THESIS_DOMAIN_COLOR = :gray70
const THESIS_INITIAL_COLOR = :seagreen
const THESIS_TARGET_COLOR = :crimson
const THESIS_TRAJ_COLOR = :steelblue4
const THESIS_START_COLOR = :black
const THESIS_END_COLOR = :crimson
const THESIS_ELLIPSOID_COLOR = :darkgoldenrod2

# Load periodic helpers once when plot.jl is used standalone.
if !isdefined(@__MODULE__, :unwrap_periodic_state_list)
    include(joinpath(@__DIR__, "periodic.jl"))
end
if !isdefined(@__MODULE__, :DEFAULT_PERIODIC_DIMS)
    const DEFAULT_PERIODIC_DIMS = SVector(3, 4)
    const DEFAULT_PERIODIC_PERIODS = SVector(2pi, 2pi)
    const DEFAULT_PERIODIC_START = SVector(-pi, -pi)
end

"""
Safely wrap a set into periodic coordinates when available.
"""
function maybe_wrap_set(set, periodic_dims, periodic_periods, periodic_start; enabled::Bool)
    set === nothing && return nothing
    enabled || return set

    try
        return UT.set_in_period(set, periodic_dims, periodic_periods, periodic_start)
    catch
        return set
    end
end

function _plot_bounds_2d(set, dims)
    set === nothing && return nothing

    if set isa UT.HyperRectangle
        d1, d2 = dims
        return (set.lb[d1], set.ub[d1], set.lb[d2], set.ub[d2])
    end

    if set isa UT.LazySetMinus
        return _plot_bounds_2d(set.A, dims)
    end

    if set isa UT.LazySetUnion
        bounds = [_plot_bounds_2d(s, dims) for s in set.sets]
        filter!(!isnothing, bounds)
        isempty(bounds) && return nothing
        return (
            minimum(b[1] for b in bounds),
            maximum(b[2] for b in bounds),
            minimum(b[3] for b in bounds),
            maximum(b[4] for b in bounds),
        )
    end

    return nothing
end

function _plot_limits_2d(set, dims; margin_ratio::Float64 = 0.06)
    bounds = _plot_bounds_2d(set, dims)
    bounds === nothing && return (; xlims = :auto, ylims = :auto)

    xmin, xmax, ymin, ymax = bounds
    width = max(float(xmax - xmin), 1.0)
    height = max(float(ymax - ymin), 1.0)
    margin = margin_ratio * max(width, height)
    return (; xlims = (xmin - margin, xmax + margin), ylims = (ymin - margin, ymax + margin))
end

function thesis_plot_kwargs(; legend = :topright, size = (640, 460))
    return (;
        title = "",
        size = size,
        dpi = 300,
        legend = legend,
        framestyle = :box,
        grid = true,
        gridalpha = 0.12,
        foreground_color_grid = :gray82,
        tickfontsize = 10,
        guidefontsize = 12,
        legendfontsize = 9,
        linewidth = 0.8,
        background_color = :white,
        foreground_color_axis = :gray25,
        foreground_color_border = :gray25,
        margin = 3 * Plots.mm,
    )
end

function _projected_ellipsoid_bounds(E::UT.Ellipsoid, dims)
    P = Matrix(E.P)
    Q = try
        inv(P)
    catch
        LA.pinv(P)
    end

    d1, d2 = dims
    r1 = sqrt(max(Float64(Q[d1, d1]), 0.0))
    r2 = sqrt(max(Float64(Q[d2, d2]), 0.0))
    return (E.c[d1] - r1, E.c[d1] + r1, E.c[d2] - r2, E.c[d2] + r2)
end

function _merge_bounds_2d(bounds)
    filter!(!isnothing, bounds)
    isempty(bounds) && return nothing
    return (
        minimum(b[1] for b in bounds),
        maximum(b[2] for b in bounds),
        minimum(b[3] for b in bounds),
        maximum(b[4] for b in bounds),
    )
end

function _expand_bounds_2d(bounds; margin_ratio::Float64 = 0.08, min_width::Float64 = 0.25, min_height::Float64 = 0.25)
    bounds === nothing && return (; xlims = :auto, ylims = :auto)

    xmin, xmax, ymin, ymax = Float64.(bounds)
    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)
    width = max(xmax - xmin, min_width)
    height = max(ymax - ymin, min_height)
    margin_x = margin_ratio * width
    margin_y = margin_ratio * height
    return (;
        xlims = (cx - 0.5 * width - margin_x, cx + 0.5 * width + margin_x),
        ylims = (cy - 0.5 * height - margin_y, cy + 0.5 * height + margin_y),
    )
end

function _trajectory_bounds_2d(xs, dims)
    isempty(xs) && return nothing
    d1, d2 = dims
    xvals = [x[d1] for x in xs]
    yvals = [x[d2] for x in xs]
    return (minimum(xvals), maximum(xvals), minimum(yvals), maximum(yvals))
end

"""
Project a full ellipsoid to a 2D ellipsoid on the selected dimensions.
"""
function project_ellipsoid_2d(E::UT.Ellipsoid; dims = (1, 2))
    i, j = dims
    P = Matrix(E.P)

    Q = try
        inv(P)
    catch
        LA.pinv(P)
    end

    Q2 = Q[[i, j], [i, j]]
    P2 = try
        inv(Q2)
    catch
        LA.pinv(Q2)
    end

    c2 = [E.c[i], E.c[j]]
    return UT.Ellipsoid(P2, c2)
end

"""
Extract certified ellipsoids from a certifier result payload.
"""
function extract_ellipsoids(cert_result; max_keep = 60)
    cert_result === nothing && return UT.Ellipsoid[]

    ells = UT.Ellipsoid[]

    if hasproperty(cert_result, :lmi_data) &&
       cert_result.lmi_data !== nothing &&
       hasproperty(cert_result.lmi_data, :ellipsoids)
        raw = cert_result.lmi_data.ellipsoids
        if !isempty(raw)
            for E in raw
                E isa UT.Ellipsoid && push!(ells, E)
            end
        end
    end

    if isempty(ells) && hasproperty(cert_result, :steps)
        for rec in cert_result.steps
            hasproperty(rec, :ellipsoid) || continue
            hasproperty(rec, :status) && rec.status != :ok && continue
            E = rec.ellipsoid
            E === nothing && continue
            E isa UT.Ellipsoid && push!(ells, E)
        end
    end

    max_keep === nothing && return ells
    length(ells) <= max_keep && return ells

    step = max(1, cld(length(ells), max_keep))
    return ells[1:step:end]
end

function _shift_ellipsoid_center(E::UT.Ellipsoid, shifts::AbstractVector{<:Real})
    c = collect(Float64, E.c) .+ collect(Float64, shifts)
    return UT.Ellipsoid(Matrix(E.P), c)
end

function _periodic_copy_offsets_for_ellipsoid(E::UT.Ellipsoid, periodic_dims, periodic_periods, periodic_start)
    isempty(periodic_dims) && return [zeros(length(E.c))]

    axis_offsets = [Float64[0.0] for _ in 1:length(E.c)]

    for i in eachindex(periodic_dims)
        d = Int(periodic_dims[i])
        p = Float64(periodic_periods[i])
        s = Float64(periodic_start[i])
        lo = s
        hi = s + p
        bd = _projected_ellipsoid_bounds(E, (d, d))
        crosses_lower = bd[1] < lo
        crosses_upper = bd[2] > hi

        crosses_lower && push!(axis_offsets[d], p)
        crosses_upper && push!(axis_offsets[d], -p)
    end

    offsets = [Float64[]]
    for d in 1:length(E.c)
        offsets = [vcat(base, off) for base in offsets for off in axis_offsets[d]]
    end
    return offsets
end

function add_periodic_ellipsoid_plot_copies(ellipsoids, periodic_dims, periodic_periods, periodic_start; enabled::Bool)
    (!enabled || isempty(ellipsoids)) && return ellipsoids

    out = UT.Ellipsoid[]
    for E in ellipsoids
        for offset in _periodic_copy_offsets_for_ellipsoid(E, periodic_dims, periodic_periods, periodic_start)
            push!(out, _shift_ellipsoid_center(E, offset))
        end
    end
    return out
end

"""
Draw domain/initial/target + trajectory in selected dimensions.
"""
function _plot_nominal_trajectory!(
    fig,
    x_traj;
    dims = [1, 2],
    label::AbstractString = "nominal trajectory",
    color = THESIS_TRAJ_COLOR,
    marker_stride = nothing,
    show_start_end::Bool = false,
    show_sparse_markers::Bool = false,
    trajectory_lw::Real = 0.9,
    trajectory_ms::Real = 2.0,
    trajectory_alpha::Real = 0.88,
)
    xs = collect(ST.enum_elems(x_traj))
    isempty(xs) && return fig

    d1, d2 = dims
    xvals = [x[d1] for x in xs]
    yvals = [x[d2] for x in xs]

    plot!(
        fig,
        xvals,
        yvals;
        color = color,
        lw = trajectory_lw,
        markershape = :circle,
        ms = trajectory_ms,
        markerstrokewidth = 0,
        alpha = trajectory_alpha,
        label = label,
    )

    if show_sparse_markers
        stride = marker_stride === nothing ? max(2, cld(length(xs), 10)) : marker_stride
        marker_idx = collect(1:stride:length(xs))
        if length(marker_idx) > 2
            marker_idx = marker_idx[2:(end - 1)]
            scatter!(
                fig,
                xvals[marker_idx],
                yvals[marker_idx];
                color = color,
                markerstrokecolor = color,
                ms = 2.0,
                label = "",
            )
        end
    end

    if show_start_end
        scatter!(
            fig,
            [first(xvals)],
            [first(yvals)];
            color = THESIS_START_COLOR,
            markershape = :circle,
            ms = 4.5,
            label = "start",
        )
        scatter!(
            fig,
            [last(xvals)],
            [last(yvals)];
            color = THESIS_END_COLOR,
            markershape = :diamond,
            ms = 5.0,
            label = "end",
        )
    end

    return fig
end

function plot_state_space_basic!(
    fig,
    domain,
    initial_set,
    target_set,
    x_traj;
    dims = [1, 2],
    trajectory_label::AbstractString = "nominal trajectory",
    show_start_end::Bool = false,
    show_sparse_markers::Bool = false,
    show_trajectory::Bool = true,
)
    domain !== nothing && plot!(
        fig,
        domain;
        dims = dims,
        color = THESIS_DOMAIN_COLOR,
        opacity = 0.08,
        hole_color = BENCH_OBSTACLE_COLOR,
        hole_alpha = 0.25,
        label = "",
    )
    initial_set !== nothing &&
        plot!(fig, initial_set; dims = dims, color = THESIS_INITIAL_COLOR, opacity = 0.16, lw = 0.8, label = L"\mathcal{X}_I")
    target_set !== nothing &&
        plot!(fig, target_set; dims = dims, color = THESIS_TARGET_COLOR, opacity = 0.16, lw = 0.8, label = L"\mathcal{X}_T")

    show_trajectory || return fig

    return _plot_nominal_trajectory!(
        fig,
        x_traj;
        dims = dims,
        label = trajectory_label,
        show_start_end = show_start_end,
        show_sparse_markers = show_sparse_markers,
    )
end

"""
Save 2D state-space plots (x,y) and (theta,phi), with optional ellipsoids.
"""
function save_state_space_plots!(
    output_dir::AbstractString,
    problem,
    candidate_traj;
    cert_result = nothing,
    show_ellipsoids::Bool = true,
    unwrap_angles::Bool = false,
    wrap_angles::Bool = true,
    periodic_dims = DEFAULT_PERIODIC_DIMS,
    periodic_periods = DEFAULT_PERIODIC_PERIODS,
    periodic_start = DEFAULT_PERIODIC_START,
    title12::AbstractString = "State Space (x,y)",
    title34::AbstractString = "State Space (theta,phi)",
)
    mkpath(output_dir)

    xs = collect(ST.enum_elems(candidate_traj.x_traj))
    if unwrap_angles
        xs = unwrap_periodic_state_list(xs, periodic_dims, periodic_periods)
    end
    if wrap_angles
        xs = wrap_periodic_state_list(xs, periodic_dims, periodic_periods, periodic_start)
    end
    x_traj = ST.Trajectory(xs)

    raw_domain =
        (hasproperty(problem, :system) && hasproperty(problem.system, :X)) ? problem.system.X : nothing
    raw_initial_set = hasproperty(problem, :initial_set) ? problem.initial_set : nothing
    raw_target_set = hasproperty(problem, :target_set) ? problem.target_set : nothing

    domain =
        maybe_wrap_set(raw_domain, periodic_dims, periodic_periods, periodic_start; enabled = wrap_angles)
    initial_set = maybe_wrap_set(
        raw_initial_set,
        periodic_dims,
        periodic_periods,
        periodic_start;
        enabled = wrap_angles,
    )
    target_set = maybe_wrap_set(
        raw_target_set,
        periodic_dims,
        periodic_periods,
        periodic_start;
        enabled = wrap_angles,
    )

    ellipsoids =
        (show_ellipsoids && cert_result !== nothing) ? extract_ellipsoids(cert_result) : UT.Ellipsoid[]

    if !isempty(ellipsoids) && unwrap_angles
        ellipsoids = unwrap_ellipsoid_centers(ellipsoids, periodic_dims, periodic_periods)
    end
    if !isempty(ellipsoids) && wrap_angles
        ellipsoids = wrap_ellipsoid_centers(
            ellipsoids,
            periodic_dims,
            periodic_periods,
            periodic_start,
        )
        ellipsoids = add_periodic_ellipsoid_plot_copies(
            ellipsoids,
            periodic_dims,
            periodic_periods,
            periodic_start;
            enabled = true,
        )
    end

    lim12 = _plot_limits_2d(domain, (1, 2))
    fig12 = plot(;
        aspect_ratio = :equal,
        legend = false,
        title = title12,
        size = (850, 650),
        xlims = lim12.xlims,
        ylims = lim12.ylims,
    )
    plot_state_space_basic!(fig12, domain, initial_set, target_set, x_traj; dims = [1, 2])
    if !isempty(ellipsoids)
        for E in ellipsoids
            E2 = project_ellipsoid_2d(E; dims = (1, 2))
            plot!(fig12, E2; color = BENCH_ELLIPSOID_COLOR, opacity = 0.20, lw = 1.0, label = "")
        end
    end
    savefig(fig12, joinpath(output_dir, "state_space_12.pdf"))
    display(fig12)

    lim34 = _plot_limits_2d(domain, (3, 4))
    fig34 = plot(;
        aspect_ratio = :equal,
        legend = false,
        title = title34,
        size = (850, 650),
        xlims = lim34.xlims,
        ylims = lim34.ylims,
    )
    plot_state_space_basic!(fig34, domain, initial_set, target_set, x_traj; dims = [3, 4])
    if !isempty(ellipsoids)
        for E in ellipsoids
            E2 = project_ellipsoid_2d(E; dims = (3, 4))
            plot!(fig34, E2; color = BENCH_ELLIPSOID_COLOR, opacity = 0.20, lw = 1.0, label = "")
        end
    end
    savefig(fig34, joinpath(output_dir, "state_space_34.pdf"))
    display(fig34)

    return nothing
end

"""
Convenience wrapper when certification payload is always present.
"""
function save_certified_state_space_plots!(
    output_dir::AbstractString,
    problem,
    candidate_traj,
    cert_result;
    show_ellipsoids::Bool = true,
    unwrap_angles::Bool = false,
    wrap_angles::Bool = true,
    periodic_dims = DEFAULT_PERIODIC_DIMS,
    periodic_periods = DEFAULT_PERIODIC_PERIODS,
    periodic_start = DEFAULT_PERIODIC_START,
    title12::AbstractString = "State Space (x,y)",
    title34::AbstractString = "State Space (theta,phi)",
)
    return save_state_space_plots!(
        output_dir,
        problem,
        candidate_traj;
        cert_result = cert_result,
        show_ellipsoids = show_ellipsoids,
        unwrap_angles = unwrap_angles,
        wrap_angles = wrap_angles,
        periodic_dims = periodic_dims,
        periodic_periods = periodic_periods,
        periodic_start = periodic_start,
        title12 = title12,
        title34 = title34,
    )
end

"""
Render and optionally save vehicle rollout animation.
"""
function plot_articulated_vehicle!(
    vehicle_module,
    concrete_system,
    params,
    x_traj,
    u_traj;
    domain = concrete_system.X,
    giffile = nothing,
    fps = 20,
    every = 1,
    dt = 0.2,
    title = nothing,
)
    gr()
    limits = _plot_limits_2d(domain, (1, 2); margin_ratio = 0.08)
    xl = limits.xlims == :auto ? (-20.0, 27.0) : limits.xlims
    yl = limits.ylims == :auto ? (-20.0, 20.0) : limits.ylims
    draw_params = vehicle_module.DrawParams(params)

    return vehicle_module.live_vehicle_progression(
        params,
        draw_params,
        x_traj,
        u_traj,
        xl,
        yl;
        domain = domain,
        every = every,
        dt = dt,
        giffile = giffile,
        fps = fps,
        title = title,
    )
end
