using Plots
using LaTeXStrings
import LinearAlgebra as LA
import StaticArrays: SVector
import Dionysos

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System

const PLOT_COLORS = Dict(
    :domain => "#F8FAFC",
    :domain_edge => "#CBD5E1",
    :obstacle => "#000000",
    :obstacle_edge => "#000000",
    :initial => "#3DB80D",
    :initial_edge => "#0B7A25",
    :target => "#DE0707",
    :target_edge => "#9F071F",
    :nominal => "#0072B2",
    :state2 => "#E69F00",
    :heading => "#009E73",
    :trailer1 => "#CC79A7",
    :trailer2 => "#D55E00",
    :ellipsoid => "#EDED13",
    :ellipsoid_edge => "#8B3A75",
    :rollout => "#94A3B8",
    :success => "#0B7A25",
    :failure => "#D55E00",
    :constraint => "#4B5563",
    :grid => "#E5E7EB",
)

const BENCH_DOMAIN_COLOR = PLOT_COLORS[:domain]
const BENCH_OBSTACLE_COLOR = PLOT_COLORS[:obstacle]
const BENCH_INITIAL_COLOR = PLOT_COLORS[:initial]
const BENCH_TARGET_COLOR = PLOT_COLORS[:target]
const BENCH_TRAJ_COLOR = PLOT_COLORS[:nominal]
const BENCH_ELLIPSOID_COLOR = PLOT_COLORS[:ellipsoid]

const THESIS_DOMAIN_COLOR = PLOT_COLORS[:domain]
const THESIS_INITIAL_COLOR = PLOT_COLORS[:initial]
const THESIS_TARGET_COLOR = PLOT_COLORS[:target]
const THESIS_TRAJ_COLOR = PLOT_COLORS[:nominal]
const THESIS_START_COLOR = PLOT_COLORS[:initial_edge]
const THESIS_END_COLOR = PLOT_COLORS[:target_edge]
const THESIS_ELLIPSOID_COLOR = PLOT_COLORS[:ellipsoid]

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
    return (;
        xlims = (xmin - margin, xmax + margin),
        ylims = (ymin - margin, ymax + margin),
    )
end

function _state_axis_label(dim::Integer)
    dim == 1 && return L"x\,[\mathrm{m}]"
    dim == 2 && return L"y\,[\mathrm{m}]"
    dim == 3 && return L"\theta_1\,[\mathrm{rad}]"
    dim == 4 && return L"\varphi\,[\mathrm{rad}]"
    return "x$(dim)"
end

wrap_to_pi(a) = mod(a + pi, 2pi) - pi

lift_angle_to_reference(a, ref) = ref + wrap_to_pi(a - ref)

function unwrap_angle_sequence(a)
    isempty(a) && return collect(Float64, a)

    out = collect(Float64, a)
    for k in 2:length(out)
        out[k] = lift_angle_to_reference(out[k], out[k - 1])
    end
    return out
end

function unwrap_angle_sequence_to_reference(a, ref)
    length(a) == length(ref) ||
        error("angle and reference sequences must have the same length")
    return [lift_angle_to_reference(a[k], ref[k]) for k in eachindex(a)]
end

function _periodic_dim_lookup(dim::Integer, periodic_dims, periodic_periods)
    periodic_dims === nothing && return nothing
    periodic_periods === nothing && return nothing

    for i in eachindex(periodic_dims)
        Int(periodic_dims[i]) == Int(dim) && return (i, Float64(periodic_periods[i]))
    end
    return nothing
end

_is_periodic_plot_dim(dim::Integer, periodic_dims, periodic_periods) =
    _periodic_dim_lookup(dim, periodic_dims, periodic_periods) !== nothing

function _lift_periodic_value_to_reference(a::Real, ref::Real, period::Real)
    period > 0 || error("period must be positive")
    return Float64(a) +
           round((Float64(ref) - Float64(a)) / Float64(period)) * Float64(period)
end

function _lift_state_to_plot_reference(x, ref, periodic_dims, periodic_periods)
    y = collect(Float64, x)
    ref === nothing && return y

    for i in eachindex(periodic_dims)
        d = Int(periodic_dims[i])
        d <= length(y) || continue
        d <= length(ref) || continue
        y[d] = _lift_periodic_value_to_reference(y[d], ref[d], periodic_periods[i])
    end
    return y
end

function _lift_state_list_to_plot_reference(
    state_list,
    ref_states,
    periodic_dims,
    periodic_periods,
)
    isempty(state_list) && return SVector{0, Float64}[]
    if ref_states === nothing ||
       isempty(ref_states) ||
       periodic_dims === nothing ||
       periodic_periods === nothing
        return [SVector{length(x), Float64}(x) for x in state_list]
    end

    nx = length(first(state_list))
    out = Vector{SVector{nx, Float64}}(undef, length(state_list))
    for k in eachindex(state_list)
        ref = ref_states[min(k, length(ref_states))]
        out[k] = SVector{nx, Float64}(
            _lift_state_to_plot_reference(
                state_list[k],
                ref,
                periodic_dims,
                periodic_periods,
            ),
        )
    end
    return out
end

function _unwrap_states_for_angle_plot(state_list, periodic_dims, periodic_periods)
    xs = [SVector{length(x), Float64}(x) for x in state_list]
    isempty(xs) && return xs
    periodic_dims === nothing && return xs
    periodic_periods === nothing && return xs
    return unwrap_periodic_state_list(xs, periodic_dims, periodic_periods)
end

function _projection_has_periodic_dim(dims, periodic_dims, periodic_periods)
    return any(d -> _is_periodic_plot_dim(d, periodic_dims, periodic_periods), dims)
end

function _reference_state_for_plot(ref_states, k::Integer)
    ref_states === nothing && return nothing
    isempty(ref_states) && return nothing
    return ref_states[clamp(Int(k), 1, length(ref_states))]
end

function _lift_ellipsoid_for_projection(
    E::UT.Ellipsoid,
    dims;
    ref_state = nothing,
    periodic_dims = nothing,
    periodic_periods = nothing,
)
    ref_state === nothing && return E
    _projection_has_periodic_dim(dims, periodic_dims, periodic_periods) || return E

    c = collect(Float64, E.c)
    for d in dims
        lookup = _periodic_dim_lookup(d, periodic_dims, periodic_periods)
        lookup === nothing && continue
        _, period = lookup
        c[Int(d)] = _lift_periodic_value_to_reference(c[Int(d)], ref_state[Int(d)], period)
    end
    return UT.Ellipsoid(Matrix(E.P), c)
end

function _lift_set_for_projection(
    set,
    dims;
    ref_state = nothing,
    periodic_dims = nothing,
    periodic_periods = nothing,
)
    set === nothing && return nothing
    ref_state === nothing && return set
    _projection_has_periodic_dim(dims, periodic_dims, periodic_periods) || return set

    if set isa UT.HyperRectangle
        lb = collect(Float64, set.lb)
        ub = collect(Float64, set.ub)
        c = 0.5 .* (lb .+ ub)
        r = 0.5 .* (ub .- lb)

        for d in dims
            lookup = _periodic_dim_lookup(d, periodic_dims, periodic_periods)
            lookup === nothing && continue
            _, period = lookup
            idx = Int(d)
            c[idx] = _lift_periodic_value_to_reference(c[idx], ref_state[idx], period)
            lb[idx] = c[idx] - r[idx]
            ub[idx] = c[idx] + r[idx]
        end

        return UT.HyperRectangle(lb, ub)
    end

    if set isa UT.LazySetUnion
        return UT.LazySetUnion([
            _lift_set_for_projection(
                s,
                dims;
                ref_state = ref_state,
                periodic_dims = periodic_dims,
                periodic_periods = periodic_periods,
            ) for s in set.sets
        ])
    end

    return set
end

function _plot_domain_with_obstacles!(fig, domain, dims)
    domain === nothing && return fig

    if domain isa UT.LazySetMinus
        plot!(
            fig,
            domain.A;
            dims = dims,
            color = THESIS_DOMAIN_COLOR,
            linecolor = PLOT_COLORS[:domain_edge],
            opacity = 0.22,
            lw = 0.8,
            label = L"\mathcal{X}",
        )

        if sort(collect(dims)) == [1, 2]
            plot!(
                fig,
                domain.B;
                dims = dims,
                color = BENCH_OBSTACLE_COLOR,
                linecolor = PLOT_COLORS[:obstacle_edge],
                opacity = 0.86,
                lw = 0.9,
                label = L"\mathcal{X}_{obs}",
            )
        end

        return fig
    end

    plot!(
        fig,
        domain;
        dims = dims,
        color = THESIS_DOMAIN_COLOR,
        linecolor = PLOT_COLORS[:domain_edge],
        opacity = 0.22,
        lw = 0.8,
        label = L"\mathcal{X}",
    )
    return fig
end

function thesis_plot_kwargs(; legend = :topright, size = (640, 460))
    return (;
        title = "",
        size = size,
        dpi = 300,
        legend = legend,
        framestyle = :box,
        grid = true,
        gridalpha = 0.18,
        foreground_color_grid = PLOT_COLORS[:grid],
        tickfontsize = 10,
        guidefontsize = 12,
        legendfontsize = 10,
        linewidth = 1.2,
        background_color = :white,
        background_color_legend = :white,
        foreground_color_legend = :gray25,
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

function _expand_bounds_2d(
    bounds;
    margin_ratio::Float64 = 0.08,
    min_width::Float64 = 0.25,
    min_height::Float64 = 0.25,
)
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
Extract certified ellipsoids from a certifier result payload, keeping the
1-based nominal state index used for plot-time angle lifting.
"""
function extract_ellipsoid_records(cert_result; max_keep = 300)
    cert_result === nothing && return NamedTuple[]

    records = NamedTuple[]

    if hasproperty(cert_result, :steps)
        for rec in cert_result.steps
            hasproperty(rec, :ellipsoid) || continue
            hasproperty(rec, :status) && rec.status != :ok && continue
            E = rec.ellipsoid
            E === nothing && continue
            if E isa UT.Ellipsoid
                state_index = hasproperty(rec, :k) ? Int(rec.k) : length(records) + 1
                push!(records, (; state_index = state_index, ellipsoid = E))
            end
        end

        if !isempty(records) &&
           hasproperty(cert_result, :lmi_data) &&
           cert_result.lmi_data !== nothing &&
           hasproperty(cert_result.lmi_data, :ellipsoids) &&
           !isempty(cert_result.lmi_data.ellipsoids) &&
           first(cert_result.lmi_data.ellipsoids) isa UT.Ellipsoid
            terminal_index = maximum(rec.state_index for rec in records) + 1
            push!(
                records,
                (;
                    state_index = terminal_index,
                    ellipsoid = first(cert_result.lmi_data.ellipsoids),
                ),
            )
        end
    end

    if isempty(records) &&
       hasproperty(cert_result, :lmi_data) &&
       cert_result.lmi_data !== nothing &&
       hasproperty(cert_result.lmi_data, :ellipsoids)
        raw = cert_result.lmi_data.ellipsoids
        if !isempty(raw)
            for (i, E) in enumerate(raw)
                E isa UT.Ellipsoid && push!(records, (; state_index = i, ellipsoid = E))
            end
        end
    end

    max_keep === nothing && return records
    length(records) <= max_keep && return records

    step = max(1, cld(length(records), max_keep))
    return records[1:step:end]
end

"""
Extract certified ellipsoids from a certifier result payload.
"""
function extract_ellipsoids(cert_result; max_keep = 300)
    return [
        rec.ellipsoid for rec in extract_ellipsoid_records(cert_result; max_keep = max_keep)
    ]
end

function _shift_ellipsoid_center(E::UT.Ellipsoid, shifts::AbstractVector{<:Real})
    c = collect(Float64, E.c) .+ collect(Float64, shifts)
    return UT.Ellipsoid(Matrix(E.P), c)
end

function _periodic_copy_offsets_for_ellipsoid(
    E::UT.Ellipsoid,
    periodic_dims,
    periodic_periods,
    periodic_start,
)
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

function add_periodic_ellipsoid_plot_copies(
    ellipsoids,
    periodic_dims,
    periodic_periods,
    periodic_start;
    enabled::Bool,
    mode::Symbol = :all,
)
    (!enabled || isempty(ellipsoids)) && return ellipsoids

    mode in (:all, :terminal_only, :none) || throw(
        ArgumentError(
            "periodic ellipsoid copy mode must be :all, :terminal_only, or :none",
        ),
    )
    mode == :none && return ellipsoids

    if mode == :terminal_only
        terminal_idx = lastindex(ellipsoids)
        terminal_copies = add_periodic_ellipsoid_plot_copies(
            ellipsoids[terminal_idx:terminal_idx],
            periodic_dims,
            periodic_periods,
            periodic_start;
            enabled = enabled,
            mode = :all,
        )
        return vcat(ellipsoids[begin:(terminal_idx - 1)], terminal_copies)
    end

    out = UT.Ellipsoid[]
    for E in ellipsoids
        for offset in _periodic_copy_offsets_for_ellipsoid(
            E,
            periodic_dims,
            periodic_periods,
            periodic_start,
        )
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
    trajectory_lw::Real = 2.0,
    trajectory_ms::Real = 2.0,
    trajectory_alpha::Real = 0.95,
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
                ms = 1.8,
                alpha = 0.75,
                label = "",
            )
        end
    end

    if show_start_end
        scatter!(
            fig,
            [first(xvals)],
            [first(yvals)];
            color = PLOT_COLORS[:initial_edge],
            markershape = :circle,
            ms = 4.5,
            label = "start",
        )
        scatter!(
            fig,
            [last(xvals)],
            [last(yvals)];
            color = PLOT_COLORS[:target_edge],
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
    _plot_domain_with_obstacles!(fig, domain, dims)
    initial_set !== nothing && plot!(
        fig,
        initial_set;
        dims = dims,
        color = THESIS_INITIAL_COLOR,
        linecolor = PLOT_COLORS[:initial_edge],
        opacity = 0.65,
        lw = 1.0,
        label = L"\mathcal{X}_I",
    )
    target_set !== nothing && plot!(
        fig,
        target_set;
        dims = dims,
        color = THESIS_TARGET_COLOR,
        linecolor = PLOT_COLORS[:target_edge],
        opacity = 0.72,
        lw = 1.0,
        label = L"\mathcal{X}_T",
    )

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

    raw_xs = collect(ST.enum_elems(candidate_traj.x_traj))
    xs =
        (unwrap_angles || wrap_angles) ?
        _unwrap_states_for_angle_plot(raw_xs, periodic_dims, periodic_periods) :
        [SVector{length(x), Float64}(x) for x in raw_xs]
    x_traj = ST.Trajectory(xs)

    raw_domain =
        (hasproperty(problem, :system) && hasproperty(problem.system, :X)) ?
        problem.system.X : nothing
    raw_initial_set = hasproperty(problem, :initial_set) ? problem.initial_set : nothing
    raw_target_set = hasproperty(problem, :target_set) ? problem.target_set : nothing

    domain = maybe_wrap_set(
        raw_domain,
        periodic_dims,
        periodic_periods,
        periodic_start;
        enabled = wrap_angles,
    )
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

    ellipsoid_records =
        (show_ellipsoids && cert_result !== nothing) ?
        extract_ellipsoid_records(cert_result) : NamedTuple[]
    ellipsoids = [rec.ellipsoid for rec in ellipsoid_records]

    lim12 = _plot_limits_2d(domain, (1, 2))
    fig12 = plot(;
        thesis_plot_kwargs(; legend = :topleft, size = (850, 650))...,
        aspect_ratio = :equal,
        xlims = lim12.xlims,
        ylims = lim12.ylims,
        xlabel = _state_axis_label(1),
        ylabel = _state_axis_label(2),
    )
    plot_state_space_basic!(fig12, domain, initial_set, target_set, x_traj; dims = [1, 2])
    if !isempty(ellipsoids)
        label_used = false
        for E in ellipsoids
            E2 = project_ellipsoid_2d(E; dims = (1, 2))
            plot!(
                fig12,
                E2;
                color = PLOT_COLORS[:ellipsoid],
                linecolor = PLOT_COLORS[:ellipsoid_edge],
                fillcolor = PLOT_COLORS[:ellipsoid],
                opacity = 0.54,
                lw = 1.0,
                label = label_used ? "" : L"\mathcal{E}_k",
            )
            label_used = true
        end
    end
    savefig(fig12, joinpath(output_dir, "state_space_12.pdf"))
    display(fig12)

    initial_set34 = _lift_set_for_projection(
        initial_set,
        (3, 4);
        ref_state = _reference_state_for_plot(xs, 1),
        periodic_dims = periodic_dims,
        periodic_periods = periodic_periods,
    )
    target_set34 = _lift_set_for_projection(
        target_set,
        (3, 4);
        ref_state = _reference_state_for_plot(xs, length(xs)),
        periodic_dims = periodic_dims,
        periodic_periods = periodic_periods,
    )
    ellipsoids34 = [
        _lift_ellipsoid_for_projection(
            rec.ellipsoid,
            (3, 4);
            ref_state = _reference_state_for_plot(xs, rec.state_index),
            periodic_dims = periodic_dims,
            periodic_periods = periodic_periods,
        ) for rec in ellipsoid_records
    ]
    lim34_bounds = Any[
        _trajectory_bounds_2d(xs, (3, 4)),
        _plot_bounds_2d(initial_set34, (3, 4)),
        _plot_bounds_2d(target_set34, (3, 4)),
    ]
    append!(lim34_bounds, [_projected_ellipsoid_bounds(E, (3, 4)) for E in ellipsoids34])
    lim34 = _expand_bounds_2d(_merge_bounds_2d(lim34_bounds))
    fig34 = plot(;
        thesis_plot_kwargs(; legend = :topleft, size = (850, 650))...,
        aspect_ratio = :equal,
        xlims = lim34.xlims,
        ylims = lim34.ylims,
        xlabel = _state_axis_label(3),
        ylabel = _state_axis_label(4),
    )
    plot_state_space_basic!(
        fig34,
        domain,
        initial_set34,
        target_set34,
        x_traj;
        dims = [3, 4],
    )
    if !isempty(ellipsoids34)
        label_used = false
        for E in ellipsoids34
            E2 = project_ellipsoid_2d(E; dims = (3, 4))
            plot!(
                fig34,
                E2;
                color = PLOT_COLORS[:ellipsoid],
                linecolor = PLOT_COLORS[:ellipsoid_edge],
                fillcolor = PLOT_COLORS[:ellipsoid],
                opacity = 0.24,
                lw = 1.0,
                label = label_used ? "" : L"\mathcal{E}_k",
            )
            label_used = true
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
function _animation_mp4_path_from_gif(giffile)
    giffile === nothing && return nothing
    root, ext = splitext(giffile)
    return lowercase(ext) == ".gif" ? root * ".mp4" : giffile * ".mp4"
end

function _animation_frame_indices(nstates::Int; every::Int = 1)
    step = max(1, Int(every))
    idxs = collect(1:step:nstates)
    isempty(idxs) || last(idxs) == nstates || push!(idxs, nstates)
    return idxs
end

function _normal_then_slow_indices(base_indices; slow_factor::Int = 3)
    factor = max(1, Int(slow_factor))
    slow = Int[]
    for idx in base_indices
        append!(slow, fill(idx, factor))
    end
    return vcat(base_indices, slow)
end

function _plot_animation_obstacles!(fig, obstacles2d)
    obstacles2d === nothing && return fig
    label_used = false
    for ob in obstacles2d
        x1l, y1l = ob.lb[1], ob.lb[2]
        x1u, y1u = ob.ub[1], ob.ub[2]
        xs = [x1l, x1u, x1u, x1l, x1l]
        ys = [y1l, y1l, y1u, y1u, y1l]
        plot!(
            fig,
            xs,
            ys;
            lw = 0.9,
            color = PLOT_COLORS[:obstacle],
            linecolor = PLOT_COLORS[:obstacle_edge],
            fill = (true, 0.86),
            label = label_used ? "" : L"\mathcal{X}_{obs}",
        )
        label_used = true
    end
    return fig
end

function _vehicle_animation_frame(
    vehicle_module,
    params,
    draw_params,
    states,
    inputs,
    k::Int,
    xl,
    yl;
    domain = nothing,
    obstacles2d = nothing,
    size = (760, 720),
)
    fig = plot(;
        thesis_plot_kwargs(; legend = :topright, size = size)...,
        aspect_ratio = :equal,
        xlims = xl,
        ylims = yl,
        xlabel = _state_axis_label(1),
        ylabel = _state_axis_label(2),
    )
    _plot_domain_with_obstacles!(fig, domain, [1, 2])
    _plot_animation_obstacles!(fig, obstacles2d)
    uk = k <= length(inputs) ? inputs[k] : inputs[end]
    vehicle_module.draw_articulated!(fig, params, draw_params, states[k], uk)
    return fig
end

function plot_articulated_vehicle!(
    vehicle_module,
    concrete_system,
    params,
    x_traj,
    u_traj;
    domain = concrete_system.X,
    obstacles2d = nothing,
    giffile = nothing,
    mp4file = _animation_mp4_path_from_gif(giffile),
    save_mp4::Bool = true,
    slow_motion_factor::Int = 3,
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
    states = collect(ST.enum_elems(x_traj))
    inputs = collect(ST.enum_elems(u_traj))
    isempty(states) && error("Cannot animate an empty state trajectory.")
    isempty(inputs) && error("Cannot animate an empty input trajectory.")

    frame_indices = _animation_frame_indices(length(states); every = every)

    anim_gif = @animate for k in frame_indices
        _vehicle_animation_frame(
            vehicle_module,
            params,
            draw_params,
            states,
            inputs,
            k,
            xl,
            yl;
            domain = domain,
            obstacles2d = obstacles2d,
        )
    end

    gif_path = nothing
    if giffile !== nothing
        mkpath(dirname(giffile))
        gif(anim_gif, giffile; fps = fps)
        gif_path = giffile
    end

    mp4_path = nothing
    if save_mp4 && mp4file !== nothing
        mkpath(dirname(mp4file))
        mp4_indices =
            _normal_then_slow_indices(frame_indices; slow_factor = slow_motion_factor)
        anim_mp4 = @animate for k in mp4_indices
            _vehicle_animation_frame(
                vehicle_module,
                params,
                draw_params,
                states,
                inputs,
                k,
                xl,
                yl;
                domain = domain,
                obstacles2d = obstacles2d,
            )
        end
        try
            mp4(anim_mp4, mp4file; fps = fps)
            mp4_path = mp4file
        catch err
            @warn "Could not save articulated vehicle MP4 animation." path = mp4file exception =
                err
        end
    end

    if gif_path === nothing && mp4_path === nothing
        for k in frame_indices
            display(
                _vehicle_animation_frame(
                    vehicle_module,
                    params,
                    draw_params,
                    states,
                    inputs,
                    k,
                    xl,
                    yl;
                    domain = domain,
                    obstacles2d = obstacles2d,
                ),
            )
            sleep(dt)
        end
    end

    return (; gif_path, mp4_path)
end
