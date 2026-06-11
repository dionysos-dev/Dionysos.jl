using Plots
import LinearAlgebra as LA
import StaticArrays: SVector
import Dionysos

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const DEFAULT_PERIODIC_DIMS = SVector(3, 4)
const DEFAULT_PERIODIC_PERIODS = SVector(2pi, 2pi)
const DEFAULT_PERIODIC_START = SVector(-pi, -pi)

function unwrap_periodic_state_list(state_list, periodic_dims, periodic_periods)
    isempty(state_list) && return state_list
    length(periodic_dims) == length(periodic_periods) ||
        error("periodic_dims et periodic_periods doivent avoir la meme longueur.")

    nx = length(state_list[1])
    xs = [collect(Float64, x) for x in state_list]

    for i in eachindex(periodic_dims)
        d = Int(periodic_dims[i])
        p = Float64(periodic_periods[i])
        1 <= d <= nx || error("Dimension periodique invalide: $d.")
        p > 0 || error("Periode invalide (<= 0): $p.")
        for k in 2:length(xs)
            Δ_raw = xs[k][d] - xs[k - 1][d]
            Δ_min = Δ_raw - round(Δ_raw / p) * p
            xs[k][d] = xs[k - 1][d] + Δ_min
        end
    end

    return [SVector{nx, Float64}(x) for x in xs]
end

function _wrap_vector_periodic!(
    x::AbstractVector{Float64},
    periodic_dims,
    periodic_periods,
    periodic_start,
)
    for i in eachindex(periodic_dims)
        d = Int(periodic_dims[i])
        p = Float64(periodic_periods[i])
        s = Float64(periodic_start[i])
        x[d] = mod(x[d] - s, p) + s
    end
    return x
end

function wrap_periodic_state_list(
    state_list,
    periodic_dims,
    periodic_periods,
    periodic_start,
)
    isempty(state_list) && return state_list
    length(periodic_dims) == length(periodic_periods) ||
        error("periodic_dims et periodic_periods doivent avoir la meme longueur.")
    length(periodic_dims) == length(periodic_start) ||
        error("periodic_dims et periodic_start doivent avoir la meme longueur.")

    nx = length(state_list[1])
    xs = [collect(Float64, x) for x in state_list]
    for x in xs
        _wrap_vector_periodic!(x, periodic_dims, periodic_periods, periodic_start)
    end
    return [SVector{nx, Float64}(x) for x in xs]
end

function _unwrap_ellipsoid_centers(ellipsoids, periodic_dims, periodic_periods)
    isempty(ellipsoids) && return ellipsoids

    cs = [collect(Float64, E.c) for E in ellipsoids]
    nx = length(cs[1])
    for i in eachindex(periodic_dims)
        d = Int(periodic_dims[i])
        p = Float64(periodic_periods[i])
        1 <= d <= nx || error("Dimension periodique invalide: $d.")
        p > 0 || error("Periode invalide (<= 0): $p.")
        for k in 2:length(cs)
            Δ_raw = cs[k][d] - cs[k - 1][d]
            Δ_min = Δ_raw - round(Δ_raw / p) * p
            cs[k][d] = cs[k - 1][d] + Δ_min
        end
    end

    out = UT.Ellipsoid[]
    for (i, E) in enumerate(ellipsoids)
        push!(out, UT.Ellipsoid(Matrix(E.P), cs[i]))
    end
    return out
end

function _wrap_ellipsoid_centers(
    ellipsoids,
    periodic_dims,
    periodic_periods,
    periodic_start,
)
    isempty(ellipsoids) && return ellipsoids
    out = UT.Ellipsoid[]
    for E in ellipsoids
        c = collect(Float64, E.c)
        _wrap_vector_periodic!(c, periodic_dims, periodic_periods, periodic_start)
        push!(out, UT.Ellipsoid(Matrix(E.P), c))
    end
    return out
end

function _maybe_wrap_set(
    set,
    periodic_dims,
    periodic_periods,
    periodic_start;
    enabled::Bool,
)
    set === nothing && return nothing
    enabled || return set
    try
        return UT.set_in_period(set, periodic_dims, periodic_periods, periodic_start)
    catch
        return set
    end
end

function projeter_ellipsoide_2d(E::UT.Ellipsoid; dims = (1, 2))
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

    length(ells) <= max_keep && return ells
    step = max(1, cld(length(ells), max_keep))
    return ells[1:step:end]
end

function plot_state_space_basic!(fig, domain, I, T, x_traj; dims = [1, 2])
    domain !== nothing &&
        plot!(fig, domain; dims = dims, color = :grey, opacity = 1.0, label = "")
    I !== nothing && plot!(fig, I; dims = dims, color = :green, opacity = 0.2, label = "I")
    T !== nothing && plot!(fig, T; dims = dims, color = :red, opacity = 0.5, label = "T")
    return plot!(
        fig,
        x_traj;
        dims = dims,
        ms = 2.0,
        arrows = false,
        color = :blue,
        label = "traj",
    )
end

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
        (hasproperty(problem, :system) && hasproperty(problem.system, :X)) ?
        problem.system.X : nothing
    raw_I = hasproperty(problem, :initial_set) ? problem.initial_set : nothing
    raw_T = hasproperty(problem, :target_set) ? problem.target_set : nothing
    domain = _maybe_wrap_set(
        raw_domain,
        periodic_dims,
        periodic_periods,
        periodic_start;
        enabled = wrap_angles,
    )
    I = _maybe_wrap_set(
        raw_I,
        periodic_dims,
        periodic_periods,
        periodic_start;
        enabled = wrap_angles,
    )
    T = _maybe_wrap_set(
        raw_T,
        periodic_dims,
        periodic_periods,
        periodic_start;
        enabled = wrap_angles,
    )

    ellipsoids =
        (show_ellipsoids && cert_result !== nothing) ? extract_ellipsoids(cert_result) :
        UT.Ellipsoid[]
    if !isempty(ellipsoids) && unwrap_angles
        ellipsoids = _unwrap_ellipsoid_centers(ellipsoids, periodic_dims, periodic_periods)
    end
    if !isempty(ellipsoids) && wrap_angles
        ellipsoids = _wrap_ellipsoid_centers(
            ellipsoids,
            periodic_dims,
            periodic_periods,
            periodic_start,
        )
    end

    dims = [1, 2]
    fig = plot(; aspect_ratio = :equal, legend = false, title = title12)
    plot_state_space_basic!(fig, domain, I, T, x_traj; dims = dims)
    if !isempty(ellipsoids)
        for E in ellipsoids
            E2 = projeter_ellipsoide_2d(E; dims = (1, 2))
            plot!(fig, E2; color = :orange, opacity = 0.22, lw = 1.0, label = "")
        end
    end
    savefig(fig, joinpath(output_dir, "state_space_12.pdf"))
    display(fig)

    dims = [3, 4]
    fig = plot(; aspect_ratio = :equal, legend = false, title = title34)
    plot_state_space_basic!(fig, domain, I, T, x_traj; dims = dims)
    if !isempty(ellipsoids)
        for E in ellipsoids
            E2 = projeter_ellipsoide_2d(E; dims = (3, 4))
            plot!(fig, E2; color = :orange, opacity = 0.22, lw = 1.0, label = "")
        end
    end
    savefig(fig, joinpath(output_dir, "state_space_34.pdf"))
    display(fig)

    return nothing
end

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
    xl = (-20.0, 20.0)
    yl = (-20.0, 20.0)
    dp = vehicle_module.DrawParams(params)

    return vehicle_module.live_vehicle_progression(
        params,
        dp,
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
