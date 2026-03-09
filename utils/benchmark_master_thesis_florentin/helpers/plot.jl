using Plots
import LinearAlgebra as LA
import StaticArrays: SVector
import Dionysos

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System

# Load periodic helpers once when plot.jl is used standalone.
if !isdefined(@__MODULE__, :unwrap_periodic_state_list)
    include(joinpath(@__DIR__, "helpers.jl"))
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

    length(ells) <= max_keep && return ells

    step = max(1, cld(length(ells), max_keep))
    return ells[1:step:end]
end

"""
Draw domain/initial/target + trajectory in selected dimensions.
"""
function plot_state_space_basic!(fig, domain, initial_set, target_set, x_traj; dims = [1, 2])
    domain !== nothing && plot!(fig, domain; dims = dims, color = :grey, opacity = 1.0, label = "")
    initial_set !== nothing &&
        plot!(fig, initial_set; dims = dims, color = :green, opacity = 0.2, label = "I")
    target_set !== nothing &&
        plot!(fig, target_set; dims = dims, color = :red, opacity = 0.5, label = "T")

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
    end

    fig12 = plot(; aspect_ratio = :equal, legend = false, title = title12)
    plot_state_space_basic!(fig12, domain, initial_set, target_set, x_traj; dims = [1, 2])
    if !isempty(ellipsoids)
        for E in ellipsoids
            E2 = project_ellipsoid_2d(E; dims = (1, 2))
            plot!(fig12, E2; color = :orange, opacity = 0.22, lw = 1.0, label = "")
        end
    end
    savefig(fig12, joinpath(output_dir, "state_space_12.pdf"))
    display(fig12)

    fig34 = plot(; aspect_ratio = :equal, legend = false, title = title34)
    plot_state_space_basic!(fig34, domain, initial_set, target_set, x_traj; dims = [3, 4])
    if !isempty(ellipsoids)
        for E in ellipsoids
            E2 = project_ellipsoid_2d(E; dims = (3, 4))
            plot!(fig34, E2; color = :orange, opacity = 0.22, lw = 1.0, label = "")
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
    xl = (-20.0, 20.0)
    yl = (-20.0, 20.0)
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
