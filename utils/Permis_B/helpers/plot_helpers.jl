"""
    plot_state_space!(optimizer, concrete_system, _I_, _T_, x_traj; ...)

Superpose sur une meme figure:
- le domaine continu,
- le domaine abstrait (grille),
- les ensembles initial/objectif,
- la trajectoire fermee.
"""
function plot_state_space!(
    optimizer,
    concrete_system,
    _I_,
    _T_,
    x_traj;
    dims = [1, 2],
    with_period = false,
    periodic_dims = SVector{0, Int}(),
    periodic_periods = SVector{0, Float64}(),
    periodic_start = SVector{0, Float64}(),
)
    abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
    abstract_system = abstract_problem.system

    X = concrete_system.X
    Xp =
        with_period ? UT.set_in_period(X, periodic_dims, periodic_periods, periodic_start) :
        X
    Ip =
        with_period ?
        UT.set_in_period(_I_, periodic_dims, periodic_periods, periodic_start) : _I_
    Tp =
        with_period ?
        UT.set_in_period(_T_, periodic_dims, periodic_periods, periodic_start) : _T_

    plot!(Xp; dims = dims, color = :grey, opacity = 1.0, label = "")
    plot!(
        abstract_system.Xdom;
        dims = dims,
        color = :blue,
        opacity = 0.2,
        efficient = false,
    )
    plot!(Ip; dims = dims, color = :green, opacity = 0.2, label = "I")
    plot!(Tp; dims = dims, color = :red, opacity = 0.5, label = "T")

    return plot!(x_traj; dims = dims, ms = 2.0, arrows = false)
end

"""
    save_state_space_plots!(output_dir, optimizer, concrete_system, _I_, _T_, x_traj; ...)

Sauvegarde les projections `(x, y)` et `(theta, phi)` en PDF.
"""
function save_state_space_plots!(
    output_dir,
    optimizer,
    concrete_system,
    _I_,
    _T_,
    x_traj;
    with_period = false,
    periodic_dims = SVector{0, Int}(),
    periodic_periods = SVector{0, Float64}(),
    periodic_start = SVector{0, Float64}(),
)
    dims = [1, 2]
    fig = plot(; aspect_ratio = :equal, legend = false)
    plot_state_space!(
        optimizer,
        concrete_system,
        _I_,
        _T_,
        x_traj;
        dims = dims,
        with_period = with_period,
        periodic_dims = periodic_dims,
        periodic_periods = periodic_periods,
        periodic_start = periodic_start,
    )
    savefig(fig, joinpath(output_dir, "state_space_12.pdf"))
    display(fig)

    dims = [3, 4]
    fig = plot(; aspect_ratio = :equal, legend = false)
    plot_state_space!(
        optimizer,
        concrete_system,
        _I_,
        _T_,
        x_traj;
        dims = dims,
        with_period = with_period,
        periodic_dims = periodic_dims,
        periodic_periods = periodic_periods,
        periodic_start = periodic_start,
    )
    savefig(fig, joinpath(output_dir, "state_space_34.pdf"))
    display(fig)
    return nothing
end

"""
    plot_articulated_vehicle!(vehicle_module, concrete_system, params, x_traj, u_traj; ...)

Animation de la trajectoire du vehicule articule.
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
