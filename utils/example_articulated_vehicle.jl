using Plots
import StaticArrays: SVector
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const OP = DI.Optim
const AB = OP.Abstraction
using JuMP
import MathOptInterface as MOI

function build_uniform_grid_abstraction(
    concrete_system,
    Δt,
    hx,
    Udom,
    jacobian_bound;
    approx_mode = AB.UniformGridAbstraction.CENTER_SIMULATION, # GROWTH, CENTER_SIMULATION
    with_period = false,
    periodic_dims = SVector{0, Int}(),
    periodic_periods = SVector{0, Float64}(),
    periodic_start = SVector{0, Float64}(),
)
    empty_problem = DI.Problem.EmptyProblem(concrete_system, concrete_system.X)

    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), empty_problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("h"), hx)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("Udom"), Udom)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), jacobian_bound)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), Δt)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("approx_mode"), approx_mode)

    MOI.set(optimizer, MOI.RawOptimizerAttribute("threaded"), true)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("efficient"), true)
    MOI.set(optimizer, MOI.Silent(), true)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)

    if with_period
        MOI.set(optimizer, MOI.RawOptimizerAttribute("use_periodic_domain"), true)
        MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_dims"), periodic_dims)
        MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_periods"), periodic_periods)
        MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_start"), periodic_start)
    end

    MOI.optimize!(optimizer)
    return optimizer
end

function build_uniform_grid_controller!(optimizer, concrete_system, _I_, _T_)
    concrete_problem = DI.Problem.OptimalControlProblem(
        concrete_system,
        _I_,
        _T_,
        nothing,
        nothing,
        DI.Problem.Infinity(),
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("early_stop"), false)

    MOI.optimize!(optimizer)

    concrete_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
    target_set = concrete_problem.target_set
    return concrete_controller, target_set
end

function simulate_closed_loop(
    concrete_system,
    concrete_controller,
    Δt,
    x0,
    target_set;
    nstep = 700,
    with_period = false,
    periodic_dims = SVector{0, Int}(),
    periodic_periods = SVector{0, Float64}(),
    periodic_start = SVector{0, Float64}(),
)
    disc = ST.discretize_continuous_system(concrete_system, Δt; num_substeps = 5)

    periodic_wrapper =
        with_period ?
        ST.get_periodic_wrapper(periodic_dims, periodic_periods; start = periodic_start) :
        (x -> x)

    reached(x) = (periodic_wrapper(x) ∈ target_set)

    traj = ST.get_closed_loop_trajectory(
        disc,
        concrete_controller,
        x0,
        nstep;
        stopping = reached,
        periodic_wrapper = periodic_wrapper,
    )
    return traj
end

function plot_state_space!(
    optimizer,
    concrete_system,
    _I_,
    _T_,
    traj;
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

    return plot!(traj; dims = dims, ms = 2.0, arrows = false)
end

function plot_articulated_vehicle!(
    concrete_system,
    params,
    traj;
    domain = concrete_system.X,
    giffile = nothing,
    fps = 20,
    every = 1,
    dt = 0.2,
)
    gr()
    xl = (-20.0, 20.0)
    yl = (-20.0, 20.0)
    dp = AV.DrawParams(params)

    return AV.live_vehicle_progression(
        params,
        dp,
        traj,
        xl,
        yl;
        domain = domain,
        every = every,
        dt = dt,
        giffile = giffile,
        fps = fps,
    )
end

function get_Udom(_U_, hu)
    u0 = zeros(SVector{2, Float64})
    input_grid = DO.GridFree(u0, hu)
    Udom = Dionysos.Domain.DomainList(input_grid)
    Dionysos.Domain.add_set!(Udom, _U_, Dionysos.Domain.CENTER)
    return Udom
end

function script()
    ### System ###
    _X_ = UT.HyperRectangle(
        SVector(-1.0, -1.0, -pi, -pi), # x1, x2, θ1, ϕ
        SVector(10.0, 9.0, pi, pi),
    )
    _X_ = AV.with_phi_limit(_X_; phi_max = 50*(pi/180.0))
    obs = [
        UT.HyperRectangle(SVector(4.0, -1.0), SVector(10.0, 4.7)),
        UT.HyperRectangle(SVector(4.0, 6.0), SVector(10.0, 9.0)),
    ]
    _X_ = AV.with_xy_obstacles(_X_; obstacles2d = obs)

    _U_ = UT.HyperRectangle(
        SVector(-2.0, -0.6), # v1, δ
        SVector(2.0, 0.6),
    )

    #CustomList([SVector(),])
    params = AV.Params(; L1 = 1.0, L2 = 1.0, Lc = 0.5)
    concrete_system = AV.system(_X_; _U_ = _U_, params = params)

    ### Abstraction params ###
    Δt = 0.2
    hx = SVector(0.4, 0.2, 5*(pi/180), 3*(pi/180))
    periodic_dims = SVector(3, 4)
    periodic_periods = SVector(2pi, 2pi)
    periodic_start = SVector(-pi, -pi)

    # Udom = get_Udom(_U_,  SVector(1.0, 0.5))

    inputs = [
        [2.0, 0.0],
        [0.0, 0.0],
        [-2.0, 0.0],
        [2.0, -0.25],
        [2.0, 0.25],
        [-2.0, 0.25],
        [-2.0, -0.25],
    ]
    Udom = Dionysos.Domain.CustomList(inputs)

    ### Control problem###
    x0 = SVector(0.0, 0.0, 0.0, 0.0)
    _I_ = UT.HyperRectangle(SVector(-1.0, -1.0, -0.4, -0.4), SVector(1.0, 1.0, 0.4, 0.4))
    _T_ = UT.HyperRectangle(
        SVector(9.0, 5.0, -5*(pi/180), -5*(pi/180)),
        SVector(10.0, 6.0, 5*(pi/180), 5*(pi/180)),
    ) # forward
    # _T_ = UT.HyperRectangle(SVector(9.0,5.0, pi-5*(pi/180), -5*(pi/180)),
    #                         SVector(10.0,6.0, pi+5*(pi/180),  5*(pi/180))) # backward

    opt = build_uniform_grid_abstraction(
        concrete_system,
        Δt,
        hx,
        Udom,
        AV.jacobian_bound(params);
        with_period = true,
        periodic_dims = periodic_dims,
        periodic_periods = periodic_periods,
        periodic_start = periodic_start,
    )

    controller, target_set = build_uniform_grid_controller!(opt, concrete_system, _I_, _T_)
    traj = simulate_closed_loop(
        concrete_system,
        controller,
        Δt,
        x0,
        target_set;
        with_period = true,
        periodic_dims = periodic_dims,
        periodic_periods = periodic_periods,
        periodic_start = periodic_start,
        nstep = 700,
    )

    dims=[1, 2]
    fig = plot(; aspect_ratio = :equal, legend = false)
    plot_state_space!(
        opt,
        concrete_system,
        _I_,
        _T_,
        traj;
        dims = dims,
        with_period = true,
        periodic_dims = periodic_dims,
        periodic_periods = periodic_periods,
        periodic_start = periodic_start,
    )
    savefig(fig, "state_space_12.pdf")
    display(fig)
    dims=[3, 4]
    fig = plot(; aspect_ratio = :equal, legend = false)
    plot_state_space!(
        opt,
        concrete_system,
        _I_,
        _T_,
        traj;
        dims = dims,
        with_period = true,
        periodic_dims = periodic_dims,
        periodic_periods = periodic_periods,
        periodic_start = periodic_start,
    )
    savefig(fig, "state_space_34.pdf")
    display(fig)
    return plot_articulated_vehicle!(concrete_system, params, traj; every = 1, dt = 0.09)
    #plot_articulated_vehicle!(concrete_system, params, traj; giffile="articulated_vehicle.gif",fps=5,every=3) 
end

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "articulated_vehicle.jl"))
AV = ArticulatedVehicle

script()

# controller = AV.get_constant_controller(SVector(1.0, 0.15))
# controller = AV.get_goal_seeking_controller(-5.0, 10; v=1.0, δmax=0.5, k=1.2)
