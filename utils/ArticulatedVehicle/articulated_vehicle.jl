using Plots
import StaticArrays: SVector
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction
using JuMP
import MathOptInterface as MOI

using JLD2

# --------------------------------------- #
# --- JLD2 Abstraction Export/Import ---- #
# --------------------------------------- #

using JLD2

function export_abstraction_jld2(opt, filename::AbstractString)
    abs_opt = opt.abstraction_solver
    abs_opt === nothing && error("No abstraction_solver in optimizer.")
    abs_sys = abs_opt.abstract_system
    abs_sys === nothing && error("No abstract_system computed yet.")

    jldopen(filename, "w") do f
        # versioning for forward compatibility
        f["format_version"] = 1
        f["abstract_system"] = abs_sys
        return f["params"] = (time_step = opt.abstraction_solver.time_step,)
    end
    return nothing
end

function import_abstraction_jld2(filename::AbstractString; opt = nothing)
    # If user didn't pass an optimizer, create one
    if opt === nothing
        opt = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    end

    # Ensure abstraction solver exists
    if opt.abstraction_solver === nothing
        opt.abstraction_solver =
            AB.UniformGridAbstraction.OptimizerAlternatingSimulationProblem()
    end

    jldopen(filename, "r") do f
        v = f["format_version"]
        v == 1 || error("Unsupported abstraction file format_version=$v")

        abs_sys = f["abstract_system"]
        return MOI.set(
            opt.abstraction_solver,
            MOI.RawOptimizerAttribute("abstract_system"),
            abs_sys,
        )
    end

    return opt
end

# --------------------------------------- #
# ------------- Builders ---------------- #
# --------------------------------------- #

function build_optimizer(
    concrete_system,
    Δt,
    hx,
    UMapping,
    jacobian_bound;
    periodic_dims = SVector{0, Int}(),
    periodic_periods = SVector{0, Float64}(),
    periodic_start = SVector{0, Float64}(),
    with_period = false,
    approx_mode = AB.UniformGridAbstraction.CENTER_SIMULATION, # GROWTH, CENTER_SIMULATION
)
    alternating_simulation_problem =
        DI.Problem.AlternatingSimulationProblem(concrete_system, concrete_system.X)

    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("concrete_problem"),
        alternating_simulation_problem,
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("h"), hx)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("UMapping"), UMapping)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), jacobian_bound)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), Δt)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("approx_mode"), approx_mode)

    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("execution_backend"),
        SY.ThreadedBackend(0.2),
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("efficient"), true)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("progress_update_interval"), Int(1e5))

    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("automaton_constructor"),
        (n, m) -> ST.NewFastIndexedAutomatonList(n, m),
    )

    if with_period
        MOI.set(optimizer, MOI.RawOptimizerAttribute("use_periodic_mapping"), true)
        MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_dims"), periodic_dims)
        MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_periods"), periodic_periods)
        MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_start"), periodic_start)
    end

    return optimizer
end

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
    abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    XMapping = SY.get_state_mapping(abstract_system)
    Xset = SY.get_state_domain(abstract_system)

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
    plot!((Xset, XMapping); dims = dims, color = :blue, opacity = 0.2, efficient = false)
    plot!(Ip; dims = dims, color = :green, opacity = 0.2, label = "I")
    plot!(Tp; dims = dims, color = :red, opacity = 0.5, label = "T")

    return plot!(x_traj; dims = dims, ms = 2.0, arrows = false)
end

function script()
    compute_abstraction = true
    save = false
    load = false
    filename = joinpath(@__DIR__, "Abstraction.jld2")

    # ------------------------------------------------------------
    # System 
    # ------------------------------------------------------------
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

    params = AV.Params(; L1 = 1.0, L2 = 1.0, Lc = 0.5)
    concrete_system = AV.system(_X_; _U_ = _U_, params = params)

    # ------------------------------------------------------------
    # Problem 
    # ------------------------------------------------------------

    _I_ = UT.HyperRectangle(SVector(-0.2, -0.2, -0.2, -0.2), SVector(0.2, 0.2, 0.2, 0.2))
    _T_ = UT.HyperRectangle(
        SVector(9.0, 5.5, -5*(pi/180), -5*(pi/180)),
        SVector(10.0, 6.0, 5*(pi/180), 5*(pi/180)),
    ) # forward
    _T_ = UT.HyperRectangle(
        SVector(9.0, 5.0, pi-5*(pi/180), -5*(pi/180)),
        SVector(10.0, 5.5, pi+5*(pi/180), 5*(pi/180)),
    ) # backward

    concrete_problem = DI.Problem.OptimalControlProblem(
        concrete_system,
        _I_,
        _T_,
        nothing,
        nothing,
        DI.Problem.Infinity(),
    )

    # ------------------------------------------------------------
    # Optimizer Parameters 
    # ------------------------------------------------------------

    periodic_dims = SVector(3, 4)
    periodic_periods = SVector(2pi, 2pi)
    periodic_start = SVector(-pi, -pi)
    with_period = true
    Δt = 0.2

    if (compute_abstraction)
        #----- Abstraction settings -----
        hx = SVector(0.4, 0.2, 5*(pi/180), 3*(pi/180))
        inputs = [
            [2.0, 0.0],
            [0.0, 0.0],
            [-2.0, 0.0],
            [2.0, -0.25],
            [2.0, 0.25],
            [-2.0, 0.25],
            [-2.0, -0.25],
        ]
        UMapping = MP.ListMapping(inputs)

        optimizer = build_optimizer(
            concrete_system,
            Δt,
            hx,
            UMapping,
            AV.jacobian_bound(params);
            periodic_dims = periodic_dims,
            periodic_periods = periodic_periods,
            periodic_start = periodic_start,
            with_period = with_period,
            approx_mode = AB.UniformGridAbstraction.CENTER_SIMULATION, # GROWTH, CENTER_SIMULATION
        )
        MOI.optimize!(optimizer)
        if (save)
            export_abstraction_jld2(optimizer, filename)
        end
    end
    if (load)
        optimizer = import_abstraction_jld2(filename)
    end

    # ------------------------------------------------------------
    # Solve concrete problem
    # ------------------------------------------------------------

    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("early_stop"), false)
    MOI.optimize!(optimizer)

    concrete_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
    target_set = concrete_problem.target_set

    # ------------------------------------------------------------
    # Simulate closed-loop trajectory
    # ------------------------------------------------------------

    x0 = SVector(0.0, 0.0, 0.0, 0.0)
    discrete_time_system =
        ST.discretize_continuous_system(concrete_system, Δt; num_substeps = 5)
    periodic_wrapper =
        with_period ?
        ST.get_periodic_wrapper(periodic_dims, periodic_periods; start = periodic_start) :
        (x -> x)
    reached(x) = (periodic_wrapper(x) ∈ target_set)
    nstep = 100

    x_traj, u_traj = ST.get_closed_loop_trajectory(
        discrete_time_system,
        concrete_controller,
        x0,
        nstep;
        stopping = reached,
        wrap = periodic_wrapper,
    )

    # ------------------------------------------------------------
    # Plot 
    # ------------------------------------------------------------

    dims = [1, 2]
    fig = plot(; aspect_ratio = :equal, legend = false)
    plot_state_space!(
        optimizer,
        concrete_system,
        _I_,
        _T_,
        x_traj;
        dims = dims,
        with_period = true,
        periodic_dims = periodic_dims,
        periodic_periods = periodic_periods,
        periodic_start = periodic_start,
    )
    # savefig(fig, "state_space_12.pdf")
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
        with_period = true,
        periodic_dims = periodic_dims,
        periodic_periods = periodic_periods,
        periodic_start = periodic_start,
    )
    # savefig(fig, "state_space_34.pdf")
    display(fig)

    # ------------------------------------------------------------
    # Animation with dashboard
    # ------------------------------------------------------------

    draw_params = AV.DrawParams(params)
    system_plot! = ArticulatedVehicle.system_plot!(;
        params = params,
        draw_params = draw_params,
        xlims = (-1.0, 10.0),
        ylims = (-1.0, 10.0),
        obstacles2d = obs,
    )

    return Dionysos.animate_trajectory_dashboard(
        system_plot!,
        x_traj,
        u_traj;
        xdims = (1, 2),        # x-y trajectory
        udims = (1,),          # velocity over time
        Δt = Δt,
        fps = 5,
        # filename = "articulated_vehicle_dashboard.mp4",
        xlabel_state = "x",
        ylabel_state = "y",
        xlabel_input = "time [s]",
        ylabel_input = "v",
        xlims_state = (-1.0, 10.0),
        ylims_state = (-1.0, 10.0),
    )
end

include("../../problems/articulated_vehicle.jl");
AV = ArticulatedVehicle

script()

# controller = AV.get_constant_controller(SVector(1.0, 0.15))
# controller = AV.get_goal_seeking_controller(-5.0, 10; v=1.0, δmax=0.5, k=1.2)
