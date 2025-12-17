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

function get_abstraction_based_controller(concrete_system)
    empty_problem = DI.Problem.EmptyProblem(concrete_system, concrete_system.X)

    # --- State grid ---
    x0 = SVector(0.0, 0.0, 0.0, 0.0)
    hx = SVector(0.4, 0.4, 10*(pi/180.0), 3*(pi/180.0))   # grid steps (θ,ϕ steps in radians) 0.4 0.15

    # --- Input grid ---
    u0 = SVector(0.0, 0.0)
    hu = SVector(1.0, 0.5)
    input_grid = DO.GridFree(u0, hu)

    with_period = true
    periodic_dims = SVector(3, 4)
    periodic_periods = SVector(2pi, 2pi)
    periodic_start   = SVector(-pi, -pi)

    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), empty_problem)
    # MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("h"), hx)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), AV.jacobian_bound(params))
    MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.2)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("approx_mode"), AB.UniformGridAbstraction.GROWTH) #GROWTH CENTER_SIMULATION

    MOI.set(optimizer, MOI.RawOptimizerAttribute("threaded"), true)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("efficient"), true)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("n_samples"), 1)

    MOI.set(optimizer, MOI.Silent(), true)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)

    if with_period
        MOI.set(optimizer, MOI.RawOptimizerAttribute("use_periodic_domain"), true)
        MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_dims"), periodic_dims)
        MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_periods"), periodic_periods)
        MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_start"), periodic_start)
    end

    MOI.optimize!(optimizer)
    println("Abstraction time = ",
        MOI.get(optimizer, MOI.RawOptimizerAttribute("abstraction_construction_time_sec"))
    )

    # Reachability
    _I_ = UT.HyperRectangle(SVector(-2.0, -2.0, -0.4, -0.4), SVector(2.0, 2.0, 0.4, 0.4))
    # _T_ = UT.HyperRectangle(SVector( 6.0,  6.0, -0.4, -0.4), SVector(9.0, 9.0, 0.4, 0.4))
    # _T_ = UT.HyperRectangle(SVector( 6.0,  12.0, pi/2-20*(pi/180), -pi), SVector(10.0, 15.0, pi/2+20*(pi/180), pi))
    _T_ = UT.HyperRectangle(SVector( 2.0,  6.0, pi-20*(pi/180), -10*(pi/180)), SVector(4.0, 8.0, pi+20*(pi/180), 10*(pi/180)))

    concrete_problem = DI.Problem.OptimalControlProblem(
        concrete_system, _I_, _T_, nothing, nothing, DI.Problem.Infinity(),
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("early_stop"), false)

    MOI.optimize!(optimizer)
    println("Abstract solve time = ",
        MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem_time_sec"))
    )

    abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))

    dims = [3, 4]
    plot_domain(concrete_problem, abstract_problem, periodic_dims, periodic_periods, periodic_start, dims)
    concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
    return concrete_controller, concrete_problem.target_set
end

function xy_obstacles()
    return [
        UT.HyperRectangle(SVector(-5.0, -2.0), SVector(-2.0,  2.0)),
        UT.HyperRectangle(SVector( 3.0,  3.0), SVector( 6.0,  6.0)),
    ]
end

function plot_domain(concrete_problem, abstract_problem, periodic_dims, periodic_periods, periodic_start, dims)
    concrete_system = concrete_problem.system
    abstract_system = abstract_problem.system
    fig = plot(; aspect_ratio = :equal, legend = false)
    state_space = concrete_system.X
    system_domain_in_periodic = UT.set_in_period(state_space, periodic_dims, periodic_periods, periodic_start)
    initial_set_in_periodic = UT.set_in_period(concrete_problem.initial_set, periodic_dims, periodic_periods, periodic_start)
    target_set_in_periodic =UT.set_in_period(concrete_problem.target_set, periodic_dims, periodic_periods, periodic_start)
    plot!(system_domain_in_periodic; dims = dims, color = :grey, opacity = 1.0, label = "")
    plot!(abstract_system.Xdom; dims=dims, color=:blue, efficient = false)
    plot!(initial_set_in_periodic; dims = dims, color = :green, opacity = 0.2, label = "Initial set")
    plot!(target_set_in_periodic; dims = dims, color = :red, opacity = 0.5, label = "Target set")
    plot!(
        Dionysos.Symbolic.get_domain_from_states(abstract_system, abstract_problem.initial_set);
        dims = dims, 
        color = :green,
        efficient = false
    )
    plot!(
        Dionysos.Symbolic.get_domain_from_states(abstract_system, abstract_problem.target_set);
        dims = dims, 
        color = :red,
        efficient = false
    )

    # Closed-loop simulate
    Δt = 0.2
    discrete_time_system = ST.discretize_continuous_system(
        concrete_system,
        Δt;
        num_substeps = 5,
    )
    nstep = 700
    x0 = SVector(0.0, 0.0, 0.0, 0.0)
    periodic_wrapper = ST.get_periodic_wrapper(SVector(3,4), SVector(2pi,2pi); start=SVector(-pi,-pi))
    reached(x) = (x ∈ target_set_in_periodic)
    traj = ST.get_closed_loop_trajectory(
        discrete_time_system,
        controller,
        x0,
        nstep;
        stopping = reached,
        periodic_wrapper = periodic_wrapper
    )
    plot!(traj; dims = dims, ms = 2.0, arrows = false)
    display(fig)
end 


function simulate(concrete_system, controller, target_set)
    Δt = 0.2
    discrete_time_system = ST.discretize_continuous_system(
        concrete_system,
        Δt;
        num_substeps = 5,
    )
    # Closed-loop simulate
    nstep = 700
    x0 = SVector(0.0, 0.0, 0.0, 0.0)
    periodic_wrapper = ST.get_periodic_wrapper(SVector(3,4), SVector(2pi,2pi); start=SVector(-pi,-pi))
    target_set_in_periodic = UT.set_in_period(target_set, SVector(3,4), SVector(2pi,2pi), SVector(-pi,-pi))
    reached(x) = (x ∈ target_set_in_periodic)
    traj = ST.get_closed_loop_trajectory(
        discrete_time_system,
        controller,
        x0,
        nstep;
        stopping = reached,
        periodic_wrapper = periodic_wrapper
    )

    gr()
    xl = (-20.0, 20.0)
    yl = (-20.0, 20.0)
    # AV.live_vehicle_progression!(params, traj, xl, yl; every=1, dt=0.008)
    dp = AV.DrawParams(params)           # derived from L1,L2,Lc
    dp2 = AV.DrawParams(params; tractor_width=1.8)  # override if you want
    domain = concrete_system.X
    AV.live_vehicle_progression_pretty(params, dp, traj, xl, yl; domain=domain, every=1, dt=0.2)
end

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "articulated_vehicle.jl"))
AV = ArticulatedVehicle

# --- Domains (example; tune to your benchmark) ---
# _X_ = UT.HyperRectangle(
#     SVector(-1.0, -1.0, -pi, -pi),   # x1, x2, θ1, ϕ
#     SVector( 15.0,  15.0,  pi,  pi),
# )

_X_ = UT.HyperRectangle(
    SVector(-1.0, -1.0, -pi, -pi),   # x1, x2, θ1, ϕ
    SVector( 10.0,  9.0,  pi,  pi),
)
_X_ = AV.with_phi_limit(_X_; phi_max=90*(pi/180.0))   # e.g. 34°
# _X_ = with_xy_obstacles(_X_; obstacles2d=xy_obstacles())

_U_ = UT.HyperRectangle(
    SVector(-2.0, -0.6),               # v1, δ
    SVector( 2.0,  0.6),
)

params = AV.Params(L1=1.0, L2=1.0, Lc=0.5)
concrete_system = AV.system(_X_; _U_=_U_, params=params)

# controller = AV.get_constant_controller(SVector(1.0, 0.15))
# controller = AV.get_goal_seeking_controller(-5.0, 10; v=1.0, δmax=0.5, k=1.2)
controller, target_set = get_abstraction_based_controller(concrete_system)

simulate(concrete_system, controller, target_set)
