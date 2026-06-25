using StaticArrays
using MathematicalSystems
using Dionysos
using JuMP
using LinearAlgebra
using JLD2
import MathOptInterface as MOI

include("./robot_vcontrol.jl")
import .RobotVelocity as RV

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction
const OPDS = OP.DiscreteSystems

#----------------------
# Save-load settings 
#----------------------
save = false
load = false
animate = true
filename = "./problems/BipedRobot/4D_model_vcontrol/optimizer.jld2"

# ------------------------------------------------------------
# 1) Define a simple 4D continuous-time system: x' = u
# ------------------------------------------------------------

robot_geometry = RV.Geometry.geometry(0.0, 0.202, 0.172)

if load
    optimizer = JLD2.load(filename, "optimizer")
    concrete_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_problem"))
    _I_ = concrete_problem.initial_set
    _T_ = concrete_problem.target_set
    abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    discrete_time_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system"))
    success = MOI.get(optimizer, MOI.RawOptimizerAttribute("success"))
else
    # ------------------------------------------------------------
    # 2) Abstraction construction (AlternatingSimulationProblem)
    # ------------------------------------------------------------

    # State space
    x_bar = 1.0
    full_X_ = UT.HyperRectangle(
        -SVector(x_bar, x_bar, x_bar, x_bar),
        SVector(x_bar, x_bar, x_bar, x_bar)
    )

    # State grid resolution (hx defines TOTAL CELL WIDTH)
    dx = 0.1
    x0 = SVector(0.0, 0.0, 0.0, 0.0)
    hx = SVector(dx, dx, dx, dx)
    state_grid = MP.GridFree(x0, hx)

    # State constraints
    _X_ = RV.RobotModel.remove_infeasible_cells(
        full_X_,
        state_grid,
        UT.HyperRectangle(SVector(-0.03, -1.0), SVector(0.03,0.03)),
        robot_geometry, true
    )
    # _X_ = full_X_

    # Input space
    u_bar = 10*dx
    _U_ = UT.HyperRectangle(
        -SVector(u_bar, u_bar, u_bar, u_bar),
        SVector(u_bar, u_bar, u_bar, u_bar)
    )

    # Input grid resolution
    du = 10*dx
    u0 = SVector(0.0, 0.0, 0.0, 0.0)
    hu = SVector(du, du, du, du)
    input_grid = MP.GridFree(u0, hu)

    # Initial set
    _I_ = UT.HyperRectangle(
        SVector(0.2, 0.0, -0.2, 0.0) - hx/2,
        SVector(0.2, 0.0, -0.2, 0.0) + hx/2
    )

    # Target set
    _T_ = UT.HyperRectangle(
        SVector(-x_bar, -x_bar, 0.2, -x_bar),
        SVector(-0.2, 0.0, x_bar, 0.0)
    )
    # _T_ = compute_target_set(_X_, state_grid, SVector(0.4, 0.0), geometry, true)

    concrete_system = RV.RobotModel.system(; _X_ = _X_, _U_ = _U_)
    jacobian_bound = RV.RobotModel.jacobian_bound()

    alternating_simulation_problem =
        DI.Problem.AlternatingSimulationProblem(concrete_system, concrete_system.X)

    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("concrete_problem"),
        alternating_simulation_problem,
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.1)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("execution_backend"), SY.ThreadedBackend(0.2))

    # choose an approx mode that exists in your setup
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("approx_mode"),
        AB.UniformGridAbstraction.CENTER_SIMULATION,
    ) # GROWTH CENTER_SIMULATION
    MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), jacobian_bound)

    MOI.set(optimizer, MOI.RawOptimizerAttribute("efficient"), true)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("use_implicit_mapping"), false)
    #MOI.set(optimizer, MOI.RawOptimizerAttribute("mapping_region"), full_X_)

    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)
    MOI.set(optimizer, MOI.Silent(), true)

    MOI.optimize!(optimizer)

    abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    discrete_time_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system"))

    println("Abstraction built.")

    # ------------------------------------------------------------
    # 3) Define co-safe LTL problem with sets labeling
    # ------------------------------------------------------------

    concrete_problem =
        DI.Problem.OptimalControlProblem(
            concrete_system,
            _I_,
            _T_,
            nothing,
            (x,u) -> 1,
            0.0
        )

    # ------------------------------------------------------------
    # 4) Solve using the SAME pipeline optimizer
    # ------------------------------------------------------------

    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    MOI.optimize!(optimizer)

    success = MOI.get(optimizer, MOI.RawOptimizerAttribute("success"))
    println("OptimalControlProblem success: $success")

    if save && success
        JLD2.jldsave(filename; optimizer = optimizer)
    end
end

if success && animate
    # ------------------------------------------------------------
    # 5) Collect results
    # ------------------------------------------------------------
    abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
    abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
    concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

    abstract_I_ = abstract_problem.initial_set
    x0_sim = SY.get_concrete_state(abstract_system, abstract_I_[1])

    dt = 0.01
    nstep = 1000
    reached(x) = x ∈ _T_

    x_traj, u_traj = ST.get_closed_loop_trajectory(
        discrete_time_system,#RV.RobotModel.dt_system(dt; _X_ = _X_, _U_ = _U_),
        concrete_controller,
        x0_sim,
        nstep;
        verbose=true
        #stopping = reached
    )

    println("Trajectory (states): ", (x_traj.seq))
    println("Trajectory (inputs): ", (u_traj.seq))
    println("Number of transitions: ", length(u_traj.seq))

    # ------------------------------------------------------------
    # 6) Plot
    # ------------------------------------------------------------
    using Plots
    # fig = plot(; aspect_ratio = :equal)
    # plot!(
    #     concrete_problem;
    #     aspect_ratio = :equal,
    # )
    # plot!(fig, x_traj; color = :blue, dims = [1, 2])
    # display(fig)

    p = RV.PostProc.animate_robot_live(
        x_traj.seq,
        0.5,
        robot_geometry;
        grounded_left_foot=true
    )

    #p = RV.PostProc.plot_configuration_interval(_T_.lb, _T_.ub, robot_geometry, true; p = p)
end