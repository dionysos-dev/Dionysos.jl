using StaticArrays
using MathematicalSystems
using Dionysos
using JuMP
using LinearAlgebra
import MathOptInterface as MOI

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction
const OPDS = OP.DiscreteSystems

# ------------------------------------------------------------
# 1) Define a simple 4D continuous-time system: x' = u
# ------------------------------------------------------------
include("./postproc.jl")
include("./geometry.jl")
include("./robot_model.jl")

_X_ = UT.HyperRectangle(
    SVector(-pi/3, -pi/3, -pi/3, -pi/3),
    SVector(pi/3, pi/3, pi/3, pi/3)
)

_U_ = UT.HyperRectangle(
    0.6*SVector(-1.0, -1.0, -1.0, -1.0),
    0.6*SVector(1.0, 1.0, 1.0, 1.0)
)

robot_geometry = geometry(0.0, 0.202, 0.172)

concrete_system = RobotModel.system(; _X_ = _X_, _U_ = _U_)
jacobian_bound = RobotModel.jacobian_bound()

# ------------------------------------------------------------
# 2) Abstraction construction (AlternatingSimulationProblem)
# ------------------------------------------------------------

alternating_simulation_problem =
    DI.Problem.AlternatingSimulationProblem(concrete_system, concrete_system.X)

# grid resolution (hx defines TOTAL CELL WIDTH)
dx = 0.25
x0 = SVector(0.0, 0.0, 0.0, 0.0)
hx = SVector(dx, dx, dx, dx)
state_grid = MP.GridFree(x0, hx)

du = 0.52
u0 = SVector(0.0, 0.0, 0.0, 0.0)
hu = SVector(du, du, du, du)
input_grid = MP.GridFree(u0, hu)

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("concrete_problem"),
    alternating_simulation_problem,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.5)

# choose an approx mode that exists in your setup
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.CENTER_SIMULATION,
) # GROWTH CENTER_SIMULATION
MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), jacobian_bound)

MOI.set(optimizer, MOI.RawOptimizerAttribute("efficient"), true)

MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)
MOI.set(optimizer, MOI.Silent(), true)

MOI.optimize!(optimizer)

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
discrete_time_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system"))

@show(discrete_time_system)

println("Abstraction built.")

# ------------------------------------------------------------
# 3) Define co-safe LTL problem with sets labeling
# ------------------------------------------------------------

_I_ = UT.HyperRectangle(
    -hx/2,
    hx/2
)
_T_ = UT.HyperRectangle(
    SVector(-3*dx, 0.0, 3*dx, 0.0) - hx/2,
    SVector(-3*dx, 0.0, 3*dx, 0.0) + hx/2
)

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

# ------------------------------------------------------------
# 5) Collect results
# ------------------------------------------------------------

abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

x0 = SVector(0.0, 0.0, 0.0, 0.0)
dt = 0.01
nstep = 1000
reached(x) = x ∈ _T_

x_traj, u_traj = ST.get_closed_loop_trajectory(
    RobotModel.dt_system(dt; _X_ = _X_, _U_ = _U_),
    concrete_controller,
    x0,
    nstep;
    stopping = reached
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

animate_robot_live(
    x_traj.seq,
    dt,
    robot_geometry;
    grounded_left_foot=true
)