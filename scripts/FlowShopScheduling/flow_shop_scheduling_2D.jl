# ============================================================================
# 2-D flowshop scheduling (2-D state + time), the executable analogue of the
# heavy `flow_shop_scheduling_3D.jl`. Three tasks (modes), each with 2-D linear
# dynamics and a clock; guarded switches reset the state and advance the clock.
# The hybrid abstraction is (x1, x2, t) per mode, small enough to build and solve.
# ============================================================================

using StaticArrays, Plots
import MathOptInterface as MOI

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const OP = DI.Optim
const AB = OP.Abstraction

include("../../problems/FlowShopScheduling/flow_shop_scheduling_2D.jl")

# -----------------------------------------------------
# Problem
# -----------------------------------------------------

concrete_problem = FlowShopScheduling2D.problem()
concrete_system = concrete_problem.system

# -----------------------------------------------------
# Per-task subsolver parameters  (dx, du, dt)
# -----------------------------------------------------

discretization_parameters = [(0.5, 0.5, 0.2), (0.5, 0.5, 0.2), (0.5, 0.5, 0.2)]

optimizer_list = [
    () -> AB.UniformGridAbstraction.Optimizer{Float64}(),
    () -> AB.UniformGridAbstraction.Optimizer{Float64}(),
    () -> AB.UniformGridAbstraction.Optimizer{Float64}(),
]

jacobian_bounds = FlowShopScheduling2D.jacobian_bounds()

make_kwargs(i) = Dict(
    "state_grid" => MP.GridFree(
        SVector(0.0, 0.0),
        SVector(discretization_parameters[i][1], discretization_parameters[i][1]),
    ),
    "input_grid" => MP.GridFree(
        SVector(0.0, 0.0),
        SVector(discretization_parameters[i][2], discretization_parameters[i][2]),
    ),
    "time_step" => discretization_parameters[i][3],
    "approx_mode" => AB.UniformGridAbstraction.GROWTH,
    "jacobian_bound" => jacobian_bounds[i],
    "print_level" => 1,
)
optimizer_kwargs_dict = [make_kwargs(1), make_kwargs(2), make_kwargs(3)]

# -----------------------------------------------------
# Solve
# -----------------------------------------------------

optimizer = MOI.instantiate(AB.HybridSystemAbstraction.Optimizer)
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("optimizer_list"), optimizer_list)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("optimizer_kwargs_dict"),
    optimizer_kwargs_dict,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 1)

MOI.optimize!(optimizer)

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
println("Abstraction states (x1, x2, t, mode): ", DI.Symbolic.get_n_state(abstract_system))

# -----------------------------------------------------
# Closed-loop trajectory  (augmented state is `([x1,x2], t, mode)`)
# -----------------------------------------------------

reached(aug_x) = AB.HybridSystemAbstraction.reached(concrete_problem, aug_x)

initial_state = concrete_problem.initial_set
tsteps = [discretization_parameters[i][3] for i in 1:3]
aug_x_traj, u_traj = AB.HybridSystemAbstraction.get_closed_loop_trajectory(
    concrete_system,
    concrete_controller,
    tsteps,
    initial_state,
    10000;
    stopping = reached,
)

println("Reached target: ", reached(aug_x_traj[end]), " ; final state: ", aug_x_traj[end])

# -----------------------------------------------------
# Dashboard animation: system view (left) + state / input / task-over-time (right)
# -----------------------------------------------------

# Switch inputs are "SWITCH ..." labels; map them to NaN so the numeric input panel skips them.
u_numeric =
    [u isa AbstractString ? SVector(NaN, NaN) : SVector{2}(Float64.(u)) for u in u_traj]

system_plot! = FlowShopScheduling2D.system_plot!(; problem = concrete_problem)
Dionysos.animate_trajectory_dashboard(
    system_plot!,
    ST.Trajectory(aug_x_traj),
    ST.Trajectory(u_numeric);
    xdims = (1, 2),
    udims = (1, 2),
    Δt = 0.2,
    fps = 4,
    state_of = s -> s[1],               # continuous state ([x1, x2]) of the augmented state
    modes = [s[3] for s in aug_x_traj], # task per step -> task-vs-time panel on the right
    ylabel_mode = "task",
    title = "2-D flowshop scheduling",
    xlabel_state = "x1",
    ylabel_state = "x2",
    xlabel_input = "u1",
    ylabel_input = "u2",
    # filename = "flow_shop_scheduling_2D.gif",
)
