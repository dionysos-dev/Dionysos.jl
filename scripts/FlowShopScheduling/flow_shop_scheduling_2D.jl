# ============================================================================
# 2-D flowshop scheduling (2-D state + time), the executable analogue of the
# heavy `flow_shop_scheduling_3D.jl`. Three tasks (modes), each with 2-D linear
# dynamics and a clock; guarded switches reset the state and advance the clock.
# The hybrid abstraction is (x1, x2, t) per mode, small enough to build and solve.
# ============================================================================

using StaticArrays, Plots, Printf
import HybridSystems as HS
import MathOptInterface as MOI
import LazySets

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
# Plot: x1(t), x2(t), mode(t)
# -----------------------------------------------------

t_traj = [s[2] for s in aug_x_traj]
x1_traj = [s[1][1] for s in aug_x_traj]
x2_traj = [s[1][2] for s in aug_x_traj]
k_traj = [s[3] for s in aug_x_traj]
switch_indices = [i for i in 2:length(k_traj) if k_traj[i] != k_traj[i - 1]]

final_spec = first(values(concrete_problem.target_set.per_mode))  # TimedSpec of the last task

fig_state = plot(;
    title = "2-D flowshop: state over time",
    xlabel = "time t",
    ylabel = "state",
    legend = :outerright,
    size = (1100, 600),
)
plot!(fig_state, t_traj, x1_traj; linewidth = 2, color = :blue, label = "x1")
plot!(fig_state, t_traj, x2_traj; linewidth = 2, color = :orange, label = "x2")
vspan!(
    fig_state,
    [final_spec.tmin, final_spec.tmax];
    color = :green,
    alpha = 0.1,
    label = "target time window",
)
if !isempty(switch_indices)
    scatter!(
        fig_state,
        t_traj[switch_indices],
        x1_traj[switch_indices];
        color = :red,
        marker = :diamond,
        markersize = 7,
        label = "switch",
    )
end

fig_mode = plot(
    t_traj,
    k_traj;
    seriestype = :steppost,
    linewidth = 2,
    title = "2-D flowshop: task (mode)",
    xlabel = "time t",
    ylabel = "task",
    yticks = ([1, 2, 3], ["1", "2", "3"]),
    ylims = (0.5, 3.5),
    legend = false,
)

fig = plot(fig_state, fig_mode; layout = (2, 1), size = (1100, 800))
display(fig)
# savefig(fig, "flow_shop_scheduling_2D.png")
