using StaticArrays, Plots, HybridSystems
## Import Dionysos sub-modules
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

## Import the 1D flowshop problem definition
include(
    joinpath(dirname(dirname(pathof(Dionysos))), "problems", "flowshopscheduling_1D.jl"),
);

## Generate the hybrid system and problem specs
HybridSystem_automaton, optimizer_factory_list, optimizer_kwargs_dict, problem_specs =
    FlowShopScheduling1D.generate_system_and_problem()

# Keep discretization parameters for compatibility with get_closed_loop_trajectory
discretization_parameters = [
    (0.25, 0.25, 0.2),
    (0.25, 0.5, 0.1),
    (0.1, 0.1, 0.2),
    (0.1, 0.1, 0.1),
    (0.2, 0.25, 0.1),
]

## Synthesize the controller
concrete_controller = AB.TimedHybridAbstraction.solve_timed_hybrid_problem(
    HybridSystem_automaton,
    optimizer_factory_list,
    optimizer_kwargs_dict,
    problem_specs,
)

## Simulate closed-loop trajectory
traj, ctrls = AB.TimedHybridAbstraction.get_closed_loop_trajectory(
    discretization_parameters,
    HybridSystem_automaton,
    problem_specs,
    concrete_controller,
    problem_specs.initial_state,
    10000;
    stopping = AB.TimedHybridAbstraction.reached,
)

# for (idx, (t, u)) in enumerate(zip(traj, ctrls))
#     println("[", idx, "] state: ", t, " - control applied: ", u)
# end
# println("Final state: ", traj[end])

## Plot the trajectory

# using Plots
# using Printf

# ## Extract trajectory data
# t_vals = [s[2] for s in traj]
# x_vals = [s[1][1] for s in traj]
# k_vals = [s[3] for s in traj]

# ## Plot limits
# t_min, t_max = minimum(t_vals), maximum(t_vals)
# x_min, x_max = minimum(x_vals), maximum(x_vals)
# t_margin = (t_max - t_min) * 0.05
# x_margin = (x_max - x_min) * 0.1

# plt = plot(;
#     title = "System state over time (1D flowshop scheduling)",
#     xlabel = "Time (s)",
#     ylabel = "System state (x)",
#     xlims = (t_min - t_margin, t_max + t_margin),
#     ylims = (x_min - x_margin, x_max + x_margin),
#     legend = :outerright,
#     legendtitle = "Legend",
#     size = (1400, 900),
#     dpi = 150,
# )

# ## Draw acceptance regions (guards)
# region_colors = palette(:tab20)
# auto = HybridSystem_automaton.automaton
# for (i, t) in enumerate(HybridSystems.transitions(auto))
#     guard = HybridSystems.guard(HybridSystem_automaton, t)
#     if guard !== nothing
#         x_accept = [guard.lb[1], guard.ub[1]]
#         t_accept = [guard.lb[2], guard.ub[2]]
#         color_task = region_colors[(i - 1) % length(region_colors) + 1]
#         plot!(
#             [t_accept[1], t_accept[2], t_accept[2], t_accept[1], t_accept[1]],
#             [x_accept[1], x_accept[1], x_accept[2], x_accept[2], x_accept[1]];
#             fill = (0, 0.15),
#             fillcolor = color_task,
#             alpha = 0.3,
#             linecolor = color_task,
#             linewidth = 2,
#             linestyle = :dash,
#             label = "Guard $(i): x∈[$(x_accept[1]),$(x_accept[2])], t∈[$(t_accept[1]),$(t_accept[2])]",
#         )
#         # Add grid for each task
#         dx, du, dt = discretization_parameters[i]
#         # Vertical grid lines
#         t_grid = t_accept[1]:dt:t_accept[2]
#         for t in t_grid
#             if t >= t_min && t <= t_max
#                 vline!(
#                     [t];
#                     color = color_task,
#                     alpha = 0.25,
#                     linewidth = 1,
#                     linestyle = :dot,
#                     label = false,
#                 )
#             end
#         end
#         # Horizontal grid lines
#         x_grid = (floor(x_accept[1] / dx) * dx):dx:(ceil(x_accept[2] / dx) * dx)
#         for x in x_grid
#             if x >= x_min && x <= x_max
#                 hline!(
#                     [x];
#                     color = color_task,
#                     alpha = 0.25,
#                     linewidth = 1,
#                     linestyle = :dot,
#                     label = false,
#                 )
#             end
#         end
#     end
# end

# ## Plot the final target region (from problem specs)
# final_x_target = problem_specs.Xs_target[1]
# final_t_target = problem_specs.Ts_target[1]
# color_final = :magenta
# if !isnothing(final_x_target) && !isnothing(final_t_target)
#     x_accept = [final_x_target.lb[1], final_x_target.ub[1]]
#     t_accept = [final_t_target.lb[1], final_t_target.ub[1]]
#     plot!(
#         [t_accept[1], t_accept[2], t_accept[2], t_accept[1], t_accept[1]],
#         [x_accept[1], x_accept[1], x_accept[2], x_accept[2], x_accept[1]];
#         fill = (0, 0.25),
#         fillcolor = color_final,
#         alpha = 0.4,
#         linecolor = color_final,
#         linewidth = 3,
#         linestyle = :solid,
#         label = "Final target (problem target)",
#     )
# end

# ## Plot main trajectory
# plot!(
#     t_vals,
#     x_vals;
#     color = :black,
#     linewidth = 2,
#     label = "Trajectory",
#     linestyle = :solid,
#     alpha = 0.8,
# )

# ## Initial point and switches
# scatter!(
#     [t_vals[1]],
#     [x_vals[1]];
#     color = :green,
#     marker = :star5,
#     markersize = 8,
#     label = "Initial state",
# )
# switch_indices = [i for i in 2:length(k_vals) if k_vals[i] != k_vals[i - 1]]
# if !isempty(switch_indices)
#     t_switches = [t_vals[i] for i in switch_indices]
#     x_switches = [x_vals[i] for i in switch_indices]
#     scatter!(
#         t_switches,
#         x_switches;
#         color = :red,
#         marker = :diamond,
#         markersize = 8,
#         label = "Switch",
#     )
#     for i in switch_indices
#         plot!(
#             [t_vals[i - 1], t_vals[i]],
#             [x_vals[i - 1], x_vals[i]];
#             color = :red,
#             linewidth = 3,
#             linestyle = :dash,
#             alpha = 0.9,
#             label = false,
#         )
#     end
# end

# ## Task starts
# if !isempty(switch_indices)
#     t_starts = [t_vals[i] for i in switch_indices]
#     x_starts = [x_vals[i] for i in switch_indices]
#     scatter!(
#         t_starts,
#         x_starts;
#         color = :blue,
#         marker = :circle,
#         markersize = 5,
#         label = "Task start",
#     )
# end

# ## Control annotations (sampled)
# n_annotations = min(15, length(ctrls))
# step = max(1, length(ctrls) ÷ n_annotations)
# for i in 1:step:length(ctrls)
#     if i <= length(t_vals) && ctrls[i] != "switch" && ctrls[i] != "SWITCH"
#         if isa(ctrls[i], Number)
#             control_str = Printf.@sprintf("u = %.2f", Float64(ctrls[i]))
#         elseif isa(ctrls[i], Vector) && length(ctrls[i]) > 0
#             control_values =
#                 join([Printf.@sprintf("%.2f", Float64(u)) for u in ctrls[i]], ", ")
#             control_str = "u = $control_values"
#         elseif isa(ctrls[i], AbstractArray) && length(ctrls[i]) == 1
#             control_str = Printf.@sprintf("u = %.2f", Float64(ctrls[i][1]))
#         else
#             control_str = "u = $(ctrls[i])"
#         end
#         t_pos = t_vals[i]
#         x_pos = x_vals[i] + x_margin * 0.3
#         annotate!(t_pos, x_pos, text(control_str, :black, 8, :left))
#     end
# end

# display(plt)
