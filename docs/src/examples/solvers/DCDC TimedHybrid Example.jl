using StaticArrays, Plots, HybridSystems

# unfixed

## Import Dionysos sub-modules
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "dcdc_hybridautomata.jl"));

# optimal control problem
hybrid_system, optimizer_factory_list, optimizer_kwargs_dict, safety_specs =
    DCDC.generate_safety_system_and_problem()
concrete_controller = AB.TimedHybridAbstraction.solve_timed_hybrid_problem(
    hybrid_system,
    optimizer_factory_list,
    optimizer_kwargs_dict,
    safety_specs,
)
# get closed loop trajectory using the concrete_controller
# Use the same discretization parameters as defined in the problem
discretization_parameters = [
    (0.0005, 0.5, 0.5),  # Mode 1: doit correspondre au module DCDC
    (0.0005, 0.5, 0.5),  # Mode 2: doit correspondre au module DCDC
]  # Same as in DCDC module
traj, ctrls = AB.TimedHybridAbstraction.get_closed_loop_trajectory(
    discretization_parameters,
    hybrid_system,
    safety_specs,
    concrete_controller,
    safety_specs.initial_state,
    1000;  # Reduced number of steps for testing
    stopping = (specs, state) -> !AB.TimedHybridAbstraction.reached(specs, state), # Stop if unsafe
)

println("Closed-loop trajectory:")
for (idx, (t, u)) in enumerate(zip(traj, ctrls))
    println("[", idx, "] state: ", t, " - control applied: ", u)
end
println("Final state: ", traj[end])

# # ========== PLOT SIMPLE DE LA TRAJECTOIRE ==========
# # Région de sécurité du système
# iL_min, iL_max = 1.15, 1.55    # Courant inductance sûr [A]
# vC_min, vC_max = 5.45, 5.85    # Tension condensateur sûre [V]
# X_region = UT.HyperRectangle(SVector(iL_min, vC_min), SVector(iL_max, vC_max))

# # Extraire les positions de la trajectoire
# positions = [SVector(state[1][1], state[1][2]) for state in traj]

# # Extraire les coordonnées x et y séparément
# x_coords = [pos[1] for pos in positions]  # iL (courant)
# y_coords = [pos[2] for pos in positions]  # vC (tension)
# modes = [state[3] for state in traj]      # Modes

# println("Premier point: (", x_coords[1], ", ", y_coords[1], ") Mode: ", modes[1])
# println("Dernier point: (", x_coords[end], ", ", y_coords[end], ") Mode: ", modes[end])

# # Séparer les points par mode
# mode1_indices = findall(m -> m == 1, modes)
# mode2_indices = findall(m -> m == 2, modes)

# # Plot classique avec couleurs par mode
# fig = plot(; aspect_ratio = :equal)
# plot!(fig, X_region; label = "Safe region", color = :grey, alpha = 0.3)

# # Mode 1 en rouge
# if !isempty(mode1_indices)
#     plot!(
#         fig,
#         x_coords[mode1_indices],
#         y_coords[mode1_indices];
#         line = false,
#         marker = :circle,
#         markersize = 3,
#         color = :red,
#         label = "Mode 1",
#     )
# end

# # Mode 2 en bleu
# if !isempty(mode2_indices)
#     plot!(
#         fig,
#         x_coords[mode2_indices],
#         y_coords[mode2_indices];
#         line = false,
#         marker = :circle,
#         markersize = 3,
#         color = :blue,
#         label = "Mode 2",
#     )
# end

# # Ajouter une ligne pour voir la trajectoire
# plot!(
#     fig,
#     x_coords,
#     y_coords;
#     line = true,
#     linewidth = 1,
#     color = :black,
#     alpha = 0.5,
#     label = "Trajectory",
#     xlabel = "iL [A]",
#     ylabel = "vC [V]",
# )

# display(fig)
