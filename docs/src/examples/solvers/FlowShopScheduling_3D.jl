# #using StaticArrays, Plots, HybridSystems

# # At this point, we import the useful Dionysos sub-modules.
# using Dionysos
# const DI = Dionysos
# const UT = DI.Utils
# const DO = DI.Domain
# const ST = DI.System
# const SY = DI.Symbolic
# const OP = DI.Optim
# const AB = OP.Abstraction

# # we can import the module containing the FlowShopScheduling3D problem like this
# include(
#     joinpath(dirname(dirname(pathof(Dionysos))), "problems", "flowshopscheduling_3D.jl"),
# );
# # generate the concrete hybrid system and the problem specifications
# HybridSystem_automaton, optimizer_factory_list, optimizer_kwargs_dict, problem_specs =
#     FlowShopScheduling3D.generate_system_and_problem()

# # Keep discretization parameters for compatibility with get_closed_loop_trajectory
# discretization_parameters = [(0.5, 0.5, 0.2), (0.5, 0.5, 0.2), (0.5, 0.5, 0.2)]

# concrete_controller = AB.TimedHybridAbstraction.solve_timed_hybrid_problem(
#     HybridSystem_automaton,
#     optimizer_factory_list,
#     optimizer_kwargs_dict,
#     problem_specs,
# )

# traj, ctrls = AB.TimedHybridAbstraction.get_closed_loop_trajectory(
#     discretization_parameters,
#     HybridSystem_automaton,
#     problem_specs,
#     concrete_controller,
#     problem_specs.initial_state,
#     1000000;
#     stopping = AB.TimedHybridAbstraction.reached,
# )

# # Display trajectory and controls 
# # for (idx, (t, u)) in enumerate(zip(traj, ctrls))
# #     println("[", idx, "] state: ", t, " - control applied: ", u)
# # end
# # println("Final state: ", traj[end])
