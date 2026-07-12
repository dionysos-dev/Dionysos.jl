using Dionysos

include(
    joinpath(dirname(dirname(pathof(Dionysos))), "problems", "DCDC", "dcdc_converter.jl"),
)
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

using JuMP, Plots
import StaticArrays: SVector

concrete_system = DCDC.system()

### Construction of the abstraction
alternating_simulation_problem =
    DI.Problem.AlternatingSimulationProblem(concrete_system, concrete_system.X)

x0 = SVector(0.0, 0.0)
hx = SVector(2.0 / 4.0e3, 2.0 / 4.0e3)
state_grid = MP.GridFree(x0, hx)
XMapping = MP.ImplicitGridMapping(state_grid, concrete_system.X; incl_mode = MP.INNER)

u0 = SVector(1)
hu = SVector(1)
input_grid = MP.GridFree(u0, hu)

Δt = 0.5

using JuMP

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

# MOI.set(
#     optimizer,
#     MOI.RawOptimizerAttribute("execution_backend"),
#     SY.SequentialBackend(),
# )

MOI.set(optimizer, MOI.RawOptimizerAttribute("execution_backend"), SY.ThreadedBackend(0.2))

# MOI.set(
#     optimizer,
#     MOI.RawOptimizerAttribute("execution_backend"),
#     SY.JuliaDistributedBackend(
#         nothing,      # use Distributed.workers()
#         nothing,      # nparts = number of workers
#         :roundrobin,
#         false,        # threaded_per_worker
#     ),
# )

# MOI.set(
#     optimizer,
#     MOI.RawOptimizerAttribute("execution_backend"),
#     SY.SlurmArrayBackend(
#         1,
#         1,
#         "out/transitions",
#         :contiguous,
#         false,                # write_only
#     ),
# )

MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("concrete_problem"),
    alternating_simulation_problem,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("XMapping"), XMapping)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("h"), hx)  # optional if you pass state_grid
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), DCDC.jacobian_bound())
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), Δt)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.CENTER_SIMULATION, # USER_DEFINED GROWTH LINEARIZED CENTER_SIMULATION RANDOM_SIMULATION
)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("transition_metadata"),
    SY.TransitionMetadata(),
) # SY.NoTransitionMetadata(), SY.TransitionMetadata()

MOI.set(optimizer, MOI.RawOptimizerAttribute("n_samples"), 1)

MOI.set(optimizer, MOI.Silent(), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("automaton_constructor"),
    (n, m) -> SY.FastIndexedAutomatonList(n, m),
)

MOI.optimize!(optimizer)

abstraction_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstraction_construction_time_sec"))
println("Time to construct the abstraction: $(abstraction_time)")

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
println(SY.has_metadata(abstract_system))
tr = first(SY.enum_transitions(abstract_system))
println(tr)
println(SY.get_metadata(abstract_system, tr))

### Solve a safety problem

# concrete_system = concrete_problem.system
_I_ = UT.box(SVector(1.19, 5.59), SVector(1.21, 5.61))
_S_ = UT.box(SVector(1.16, 5.46), SVector(1.53, 5.82))
concrete_problem_safety =
    Dionysos.Problem.SafetyProblem(concrete_system, _I_, _S_, Dionysos.Problem.Infinity())
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem_safety)

MOI.optimize!(optimizer)
abstract_problem_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem_time_sec"))
println("Time to solve the abstract problem: $(abstract_problem_time)")

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
invariant_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("invariant_set"))
invariant_set_complement =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("invariant_set_complement"))

automaton = SY.get_automaton(abstract_system)
fig = histogram(SY.nondeterminism_counts(automaton); legend = false)
display(fig)
println("Number of self-loops: $(SY.count_self_loops(automaton))")

nstep = 300
x0 = SVector(1.2, 5.6)
traj = ST.get_closed_loop_trajectory(
    MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system")),
    concrete_controller,
    x0,
    nstep,
);

XMapping = SY.get_state_mapping(abstract_system)
Xset = SY.get_state_set(abstract_system)

fig = plot(; aspect_ratio = :equal);
plot!(concrete_problem_safety; opacity = 1.0);
plot!(XMapping; efficient = true, color = :grey)
plot!((Xset, XMapping); efficient = true, color = :grey)
plot!((invariant_set, XMapping); color = :blue, linecolor = :blue)
plot!((invariant_set_complement, XMapping); color = :red, linecolor = :red)
plot!(traj)
display(fig)

# Export in csv file the controller, and reload it
# using CSV, DataFrames
# filename = "concrete_controller"
# Dionysos.export_controller_csv(optimizer, filename)
# Dionysos.import_controller_csv(filename)

## Solve a reachability problem
_T_ = UT.box(SVector(1.20, 5.75), SVector(1.25, 5.80))

# _T_ = UT.box(SVector(1.20, 5.75), SVector(1.25, 5.80))
concrete_problem_reachability = Dionysos.Problem.OptimalControlProblem(
    concrete_system,
    _I_,
    _T_,
    nothing,
    nothing,
    Dionysos.Problem.Infinity(),
)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("concrete_problem"),
    concrete_problem_reachability,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("early_stop"), false)
MOI.optimize!(optimizer)
abstract_problem_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem_time_sec"))
println("Time to solve the abstract problem: $(abstract_problem_time)")

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
controllable_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("controllable_set"))
uncontrollable_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("uncontrollable_set"))

nstep = 300
x0 = SVector(1.2, 5.6)
reached(x) = x ∈ concrete_problem_reachability.target_set

traj = ST.get_closed_loop_trajectory(
    MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system")),
    concrete_controller,
    x0,
    nstep;
    stopping = reached,
);

XMapping = SY.get_state_mapping(abstract_system)
Xset = SY.get_state_set(abstract_system)

fig = plot(; aspect_ratio = :equal);
plot!(concrete_problem_reachability);
# plot!((Xset, XMapping); color = :grey, linecolor = :grey, label = "Domain")
plot!(
    (controllable_set, XMapping);
    color = :yellow,
    linecolor = :yellow,
    label = "Controllable set",
)
plot!(
    (uncontrollable_set, XMapping);
    color = :black,
    linecolor = :black,
    label = "Uncontrollable set",
)
plot!(traj)
display(fig)

# ------------------------------------------------------------
# Animation with dashboard
# ------------------------------------------------------------

system_plot! = DCDC.system_plot!()
Dionysos.animate_trajectory_dashboard(
    system_plot!,
    traj;
    xdims = (1, 2),
    udims = (1,),
    Δt = Δt,
    fps = 5,
    # filename = "dcdc_converter_dashboard.mp4",
    xlabel_state = "iL",
    ylabel_state = "vC",
    xlabel_input = "time [s]",
    ylabel_input = "mode",
    xlims_state = (1.15, 1.55),
    ylims_state = (5.45, 5.85),
    ylims_input = (0.5, 2.5),
)
