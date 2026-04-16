using Distributed

# ----------------------------------------------------------------------
# Worker configuration
# ----------------------------------------------------------------------
# const NWORKERS = 3
# length(workers()) < NWORKERS && addprocs(NWORKERS - length(workers()))

@everywhere using Dionysos

dcdc_path = joinpath(dirname(dirname(pathof(Dionysos))), "problems", "dc_dc.jl")
@everywhere include($dcdc_path)

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

using JuMP
using MathOptInterface
using Plots
@everywhere using StaticArrays

const MOI = MathOptInterface

# ----------------------------------------------------------------------
# Concrete system
# ----------------------------------------------------------------------
concrete_system = DCDC.system()

# ----------------------------------------------------------------------
# Construction of the abstraction
# ----------------------------------------------------------------------
alternating_simulation_problem =
    DI.Problem.AlternatingSimulationProblem(concrete_system, concrete_system.X)

x0_grid = SVector(0.0, 0.0)
hx = SVector(2.0 / 4.0e3, 2.0 / 4.0e3)
state_grid = MP.GridFree(x0_grid, hx)
XMapping = MP.ImplicitGridMapping(state_grid, concrete_system.X; incl_mode = MP.INNER)

u0 = SVector(1)
hu = SVector(1)
input_grid = MP.GridFree(u0, hu)

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("concrete_problem"),
    alternating_simulation_problem,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("XMapping"), XMapping)
# MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
# MOI.set(optimizer, MOI.RawOptimizerAttribute("h"), hx) # optional if state_grid is provided
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), DCDC.jacobian_bound())
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.5)

MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.GROWTH,
    # USER_DEFINED / GROWTH / LINEARIZED / CENTER_SIMULATION / RANDOM_SIMULATION
)

# ----------------------------------------------------------------------
# Parallelism settings
# ----------------------------------------------------------------------
USE_DISTRIBUTED = !isempty(workers())
DISTRIBUTED_NPARTS = max(length(workers()), 1)
DISTRIBUTED_PARTITION_STRATEGY = :roundrobin
USE_THREADED_PER_WORKER = false

MOI.set(optimizer, MOI.RawOptimizerAttribute("distributed"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("distributed_nparts"), DISTRIBUTED_NPARTS)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("distributed_partition_strategy"),
    DISTRIBUTED_PARTITION_STRATEGY,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("threaded"), USE_THREADED_PER_WORKER)

MOI.set(optimizer, MOI.RawOptimizerAttribute("efficient"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("n_samples"), 1)

MOI.set(optimizer, MOI.Silent(), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("automaton_constructor"),
    (n, m) -> ST.NewIndexedAutomatonList(n, m),
)

println("distributed = ", USE_DISTRIBUTED)
println("nworkers = ", length(workers()))
println("distributed_nparts = ", DISTRIBUTED_NPARTS)

if USE_DISTRIBUTED
    @info "Warming up distributed workers..."
    for p in workers()
        remotecall_wait(p) do
            sys = DCDC.system()
            x = SVector(0.0, 0.0)
            u = SVector(1)
            return sys.f(x, u)
        end
    end
end

MOI.optimize!(optimizer)

abstraction_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstraction_construction_time_sec"))
println("Time to construct the abstraction: $(abstraction_time) sec")

# ----------------------------------------------------------------------
# Solve a safety problem
# ----------------------------------------------------------------------
_I_ = UT.HyperRectangle(SVector(1.19, 5.59), SVector(1.21, 5.61))
_S_ = UT.HyperRectangle(SVector(1.16, 5.46), SVector(1.53, 5.82))

concrete_problem_safety =
    DI.Problem.SafetyProblem(concrete_system, _I_, _S_, DI.Problem.Infinity())

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem_safety)

MOI.optimize!(optimizer)

abstract_problem_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem_time_sec"))
println("Time to solve the abstract problem: $(abstract_problem_time) sec")

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
invariant_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("invariant_set"))
invariant_set_complement =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("invariant_set_complement"))

# ----------------------------------------------------------------------
# Closed-loop simulation
# ----------------------------------------------------------------------
nstep = 300
x0 = SVector(1.2, 5.6)

x_traj, u_traj = ST.get_closed_loop_trajectory(
    MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system")),
    concrete_controller,
    x0,
    nstep,
)

# ----------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------
XMapping_plot = SY.get_state_mapping(abstract_system)
Xset = SY.get_state_domain(abstract_system)

fig = plot(; aspect_ratio = :equal)
plot!(concrete_problem_safety; opacity = 1.0)
plot!(XMapping_plot; efficient = true, color = :grey)
plot!((Xset, XMapping_plot); efficient = true, color = :grey)
plot!((invariant_set, XMapping_plot); color = :blue, linecolor = :blue)
plot!((invariant_set_complement, XMapping_plot); color = :red, linecolor = :red)
plot!(x_traj)

display(fig)
