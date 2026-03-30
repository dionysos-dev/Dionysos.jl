# ---------------------------------------------------------------------------
#  Configuration from environment variables
# ---------------------------------------------------------------------------
const USE_DISTRIBUTED = lowercase(get(ENV, "DIONYSOS_DISTRIBUTED", "false")) == "true"
const USE_THREADED = lowercase(get(ENV, "DIONYSOS_THREADED", "false")) == "true"
const N_PARTS = parse(Int, get(ENV, "DIONYSOS_NPARTS", "8"))

using Distributed
if USE_DISTRIBUTED && length(workers()) < 2
    addprocs(max(N_PARTS, 2) - length(workers()))
end

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

using JuMP, Plots
import StaticArrays: SVector

concrete_system = DCDC.system()

### Construction of the abstraction
empty_problem = DI.Problem.EmptyProblem(concrete_system, concrete_system.X)

x0 = SVector(0.0, 0.0)
hx = SVector(2.0 / 4.0e3, 2.0 / 4.0e3) * 0.5
state_grid = MP.GridFree(x0, hx)
XMapping = MP.ImplicitGridMapping(state_grid, concrete_system.X; incl_mode = MP.INNER)

u0 = SVector(1)
hu = SVector(1)
input_grid = MP.GridFree(u0, hu)
optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), empty_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("XMapping"), XMapping)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("h"), hx)  # optional if you pass state_grid
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), DCDC.jacobian_bound())
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.5)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.GROWTH, # USER_DEFINED GROWTH LINEARIZED CENTER_SIMULATION RANDOM_SIMULATION
)

MOI.set(optimizer, MOI.RawOptimizerAttribute("distributed"), USE_DISTRIBUTED)
MOI.set(optimizer, MOI.RawOptimizerAttribute("distributed_nparts"), N_PARTS)
MOI.set(optimizer, MOI.RawOptimizerAttribute("distributed_partition_strategy"), :roundrobin)
MOI.set(optimizer, MOI.RawOptimizerAttribute("threaded"), USE_THREADED)

MOI.set(optimizer, MOI.RawOptimizerAttribute("efficient"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("n_samples"), 1)

MOI.set(optimizer, MOI.Silent(), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("automaton_constructor"),
    (n, m) -> SY.NewIndexedAutomatonList(n, m),
)

MOI.optimize!(optimizer)
abstraction_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstraction_construction_time_sec"))
println("Time to construct the abstraction: $(abstraction_time)")
