# The same reachability task, solved two ways on the *same* abstraction:
#
#   A. as an OptimalControlProblem   -> the reachability fixed point runs on the
#      abstraction itself;
#   B. as a CoSafeLTLProblem whose monitor is the two-state automaton of `F(goal)`
#      -> the solver builds the synchronous product of the abstraction with that
#      monitor and runs the identical fixed point on the product.
#
# Both return the same controllable set. What differs is the graph searched.
# Reusing one automaton for both makes the comparison exact: only the solver
# changes.
#
# Output: printed table (numbers go on the slide by hand).
#
# Run: julia --project=test presentation/asset_spec_vs_ltl.jl

ENV["GKSwstype"] = "100"

using StaticArrays, JuMP
import LazySets
import MathOptInterface as MOI
import HybridSystems
using Symbolics, MathOptSymbolicAD
using Dionysos
const DI = Dionysos
const PR = DI.Problem
const MP = DI.Mapping
const SY = DI.Symbolic
const OPDS = DI.Optim.DiscreteSystems

const GOAL = LazySets.Hyperrectangle(; low = [3.0, 3.0], high = [4.0, 4.0])

# ------------------------------------------------------------
# A. the specialized solver, through the front-end
# ------------------------------------------------------------
model = Model(Dionysos.Optimizer)
@variable(model, -5 <= x[1:2] <= 5)
@variable(model, -1 <= u[1:2] <= 1)
@constraint(model, ∂(x[1]) == u[1])
@constraint(model, ∂(x[2]) == u[2])
@constraint(model, x in Final(GOAL))

set_attribute(model, "time_step", 0.3)
set_attribute(model, "state_grid", MP.GridFree(SVector(0.0, 0.0), SVector(0.2, 0.2)))
set_attribute(model, "input_grid", MP.GridFree(SVector(0.0, 0.0), SVector(0.5, 0.5)))
set_attribute(model, "jacobian_bound", u -> SMatrix{2, 2}(0.0, 0.0, 0.0, 0.0))
set_attribute(model, "print_level", 0)

println("— A: OptimalControlProblem, Final(goal) —")
optimize!(model)
println("  ", termination_status(model))

abs_sys = get_attribute(model, "abstract_system")
ctrl_solver = get_attribute(model, "control_solver")
abs_prob = ctrl_solver.abstract_problem
autom = abs_prob.system

n_abs = SY.get_n_state(autom)
n_trans = HybridSystems.ntransitions(autom)
time_a = ctrl_solver.abstract_problem_time_sec
win_a = length(ctrl_solver.abstract_optimizer.controllable_set)

println("  $n_abs states, $n_trans transitions")
println("  synthesis $(round(time_a; digits = 2)) s, controllable $win_a")

# ------------------------------------------------------------
# B. the same task as co-safe LTL, on the same automaton
# ------------------------------------------------------------
# The monitor of `F(goal)`: state 1 until the goal label is seen, state 2 for
# ever after, and state 2 accepts. Written directly rather than through Spot so
# that the only difference with A is the product construction.
target_states = collect(abs_prob.target_set)
monitor = OPDS.FunctionMonitor(1, Set([2]), (qa, ap) -> (qa == 2 || :g in ap) ? 2 : 1)

ltl_problem = PR.CoSafeLTLProblem(
    autom,
    collect(abs_prob.initial_set),
    monitor,
    Dict(:g => target_states),
    Dict{Symbol, Any}(),
)

println("— B: CoSafeLTLProblem, monitor of F(goal) —")
opt_b = MOI.instantiate(OPDS.OptimizerCoSafeLTLProblem)
MOI.set(opt_b, MOI.RawOptimizerAttribute("problem"), ltl_problem)
MOI.set(opt_b, MOI.RawOptimizerAttribute("print_level"), 0)
time_b = @elapsed MOI.optimize!(opt_b)

n_prod = SY.get_n_state(opt_b.product_automaton_optimizer.problem.system)
n_prod_trans = HybridSystems.ntransitions(opt_b.product_automaton_optimizer.problem.system)
win_b = length(opt_b.controllable_set)

println("  product $n_prod states, $n_prod_trans transitions")
println("  synthesis $(round(time_b; digits = 2)) s, controllable $win_b")

println()
println("=== same task, same abstraction, two formulations ===")
println("abstraction            : $n_abs states, $n_trans transitions")
println("graph searched      A  : $n_abs")
println("graph searched      B  : $n_prod  (×$(round(n_prod / n_abs; digits = 2)))")
println("synthesis time      A  : $(round(time_a; digits = 3)) s")
println(
    "synthesis time      B  : $(round(time_b; digits = 3)) s  (×$(round(time_b / time_a; digits = 2)))",
)
println("controllable states A  : $win_a")
println("controllable states B  : $win_b   (must match A)")
