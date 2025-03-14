using Test, Plots     #src
import StaticArrays: SVector

# At this point, we import the useful Dionysos sub-modules.
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "dc_dc.jl"))
concrete_system = DCDC.system()

### Construction of the abstraction
empty_problem = DI.Problem.EmptyProblem(concrete_system, concrete_system.X)

x0 = SVector(0.0, 0.0)
hx = SVector(2.0 / 4.0e3, 2.0 / 4.0e3)
state_grid = DO.GridFree(x0, hx)
u0 = SVector(1)
hu = SVector(1)
input_grid = DO.GridFree(u0, hu)
using JuMP

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), empty_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), DCDC.jacobian_bound())
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.5)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.RANDOM_SIMULATION,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("efficient"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("n_samples"), 1)
MOI.set(optimizer, MOI.Silent(), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)

# USER_DEFINED GROWTH LINEARIZED CENTER_SIMULATION RANDOM_SIMULATION

MOI.optimize!(optimizer)
abstraction_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstraction_construction_time_sec"))
println("Time to construct the abstraction: $(abstraction_time)")

### Solve a safety problem

# concrete_system = concrete_problem.system
_I_ = UT.HyperRectangle(SVector(1.19, 5.59), SVector(1.21, 5.61))
_S_ = UT.HyperRectangle(SVector(1.16, 5.46), SVector(1.53, 5.82))
concrete_problem_safety =
    Dionysos.Problem.SafetyProblem(concrete_system, _I_, _S_, Dionysos.Problem.Infinity())
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem_safety)
MOI.optimize!(optimizer)
abstract_problem_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem_time_sec"))
println("Time to solve the abstract problem: $(abstract_problem_time)")

abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
# @test length(abstract_controller.data) == 893803 #src
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
invariant_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("invariant_set"))
uninvariant_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("uninvariant_set"))

nstep = 300
x0 = SVector(1.2, 5.6)
control_trajectory = ST.get_closed_loop_trajectory(
    MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system")),
    concrete_controller,
    x0,
    nstep,
);

fig = plot(; aspect_ratio = :equal);
plot!(concrete_problem_safety; opacity = 1.0);
plot!(invariant_set; color = :blue, linecolor = :blue)
plot!(uninvariant_set; color = :red, linecolor = :red)
plot!(control_trajectory)
display(fig)

### Solve a reachability problem
_T_ = UT.HyperRectangle(SVector(1.20, 5.75), SVector(1.25, 5.80)) # T = Target, I = Init

# _T_ = UT.HyperRectangle(SVector(1.20, 5.75), SVector(1.25, 5.80))
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

abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
# @test length(abstract_controller.data) == 893803 #src
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
controllable_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("controllable_set"))
uncontrollable_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("uncontrollable_set"))

nstep = 300
x0 = SVector(1.2, 5.6)
function reached(x)
    if x ∈ concrete_problem_reachability.target_set
        return true
    else
        return false
    end
end

control_trajectory = ST.get_closed_loop_trajectory(
    MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system")),
    concrete_controller,
    x0,
    nstep;
    stopping = reached,
);

fig = plot(; aspect_ratio = :equal);
plot!(concrete_problem_reachability);
plot!(controllable_set; color = :yellow, linecolor = :yellow, label = "Controllable set")
plot!(uncontrollable_set; color = :black, linecolor = :black, label = "Uncontrollable set")
plot!(control_trajectory)
display(fig)
