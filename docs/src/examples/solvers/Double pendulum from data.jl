using Test     #src
using StaticArrays, Plots
# At this point, we import the useful Dionysos sub-modules.
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "double_pendulum.jl"))
concrete_problem = DoublePendulum.problem(; approx_mode = "growth")
concrete_system = concrete_problem.system
x0 = SVector(0.0, 0.0, 0.0, 0.0)

hx_param = 0.3
hx = SVector(hx_param, hx_param, hx_param, hx_param)

state_grid = DO.GridFree(x0, hx)
u0 = SVector(0.0);
h = SVector(1.); # discretization step of the input space
input_grid = DO.GridFree(u0, h);

using JuMP
optimizer = MOI.instantiate(AB.SampleBasedAbstraction.Optimizer)
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.optimize!(optimizer)

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
automaton = abstract_system.autom
# UT.analyze_non_determinism(automaton)
# n_sl = UT.analyze_self_loops(automaton)
# println("Number of self loops: $n_sl")

abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
println()

concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
# ### Trajectory display
# We choose the number of steps `nsteps` for the sampled system, i.e. the total elapsed time: `nstep`*`tstep`
# as well as the true initial state `x0` which is contained in the initial state-space defined previously.
nstep = 100
x0 = SVector(0.15,0.0,0.0,0.0) # SVector(0.15,0.0) #
control_trajectory =
    ST.get_closed_loop_trajectory(concrete_system.f, concrete_controller, x0, nstep)
fig = plot(; aspect_ratio = :equal);
plot!(concrete_system.X);
plot!(control_trajectory; markersize=1,arrows=false)