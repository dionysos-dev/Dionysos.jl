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
include(joinpath(dirname(dirname(pathof(Dionysos))), "problems/simple_pendulum.jl")) 
# and we can instantiate the DC system with the provided system
concrete_problem = Pendulum.problem(; approx_mode = "growth")
concrete_system = concrete_problem.system
x0 = SVector(0.0, 0.0)

hx_param = 0.2

hx = SVector(hx_param, hx_param)
state_grid = DO.GridFree(x0, hx)
u0 = SVector(0.0);
h = SVector(0.3);
input_grid = DO.GridFree(u0, h);

using JuMP
optimizer = MOI.instantiate(AB.SampleBasedAbstraction.Optimizer)

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.optimize!(optimizer)
abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

automaton = abstract_system.autom
UT.analyze_non_determinism(automaton, abstract_system)
n_sl = UT.analyze_self_loops(automaton)
println("Number of self loops: $n_sl")

# ### Trajectory display
# We choose the number of steps `nsteps` for the sampled system, i.e. the total elapsed time: `nstep`*`tstep`
# as well as the true initial state `x0` which is contained in the initial state-space defined previously.
nstep = 10000
function reached(x)
    if x ∈ concrete_problem.target_set    #target_set
        return true
    else
        return false
    end
end

x0 = SVector(0.0,0.0)
#x0 = SVector(4.5*pi/180, 0.75) # SVector(pi+0.15,0.5)
control_trajectory =
    ST.get_closed_loop_trajectory(concrete_system.f, concrete_controller, x0, nstep;
    stopping = reached)
fig = plot(; aspect_ratio = :equal);
# plot!(concrete_system.X);
plot!(abstract_system.Xdom; color = :blue, opacity = 0.1);
plot!(
    SY.get_domain_from_symbols(abstract_system, abstract_problem.initial_set);
    color = :green,
    opacity= 1.0
);
plot!(
    SY.get_domain_from_symbols(abstract_system, abstract_problem.target_set); #target_set
    color = :red,
    opacity= 1.0
);
plot!(control_trajectory; markersize=1,arrows=false)