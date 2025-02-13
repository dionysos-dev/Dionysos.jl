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
    MOI.get(optimizer, MOI.RawOptimizerAttribute("discretized_system")),
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
_T_ = UT.HyperRectangle(SVector(1.20, 5.75), SVector(1.25, 5.80))

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
    MOI.get(optimizer, MOI.RawOptimizerAttribute("discretized_system")),
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

# # # Example: DC-DC converter solved by [Uniform grid abstraction] (https://github.com/dionysos-dev/Dionysos.jl/blob/master/docs/src/manual/manual.md#solvers) by exploiting the incremental stability of the system.
# # ### Definition of the system
# # we can import the module containing the DCDC problem like this 

# ### Construction of the abstraction

# origin = SVector(0.0, 0.0)
# η = (2 / 4.0) * 10^(-3)

# # Note: In the following, `P` and `ϵ` are computed by hand, but their computation is not crucial since they only affect the visualization of the abstraction. See https://github.com/dionysos-dev/Dionysos.jl/issues/345
# ϵ = 0.1 * 0.01
# P = SMatrix{2, 2}(1.0224, 0.0084, 0.0084, 1.0031)
# state_grid = DO.GridEllipsoidalRectangular(origin, SVector(η, η), P / ϵ, concrete_system.X)

# u0 = SVector(1)
# hu = SVector(1)
# input_grid = DO.GridFree(u0, hu)

# optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
# MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem_safety)
# MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
# MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
# MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), DCDC.jacobian_bound())
# MOI.set(
#     optimizer,
#     MOI.RawOptimizerAttribute("approx_mode"),
#     Dionysos.Optim.Abstraction.UniformGridAbstraction.DELTA_GAS,
# )
# MOI.set(optimizer, MOI.RawOptimizerAttribute("δGAS"), true)
# MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.5)
# MOI.optimize!(optimizer)

# abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
# concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

# # ### Trajectory display
# # We choose the number of steps `nsteps` for the sampled system, i.e. the total elapsed time: `nstep`*`tstep`
# # as well as the true initial state `x0` which is contained in the initial state-space defined previously.
# nstep = 300
# x0 = SVector(1.2, 5.6)
# control_trajectory = ST.get_closed_loop_trajectory(
#     MOI.get(optimizer, MOI.RawOptimizerAttribute("discretized_system")),
#     concrete_controller,
#     x0,
#     nstep,
# )

# fig = plot(; aspect_ratio = :equal);
# plot!(concrete_problem_safety; opacity = 1.0);
# plot!(invariant_set, color = :blue, linecolor = :blue)
# plot!(uninvariant_set; color = :red, linecolor = :red)
# plot!(concrete_problem_safety.initial_set; color = :green, opacity = 1.0, label = "");
# plot!(control_trajectory)

# ### References
# 1. A. Girard, G. Pola and P. Tabuada, "Approximately Bisimilar Symbolic Models for Incrementally Stable Switched Systems," in IEEE Transactions on Automatic Control, vol. 55, no. 1, pp. 116-126, Jan. 2010.
# 2. S. Mouelhi, A. Girard, and G. Gössler. “CoSyMA: a tool for controller synthesis using multi-scale abstractions”. In: HSCC. ACM. 2013, pp. 83–88.
# 3. A. Girard. “Controller synthesis for safety and reachability via approximate bisimulation”. In: Automatica 48.5 (2012), pp. 947–953.
# 4. G. Reissig, A. Weber and M. Rungger, "Feedback Refinement Relations for the Synthesis of Symbolic Controllers," in IEEE Transactions on Automatic Control, vol. 62, no. 4, pp. 1781-1796.
