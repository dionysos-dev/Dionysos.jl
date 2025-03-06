using Plots     #src
import StaticArrays: SVector, SMatrix
import MathematicalSystems

# At this point, we import the useful Dionysos sub-modules.
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const PB = DI.Problem
const OP = DI.Optim
const AB = OP.Abstraction

function f(x, u)
    return SVector{4}(x[3], x[4], u[1], u[2])
end

function jacobian_bound(u)
    return SMatrix{4, 4}(
        0.0,
        0.0,
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    )
end

_X_ = UT.HyperRectangle(SVector(-10.0, -10.0, -2.0, -2.0), SVector(10.0, 10.0, 2.0, 2.0))
obs1 = UT.HyperRectangle(
    SVector(-2.0, -2.0, _X_.lb[3], _X_.lb[4]),
    SVector(2.0, 2.0, _X_.ub[3], _X_.ub[4]),
)
obs2 = UT.HyperRectangle(
    SVector(7.0, -5.0, _X_.lb[3], _X_.lb[4]),
    SVector(8.0, 5.0, _X_.ub[3], _X_.ub[4]),
)
obstacles = UT.LazyUnionSetArray([obs1, obs2])
_X_ = UT.LazySetMinus(_X_, obstacles)

_U_ = UT.HyperRectangle(SVector(-2.0, -2.0), SVector(2.0, 2.0))

## Concrete system
concrete_system =
    MathematicalSystems.MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        f,
        UT.get_dims(_X_),
        UT.get_dims(_U_),
        _X_,
        _U_,
    )

## Concrete problem
_I_ = UT.HyperRectangle(SVector(0.0, -9.0, 0.0, 0.0), SVector(1.0, -8.0, 0.0, 0.0))
_T_ = UT.HyperRectangle(SVector(-2.0, 6.0, -1.0, -1.0), SVector(2.0, 10.0, 1.0, 1.0))
concrete_problem =
    PB.OptimalControlProblem(concrete_system, _I_, _T_, nothing, nothing, PB.Infinity())

fig = plot(; aspect_ratio = :equal);
plot!(concrete_problem; minuscolor = :black);
display(fig)

## Solver parameters
x0 = SVector(0.0, 0.0, 0.0, 0.0)
hx = SVector(0.4, 0.4, 0.2, 0.2)
state_grid = DO.GridFree(x0, hx)
u0 = SVector(0.0, 0.0)
hu = SVector(0.4, 0.4)
input_grid = DO.GridFree(u0, hu)
using JuMP

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.CENTER_SIMULATION,
) # USER_DEFINED GROWTH LINEARIZED CENTER_SIMULATION RANDOM_SIMULATION
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), jacobian_bound)
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.5)
MOI.set(optimizer, MOI.RawOptimizerAttribute("efficient"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("early_stop"), false)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)

## Solve the problem
MOI.optimize!(optimizer)

## Get the results
abstraction_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstraction_construction_time_sec"))
println("Time to construct the abstraction: $(abstraction_time)")
abstract_problem_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem_time_sec"))
println("Time to solve the abstract problem: $(abstract_problem_time)")

success = MOI.get(optimizer, MOI.RawOptimizerAttribute("success"))
println("Success?: $(success)")
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
controllable_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("controllable_set"))
uncontrollable_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("uncontrollable_set"))

nstep = 300
x0 = SVector(1.0, -8.0, 0.0, 0.0)
function reached(x)
    if x ∈ concrete_problem.target_set
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
plot!(concrete_problem);
# plot!(controllable_set; color = :yellow, linecolor = :yellow, label = "Controllable set", efficient = true)
plot!(control_trajectory; arrows = false, ms = 2.0, color = :blue)
display(fig)
