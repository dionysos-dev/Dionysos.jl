using StaticArrays
using MathematicalSystems
using Dionysos
using Spot
using JuMP
import MathOptInterface as MOI

const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

# ------------------------------------------------------------
# 1) Define a simple 2D continuous-time system: x' = u
# ------------------------------------------------------------

# Domain (like your examples)
_X_ = UT.HyperRectangle(SVector(-2.0, -2.0), SVector(2.0, 2.0))
_U_ = UT.HyperRectangle(SVector(-1.0, -1.0), SVector(1.0, 1.0))

dynamic = (x, u) -> SVector(u[1], u[2])

concrete_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
    dynamic,
    2,  # nx
    2,  # nu
    _X_,
    _U_,
)

# ------------------------------------------------------------
# 2) Abstraction construction (EmptyProblem), same as you do
# ------------------------------------------------------------

empty_problem = DI.Problem.EmptyProblem(concrete_system, concrete_system.X)

# grid resolution
x0 = SVector(-2.0, -2.0)
hx = SVector(0.1, 0.1)                 # coarse so it runs fast
state_grid = DO.GridFree(x0, hx)

u0 = SVector(-1.0, -1.0)
hu = SVector(0.5, 0.5)
input_grid = DO.GridFree(u0, hu)

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), empty_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.5)

# choose an approx mode that exists in your setup
MOI.set(optimizer, MOI.RawOptimizerAttribute("approx_mode"), AB.UniformGridAbstraction.CENTER_SIMULATION)
MOI.set(optimizer, MOI.RawOptimizerAttribute("n_samples"), 1)
MOI.set(optimizer, MOI.RawOptimizerAttribute("threaded"), false)
MOI.set(optimizer, MOI.RawOptimizerAttribute("efficient"), true)

MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)
MOI.set(optimizer, MOI.Silent(), true)

MOI.optimize!(optimizer)

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
discrete_time_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system"))

println("Abstraction built.")

# ------------------------------------------------------------
# 3) Define co-safe LTL problem with LazySets-style labeling
# ------------------------------------------------------------

_I_  = UT.HyperRectangle(SVector(-1.7, -1.7), SVector(-1.6, -1.6))

g1   = UT.HyperRectangle(SVector(-0.6, -0.6), SVector(-0.2, -0.2))
g2   = UT.HyperRectangle(SVector( 1.0,  1.0), SVector( 1.7,  1.7))
obs  = UT.HyperRectangle(SVector(-0.1, -0.1), SVector( 0.1,  0.1))

# co-safe formula
# φ = ltl"G(!obs) & F(g1 & F(g2))"

struct MonitorG1ThenG2NoObs end

@inline function mon_next(::MonitorG1ThenG2NoObs, q::Int, ap::Tuple{Vararg{Symbol}})
    obs = (:obs in ap)
    g1  = (:g1 in ap)
    g2  = (:g2 in ap)

    obs && return 0
    q == 0 && return 0
    q == 3 && return 3

    if q == 1
        return g1 ? (g2 ? 3 : 2) : 1
    else
        # q == 2
        return g2 ? 3 : 2
    end
end

mon = AB.UniformGridAbstraction.FunctionMonitor(
    1,         # initial
    Set([3]),  # accepting
    (qa, ap) -> mon_next(MonitorG1ThenG2NoObs(), qa, ap),
)

# labeling dictionary: AP => concrete set (LazySet / HyperRectangle)
labeling = Dict{Symbol, Any}(
    :g1  => g1,
    :g2  => g2,
    :obs => obs,
)

# semantics per AP
ap_semantics = Dict{Symbol, Any}(
    :g1  => DO.INNER,
    :g2  => DO.INNER,
    :obs => DO.OUTER,
)

# This assumes your PR.CoSafeLTLProblem signature:
# CoSafeLTLProblem(system, initial_set, spec; labeling=..., ap_semantics=..., strict_spot=...)
concrete_problem_ltl = DI.Problem.CoSafeLTLProblem(
    concrete_system,
    _I_,
    mon, # φ, mon
    labeling,
    ap_semantics,
    false,
)

# ------------------------------------------------------------
# 4) Solve using the SAME pipeline optimizer
# ------------------------------------------------------------

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem_ltl)
MOI.optimize!(optimizer)

success = MOI.get(optimizer, MOI.RawOptimizerAttribute("success"))
println("Co-safe LTL success: $success")

fm_abs_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
# fm_abs_controller should be your ST.FiniteMemorySymbolicController

# ------------------------------------------------------------
# 5) Build concrete controller + memory updater and simulate
# ------------------------------------------------------------

# Requires the helper we discussed:
# ctrl_concrete, qa_ref, reset_memory!, update_memory! = solve_concrete_problem_fm(abstract_system, fm_abs_controller)
ctrl_concrete, qa_ref, reset_memory!, update_memory! =
    AB.UniformGridAbstraction.solve_concrete_problem_fm(abstract_system, fm_abs_controller; randomize=false)

x0 = SVector(-1.65, -1.65)
nstep = 60

traj = AB.UniformGridAbstraction.get_closed_loop_trajectory_fm(
    discrete_time_system,
    abstract_system,
    ctrl_concrete,
    update_memory!,
    x0,
    nstep;
    stopping = x -> false,
)

println("Trajectory length: ", length(traj.states.seq))

# ------------------------------------------------------------
# 6) Plot (optional)
# ------------------------------------------------------------
using Plots

fig = plot(; aspect_ratio=:equal)
plot!(fig, _X_; opacity=0.05)
plot!(fig, g1; color=:green, opacity=0.3)
plot!(fig, g2; color=:green, opacity=0.3)
plot!(fig, obs; color=:red, opacity=0.5)
plot!(fig, traj; dims=[1,2])
display(fig)
