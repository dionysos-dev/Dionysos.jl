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

# For x' = u, ∂f/∂x = 0, so bound is the zero matrix (independent of u)
jacobian_bound = u -> @SMatrix [
    0.0 0.0;
    0.0 0.0
]

# ------------------------------------------------------------
# 2) Abstraction construction (EmptyProblem)
# ------------------------------------------------------------

empty_problem = DI.Problem.EmptyProblem(concrete_system, concrete_system.X)

# grid resolution
x0 = SVector(-2.0, -2.0)
hx = SVector(0.2, 0.2)
state_grid = DO.GridFree(x0, hx)

u0 = SVector(-1.0, -1.0)
hu = SVector(0.5, 0.5)
input_grid = DO.GridFree(u0, hu)

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), empty_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.3)

# choose an approx mode that exists in your setup
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.GROWTH,
) # GROWTH CENTER_SIMULATION
MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), jacobian_bound)

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

_I_ = UT.HyperRectangle(SVector(-1.7, -1.7), SVector(-1.6, -1.6))

g11 = UT.HyperRectangle(SVector(-1.0, 1.0), SVector(-0.3, 1.7))
g12 = UT.HyperRectangle(SVector(1.0, 1.0), SVector(1.7, 1.7))
g1 = UT.LazyUnionSetArray([g11, g12])

g2_big = UT.HyperRectangle(SVector(-1.5, -1.2), SVector(-0.6, -0.2))
g2_hole = UT.HyperRectangle(SVector(-1.2, -1.0), SVector(-0.9, -0.8))
g2 = UT.LazySetMinus(g2_big, g2_hole)

g3 = UT.HyperRectangle(SVector(1.0, -1.8), SVector(1.5, -1.1))

obs1 = UT.HyperRectangle(SVector(-0.5, -0.5), SVector(0.5, 0.5))
obs2 = UT.HyperRectangle(SVector(1.3, -0.5), SVector(2.0, 0.5))
obs = UT.LazyUnionSetArray([obs1, obs2])

# co-safe formula
φ = ltl"G(!obs) & F(g1 & F(g2 & F(g3  & F(g1))))"

struct MonitorG1G2G3G1NoObs end

@inline function mon_next(::MonitorG1G2G3G1NoObs, q::Int, ap::Tuple{Vararg{Symbol}})
    obs = (:obs in ap)
    g1 = (:g1 in ap)
    g2 = (:g2 in ap)
    g3 = (:g3 in ap)

    obs && return 0        # safety violation -> dead
    q == 0 && return 0     # dead sink
    q == 5 && return 5     # done sink

    if q == 1
        # waiting for first g1
        return g1 ? 2 : 1
    elseif q == 2
        # have first g1, waiting for g2
        return g2 ? 3 : 2
    elseif q == 3
        # have g2, waiting for g3
        return g3 ? 4 : 3
    else
        @assert q == 4
        # have g3, waiting for final g1
        return g1 ? 5 : 4
    end
end

mon = AB.UniformGridAbstraction.FunctionMonitor(
    1,         # initial
    Set([5]),  # accepting
    (qa, ap) -> mon_next(MonitorG1G2G3G1NoObs(), qa, ap),
)

# labeling dictionary: AP => concrete set (LazySet / HyperRectangle)
labeling = Dict{Symbol, Any}(:g1 => g1, :g2 => g2, :g3 => g3, :obs => obs)

# semantics per AP
ap_semantics =
    Dict{Symbol, Any}(:g1 => DO.INNER, :g2 => DO.INNER, :g3 => DO.INNER, :obs => DO.OUTER)

concrete_problem = DI.Problem.CoSafeLTLProblem(
    concrete_system,
    _I_,
    φ, # φ, mon
    labeling,
    ap_semantics,
    false,
)

# ------------------------------------------------------------
# 4) Solve using the SAME pipeline optimizer
# ------------------------------------------------------------

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.optimize!(optimizer)

success = MOI.get(optimizer, MOI.RawOptimizerAttribute("success"))
println("Co-safe LTL success: $success")

# ------------------------------------------------------------
# 5) Collect results
# ------------------------------------------------------------

abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
q0 = MOI.get(optimizer, MOI.RawOptimizerAttribute("qa0"))

x0 = SVector(-1.65, -1.65)
nstep = 60

x_traj, u_traj, q_traj = ST.get_closed_loop_trajectory(
    discrete_time_system,
    concrete_controller,
    x0,
    q0,
    nstep;
    update_on_next = true,
    stopping = x -> false,
)

println("Trajectory length: ", length(x_traj.seq))

# ------------------------------------------------------------
# 6) Plot
# ------------------------------------------------------------
using Plots
φ_str = string(φ)
fig = plot(; aspect_ratio = :equal, title = "$φ_str")
plot!(
    concrete_problem;
    ap_colors = Dict(:g1 => :red, :g2 => :cyan, :g3 => :orange, :obs => :black),
    aspect_ratio = :equal,
)
plot!(fig, x_traj; color = :blue, dims = [1, 2])
display(fig)
