using StaticArrays
using MathematicalSystems
using Dionysos
using JuMP
import MathOptInterface as MOI

# using Spot

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction
const OPDS = OP.DiscreteSystems

# ------------------------------------------------------------
# 1) Define a simple 2D continuous-time system: x' = u
# ------------------------------------------------------------

include(
    joinpath(dirname(dirname(pathof(Dionysos))), "problems", "Integrator", "integrator.jl"),
)

_X_ = UT.box(SVector(-2.0, -2.0), SVector(2.0, 2.0))
_U_ = UT.box(SVector(-1.0, -1.0), SVector(1.0, 1.0))

concrete_system = Integrator.system(; _X_ = _X_, _U_ = _U_)
jacobian_bound = Integrator.jacobian_bound()

# ------------------------------------------------------------
# 2) Abstraction construction (AlternatingSimulationProblem)
# ------------------------------------------------------------

alternating_simulation_problem =
    DI.Problem.AlternatingSimulationProblem(concrete_system, concrete_system.X)

# grid resolution
x0 = SVector(-2.0, -2.0)
hx = SVector(0.2, 0.2)
state_grid = MP.GridFree(x0, hx)

u0 = SVector(-1.0, -1.0)
hu = SVector(0.5, 0.5)
input_grid = MP.GridFree(u0, hu)

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("concrete_problem"),
    alternating_simulation_problem,
)
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

MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)
MOI.set(optimizer, MOI.Silent(), true)

MOI.optimize!(optimizer)

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
discrete_time_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system"))

println("Abstraction built.")

# ------------------------------------------------------------
# 3) Define co-safe LTL problem with sets labeling
# ------------------------------------------------------------

_I_ = UT.box(SVector(-1.7, -1.7), SVector(-1.6, -1.6))

g1 = UT.box(SVector(1.0, 1.0), SVector(1.7, 1.7))
g2 = UT.box(SVector(-1.5, -1.2), SVector(-0.6, -0.2))

obs = UT.box(SVector(-1.8, 0.0), SVector(-0.6, 1.0))

danger1 = UT.box(SVector(-0.5, -0.5), SVector(0.5, 0.5))
danger2 = UT.box(SVector(1.3, -0.5), SVector(2.0, 0.5))
danger = UT.set_union([danger1, danger2])

# co-safe formula
# φ = ltl"G(!obs) & F(g1 & ((!danger) U g2))"
φ = "G(!obs) & F(g1 & ((!danger) U g2))"

struct MonitorG1NoDangerUntilG2 end

@inline function mon_next(::MonitorG1NoDangerUntilG2, q::Int, ap::Tuple{Vararg{Symbol}})
    obs = (:obs in ap)
    g1 = (:g1 in ap)
    g2 = (:g2 in ap)
    danger = (:danger in ap)

    obs && return 0                 # always forbidden

    q == 0 && return 0              # dead absorbs
    q == 3 && return 3              # done absorbs

    if q == 1
        # before g1
        if g1
            # reaching g1 starts phase2; if g2 is already true, we can finish immediately
            return g2 ? 3 : 2
        else
            return 1
        end
    end

    @assert q == 2
    # after g1: danger forbidden until g2
    if g2
        return 3
    elseif danger
        return 0
    else
        return 2
    end
end

mon = OPDS.FunctionMonitor(
    1,         # initial
    Set([3]),  # accepting
    (qa, ap) -> mon_next(MonitorG1NoDangerUntilG2(), qa, ap),
)

labeling = Dict{Symbol, Any}(:g1 => g1, :g2 => g2, :danger => danger, :obs => obs)

# semantics per AP
ap_semantics = Dict{Symbol, Any}(
    :g1 => MP.INNER,
    :g2 => MP.INNER,
    :danger => MP.OUTER,
    :obs => MP.OUTER,
)

concrete_problem = DI.Problem.CoSafeLTLProblem(
    concrete_system,
    _I_,
    mon, # φ, mon
    labeling,
    ap_semantics,
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
x0 = SVector(-1.65, -1.65)
nstep = 60

traj = ST.get_closed_loop_trajectory(
    discrete_time_system,
    concrete_controller,
    x0,
    nstep;
    update_on_next = true,
    stopping = x -> false,
)

println("Trajectory length: ", length(ST.states(traj)))

# ------------------------------------------------------------
# 6) Plot
# ------------------------------------------------------------
using Plots
φ_str = string(φ)
fig = plot(; aspect_ratio = :equal, title = "$φ_str")
plot!(
    concrete_problem;
    ap_colors = Dict(:g1 => :red, :g2 => :cyan, :danger => :orange, :obs => :black),
    aspect_ratio = :equal,
)
plot!(fig, traj; color = :blue, dims = [1, 2])
display(fig)
