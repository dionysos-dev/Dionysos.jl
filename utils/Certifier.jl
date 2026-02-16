using StaticArrays
using JuMP
using Plots
import MathOptInterface as MOI

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const PR = DI.Problem
const OP = DI.Optim
const AB = OP.Abstraction

# Certifier module path inside Dionysos
const SC = OP.SymbolicCertifier

# --- Toy problem definition ---
include("../problems/toy_problem.jl")

# ------------------------------------------------------------
# 1) Build the optimizer
# ------------------------------------------------------------
concrete_problem = ToyProblem.optimal_control_problem()

# State discretization
hx = SVector(0.1, 0.1)
# Input discretization
u0 = SVector(0.0, 0.0)
hu = SVector(0.2, 0.2)
# Time discretization
tstep = 0.2

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("h"), hx)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), DO.GridFree(u0, hu))
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("jacobian_bound"),
    ToyProblem.jacobian_bound(),
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), tstep)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.GROWTH, # GROWTH, CENTER_SIMULATION
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("early_stop"), true)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("automaton_constructor"),
    (n, m) -> SY.NewIndexedAutomatonList(n, m),
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)

# ------------------------------------------------------------
# 2) Build a *candidate* trajectory (simple heuristic)
# ------------------------------------------------------------
center(rect::UT.HyperRectangle) = (rect.lb + rect.ub) / 2

# Sample initial point
x0 = SVector(UT.get_center(concrete_problem.initial_set)...)
target_set = concrete_problem.target_set
xT = center(target_set)

# Build a discrete-time system for simulation via growth-bound discretization (same idea as pendulum)
ctsys = concrete_problem.system
ctapprox = ST.ContinuousTimeGrowthBound_from_jacobian_bound(ctsys, ToyProblem.jacobian_bound())
dtapprox = ST.discretize(ctapprox, tstep)
dtsys = ST.get_system(dtapprox)

# Input constraints (assume HyperRectangle for toy; if not, adapt)
Uset = concrete_problem.system.U

# Clamp utility for HyperRectangle inputs
function clamp_to_U(u::SVector{2,Float64}, U::UT.HyperRectangle)
    return SVector(
        clamp(u[1], U.lb[1], U.ub[1]),
        clamp(u[2], U.lb[2], U.ub[2]),
    )
end

nstep = 60

# roll out open-loop
x_seq = Vector{SVector{2,Float64}}(undef, nstep + 1)
u_seq = Vector{SVector{2,Float64}}(undef, nstep)

x_seq[1] = SVector{2,Float64}(x0)

reached(x) = (x ∈ target_set)

function rollout_candidate(x0, nstep, xT, tstep, dtsys, Uset, target_set)
    x_seq = Vector{SVector{2,Float64}}(undef, nstep + 1)
    u_seq = Vector{SVector{2,Float64}}(undef, nstep)
    x_seq[1] = SVector{2,Float64}(x0)

    reached(x) = (x ∈ target_set)

    for k in 1:nstep
        xk = x_seq[k]
        u_des = (xT - xk) / tstep
        uk = clamp_to_U(SVector{2,Float64}(u_des), Uset)
        u_seq[k] = uk

        xnext = dtsys.f(xk, uk)
        x_seq[k + 1] = xnext

        if reached(xnext)
            return x_seq[1:(k + 1)], u_seq[1:k]
        end
    end

    return x_seq, u_seq
end

x_seq, u_seq = rollout_candidate(x0, nstep, xT, tstep, dtsys, Uset, target_set)
ts = collect(0.0:tstep:((length(x_seq) - 1) * tstep))
traj = (; X = x_seq, U = u_seq, t = ts)
x_traj = ST.Trajectory(traj.X)

println("Candidate trajectory length: ", length(traj.X))
println("Candidate final state: ", traj.X[end])
println("Reached target (candidate): ", reached(traj.X[end]))

# ------------------------------------------------------------
# 3) Call your local tube certifier
# ------------------------------------------------------------
cert = SC.UniformGridLocalTubeCertifier()

SC.set_optimizer!(cert, optimizer)
SC.set_trajectory!(cert, x_traj)
cert.radius = SVector(0.25, 0.25)
cert.margin = 0.0
cert.incl_mode = DO.OUTER

SC.certify!(cert)

println("\n=== Local Certification Result ===")
println("success:    ", SC.get_success(cert))
println("time (sec): ", SC.get_solve_time(cert))
# controller = SC.get_controller(cert)
# println("controller: ", controller === nothing ? "nothing" : string(typeof(controller)))

# ------------------------------------------------------------
# 4) Plot 
# ------------------------------------------------------------
fig = plot(; aspect_ratio = :equal, title = "ToyExample: candidate traj + sets")

plot!(
    concrete_problem.system.X;
    color = :grey,
    opacity = 0.15,
    label = "X",
)

plot!(
    concrete_problem.initial_set;
    color = :green,
    opacity = 0.25,
    label = "Initial set",
)

plot!(
    concrete_problem.target_set;
    color = :red,
    opacity = 0.35,
    label = "Target set",
)

# plot!(
#     cert.optimizer.abstraction_solver.abstraction_region;
#     color = :blue,
#     opacity = 0.35,
#     label = "Tube",
# )

plot!(
    cert.optimizer.abstraction_solver.abstract_system;
    color = :blue,
    opacity = 0.35,
    label = "Tube",
    efficient = false,
)

plot!(x_traj; ms = 2.0, arrows = false, label = "Candidate")

display(fig)
