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
const SC = OP.Symbolic_certifier

# --- Toy problem definition ---
include("../problems/toy_problem.jl")

# ------------------------------------------------------------
# 1) Build the concrete problem
# ------------------------------------------------------------
concrete_problem = ToyProblem.optimal_control_problem()

# state grid resolution (2D)
hx = SVector(0.1, 0.1)

# input grid (2D)
u0 = SVector(0.0, 0.0)
hu = SVector(0.2, 0.2)

tstep = 0.2

target_set = concrete_problem.target_set

# ------------------------------------------------------------
# 2) Build a *candidate* trajectory (simple heuristic)
# ------------------------------------------------------------
# We'll generate a candidate by simulating an open-loop "go-to-target center" input:
# u_k = clamp((x_target_center - x_k)/tstep, U bounds)
#
# This is not optimal, but usually produces a reaching trajectory for ẋ=u.

# Helper: center of a HyperRectangle
center(rect::UT.HyperRectangle) = (rect.lb + rect.ub) / 2

xT = center(target_set)

# Sample initial point
x0 = SVector(UT.sample(concrete_problem.initial_set)...)

# Build a discrete-time system for simulation via growth-bound discretization (same idea as pendulum)
ctsys = concrete_problem.system
ctapprox = ST.ContinuousTimeGrowthBound_from_jacobian_bound(ctsys, TE.jacobian_bound())
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
X = Vector{SVector{2,Float64}}(undef, nstep + 1)
U = Vector{SVector{2,Float64}}(undef, nstep)

X[1] = SVector{2,Float64}(x0)

reached(x) = (x ∈ target_set)

for k in 1:nstep
    xk = X[k]
    # naive steering
    u_des = (xT - xk) / tstep
    uk = clamp_to_U(SVector{2,Float64}(u_des), Uset)
    U[k] = uk

    # one step using the discrete-time system
    xnext = dtsys.f(xk, uk)
    X[k + 1] = xnext

    if reached(xnext)
        # truncate
        X = X[1:(k + 1)]
        U = U[1:k]
        break
    end
end

ts = collect(0.0:tstep:((length(X) - 1) * tstep))
traj = (; X = X, U = U, t = ts)

println("Candidate trajectory length: ", length(traj.X))
println("Candidate final state: ", traj.X[end])
println("Reached target (candidate): ", reached(traj.X[end]))

# ------------------------------------------------------------
# 3) Call your local tube certifier
# ------------------------------------------------------------
cert = SC.CertifierUniformGridLocalTube()

SC.set_problem!(cert, concrete_problem)
SC.set_trajectory!(cert, traj)

# Tube radius (2D)
cert.radius = SVector(0.25, 0.25)
cert.margin = 0.0

# Local abstraction settings
cert.h = hx
cert.input_grid = DO.GridFree(u0, hu)
cert.time_step = tstep
cert.approx_mode = AB.UniformGridAbstraction.GROWTH
cert.threaded = false
cert.print_level = 2

# No periodic settings for the toy system

SC.certify!(cert)

println("\n=== Local Certification Result ===")
println("success:    ", SC.get_success(cert))
println("time (sec): ", SC.get_solve_time(cert))
controller = SC.get_controller(cert)
println("controller: ", controller === nothing ? "nothing" : string(typeof(controller)))

# ------------------------------------------------------------
# 4) Plot (toy: 2D sets + trajectory)
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

plot!(traj.X; ms = 2.0, arrows = false, label = "Candidate")

display(fig)
