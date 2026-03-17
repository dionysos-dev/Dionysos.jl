using StaticArrays
using LinearAlgebra
using JuMP
using Clarabel
using JuMP

import HybridSystems
using LazySets
using Plots

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const OP = DI.Optim
const AB = OP.Abstraction

const PCLF = UT.PathCompleteFramework

# ---------------------------------------------------------
# Define a stable switched system
# ---------------------------------------------------------
A1 = @SMatrix [
    0.70 0.10;
    0.00 0.65
]

A2 = @SMatrix [
    0.60 -0.15;
    0.10 0.55
]

f = HybridSystems.discreteswitchedsystem([Matrix(A1), Matrix(A2)])

# ---------------------------------------------------------
# Define problem:
#  Define region X, terminal region D, and observation regions
# ---------------------------------------------------------

# X = [-2,2]^2
X = Hyperrectangle(; low = [-2.0, -2.0], high = [2.0, 2.0])

# D = [-0.2,0.2]^2
D = Hyperrectangle(; low = [-0.2, -0.2], high = [0.2, 0.2])

# Two observation regions
R1 = Hyperrectangle(; low = [0.8, 0.8], high = [1.5, 1.5])
R2 = Hyperrectangle(; low = [-1.5, 0.8], high = [-0.8, 1.5])

observation_regions = [R1, R2]

problem = PR.BisimulationQuotientProblem(f, X, D, observation_regions)

# ---------------------------------------------------------
# Compute a common polyhedral PCLF (k = 0 gives 1-node graph)
# ---------------------------------------------------------
G = PCLF.generate_DeBruijn_edges(2, 1; dual = false)
sdp_optimizer = JuMP.optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => 1000)

θ = π / 6
R = [
    cos(θ) -sin(θ)
    sin(θ) cos(θ)
]

# Gmats = :identity
Gmats = Dict(
    (1,) => R,   # rotation
    (2,) => R,   # rotation
)

pclf = PCLF.compute_polyhedral_pieces_pclf(
    f,
    G,
    sdp_optimizer;
    Gmats = Gmats,
    MLF = true,
    verbose = false,
)

println("Computed JSR upper bound / contraction rate = ", pclf.JSRapprox)

# ---------------------------------------------------------
# Instantiate abstraction optimizer
# ---------------------------------------------------------
optimizer = MOI.instantiate(AB.PCLFBisimulationQuotient.OptimizerBisimulationQuotient)

MOI.set(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient_problem"), problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("pclf"), pclf)
MOI.set(optimizer, MOI.RawOptimizerAttribute("verbose"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("atol"), 1e-3)
MOI.set(optimizer, MOI.RawOptimizerAttribute("num_levels"), 8)
MOI.set(optimizer, MOI.RawOptimizerAttribute("max_slices"), 4) # nothing

# ---------------------------------------------------------
# Solve
# ---------------------------------------------------------

using StatProfilerHTML
@profilehtml MOI.optimize!(optimizer)
# MOI.optimize!(optimizer)

# MOI.optimize!(optimizer)
construction_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstraction_construction_time_sec"))
println("Construction time = ", construction_time)

const FILENAME = joinpath(@__DIR__, "example1.jld2")
AB.PCLFBisimulationQuotient.export_optimizer_jld2(optimizer, FILENAME)
# optimizer = AB.PCLFBisimulationQuotient.import_optimizer_jld2(FILENAME)

bisimulation = MOI.get(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient"))
obs_partition = MOI.get(optimizer, MOI.RawOptimizerAttribute("obs_partition"))

AB.PCLFBisimulationQuotient.print_bisimulation_stats(bisimulation)

fig = plot(; aspect_ratio = :equal);
# plot!(problem)                # problem geometry
# plot!(obs_partition)          # observation partition
# plot!(bisimulation; what = :slices, node = (1,), show_contours=false)
# plot(T, what=:states, by=:obs, show_contours=false)
# plot!(bisimulation; what = :states, node = 1, by = :obs)
# println(bisimulation.states[1].node)
plot!(bisimulation; what = :states, node = (1,), show_contours = false)
display(fig)

# # ---------------------------------------------------------
# # CoSafe LTL control synthesis on the quotient
# # ---------------------------------------------------------

# using Spot

# # φ = Spot.parse_formula("F(p)")
# # φ = ltl"F(p)"
# φ = ltl"F(p & F(q))"

# # qa = 0  -> not yet reached D
# # qa = 1  -> reached D, stay there
# M = DI.Symbolic.FunctionMonitor(0, Set([1]), (qa, ap) -> (qa == 1 || :D in ap) ? 1 : 0)

# _I_ = Hyperrectangle(; low = [0.39, 0.39], high = [0.41, 0.41])
# prob = PR.CoSafeLTLProblem(
#     f,
#     _I_,
#     M,
#     Dict(:D => D), # no really useful since we have ap_to_obs, but let's be explicit
#     Dict{Symbol, Any}(:D => MP.INNER), # no really useful, but let's be explicit
#     true,
# )

# opt = MOI.instantiate(AB.PCLFBisimulationQuotient.OptimizerCoSafeLTLOnQuotient)
# MOI.set(opt, MOI.RawOptimizerAttribute("concrete_problem"), prob)
# MOI.set(opt, MOI.RawOptimizerAttribute("bisimulation_quotient"), bisimulation)
# MOI.set(opt, MOI.RawOptimizerAttribute("ap_to_obs"), Dict(:D => -1))
# MOI.set(opt, MOI.RawOptimizerAttribute("print_level"), 1)
# MOI.optimize!(opt)

# concrete_controller = AB.PCLFBisimulationQuotient.solve_concrete_problem(opt)

# x0 = SVector(0.4, 0.4)
# mem0 = AB.PCLFBisimulationQuotient.initial_controller_memory(opt, x0)

# (X_seq, U_seq, M_seq) = AB.PCLFBisimulationQuotient.simulate_closed_loop(
#     f,
#     concrete_controller,
#     x0,
#     mem0;
#     N = 50,
# )
# println(X_seq)
# println(U_seq)
# println(M_seq)
# plot!(ST.Trajectory(X_seq); label = "Trajectory")
