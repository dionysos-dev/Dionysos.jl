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
# ---------------------------------------------------------

X = Hyperrectangle(; low = [-2.0, -2.0], high = [2.0, 2.0])
R1 = Hyperrectangle(; low = [0.5, 0.5], high = [0.9, 0.9])
R2 = Hyperrectangle(; low = [0.5, -0.5], high = [0.9, -0.1])
observation_regions = [R1, R2]
problem = PR.BisimulationQuotientProblem(f, X, observation_regions)

# ---------------------------------------------------------
# Compute a common polyhedral PCLF (k = 0 gives 1-node graph)
# ---------------------------------------------------------

G = PCLF.generate_DeBruijn_edges(2, 1; dual = false)
sdp_optimizer = JuMP.optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => 1000)

Gmats = :identity
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
MOI.set(optimizer, MOI.RawOptimizerAttribute("max_levels"), 100)
MOI.set(optimizer, MOI.RawOptimizerAttribute("max_slices"), 8)

# ---------------------------------------------------------
# Solve
# ---------------------------------------------------------
MOI.optimize!(optimizer)
construction_time = MOI.get(optimizer, MOI.RawOptimizerAttribute("construction_time_sec"))
println("Construction time = ", construction_time)

# const FILENAME = joinpath(@__DIR__, "example2.jld2")
# AB.PCLFBisimulationQuotient.export_optimizer_jld2(optimizer, FILENAME)
# optimizer = AB.PCLFBisimulationQuotient.import_optimizer_jld2(FILENAME)

bisimulation = MOI.get(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient"))

AB.PCLFBisimulationQuotient.print_bisimulation_stats(bisimulation)

fig = plot(; aspect_ratio = :equal);
plot!(problem)
plot!(bisimulation; what = :states, node = (1,), by = :obs)
display(fig)

# ---------------------------------------------------------
# CoSafe LTL control synthesis on the quotient
# ---------------------------------------------------------
φ = ltl"(F(D)) & (G(!R1))"
φ = ltl"F(R1 & F(D))"

M = DI.Symbolic.FunctionMonitor(0, Set([2]), (qa, ap) -> begin
    if qa == 2
        return 2
    elseif qa == 0
        return (:R1 in ap) ? 1 : 0
    elseif qa == 1
        return (:D in ap) ? 2 : 1
    else
        return qa
    end
end)

_I_ = Hyperrectangle(; low = [1.3, 1.3], high = [1.5, 1.5])
prob = PR.CoSafeLTLProblem(
    f,
    _I_,
    φ,
    Dict(:R1 => R1, :D => D),
    Dict{Symbol, Any}(:R1 => MP.INNER, :D => MP.INNER),
    true,
)

opt = MOI.instantiate(AB.PCLFBisimulationQuotient.OptimizerCoSafeLTLOnQuotient)
MOI.set(opt, MOI.RawOptimizerAttribute("concrete_problem"), prob)
MOI.set(opt, MOI.RawOptimizerAttribute("bisimulation_quotient"), bisimulation)
MOI.set(opt, MOI.RawOptimizerAttribute("ap_to_obs"), Dict(:R1 => 1, :D => -1))
MOI.set(opt, MOI.RawOptimizerAttribute("print_level"), 1)
MOI.optimize!(opt)

concrete_controller = AB.PCLFBisimulationQuotient.solve_concrete_problem(opt)
mem0 = AB.PCLFBisimulationQuotient.initial_controller_memory(opt, x0)

x0 = SVector(1.4, 1.4)
(X_seq, U_seq, M_seq) = AB.PCLFBisimulationQuotient.simulate_closed_loop(
    f,
    concrete_controller,
    x0,
    mem0;
    N = 50,
)
println(X_seq)
println(U_seq)
println(M_seq)
φ_str = string(φ)
plot!(fig; title = "$φ_str")
plot!(ST.Trajectory(X_seq); label = "Trajectory")
