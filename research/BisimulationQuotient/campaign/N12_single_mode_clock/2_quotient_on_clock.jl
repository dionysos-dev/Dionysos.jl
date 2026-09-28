# N12 step 1: build the bisimulation, not just the certificate.
#
# `1_certificate_separation.jl` settles the Lyapunov side -- one node at four facets certifies nothing for
# x⁺ = ρ R(2π/q) x, while the q-cycle attains ρ exactly. That is a statement about certificates. The
# claim this folder actually makes is about ABSTRACTIONS, so the quotient has to be built and a
# specification answered on it, and the single-node arm has to be shown to have nothing to build FROM.
#
# Two properties of this setting are worth stating before the numbers, because they are what make it
# a test rather than a demonstration:
#
#   ∃ and ∀ coincide.  There is one mode, so there is no switching signal for anyone to choose. The
#                      controller has no decision and the adversary has no move. Running both and
#                      getting different answers would mean the semantics are wrong somewhere; this
#                      script asserts they agree, which is a free correctness check the switched
#                      experiments cannot perform.
#   the node is a clock. With M = 1 the memory is the phase k mod q. The certificate decreases every
#                      q steps, while the abstraction still respects PER-STEP transitions and
#                      per-step observations -- which is exactly why the q steps cannot be lumped
#                      into one transition and the result is not a restatement of a periodic
#                      Lyapunov function.
#
# A note on what this graph does NOT exercise. The q-cycle is over a ONE-letter alphabet, so every
# node enables the whole alphabet and `PCLF.restricts_future` is false: the completion defect fixed
# in `notes/constrained-language-review.md` never touched this case. Its periodic-LTV sibling
# (`3_periodic_ltv.jl`) is the opposite -- node k enables only letter k -- and that one would have
# returned an empty universal answer before the fix.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))

using Spot
using LinearAlgebra

include(joinpath(@__DIR__, "1_certificate_separation.jl"))

const Q = 5
const RHO = 0.9

f = ST.with_switching(one_mode_system(RHO, Q), HybridSystems.ControlledSwitching())

# The working set and two regions to talk about. Kept well inside the sublevel family so the
# observation refinement has something to cut, and away from the origin so reaching a region is a
# real obligation rather than a consequence of contraction.
X = LazySets.Hyperrectangle(; low = [-6.0, -6.0], high = [6.0, 6.0])
R1 = LazySets.Hyperrectangle(; low = [2.0, 0.5], high = [4.0, 2.5])
R2 = LazySets.Hyperrectangle(; low = [-4.0, -2.5], high = [-2.0, -0.5])
problem = PR.BisimulationQuotientProblem(f, X, [R1, R2])

println("x⁺ = ρ R(2π/q) x with q = ", Q, ", ρ = ", RHO, " -- ONE mode, no switching\n")

# ---------------------------------------------------------
# The single-node arm has nothing to build from
# ---------------------------------------------------------

one_node = single_node_rate(RHO, Q)
println(
    "single node, four facets, best of 41 orientations: rate = ",
    isfinite(one_node) ? round(one_node; digits = 5) : "none",
)
one_node >= 1.0 || error(
    "the single node certified at $one_node < 1; the separation this experiment rests on is gone",
)
println(
    "  >= 1 at a four-facet budget: GLB's single node has no certificate and nothing to build.\n" *
    "  The INDUCED common is a DIFFERENT baseline and does start, at 4q facets -- and on this\n" *
    "  complete graph it is the faster of the two. See 3_route1_vs_route2.jl; N12 is a feasibility\n" *
    "  result, not a cost result.\n",
)

# ---------------------------------------------------------
# The cycle arm builds
# ---------------------------------------------------------

pclf = cycle_pclf(RHO, Q)
println("q-cycle certificate rate = ", pclf.JSRapprox, " (analytic, exactly ρ)")
println(
    "graph: ",
    length(pclf.graph.verts),
    " nodes, restricts_future = ",
    PCLF.restricts_future(pclf.graph),
    " (one letter, so every node enables all of it)\n",
)

const CACHE = joinpath(@__DIR__, "clock_q$(Q)_rho$(RHO).jld2")

quotient, D = if isfile(CACHE)
    println("reusing the cached quotient from ", basename(CACHE))
    opt = import_optimizer_jld2(CACHE)
    MOI.get(opt, MOI.RawOptimizerAttribute("bisimulation_quotient")),
    MOI.get(opt, MOI.RawOptimizerAttribute("D"))
else
    r = build_quotient(
        problem,
        pclf;
        atol = 1e-4,
        level_tol = 1e-2,
        max_slices = 30,
        print_level = 1,
    )
    export_optimizer_jld2(r.optimizer, CACHE)
    r.quotient, r.D
end

println("\nquotient: ", length(quotient.states), " cells on ", Q, " clock phases")
println(
    "slices: ",
    PCQ.num_slices(quotient),
    "   out-degree: ",
    PCQ.outgoing_degree_stats(quotient),
)

bynode = Dict{Any, Int}()
for s in values(quotient.states)
    bynode[s.node] = get(bynode, s.node, 0) + 1
end
println(
    "cells per phase: ",
    join(["$k => $v" for (k, v) in sort(collect(bynode); by = string ∘ first)], " · "),
)

# ---------------------------------------------------------
# A specification, answered both ways
# ---------------------------------------------------------

# Reach R1 while avoiding R2 until then. Co-safe, and it constrains the ORDER of events, which is the
# class N14 showed an over-approximating grid handles badly.
φ = ltl"(!R2 U R1)"

function certified_set(system)
    p = PR.CoSafeLTLProblem(
        system,
        LazySets.Hyperrectangle(; low = [3.0, 1.0], high = [3.0, 1.0]),
        Dionysos.spot_stepper(φ),
        Dict(:D => D, :R1 => R1, :R2 => R2),
        Dict{Symbol, Any}(ap => MP.INNER for ap in (:D, :R1, :R2)),
    )
    o = MOI.instantiate(PCQ.OptimizerCoSafeLTLOnQuotient)
    MOI.set(o, MOI.RawOptimizerAttribute("concrete_problem"), p)
    MOI.set(o, MOI.RawOptimizerAttribute("bisimulation_quotient"), quotient)
    MOI.set(o, MOI.RawOptimizerAttribute("ap_to_obs"), Dict(:D => -1, :R1 => 1, :R2 => 2))
    MOI.set(o, MOI.RawOptimizerAttribute("early_stop"), false)
    MOI.set(o, MOI.RawOptimizerAttribute("switching_semantics"), :language)
    MOI.set(o, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.optimize!(o)
    return (
        set = MOI.get(o, MOI.RawOptimizerAttribute("controllable_set")),
        completions = MOI.get(o, MOI.RawOptimizerAttribute("num_completions")),
    )
end

println("\nsolving (!R2 U R1) ...")
exists = certified_set(f)
forall = certified_set(ST.with_switching(f, HybridSystems.AutonomousSwitching()))

vol(set) = sum(PCQ.get_volume(quotient, set; backend = CDDLib.Library()))
v_exists, v_forall = vol(exists.set), vol(forall.set)
v_total = vol(collect(keys(quotient.states)))

println(
    "∃ (synthesis)    ",
    rpad(length(exists.set), 8),
    "cells   volume ",
    round(v_exists; digits = 3),
)
println(
    "∀ (verification) ",
    rpad(length(forall.set), 8),
    "cells   volume ",
    round(v_forall; digits = 3),
    "   completions ",
    forall.completions,
)
println(
    "covered          ",
    rpad(length(quotient.states), 8),
    "cells   volume ",
    round(v_total; digits = 3),
)

# THE check this setting affords. One mode means no choice for anyone, so the two semantics are the
# same question asked twice. A discrepancy would be a bug in the fold, not a property of the system.
if Set(exists.set) == Set(forall.set)
    println("\n∃ == ∀ exactly, as a one-mode system requires. The fold is consistent here.")
else
    only_e = length(setdiff(Set(exists.set), Set(forall.set)))
    only_f = length(setdiff(Set(forall.set), Set(exists.set)))
    error(
        "∃ and ∀ DIFFER on a single-mode system ($only_e only-∃, $only_f only-∀). " *
        "With one mode neither player has a decision, so this is a defect in the fold.",
    )
end

println(
    """

Route 1 returns nothing here and route 2 returns a finite bisimulation answering a co-safe LTL
formula, on a system with no switching, no input and no adversary -- a single linear map. The graph
node is the phase k mod $Q, and it is the only thing that makes the certificate exist.""",
)
