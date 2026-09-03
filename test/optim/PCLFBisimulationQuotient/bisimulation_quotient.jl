module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

using LinearAlgebra
using JuMP
import MathOptInterface as MOI
import Clarabel
import HybridSystems
import CDDLib
import Spot
using LazySets

const PCLF = UT.PathCompleteFramework

# Contract test for the PCLF bisimulation-quotient optimizer. Small discrete switched
# linear system with a polyhedral path-complete Lyapunov function; adapted from the
# BisimulationQuotient case-study script.
@testset "PCLF bisimulation quotient" begin
    A1 = (1.0 / 10.0) * [1.5519 0.4474; 7.6412 7.4716]
    A2 = (1.0 / 10.0) * [0.4750 9.1755; 1.8955 0.1850]
    # The synthesis below chooses the switching, so the modes are declared the controller's;
    # `discreteswitchedsystem` alone would declare them autonomous and the solve would answer
    # the universal question instead.
    f = ST.with_switching(
        HybridSystems.discreteswitchedsystem([A1, A2]),
        HybridSystems.ControlledSwitching(),
    )

    p = 2.5
    θ = deg2rad(10.0)
    R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    normals = [R * [1.0, 0.0], R * [-1.0, 0.0], R * [0.0, 1.0], R * [0.0, -1.0]]
    X = HPolytope([HalfSpace(normals[i], p) for i in 1:4])

    R1 = HPolytope([
        HalfSpace([1.0, 0.0], 1.5),
        HalfSpace([-1.0, 0.0], -0.8),
        HalfSpace([0.0, 1.0], 1.5),
        HalfSpace([0.0, -1.0], -0.8),
    ])
    problem = PR.BisimulationQuotientProblem(f, X, [R1])

    graph = PCLF.edgeList_to_LabDigraph([
        (1, 2, 1),
        (2, 1, 1),
        (2, 4, 1),
        (2, 3, 1),
        (3, 4, 1),
        (4, 3, 2),
        (4, 4, 2),
        (4, 1, 2),
    ])
    v1, v2, v3, v4, v5 = [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [-1.0, 1.0], [-1.0, 0.0]
    pieces = [hcat(v1, v2), hcat(v2, v3), hcat(v3, v4), hcat(v4, v5)]
    partitions = Dict(k => pieces for k in 1:4)

    lp_optimizer = JuMP.optimizer_with_attributes(
        Clarabel.Optimizer,
        "max_iter" => 1000,
        MOI.Silent() => true,
    )
    pclf_poly =
        PCLF.compute_polyhedral_pieces_pclf(f, graph, lp_optimizer, partitions; MLF = true)
    @test isfinite(pclf_poly.JSRapprox)

    optimizer = MOI.instantiate(AB.PCLFBisimulationQuotient.OptimizerBisimulationQuotient)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient_problem"), problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("pclf"), pclf_poly)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("atol"), 1e-4)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("level_tol"), 1e-2)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("max_slices"), 6)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("max_levels"), 20)

    MOI.optimize!(optimizer)

    bisimulation = MOI.get(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient"))
    D = MOI.get(optimizer, MOI.RawOptimizerAttribute("D"))
    construction_time =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("construction_time_sec"))

    @test bisimulation !== nothing
    @test D !== nothing
    @test construction_time >= 0.0

    # Co-safe LTL control synthesis on the quotient.
    φ = Spot.@ltl_str "F(R1 & F(D))"
    spec = DI.spot_stepper(φ)
    x0 = SVector(1.0, 1.0)               # inside R1
    _I_ = LazySets.Hyperrectangle(; low = [x0[1], x0[2]], high = [x0[1], x0[2]])
    cosafe_problem = PR.CoSafeLTLProblem(
        f,
        _I_,
        spec,
        Dict(:D => D, :R1 => R1),
        Dict{Symbol, Any}(:D => MP.INNER, :R1 => MP.INNER),
    )

    cosafe_optimizer =
        MOI.instantiate(AB.PCLFBisimulationQuotient.OptimizerCoSafeLTLOnQuotient)
    MOI.set(cosafe_optimizer, MOI.RawOptimizerAttribute("concrete_problem"), cosafe_problem)
    MOI.set(
        cosafe_optimizer,
        MOI.RawOptimizerAttribute("bisimulation_quotient"),
        bisimulation,
    )
    MOI.set(
        cosafe_optimizer,
        MOI.RawOptimizerAttribute("ap_to_obs"),
        Dict(:D => -1, :R1 => 1),
    )
    MOI.set(cosafe_optimizer, MOI.RawOptimizerAttribute("early_stop"), false)
    MOI.set(cosafe_optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.optimize!(cosafe_optimizer)

    concrete_controller =
        AB.PCLFBisimulationQuotient.solve_concrete_problem(cosafe_optimizer)
    controllable_set =
        MOI.get(cosafe_optimizer, MOI.RawOptimizerAttribute("controllable_set"))
    @test concrete_controller !== nothing
    @test controllable_set !== nothing

    # `success` must answer whether the problem's own initial set is controllable, not whether
    # the whole domain is. Under `early_stop = false` the controller is built over every
    # state, so forwarding the sub-solver's flag reported failure as soon as any single state
    # could not satisfy the specification -- whatever was asked about.
    @test cosafe_optimizer.success

    q0 = [q.id for q in values(bisimulation.states) if x0 in q.set]
    @test !isempty(q0)
    @test all(q -> q in controllable_set, q0)

    # ... and it must not depend on `early_stop`, which only changes how much of the domain
    # the controller is constructed over.
    early_stop_optimizer =
        MOI.instantiate(AB.PCLFBisimulationQuotient.OptimizerCoSafeLTLOnQuotient)
    MOI.set(
        early_stop_optimizer,
        MOI.RawOptimizerAttribute("concrete_problem"),
        cosafe_problem,
    )
    MOI.set(
        early_stop_optimizer,
        MOI.RawOptimizerAttribute("bisimulation_quotient"),
        bisimulation,
    )
    MOI.set(
        early_stop_optimizer,
        MOI.RawOptimizerAttribute("ap_to_obs"),
        Dict(:D => -1, :R1 => 1),
    )
    MOI.set(early_stop_optimizer, MOI.RawOptimizerAttribute("early_stop"), true)
    MOI.set(early_stop_optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.optimize!(early_stop_optimizer)
    @test early_stop_optimizer.success == cosafe_optimizer.success

    # ------------------------------------------------------------------
    # The same problem with the switching declared autonomous: the modes become the
    # environment's, the quotient is folded, and the same optimizer answers the universal
    # question — from which states does EVERY switching sequence satisfy the specification?
    # ------------------------------------------------------------------
    verification_problem = PR.remake(
        cosafe_problem;
        system = ST.with_switching(f, HybridSystems.AutonomousSwitching()),
    )

    verification_optimizer =
        MOI.instantiate(AB.PCLFBisimulationQuotient.OptimizerCoSafeLTLOnQuotient)
    MOI.set(
        verification_optimizer,
        MOI.RawOptimizerAttribute("concrete_problem"),
        verification_problem,
    )
    MOI.set(
        verification_optimizer,
        MOI.RawOptimizerAttribute("bisimulation_quotient"),
        bisimulation,
    )
    MOI.set(
        verification_optimizer,
        MOI.RawOptimizerAttribute("ap_to_obs"),
        Dict(:D => -1, :R1 => 1),
    )
    MOI.set(verification_optimizer, MOI.RawOptimizerAttribute("early_stop"), false)
    MOI.set(verification_optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.set(
        verification_optimizer,
        MOI.RawOptimizerAttribute("coverage_backend"),
        CDDLib.Library(),
    )
    MOI.optimize!(verification_optimizer)

    # The atol-erosion caveat is measured, not assumed: a small sliver of the slice family is
    # uncovered by the cells, and the verified set says nothing there.
    uncovered = verification_optimizer.uncovered_fraction
    @test uncovered !== nothing
    @test 0.0 <= uncovered < 0.1

    verified_set =
        MOI.get(verification_optimizer, MOI.RawOptimizerAttribute("controllable_set"))

    # A run every environment survives is in particular a run some controller wins, so the
    # verified set is contained in the synthesised one — and on this system strictly, since
    # the synthesis had to steer.
    @test verified_set ⊆ controllable_set
    @test length(verified_set) < length(controllable_set)

    # Verification returns a set, never a controller.
    @test verification_optimizer.environment_folded
    @test_throws ErrorException AB.PCLFBisimulationQuotient.solve_concrete_problem(
        verification_optimizer,
    )

    # ------------------------------------------------------------------
    # A failed verification owes evidence: from a state that synthesis wins but verification
    # does not, there is a switching word no controller survives.
    # ------------------------------------------------------------------
    gap = setdiff(Set(controllable_set), Set(verified_set))
    @test !isempty(gap)

    find_point(qid_set) = begin
        for xx in range(-2.0, 2.0; length = 41), yy in range(-2.0, 2.0; length = 41)
            pt = [xx, yy]
            for qid in qid_set
                pt ∈ bisimulation.states[qid].set && return pt
            end
        end
        return nothing
    end

    x0_gap = find_point(gap)
    @test x0_gap !== nothing

    cex = AB.PCLFBisimulationQuotient.verification_counterexample(
        verification_optimizer,
        x0_gap,
    )
    @test !isempty(cex.modes)
    @test all(m -> m in 1:2, cex.modes)
    # Either the environment loops forever outside the verified set, or it drives the run to
    # where the abstraction's coverage ends. Both are evidence; neither is silence.
    @test cex.entered_sink || cex.lasso_start >= 1
    @test x0_gap ∈ bisimulation.states[first(cex.qids)].set
    @test length(cex.X) == length(cex.modes) + 1

    # No counterexample exists from a verified state — when this small system verifies anything
    # at all: a two-mode adversary on a tight domain can legitimately leave the verified set
    # empty, in which case there is no such state to probe.
    x0_ok = isempty(verified_set) ? nothing : find_point(Set(verified_set))
    if x0_ok !== nothing
        @test_throws ErrorException AB.PCLFBisimulationQuotient.verification_counterexample(
            verification_optimizer,
            x0_ok,
        )
    end

    # And none from a synthesis run.
    @test_throws ErrorException AB.PCLFBisimulationQuotient.verification_counterexample(
        cosafe_optimizer,
        x0_gap,
    )
end

end # module TestMain
