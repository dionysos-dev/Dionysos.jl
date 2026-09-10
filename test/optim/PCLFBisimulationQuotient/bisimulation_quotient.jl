module TestMain

import Dionysos


include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

const PCQ = AB.PCLFBisimulationQuotient

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
#@testset "PCLF bisimulation quotient" begin
#    A1 = (1.0 / 10.0) * [1.5519 0.4474; 7.6412 7.4716]
#    A2 = (1.0 / 10.0) * [0.4750 9.1755; 1.8955 0.1850]
#    # The synthesis below chooses the switching, so the modes are declared the controller's;
#    # `discreteswitchedsystem` alone would declare them autonomous and the solve would answer
#    # the universal question instead.
#    f = ST.with_switching(
#        HybridSystems.discreteswitchedsystem([A1, A2]),
#        HybridSystems.ControlledSwitching(),
#    )
#
#    p = 2.5
#    θ = deg2rad(10.0)
#    R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
#    normals = [R * [1.0, 0.0], R * [-1.0, 0.0], R * [0.0, 1.0], R * [0.0, -1.0]]
#    X = HPolytope([HalfSpace(normals[i], p) for i in 1:4])
#
#    R1 = HPolytope([
#        HalfSpace([1.0, 0.0], 1.5),
#        HalfSpace([-1.0, 0.0], -0.8),
#        HalfSpace([0.0, 1.0], 1.5),
#        HalfSpace([0.0, -1.0], -0.8),
#    ])
#    problem = PR.BisimulationQuotientProblem(f, X, [R1])
#
#    graph = PCLF.edgeList_to_LabDigraph([
#        (1, 2, 1),
#        (2, 1, 1),
#        (2, 4, 1),
#        (2, 3, 1),
#        (3, 4, 1),
#        (4, 3, 2),
#        (4, 4, 2),
#        (4, 1, 2),
#    ])
#    v1, v2, v3, v4, v5 = [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [-1.0, 1.0], [-1.0, 0.0]
#    pieces = [hcat(v1, v2), hcat(v2, v3), hcat(v3, v4), hcat(v4, v5)]
#    partitions = Dict(k => pieces for k in 1:4)
#
#    lp_optimizer = JuMP.optimizer_with_attributes(
#        Clarabel.Optimizer,
#        "max_iter" => 1000,
#        MOI.Silent() => true,
#    )
#    pclf_poly =
#        PCLF.compute_polyhedral_pieces_pclf(f, graph, lp_optimizer, partitions; MLF = true)
#    @test isfinite(pclf_poly.JSRapprox)
#
#    optimizer = MOI.instantiate(PCQ.OptimizerBisimulationQuotient)
#    MOI.set(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient_problem"), problem)
#    MOI.set(optimizer, MOI.RawOptimizerAttribute("pclf"), pclf_poly)
#    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
#    MOI.set(optimizer, MOI.RawOptimizerAttribute("atol"), 1e-4)
#    MOI.set(optimizer, MOI.RawOptimizerAttribute("level_tol"), 1e-2)
#    MOI.set(optimizer, MOI.RawOptimizerAttribute("max_slices"), 6)
#    MOI.set(optimizer, MOI.RawOptimizerAttribute("max_levels"), 20)
#
#    MOI.optimize!(optimizer)
#
#    bisimulation = MOI.get(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient"))
#    D = MOI.get(optimizer, MOI.RawOptimizerAttribute("D"))
#    construction_time =
#        MOI.get(optimizer, MOI.RawOptimizerAttribute("construction_time_sec"))
#
#    @test bisimulation !== nothing
#    @test D !== nothing
#    @test construction_time >= 0.0
#
#    # Co-safe LTL control synthesis on the quotient.
#    φ = Spot.@ltl_str "F(R1 & F(D))"
#    spec = DI.spot_stepper(φ)
#    x0 = SVector(1.0, 1.0)               # inside R1
#    _I_ = LazySets.Hyperrectangle(; low = [x0[1], x0[2]], high = [x0[1], x0[2]])
#    cosafe_problem = PR.CoSafeLTLProblem(
#        f,
#        _I_,
#        spec,
#        Dict(:D => D, :R1 => R1),
#        Dict{Symbol, Any}(:D => MP.INNER, :R1 => MP.INNER),
#    )
#
#    cosafe_optimizer =
#        MOI.instantiate(PCQ.OptimizerCoSafeLTLOnQuotient)
#    MOI.set(cosafe_optimizer, MOI.RawOptimizerAttribute("concrete_problem"), cosafe_problem)
#    MOI.set(
#        cosafe_optimizer,
#        MOI.RawOptimizerAttribute("bisimulation_quotient"),
#        bisimulation,
#    )
#    MOI.set(
#        cosafe_optimizer,
#        MOI.RawOptimizerAttribute("ap_to_obs"),
#        Dict(:D => -1, :R1 => 1),
#    )
#    MOI.set(cosafe_optimizer, MOI.RawOptimizerAttribute("early_stop"), false)
#    MOI.set(cosafe_optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
#    MOI.optimize!(cosafe_optimizer)
#
#    concrete_controller =
#        PCQ.solve_concrete_problem(cosafe_optimizer)
#    controllable_set =
#        MOI.get(cosafe_optimizer, MOI.RawOptimizerAttribute("controllable_set"))
#    @test concrete_controller !== nothing
#    @test controllable_set !== nothing
#
#    # `success` must answer whether the problem's own initial set is controllable, not whether
#    # the whole domain is. Under `early_stop = false` the controller is built over every
#    # state, so forwarding the sub-solver's flag reported failure as soon as any single state
#    # could not satisfy the specification -- whatever was asked about.
#    @test cosafe_optimizer.success
#
#    q0 = [q.id for q in values(bisimulation.states) if x0 in q.set]
#    @test !isempty(q0)
#    @test all(q -> q in controllable_set, q0)
#
#    # ... and it must not depend on `early_stop`, which only changes how much of the domain
#    # the controller is constructed over.
#    early_stop_optimizer =
#        MOI.instantiate(PCQ.OptimizerCoSafeLTLOnQuotient)
#    MOI.set(
#        early_stop_optimizer,
#        MOI.RawOptimizerAttribute("concrete_problem"),
#        cosafe_problem,
#    )
#    MOI.set(
#        early_stop_optimizer,
#        MOI.RawOptimizerAttribute("bisimulation_quotient"),
#        bisimulation,
#    )
#    MOI.set(
#        early_stop_optimizer,
#        MOI.RawOptimizerAttribute("ap_to_obs"),
#        Dict(:D => -1, :R1 => 1),
#    )
#    MOI.set(early_stop_optimizer, MOI.RawOptimizerAttribute("early_stop"), true)
#    MOI.set(early_stop_optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
#    MOI.optimize!(early_stop_optimizer)
#    @test early_stop_optimizer.success == cosafe_optimizer.success
#
#    # ------------------------------------------------------------------
#    # The same problem with the switching declared autonomous: the modes become the
#    # environment's, the quotient is folded, and the same optimizer answers the universal
#    # question — from which states does EVERY switching sequence satisfy the specification?
#    # ------------------------------------------------------------------
#    verification_problem = PR.remake(
#        cosafe_problem;
#        system = ST.with_switching(f, HybridSystems.AutonomousSwitching()),
#    )
#
#    verification_optimizer =
#        MOI.instantiate(PCQ.OptimizerCoSafeLTLOnQuotient)
#    MOI.set(
#        verification_optimizer,
#        MOI.RawOptimizerAttribute("concrete_problem"),
#        verification_problem,
#    )
#    MOI.set(
#        verification_optimizer,
#        MOI.RawOptimizerAttribute("bisimulation_quotient"),
#        bisimulation,
#    )
#    MOI.set(
#        verification_optimizer,
#        MOI.RawOptimizerAttribute("ap_to_obs"),
#        Dict(:D => -1, :R1 => 1),
#    )
#    MOI.set(verification_optimizer, MOI.RawOptimizerAttribute("early_stop"), false)
#    MOI.set(verification_optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
#    MOI.set(
#        verification_optimizer,
#        MOI.RawOptimizerAttribute("coverage_backend"),
#        CDDLib.Library(),
#    )
#    MOI.optimize!(verification_optimizer)
#
#    # The atol-erosion caveat is measured, not assumed: a small sliver of the slice family is
#    # uncovered by the cells, and the verified set says nothing there.
#    uncovered = verification_optimizer.uncovered_fraction
#    @test uncovered !== nothing
#    @test 0.0 <= uncovered < 0.1
#
#    verified_set =
#        MOI.get(verification_optimizer, MOI.RawOptimizerAttribute("controllable_set"))
#
#    # A run every environment survives is in particular a run some controller wins, so the
#    # verified set is contained in the synthesised one — and on this system strictly, since
#    # the synthesis had to steer.
#    @test verified_set ⊆ controllable_set
#    @test length(verified_set) < length(controllable_set)
#
#    # Verification returns a set, never a controller.
#    @test verification_optimizer.environment_folded
#    @test_throws ErrorException PCQ.solve_concrete_problem(
#        verification_optimizer,
#    )
#
#    # ------------------------------------------------------------------
#    # A failed verification owes evidence: from a state that synthesis wins but verification
#    # does not, there is a switching word no controller survives.
#    # ------------------------------------------------------------------
#    gap = setdiff(Set(controllable_set), Set(verified_set))
#    @test !isempty(gap)
#
#    find_point(qid_set) = begin
#        for xx in range(-2.0, 2.0; length = 41), yy in range(-2.0, 2.0; length = 41)
#            pt = [xx, yy]
#            for qid in qid_set
#                pt ∈ bisimulation.states[qid].set && return pt
#            end
#        end
#        return nothing
#    end
#
#    x0_gap = find_point(gap)
#    @test x0_gap !== nothing
#
#    cex = PCQ.verification_counterexample(
#        verification_optimizer,
#        x0_gap,
#    )
#    @test !isempty(cex.modes)
#    @test all(m -> m in 1:2, cex.modes)
#    # Either the environment loops forever outside the verified set, or it drives the run to
#    # where the abstraction's coverage ends. Both are evidence; neither is silence.
#    @test cex.entered_sink || cex.lasso_start >= 1
#    @test x0_gap ∈ bisimulation.states[first(cex.qids)].set
#    @test length(cex.X) == length(cex.modes) + 1
#
#    # No counterexample exists from a verified state — when this small system verifies anything
#    # at all: a two-mode adversary on a tight domain can legitimately leave the verified set
#    # empty, in which case there is no such state to probe.
#    x0_ok = isempty(verified_set) ? nothing : find_point(Set(verified_set))
#    if x0_ok !== nothing
#        @test_throws ErrorException PCQ.verification_counterexample(
#            verification_optimizer,
#            x0_ok,
#        )
#    end
#
#    # And none from a synthesis run.
#    @test_throws ErrorException PCQ.verification_counterexample(
#        cosafe_optimizer,
#        x0_gap,
#    )
#end

@testset "num_slices" begin
    # An empty quotient has no slices.
    empty_slices = Dict{Int, Vector{Int}}()
    T = PCQ.PCBisimulationQuotient{Int, Int}(empty_slices)

    @test PCQ.num_slices(T) == 0

    # For a non-empty quotient, `num_slices` returns the number of slices
    # stored for the first node.
    slices = Dict(
        1 => [1, 2, 3],
        2 => [4, 5],
    )
    T = PCQ.PCBisimulationQuotient{Int, Int}(slices)

    @test PCQ.num_slices(T) == 2
end
@testset "states_by_obs" begin
    # ------------------------------------------------------------------
    # An empty quotient has no states and therefore no observations.
    # ------------------------------------------------------------------
    T = AB.PCLFBisimulationQuotient.PCBisimulationQuotient{Int, Int}(
        Dict{Int, Vector{Int}}(),
    )

    @test AB.PCLFBisimulationQuotient.states_by_obs(T) == Dict{Int, Int}()

    # ------------------------------------------------------------------
    # States are grouped by observation and repeated observations are
    # counted.
    # ------------------------------------------------------------------
    T = AB.PCLFBisimulationQuotient.PCBisimulationQuotient{Int, Int}(
        Dict{Int, Vector{Int}}(),
    )

    T.states[1] = AB.PCLFBisimulationQuotient.PCAbstractState(
        1, 1, 1, 2, 1, Tuple{Int, Int}[],
    )
    T.states[2] = AB.PCLFBisimulationQuotient.PCAbstractState(
        2, 1, 2, 1, 1, Tuple{Int, Int}[],
    )
    T.states[3] = AB.PCLFBisimulationQuotient.PCAbstractState(
        3, 2, 3, 2, 1, Tuple{Int, Int}[],
    )
    T.states[4] = AB.PCLFBisimulationQuotient.PCAbstractState(
        4, 2, 4, 3, 1, Tuple{Int, Int}[],
    )

    result = AB.PCLFBisimulationQuotient.states_by_obs(T)

    @test result == Dict(
        1 => 1,
        2 => 2,
        3 => 1,
    )
end
@testset "states_by_slice" begin
    # ------------------------------------------------------------------
    # An empty quotient has no states and therefore no slices.
    # ------------------------------------------------------------------
    T = AB.PCLFBisimulationQuotient.PCBisimulationQuotient{Int, Int}(
        Dict{Int, Vector{Int}}(),
    )

    @test AB.PCLFBisimulationQuotient.states_by_slice(T) == Dict{Int, Int}()

    # ------------------------------------------------------------------
    # States are grouped by slice and repeated slices are counted.
    # ------------------------------------------------------------------
    T = AB.PCLFBisimulationQuotient.PCBisimulationQuotient{Int, Int}(
        Dict{Int, Vector{Int}}(),
    )

    T.states[1] = AB.PCLFBisimulationQuotient.PCAbstractState(
        1, 1, 1, 1, 1, Tuple{Int, Int}[],
    )
    T.states[2] = AB.PCLFBisimulationQuotient.PCAbstractState(
        2, 1, 2, 2, 1, Tuple{Int, Int}[],
    )
    T.states[3] = AB.PCLFBisimulationQuotient.PCAbstractState(
        3, 2, 3, 1, 2, Tuple{Int, Int}[],
    )
    T.states[4] = AB.PCLFBisimulationQuotient.PCAbstractState(
        4, 2, 4, 3, 2, Tuple{Int, Int}[],
    )

    result = AB.PCLFBisimulationQuotient.states_by_slice(T)

    @test result == Dict(
        1 => 2,
        2 => 2,
    )
end
@testset "states_by_node" begin
    # ------------------------------------------------------------------
    # States are grouped by node and repeated nodes are counted.
    # ------------------------------------------------------------------
    T = AB.PCLFBisimulationQuotient.PCBisimulationQuotient{Int, Int}(
        Dict{Int, Vector{Int}}(),
    )

    T.states[1] = AB.PCLFBisimulationQuotient.PCAbstractState(
        1, 1, 1, 1, 1, Tuple{Int, Int}[],
    )
    T.states[2] = AB.PCLFBisimulationQuotient.PCAbstractState(
        2, 1, 2, 2, 1, Tuple{Int, Int}[],
    )
    T.states[3] = AB.PCLFBisimulationQuotient.PCAbstractState(
        3, 2, 3, 1, 2, Tuple{Int, Int}[],
    )
    T.states[4] = AB.PCLFBisimulationQuotient.PCAbstractState(
        4, 2, 4, 3, 2, Tuple{Int, Int}[],
    )

    # `part_ids` is used to determine the node key type.
    T.part_ids[1] = [1, 2]
    T.part_ids[2] = [3, 4]

    result = AB.PCLFBisimulationQuotient.states_by_node(T)

    @test result == Dict(
        1 => 2,
        2 => 2,
    )
end
@testset "transitions_by_mode" begin
    # ------------------------------------------------------------------
    # An empty quotient has no transitions and therefore no modes.
    # ------------------------------------------------------------------
    T = AB.PCLFBisimulationQuotient.PCBisimulationQuotient{Int, Int}(
        Dict{Int, Vector{Int}}(),
    )

    @test AB.PCLFBisimulationQuotient.transitions_by_mode(T) == Dict{Int, Int}()

    # ------------------------------------------------------------------
    # Transitions are grouped by mode and repeated modes are counted.
    # ------------------------------------------------------------------
    T = AB.PCLFBisimulationQuotient.PCBisimulationQuotient{Int, Int}(
        Dict{Int, Vector{Int}}(),
    )

    T.states[1] = AB.PCLFBisimulationQuotient.PCAbstractState(
        1, 1, 1, 1, 1, [(1, 2), (2, 3)],
    )
    T.states[2] = AB.PCLFBisimulationQuotient.PCAbstractState(
        2, 1, 2, 2, 1, [(1, 3), (2, 4)],
    )
    T.states[3] = AB.PCLFBisimulationQuotient.PCAbstractState(
        3, 2, 3, 1, 2, [(1, 4)],
    )

    result = AB.PCLFBisimulationQuotient.transitions_by_mode(T)

    @test result == Dict(
        1 => 3,
        2 => 2,
    )
end
@testset "outgoing_degree_stats" begin
    # ------------------------------------------------------------------
    # An empty quotient has zero outgoing degree statistics.
    # ------------------------------------------------------------------
    T = AB.PCLFBisimulationQuotient.PCBisimulationQuotient{Int, Int}(
        Dict{Int, Vector{Int}}(),
    )

    @test AB.PCLFBisimulationQuotient.outgoing_degree_stats(T) == Dict(
        :min => 0,
        :max => 0,
        :mean => 0.0,
        :median => 0.0,
    )

    # ------------------------------------------------------------------
    # Outgoing degree statistics are computed from the transitions of
    # each state.
    # ------------------------------------------------------------------
    T = AB.PCLFBisimulationQuotient.PCBisimulationQuotient{Int, Int}(
        Dict{Int, Vector{Int}}(),
    )

    T.states[1] = AB.PCLFBisimulationQuotient.PCAbstractState(
        1, 1, 1, 1, 1, Tuple{Int, Int}[],
    )
    T.states[2] = AB.PCLFBisimulationQuotient.PCAbstractState(
        2, 1, 2, 1, 1, [(1, 3)],
    )
    T.states[3] = AB.PCLFBisimulationQuotient.PCAbstractState(
        3, 2, 3, 1, 1, [(1, 4), (2, 4)],
    )
    T.states[4] = AB.PCLFBisimulationQuotient.PCAbstractState(
        4, 2, 4, 1, 1, [(1, 1), (2, 2), (3, 3)],
    )

    result = AB.PCLFBisimulationQuotient.outgoing_degree_stats(T)

    # Degrees are [0, 1, 2, 3].
    @test result[:min] == 0
    @test result[:max] == 3
    @test result[:mean] == 1.5
    @test result[:median] == 1
end
@testset "state_ids_in_node" begin
    # ------------------------------------------------------------------
    # Return the state IDs belonging to a given node.
    # ------------------------------------------------------------------
    T = AB.PCLFBisimulationQuotient.PCBisimulationQuotient{Int, Int}(
        Dict{Int, Vector{Int}}(),
    )

    T.states[1] = AB.PCLFBisimulationQuotient.PCAbstractState(
        1, 1, 1, 1, 1, Tuple{Int, Int}[],
    )
    T.states[2] = AB.PCLFBisimulationQuotient.PCAbstractState(
        2, 1, 2, 1, 1, Tuple{Int, Int}[],
    )
    T.states[3] = AB.PCLFBisimulationQuotient.PCAbstractState(
        3, 2, 3, 1, 1, Tuple{Int, Int}[],
    )
    T.states[4] = AB.PCLFBisimulationQuotient.PCAbstractState(
        4, 2, 4, 1, 1, Tuple{Int, Int}[],
    )

    @test Set(AB.PCLFBisimulationQuotient.state_ids_in_node(T, 1)) == Set([1, 2])
    @test Set(AB.PCLFBisimulationQuotient.state_ids_in_node(T, 2)) == Set([3, 4])
    @test AB.PCLFBisimulationQuotient.state_ids_in_node(T, 3) == Int[]

    # ------------------------------------------------------------------
    # An explicit collection of state IDs can be supplied.
    # ------------------------------------------------------------------
    @test Set(
        AB.PCLFBisimulationQuotient.state_ids_in_node(
            T,
            1;
            state_ids = [1, 3, 4],
        ),
    ) == Set([1])

    # ------------------------------------------------------------------
    # State IDs that no longer exist are ignored.
    # ------------------------------------------------------------------
    @test Set(
        AB.PCLFBisimulationQuotient.state_ids_in_node(
            T,
            1;
            state_ids = [1, 99],
        ),
    ) == Set([1])
end
@testset "deadend_states" begin
    # ------------------------------------------------------------------
    # An empty quotient has no deadend states.
    # ------------------------------------------------------------------
    T = AB.PCLFBisimulationQuotient.PCBisimulationQuotient{Int, Int}(
        Dict{Int, Vector{Int}}(),
    )

    @test AB.PCLFBisimulationQuotient.deadend_states(T) == Int[]

    # ------------------------------------------------------------------
    # States with no outgoing transitions are identified as deadends.
    # ------------------------------------------------------------------
    T.states[1] = AB.PCLFBisimulationQuotient.PCAbstractState(
        1, 1, 1, 1, 1, Tuple{Int, Int}[],
    )
    T.states[2] = AB.PCLFBisimulationQuotient.PCAbstractState(
        2, 1, 2, 1, 1, [(1, 3)],
    )
    T.states[3] = AB.PCLFBisimulationQuotient.PCAbstractState(
        3, 2, 3, 1, 1, Tuple{Int, Int}[],
    )
    T.states[4] = AB.PCLFBisimulationQuotient.PCAbstractState(
        4, 2, 4, 1, 1, [(1, 1), (2, 2)],
    )

    result = AB.PCLFBisimulationQuotient.deadend_states(T)

    @test Set(result) == Set([1, 3])
end
@testset "self_loop_count" begin
    # ------------------------------------------------------------------
    # An empty quotient has no self-loops.
    # ------------------------------------------------------------------
    T = AB.PCLFBisimulationQuotient.PCBisimulationQuotient{Int, Int}(
        Dict{Int, Vector{Int}}(),
    )

    @test AB.PCLFBisimulationQuotient.self_loop_count(T) == 0

    # ------------------------------------------------------------------
    # Self-loops are counted when a transition targets its own state.
    # ------------------------------------------------------------------
    T.states[1] = AB.PCLFBisimulationQuotient.PCAbstractState(
        1, 1, 1, 1, 1, [(1, 1), (2, 2)],
    )
    T.states[2] = AB.PCLFBisimulationQuotient.PCAbstractState(
        2, 1, 2, 1, 1, [(1, 3)],
    )
    T.states[3] = AB.PCLFBisimulationQuotient.PCAbstractState(
        3, 2, 3, 1, 1, [(1, 3), (2, 1)],
    )

    @test AB.PCLFBisimulationQuotient.self_loop_count(T) == 2

    # ------------------------------------------------------------------
    # Adding a state with two self-loops increases the count by two.
    # ------------------------------------------------------------------
    T.states[4] = AB.PCLFBisimulationQuotient.PCAbstractState(
        4, 2, 4, 1, 1, [(1, 4), (2, 4), (3, 2)],
    )

    @test AB.PCLFBisimulationQuotient.self_loop_count(T) == 4
end
@testset "num_parts" begin
    # ------------------------------------------------------------------
    # A semilinear set with one part has one part.
    # ------------------------------------------------------------------
    S = UT.semilinear_set([
        LazySets.Hyperrectangle(
            low = [0.0, 0.0],
            high = [1.0, 1.0],
        ),
    ])

    @test AB.PCLFBisimulationQuotient.num_parts(S) == 1

    # ------------------------------------------------------------------
    # The number of parts is the number of sets in the semilinear set.
    # ------------------------------------------------------------------
    S = UT.semilinear_set([
        LazySets.Hyperrectangle(
            low = [0.0, 0.0],
            high = [1.0, 1.0],
        ),
        LazySets.Hyperrectangle(
            low = [2.0, 2.0],
            high = [3.0, 3.0],
        ),
        LazySets.Hyperrectangle(
            low = [4.0, 4.0],
            high = [5.0, 5.0],
        ),
    ])

    @test AB.PCLFBisimulationQuotient.num_parts(S) == 3
end
@testset "num_faces" begin
    # ------------------------------------------------------------------
    # A semilinear set with one rectangular part has four faces.
    # ------------------------------------------------------------------
    S = UT.semilinear_set([
        LazySets.Hyperrectangle(
            low = [0.0, 0.0],
            high = [1.0, 1.0],
        ),
    ])

    @test AB.PCLFBisimulationQuotient.num_faces(S) == 4

    # ------------------------------------------------------------------
    # The total number of faces is the sum over all parts.
    # ------------------------------------------------------------------
    S = UT.semilinear_set([
        LazySets.Hyperrectangle(
            low = [0.0, 0.0],
            high = [1.0, 1.0],
        ),
        LazySets.Hyperrectangle(
            low = [2.0, 2.0],
            high = [3.0, 3.0],
        ),
    ])

    @test AB.PCLFBisimulationQuotient.num_faces(S) == 8
end
@testset "cell_complexities" begin
    # ------------------------------------------------------------------
    # An empty quotient has no cell complexities.
    # ------------------------------------------------------------------
    S = UT.semilinear_set([
        LazySets.Hyperrectangle(
            low = [0.0, 0.0],
            high = [1.0, 1.0],
        ),
    ])

    T = AB.PCLFBisimulationQuotient.PCBisimulationQuotient{typeof(S), Int}(
        Dict{Int, Vector{typeof(S)}}(),
    )

    @test AB.PCLFBisimulationQuotient.cell_complexities(T) == (Int[], Int[])

    # ------------------------------------------------------------------
    # Cell complexities contain the number of parts and faces for each
    # state.
    # ------------------------------------------------------------------
    S1 = UT.semilinear_set([
        LazySets.Hyperrectangle(
            low = [0.0, 0.0],
            high = [1.0, 1.0],
        ),
    ])

    S2 = UT.semilinear_set([
        LazySets.Hyperrectangle(
            low = [0.0, 0.0],
            high = [1.0, 1.0],
        ),
        LazySets.Hyperrectangle(
            low = [2.0, 2.0],
            high = [3.0, 3.0],
        ),
    ])

    T = AB.PCLFBisimulationQuotient.PCBisimulationQuotient{typeof(S1), Int}(
        Dict{Int, Vector{typeof(S1)}}(),
    )

    T.states[1] = AB.PCLFBisimulationQuotient.PCAbstractState(
        1, 1, S1, 1, 1, Tuple{Int, Int}[],
    )
    T.states[2] = AB.PCLFBisimulationQuotient.PCAbstractState(
        2, 2, S2, 1, 1, Tuple{Int, Int}[],
    )

    n_parts, n_faces = AB.PCLFBisimulationQuotient.cell_complexities(T)

    @test Set(n_parts) == Set([1, 2])
    @test Set(n_faces) == Set([4, 8])
end
@testset "bisimulation_stats" begin
    # ------------------------------------------------------------------
    # A populated quotient combines all individual statistics.
    # ------------------------------------------------------------------
    T = AB.PCLFBisimulationQuotient.PCBisimulationQuotient{Int, Int}(
        Dict(
            1 => [1, 2],
            2 => [3, 4],
        ),
    )

    T.states[1] = AB.PCLFBisimulationQuotient.PCAbstractState(
        1, 1, 1, 1, 1, [(1, 1), (2, 2)],
    )
    T.states[2] = AB.PCLFBisimulationQuotient.PCAbstractState(
        2, 1, 2, 2, 1, [(1, 3)],
    )
    T.states[3] = AB.PCLFBisimulationQuotient.PCAbstractState(
        3, 2, 3, 1, 2, Tuple{Int, Int}[],
    )
    T.states[4] = AB.PCLFBisimulationQuotient.PCAbstractState(
        4, 2, 4, 3, 2, [(2, 4)],
    )

    T.part_ids[1] = [1, 2]
    T.part_ids[2] = [3, 4]

    stats = AB.PCLFBisimulationQuotient.bisimulation_stats(T)

    @test stats[:num_nodes] == 2
    @test stats[:num_slices] == 2
    @test stats[:num_states] == 4
    @test stats[:num_transitions] == 4

    @test stats[:states_by_obs] == Dict(
        1 => 2,
        2 => 1,
        3 => 1,
    )

    @test stats[:states_by_slice] == Dict(
        1 => 2,
        2 => 2,
    )

    @test stats[:states_by_node] == Dict(
        1 => 2,
        2 => 2,
    )

    @test stats[:transitions_by_mode] == Dict(
        1 => 2,
        2 => 2,
    )

    @test stats[:outgoing_degree_stats] == Dict(
        :min => 0,
        :max => 2,
        :mean => 1.0,
        :median => 1,
    )

    @test Set(stats[:deadend_states]) == Set([3])
    @test stats[:num_deadend_states] == 1
    @test stats[:self_loop_count] == 2
end
@testset "print_bisimulation_stats" begin
    # ------------------------------------------------------------------
    # Create a populated quotient with known statistics.
    # ------------------------------------------------------------------
    T = AB.PCLFBisimulationQuotient.PCBisimulationQuotient{Int, Int}(
        Dict(
            1 => [1, 2],
            2 => [3, 4],
        ),
    )

    T.states[1] = AB.PCLFBisimulationQuotient.PCAbstractState(
        1, 1, 1, 1, 1, [(1, 1), (2, 2)],
    )
    T.states[2] = AB.PCLFBisimulationQuotient.PCAbstractState(
        2, 1, 2, 2, 1, [(1, 3)],
    )
    T.states[3] = AB.PCLFBisimulationQuotient.PCAbstractState(
        3, 2, 3, 1, 2, Tuple{Int, Int}[],
    )
    T.states[4] = AB.PCLFBisimulationQuotient.PCAbstractState(
        4, 2, 4, 3, 2, [(2, 4)],
    )

    T.part_ids[1] = [1, 2]
    T.part_ids[2] = [3, 4]

    # ------------------------------------------------------------------
    # The function prints the statistics and returns nothing.
    # ------------------------------------------------------------------
    @test AB.PCLFBisimulationQuotient.print_bisimulation_stats(T) === nothing
end
end # module TestMain
