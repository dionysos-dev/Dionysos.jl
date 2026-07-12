module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

using LinearAlgebra
using JuMP
import MathOptInterface as MOI
import Clarabel
import HybridSystems
import Spot
using LazySets

const PCLF = UT.PathCompleteFramework

# Contract test for the PCLF bisimulation-quotient optimizer. Small discrete switched
# linear system with a polyhedral path-complete Lyapunov function; adapted from the
# BisimulationQuotient case-study script.
@testset "PCLF bisimulation quotient" begin
    A1 = (1.0 / 10.0) * [1.5519 0.4474; 7.6412 7.4716]
    A2 = (1.0 / 10.0) * [0.4750 9.1755; 1.8955 0.1850]
    f = HybridSystems.discreteswitchedsystem([A1, A2])

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
end

end # module TestMain
