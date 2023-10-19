module TestMain
include("solvers.jl")
include("../../problems/gol_lazar_belta.jl")

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const OP = DI.Optim

using LinearAlgebra, Test
import CDDLib
using Polyhedra
using HybridSystems

_name(o::MOI.OptimizerWithAttributes) = split(string(o.optimizer_constructor), ".")[2]

function _prob(N, q_0, x_0::Vector{T}, zero_cost::Bool) where {T}
    return GolLazarBelta.problem(CDDLib.Library(), T; N, q_0, x_0, zero_cost)
end
function _test(
    algo,
    N,
    q0,
    x0,
    x_expected,
    u_expected,
    obj_expected,
    zero_cost::Bool,
    mi::Bool;
    kws...,
)
    problem = _prob(N, q0, x0, zero_cost)
    optimizer = MOI.instantiate(algo)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("problem"), problem)
    @info("Solving... depth: $N")
    @time MOI.optimize!(optimizer)
    @info("Solved.")
    if x_expected === nothing
        @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.INFEASIBLE
    else
        @test MOI.get(optimizer, MOI.TerminationStatus()) in
              [MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED]
        xu = MOI.get(optimizer, ST.ContinuousTrajectoryAttribute())
        @test typeof(x_expected) == typeof(xu.x)
        @test typeof(u_expected) == typeof(xu.u)
        if isempty(x_expected)
            @test isempty(u_expected)
            @test isempty(xu.x)
            @test isempty(xu.u)
        else
            @test xu.x ≈ x_expected atol = 1e-2
            @test xu.u ≈ u_expected atol = 1e-2
        end
        @test MOI.get(optimizer, MOI.ObjectiveValue()) ≈ obj_expected atol = 1e-2
    end
    if optimizer isa OP.BranchAndBound.Optimizer
        @test optimizer.num_done + optimizer.num_pruned_bound + optimizer.num_pruned_inf ==
              optimizer.num_total
        if x_expected !== nothing
            return optimizer.Q_function
        end
    end
    return
end
function learn_test(qp_solver, x0 = [-1.645833614657878, 1.7916672467705592])
    prob = _prob(1, 15, x0, false)
    t(i, j) = first(transitions(prob.system, i, j))
    dtraj = ST.DiscreteTrajectory(15, [t(15, 20)])
    ctraj = ST.ContinuousTrajectory([[-0.5, 0.5]], [[-1.2916666674915085]])
    algo = OP.HybridDualDynamicProgrammingAlgo(qp_solver, CDDLib.Library(), 1e-5, 1e-4, 1)
    Q_function = OP.instantiate(prob, algo)
    OP.learn(Q_function, prob, dtraj, ctraj, algo)
    @test isempty(Q_function.cuts[0, 15])
    @test !hasallhalfspaces(Q_function.domains[0, 15])
    @test length(Q_function.cuts[1, 15]) == 1
    @test first(Q_function.cuts[1, 15]) ≈
          UT.AffineFunction([0.0, 2.583334480953581], -2.960071516682004) rtol = 1e-6
    @test !hashyperplanes(Q_function.domains[1, 15]) == 1
    @test nhalfspaces(Q_function.domains[1, 15]) == 1
    a = normalize(-[2, 1])
    @test first(halfspaces(Q_function.domains[1, 15])) ≈ HalfSpace(a, a ⋅ x0)
    @test isempty(Q_function.cuts[0, 20])
    @test !hasallhalfspaces(Q_function.domains[0, 20])
    @test isempty(Q_function.cuts[1, 20])
    @test !hasallhalfspaces(Q_function.domains[1, 20])
end
const x_expected_9 = [
    [-0.925, -2.35],
    [-3.0625, -1.925],
    [-4.6375, -1.225],
    [-5.375, -0.25],
    [-5.0, 1.0],
    [-3.8375, 1.325],
    [-2.5, 1.35],
    [-1.2875, 1.075],
    [-0.5, 0.5],
]
const u_expected_9 =
    [[0.15], [0.425], [0.7], [0.975], [1.25], [0.325], [0.025], [-0.275], [-0.575]]
const x_expected_11 = [
    [-4.02204, -4.04409],
    [-7.23015, -2.37212],
    [-8.90827, -0.984118],
    [-9.34036, 0.119934],
    [-8.81038, 0.940033],
    [-7.60227, 1.47618],
    [-6.0, 1.72837],
    [-4.26869, 1.73426],
    [-2.63582, 1.53149],
    [-1.31004, 1.12007],
    [-0.5, 0.5],
]
const u_expected_11 = [
    [1.9559145673603502],
    [1.6719605695509308],
    [1.3880065717415113],
    [1.1040525739320919],
    [0.8200985761226725],
    [0.5361445783132529],
    [0.2521905805038335],
    [0.0058871851040525],
    [-0.2027656078860899],
    [-0.4114184008762322],
    [-0.6200711938663745],
]
function _test9(algo, T::Type; kws...)
    return _test(
        algo,
        9,
        8,
        T[1.5, -2.5],
        x_expected_9,
        u_expected_9,
        11.71875,
        false,
        true;
        kws...,
    )
end
function _test11(algo, T::Type; kws...)
    return _test(
        algo,
        11,
        3,
        T[1.0, -6.0],
        x_expected_11,
        u_expected_11,
        21.385062979189478,
        false,
        true;
        kws...,
    )
end
@testset "Gol-Lazar-Belta" begin
    function tests(qp_solver, miqp_solver)
        # Pavito does not support indicator constraints yet so we use `false` here
        @testset "$(_name(algo))" for algo in [
            optimizer_with_attributes(
                OP.BemporadMorari.Optimizer{Float64},
                "continuous_solver" => qp_solver,
                "mixed_integer_solver" => miqp_solver,
                "indicator" => false,
                "log_level" => 0,
            ),
        ]
            @testset "Depth: 0" begin
                _test(algo, 0, 18, [0.0, 1.0], nothing, nothing, nothing, true, false)
                _test(
                    algo,
                    0,
                    20,
                    [0.5, 0.0],
                    Vector{Float64}[],
                    Vector{Float64}[],
                    0.0,
                    true,
                    false,
                )
            end
            @testset "Depth: 1" begin
                _test(algo, 1, 18, [0.0, 1.0], [[0.5, 0.0]], [[-1.0]], 1.0, true, false)
                _test(algo, 1, 7, [0.0, 1.0], nothing, nothing, nothing, true, false)
            end
            @testset "Depth: 2" begin
                _test(
                    algo,
                    2,
                    18,
                    [0.0, 1.0],
                    [[0.55, 0.1], [0.5, -0.2]],
                    [[-0.9], [-0.3]],
                    0.9,
                    true,
                    true,
                )
                _test(algo, 2, 7, [0.0, 1.0], nothing, nothing, nothing, true, true)
            end
            @testset "Depth: 9" begin
                _test9(algo, Float64)
            end
            @testset "Depth: 11" begin
                _test11(algo, Float64)
            end
        end
    end
    function test_Q_reuse(qp_solver, miqp_solver)
        # Pavito does not support indicator constraints yet so we use `false` here

        algo(max_iter, Q_function_init) = optimizer_with_attributes(
            OP.BranchAndBound.Optimizer{Float64},
            "continuous_solver" => qp_solver,
            "mixed_integer_solver" => miqp_solver,
            "max_iter" => max_iter,
            "Q_function_init" => Q_function_init,
        )
        qalgo(max_iter) = optimizer_with_attributes(
            OP.BranchAndBound.Optimizer{Float64},
            "continuous_solver" => qp_solver,
            "mixed_integer_solver" => miqp_solver,
            "max_iter" => max_iter,
            "max_time" => 300.0,
            "lower_bound" => OP.HybridDualDynamicProgrammingAlgo(
                qp_solver,
                CDDLib.Library(),
                1e-5,
                1e-4,
                1,
            ),
        )
        # Gurobi | OSQP
        Q9 = _test9(qalgo(789), Float64)   #    746 | 789
        @show sum(length.(Q9.cuts))
        _test9(algo(821, Q9), Float64)     #    821 | 782
        _test11(algo(74, Q9), Float64)     #     74 | 74
        Q11 = _test11(qalgo(75), Float64)  #     75 | 75
        @show sum(length.(Q11.cuts))
        _test9(algo(800, Q11), Float64)    #    785 | 782
        _test11(algo(74, Q11), Float64)    #     74 | 74
        Q = OP.q_merge(Q9, Q11)
        _test9(algo(818, Q), Float64)      #    818 | 775
        return _test11(algo(74, Q), Float64)      #     74 | 74
    end
    tests(qp_solver, miqp_solver)
    @testset "Q_reuse" begin
        test_Q_reuse(qp_solver, miqp_solver)
    end
    @testset "Learn test" begin
        learn_test(qp_solver)
    end
end
end
