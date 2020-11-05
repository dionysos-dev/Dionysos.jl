include("../examples/gol_lazar_belta.jl")
include("solvers.jl")
using Test
import CDDLib
using Dionysos

@testset "Gol-Lazar-Belta" begin
    function _test(algo, N, q0, x0, x_expected, u_expected, obj_expected, zero_cost::Bool, mi::Bool)
        system = gol_lazar_belta(CDDLib.Library())
        if zero_cost
            state_cost = Fill(ZeroCost(), nmodes(system))
        else
            state_cost = [mode == system.ext[:q_T] ? ConstantCost(0.0) : ConstantCost(1.0)
                          for mode in modes(system)]
        end
        problem = OptimalControlProblem(
            system,
            q0, x0,
            Fill(state_cost, N),
            Fill(Fill(QuadraticControlCost(ones(1, 1)), ntransitions(system)), N),
            system.ext[:q_T],
            N
        )
        @info("Solving... depth: $N")
        @time xu = optimal_control(problem, algo)
        @info("Solved.")
        if x_expected === nothing
            @test xu === nothing
        else
            @test xu !== nothing
            @test typeof(x_expected) == typeof(xu[1].x)
            @test typeof(u_expected) == typeof(xu[1].u)
            if isempty(x_expected)
                @test isempty(u_expected)
                @test isempty(xu[1].x)
                @test isempty(xu[1].u)
            else
                @test xu[1].x ≈ x_expected atol=1e-2
                @test xu[1].u ≈ u_expected atol=1e-2
            end
            @test xu[2] ≈ obj_expected atol=1e-2
        end
    end
    function tests(qp_solver, miqp_solver)
        # Pavito does not support indicator constraints yet so we use `false` here
        @testset "$(typeof(algo))" for algo in [
            BemporadMorari(qp_solver, miqp_solver, false, 0),
            BranchAndBound(qp_solver, miqp_solver, max_iter = 4000)
        ]
            @testset "Depth: 0" begin
            _test(algo, 0, 18, [0.0, 1.0], nothing, nothing, nothing, true, false)
            _test(algo, 0, 20, [0.5, 0.0], Vector{Float64}[], Vector{Float64}[], 0.0, true, false)
            end
            @testset "Depth: 1" begin
            _test(algo, 1, 18, [0.0, 1.0], [[0.5, 0.0]], [[-1.0]], 1.0, true, false)
            _test(algo, 1, 7, [0.0, 1.0], nothing, nothing, nothing, true, false)
            end
            @testset "Depth: 2" begin
            _test(algo, 2, 18, [0.0, 1.0], [[0.55, 0.1], [0.5, -0.2]], [[-0.9], [-0.3]], 0.9, true, true)
            _test(algo, 2, 7, [0.0, 1.0], nothing, nothing, nothing, true, true)
            end
            @testset "Depth: 9" begin
            _test(algo, 9, 8, [1.5, -2.5],
                [[-0.925 , -2.35],
                 [-3.0625, -1.925],
                 [-4.6375, -1.225],
                 [-5.375 , -0.25],
                 [-5.0   ,  1.0],
                 [-3.8375,  1.325],
                 [-2.5   ,  1.35],
                 [-1.2875,  1.075],
                 [-0.5   ,  0.5]],
                [[0.15], [0.425], [0.7], [0.975], [1.25], [0.325], [0.025], [-0.275], [-0.575]],
                11.71875, false, true)
            end
            @testset "Depth: 11" begin
            _test(algo, 11, 3, [1.0, -6.0],
                [
     [-4.02204,  -4.04409],
     [-7.23015,  -2.37212],
     [-8.90827,  -0.984118],
     [-9.34036,   0.119934],
     [-8.81038,   0.940033],
     [-7.60227,   1.47618],
     [-6.0    ,   1.72837],
     [-4.26869,   1.73426],
     [-2.63582,   1.53149],
     [-1.31004,   1.12007],
     [-0.5    ,   0.5]],
                [
     [ 1.9559145673603502],
     [ 1.6719605695509308],
     [ 1.3880065717415113],
     [ 1.1040525739320919],
     [ 0.8200985761226725],
     [ 0.5361445783132529],
     [ 0.2521905805038335],
     [ 0.0058871851040525],
     [-0.2027656078860899],
     [-0.4114184008762322],
     [-0.6200711938663745]],
                21.385062979189478, false, true)
            end
        end
    end
    tests(qp_solver, miqp_solver)
end
