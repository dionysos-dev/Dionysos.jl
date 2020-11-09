include("../examples/gol_lazar_belta.jl")
include("solvers.jl")
using Test
import CDDLib
using Dionysos

@testset "Gol-Lazar-Belta" begin
    function _prob( N, q0, x0, zero_cost::Bool)
        system = gol_lazar_belta(CDDLib.Library())
        if zero_cost
            state_cost = Fill(ZeroFunction(), nmodes(system))
        else
            state_cost = [mode == system.ext[:q_T] ? ConstantFunction(0.0) : ConstantFunction(1.0)
                          for mode in modes(system)]
        end
        return OptimalControlProblem(
            system,
            q0, x0,
            Fill(state_cost, N),
            Fill(Fill(QuadraticControlFunction(ones(1, 1)), ntransitions(system)), N),
            system.ext[:q_T],
            N
        )
    end
    function _test(algo, N, q0, x0, x_expected, u_expected, obj_expected, zero_cost::Bool, mi::Bool; kws...)
        problem = _prob(N, q0, x0, zero_cost)
        @info("Solving... depth: $N")
        @time xu = optimal_control(problem, algo; kws...)
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
        xu !== nothing && length(xu) >= 3 && return xu[3]
        return
    end
#    function learn_test(qp_solver)
#        prob = _prob(11, 3, [1.0, -6.0], false)
#        t(i, j) = first(transitions(prob.system, i, j))
#        dtraj = Dionysos.DiscreteTrajectory(3, [t(3, 3), t(3, 2), t(2, 2), t(2, 2), t(2, 2), t(2, 2), t(2, 2), t(2, 10), t(10, 15), t(15, 15), t(15, 20)])
#        ctraj = Dionysos.ContinuousTrajectory(
#            [[-4.08326083652693, -4.16652157792225], [-7.452103692575193, -2.571163630780764], [-9.344647629427765, -1.2139233592305483], [-9.999004371331917, -0.09478891152710031], [-9.653269015312578, 0.7862610956032603], [-8.545511305383734, 1.4292559716987285], [-6.913769471419556, 1.8342294340939918], [-4.996048425321604, 2.00121442419686], [-3.0616913911649153, 1.867491079841009], [-1.4429744774058832, 1.3699304105343828], [-0.5037481670764479, 0.5085061781102557]],
#            [[1.833479150334962], [1.5953582598989637], [1.3572401773353209], [1.119133970839774], [0.8810491889863714], [0.642993774334867], [0.40497214759387024], [0.1669835371039382], [-0.13372239101859426], [-0.49755590489234236], [-0.8614157126243801]]
#        )
#        algo = HybridDualDynamicProgrammingAlgo(qp_solver)
#        Q_function = Dionysos.instantiate(prob, algo)
#        model = Model(qp_solver)
#        params = ctraj.x[end - 1]
#        t = dtraj.transitions[end]
#        r = resetmap(prob.system, t)
#        #@variable(model, u[1:inputdim(r)] in inputset(r))
#        @show hrep(stateset(prob.system, 20))
#        u = ctraj.u[end]
#        verts = [[0.5, 0.5, 0.0], [0.5, -0.5, 0.0], [-0.5, -0.5, 0.0], [-0.5, 0.5, 0.0]]
#        @variable(model, λ[eachindex(verts)] ≥ 0)
#        @constraint(model, sum(λ) == 1)
#        @variable(model, θ)
#        epi(i) = sum(λ[v] * verts[v][i] for v in eachindex(verts))
#        x_next = r.A * params + r.B * u
#        for j in eachindex(x_next)
#            @constraint(model, x_next[j] == epi(j))
#        end
#        @constraint(model, θ >= epi(length(x_next) + 1))
#        println(model)
#        optimize!(model)
#        @show termination_status(model)
#        #Dionysos.learn(Q_function, prob, dtraj, ctraj, algo)
#    end
    function _test9(algo; kws...)
        _test(algo, 9, 8, [1.5, -2.5],
              x_expected_9,
              u_expected_9,
              11.71875, false, true; kws...)
    end
    function _test11(algo; kws...)
        _test(algo, 11, 3, [1.0, -6.0],
              x_expected_11,
              u_expected_11,
              21.385062979189478, false, true; kws...)
    end
    x_expected_9 = [
        [-0.925 , -2.35],
        [-3.0625, -1.925],
        [-4.6375, -1.225],
        [-5.375 , -0.25],
        [-5.0   ,  1.0],
        [-3.8375,  1.325],
        [-2.5   ,  1.35],
        [-1.2875,  1.075],
        [-0.5   ,  0.5]]
    u_expected_9 = [[0.15], [0.425], [0.7], [0.975], [1.25], [0.325], [0.025], [-0.275], [-0.575]]
    x_expected_11 = [
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
        [-0.5    ,   0.5]]
    u_expected_11 = [
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
        [-0.6200711938663745]]
    function tests(qp_solver, miqp_solver)
        # Pavito does not support indicator constraints yet so we use `false` here
        @testset "$(split(string(typeof(algo)), "{")[1])" for algo in [
#            BemporadMorari(qp_solver, miqp_solver, false, 0),
#            BranchAndBound(qp_solver, miqp_solver, DiscreteLowerBoundAlgo(qp_solver), max_iter = 1111),
            BranchAndBound(qp_solver, miqp_solver, HybridDualDynamicProgrammingAlgo(qp_solver), max_iter = 871)
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
                _test9(algo)
            end
            @testset "Depth: 11" begin
                Q = _test11(algo)
            end
        end
    end
    function test_Q_reuse(qp_solver, miqp_solver)
        # Pavito does not support indicator constraints yet so we use `false` here
        algo(max_iter) = BranchAndBound(qp_solver, miqp_solver, DiscreteLowerBoundAlgo(qp_solver), max_iter = max_iter)
        qalgo(max_iter) = BranchAndBound(qp_solver, miqp_solver, HybridDualDynamicProgrammingAlgo(qp_solver), max_iter = max_iter)
        Q9 = _test9(qalgo(871))
        @show sum(length.(Q9.cuts))
        _test9(algo(761), Q_function_init = Q9)
        _test11(algo(85), Q_function_init = Q9)
        Q11 = _test11(qalgo(85))
        @show sum(length.(Q11.cuts))
        _test9(algo(880), Q_function_init = Q11)
        _test11(algo(85), Q_function_init = Q11)
        Q = Dionysos.q_merge(Q9, Q11)
        _test9(algo(747), Q_function_init = Q)
        _test11(algo(85), Q_function_init = Q)
    end
    tests(qp_solver, miqp_solver)
    test_Q_reuse(qp_solver, miqp_solver)
    #learn_test(qp_solver)
end
