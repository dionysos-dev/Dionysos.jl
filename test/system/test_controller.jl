module TestMain
using Test

using StaticArrays, LinearAlgebra, IntervalArithmetic, Random
using JuMP
using Clarabel
using Plots

using Dionysos
const DI = Dionysos
const ST = DI.System
const UT = DI.Utils
const SY = DI.Symbolic
Random.seed!(0)

@testset "ConstantController" begin
    c = [0.1; 0.2]
    cont = ST.ConstantController(c)
    @test cont.c === c
    @test (cont.c_eval)(1.0) === c
end

include("../../problems/non_linear.jl")

function trial(E2, c, μ, U, W, λ)
    problem = NonLinear.problem(; U = U, W = W, noise = true, μ = μ)
    sys = problem.system

    # Construct the linear approximation
    unew = zeros(sys.nu)
    wnew = zeros(sys.nw)
    X̄ = IntervalBox(c .+ sys.ΔX)
    Ū = IntervalBox(unew .+ sys.ΔU)
    W̄ = IntervalBox(wnew .+ sys.ΔW)
    (affineSys, L) = ST.buildAffineApproximation(
        sys.fsymbolic,
        sys.x,
        sys.u,
        sys.w,
        c,
        unew,
        wnew,
        X̄,
        Ū,
        W̄,
    )

    # Solve the control problem
    S = UT.get_full_psd_matrix(problem.transition_cost)
    sdp_opt = optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => true)

    maxδx = 100.0
    maxδu = 100.0
    E1, cont, max_cost = SY.transition_backward(
        affineSys,
        E2,
        c,
        unew,
        sys.Uformat,
        sys.Wformat,
        S,
        L,
        sdp_opt;
        λ = λ,
        maxδx = maxδx,
        maxδu = maxδu,
    )

    if cont === nothing
        success = false
        init_set_volume = 0.0
        ETilde = nothing
        U_used = nothing
        input_set_volume = 0.0
    else
        success = ST.check_feasibility(
            E1,
            E2,
            sys.f_eval,
            cont.c_eval,
            sys.U,
            sys.W;
            N = 500,
            input_check = true,
            noise_check = true,
        )
        init_set_volume = UT.get_volume(E1)
        ETilde = UT.affine_transformation(
            E1,
            affineSys.A + affineSys.B * cont.K,
            affineSys.B * (cont.ℓ - cont.K * cont.c) + affineSys.c,
        )
        U_used = UT.affine_transformation(E1, cont.K, cont.ℓ - cont.K * cont.c)
        input_set_volume = UT.get_volume(U_used)
    end

    fig1 = plot(; aspect_ratio = :equal)
    ST.plot_check_feasibility!(
        E1,
        E2,
        sys.f_eval,
        cont.c_eval,
        sys.W;
        dims = [1, 2],
        N = 500,
    )
    @test isa(fig1, Plots.Plot{Plots.GRBackend})

    fig2 = plot(; aspect_ratio = :equal)
    cost_eval(x, u) = UT.function_value(problem.transition_cost, x, u)
    ST.plot_controller_cost!(
        E1,
        cont.c_eval,
        cost_eval;
        N = 3000,
        scale = 0.01,
        dims = [1, 2],
        color = :white,
        linewidth = 7,
    )
    @test isa(fig2, Plots.Plot{Plots.GRBackend})

    return (success, max_cost, init_set_volume, input_set_volume)
end

@testset "Backward_transition" begin
    E2 = UT.Ellipsoid([2.0 0.2; 0.2 0.5], [3.0; 3.0])
    c = SVector{2, Float64}([1.0; 1.0])
    U = UT.IntersectionSet([
        UT.Ellipsoid([1/25.0 0.0; 0.0 1/25.0], [0.0; 0.0]),
        UT.Ellipsoid([1/20.0 0.0; 0.0 1/30.0], [0.0; 0.0]),
    ])
    Wbound = 0.01
    W = UT.HyperRectangle(SVector(-Wbound, -Wbound), SVector(Wbound, Wbound))
    μ = 0.0008
    λ = 0.01

    (success, max_cost, init_set_volume, input_set_volume) = trial(E2, c, μ, U, W, λ)
    @test success === true
    @test max_cost ≈ 33.98 atol = 1
    @test init_set_volume ≈ 22.17 atol = 1
    @test input_set_volume ≈ 22.21 atol = 1
end

end
