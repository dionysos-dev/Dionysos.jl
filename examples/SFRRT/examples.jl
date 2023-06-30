
using StaticArrays, LinearAlgebra, Random, IntervalArithmetic
using MathematicalSystems, HybridSystems
using JuMP, Mosek, MosekTools
using Plots, Colors
using Test
Random.seed!(0)

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const CO = DI.Control
const PR = DI.Problem
const OP = DI.Optim
const AB = OP.Abstraction

function example_box_ellipsoid()
    c = [-10.0; -10.0]
    P = [2.0 6.0; 6.0 20.0]
    E = UT.Ellipsoid(P, c)
    box = UT.get_min_bounding_box(E)
    fig = plot(; aspect_ratio = :equal)
    plot!(fig, box)
    plot!(fig, E)
    display(fig)
    return
end

include("../../problems/non_linear.jl")

function test_backward_transition()
    E2 = UT.Ellipsoid([2.0 0.2; 0.2 0.5], [3.0; 3.0])
    Ubound = 5.0
    # U = UT.HyperRectangle(SVector(-Ubound, -Ubound), SVector(Ubound, Ubound))
    # U = UT.Ellipsoid([1/25.0 0.0; 0.0 1/25.0], [0.0; 0.0])
    # U = UT.IntersectionSet([UT.Ellipsoid([1/25.0 0.0; 0.0 1/25.0], [0.0; 0.0]), UT.Ellipsoid([1/20.0 0.0; 0.0 1/30.0], [0.0; 0.0])])
    U = UT.IntersectionSet([
        UT.Ellipsoid([1/25.0 0.0; 0.0 1/25.0], [0.0; 0.0]),
        UT.Ellipsoid([1/20.0 0.0; 0.0 1/30.0], [0.0; 0.0]),
        UT.HyperRectangle(SVector(-4.0, -5.0), SVector(4.0, 5.0)),
    ])
    problem = NonLinear.problem(; U = U)
    sys = problem.system
    # Construct the linear approximation
    xnew = SVector{2, Float64}([0.0; 0.0])
    unew = [0.0; 0.0]
    nw = sys.nw
    wnew = zeros(nw)
    X̄ = IntervalBox(xnew .+ sys.ΔX)
    Ū = IntervalBox(unew .+ sys.ΔU)
    W̄ = IntervalBox(wnew .+ sys.ΔW)
    (affineSys, L) = ST.buildAffineApproximation(
        sys.fsymbolic,
        sys.x,
        sys.u,
        sys.w,
        xnew,
        unew,
        wnew,
        X̄,
        Ū,
        W̄,
    )
    # Solve the control problem
    S = UT.get_full_psd_matrix(problem.transition_cost)
    sdp_opt = optimizer_with_attributes(Mosek.Optimizer, MOI.Silent() => true)
    maxδx = 100
    maxδu = 10 * 2
    λ = 0.01
    E1, cont, cost = SY.transition_backward(
        affineSys,
        E2,
        xnew,
        unew,
        sys.Uformat,
        S,
        L,
        sdp_opt;
        λ = λ,
        maxδx = maxδx,
        maxδu = maxδu,
    )

    println("Max cost : ", cost)
    println(
        "Controller feasible : ",
        ST.check_feasibility(
            E1,
            E2,
            sys.f_eval,
            cont.c_eval,
            nw,
            sys.U;
            N = 500,
            input_check = true,
        ),
    )

    fig1 = plot(; aspect_ratio = :equal)
    ST.plot_check_feasibility!(E1, E2, sys.f_eval, cont.c_eval, nw; dims = [1, 2], N = 500)
    display(fig1)

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
    plot!(E2; color = :red)
    display(fig2)

    U_used = UT.affine_transformation(E1, cont.K, cont.ℓ - cont.K * cont.c)
    fig3 = plot(; aspect_ratio = :equal)
    plot!(
        fig3,
        sys.U;
        color = :green,
        label = "Feasible input set",
        fillalpha = 0.4,
        linealpha = 1.0,
        linewidth = 2,
    )
    plot!(fig3, U_used; color = :red, label = "Input set used")
    return display(fig3)
end

test_backward_transition()
