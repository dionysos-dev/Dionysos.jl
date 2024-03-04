using StaticArrays, LinearAlgebra, IntervalArithmetic, Random
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
const FALLBACK_URL = "mosek://solve.mosek.com:30080"

function test_backward_transition(Wbound, E2, xnew, U, λ, ρ)
    W = UT.HyperRectangle(SVector(-Wbound, -Wbound), SVector(Wbound, Wbound))
    problem = NonLinear.problem(; U = U, W = W, noise = true, μ = ρ)
    sys = problem.system
    # Construct the linear approximation
    unew = zeros(sys.nu)
    wnew = zeros(sys.nw)
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
    MOI.set(sdp_opt, MOI.RawOptimizerAttribute("fallback"), FALLBACK_URL)
    maxδx = 100.0
    maxδu = 100.0

    E1, cont, cost = SY.transition_backward(
        affineSys,
        E2,
        xnew,
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

    # Get results
    cost_eval(x, u) = UT.function_value(problem.transition_cost, x, u)
    ETilde = UT.affine_transformation(
        E1,
        affineSys.A + affineSys.B * cont.K,
        affineSys.B * (cont.ℓ - cont.K * cont.c) + affineSys.c,
    )
    U_used = UT.affine_transformation(E1, cont.K, cont.ℓ - cont.K * cont.c)
    # Display results
    println()
    println("Max cost : ", cost)
    println("Volume of initial ellipsoid : ", UT.get_volume(E1))
    println("Input set volume : ", UT.get_volume(U_used))
    println(
        "Controller feasible : ",
        ST.check_feasibility(
            E1,
            E2,
            sys.f_eval,
            cont.c_eval,
            sys.U,
            sys.W;
            N = 500,
            input_check = true,
            noise_check = true,
        ),
    )

    # Display the initial set, target set and the image of the initial ellipsoid under the linear model approximation
    fig1 = plot(; aspect_ratio = :equal)
    plot!(fig1, E1; color = :green)
    plot!(fig1, E2; color = :red)

    ST.plot_transitions!(E1, sys.f_eval, cont.c_eval, sys.W; N = 300)
    plot!(fig1, ETilde; color = :blue)
    display(fig1)

    # Display the data-driven test of the controller
    fig2 = plot(; aspect_ratio = :equal)
    ST.plot_check_feasibility!(
        E1,
        E2,
        sys.f_eval,
        cont.c_eval,
        sys.W;
        dims = [1, 2],
        N = 500,
    )
    display(fig2)

    # Display the cost of the controller
    fig3 = plot(;
        aspect_ratio = :equal,
        xtickfontsize = 12,
        ytickfontsize = 12,
        guidefontsize = 18,
        titlefontsize = 16,
        legend = false,
    )
    xlims!(-1, 6)
    ylims!(-1, 6)
    xticks!(-1:1:6)
    yticks!(-1:1:6)
    xlabel!("\$x_1\$")
    ylabel!("\$x_2\$")
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
    col1 = parse(ColorTypes.RGB, "#ff756f")
    col2 = parse(ColorTypes.RGB, "#64bc60")
    plot!(E2; color = col1)
    plot!(ETilde; color = col2)
    display(fig3)

    #Display the feasible input set and the input set effectively used
    col3 = RGB(0.2, 0.4, 0.8)
    fig4 = plot(;
        aspect_ratio = :equal,
        xtickfontsize = 12,
        ytickfontsize = 12,
        guidefontsize = 18,
        titlefontsize = 16,
        legend = false,
    )
    xlims!(-6, 6)
    ylims!(-6, 6)
    xticks!(-6:2:6)
    yticks!(-6:2:6)
    xlabel!("\$u_1\$")
    ylabel!("\$u_2\$")

    plot!(
        fig4,
        sys.U;
        color = col3,
        label = "",
        fillalpha = 0.4,
        linealpha = 1.0,
        linewidth = 2,
    )
    plot!(fig4, U_used; color = :green, label = "Input set used")
    display(fig4)

    savefig(fig1, "fig1.png")
    savefig(fig2, "fig2.png")
    savefig(fig3, "fig3.png")
    savefig(fig4, "fig4.png")
    return
end

E2 = UT.Ellipsoid([2.0 0.2; 0.2 0.5], [4.0; 4.0])
U = UT.IntersectionSet([
    UT.Ellipsoid([1/25.0 0.0; 0.0 1/25.0], [0.0; 0.0]),
    UT.Ellipsoid([1/20.0 0.0; 0.0 1/30.0], [0.0; 0.0]),
    UT.HyperRectangle(SVector(-4.0, -5.0), SVector(4.0, 5.0)),
])
xnew = SVector{2, Float64}([1.0; 1.0])
#fig 1
ρ = 0.0005
Wbound = 0.1
λ = 0.01
test_backward_transition(Wbound, E2, xnew, U, λ, ρ)

#fig 2
ρ = 0.0005
Wbound = 0.1
λ = 0.3
test_backward_transition(Wbound, E2, xnew, U, λ, ρ)

# #fig 3
ρ = 0.001
Wbound = 0.15
λ = 0.01
test_backward_transition(Wbound, E2, xnew, U, λ, ρ)
