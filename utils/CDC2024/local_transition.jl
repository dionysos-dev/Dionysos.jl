using StaticArrays, Random
using MathematicalSystems, HybridSystems
using JuMP, Clarabel
import LazySets
using Plots, Colors
using Test
Random.seed!(0)

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System

function example_box_ellipsoid()
    c = [-10.0; -10.0]
    P = [2.0 6.0; 6.0 20.0]
    E = LazySets.Ellipsoid(c, inv(P))
    box = UT.get_min_bounding_box(E)
    fig = plot(; aspect_ratio = :equal)
    plot!(fig, box)
    plot!(fig, E)
    display(fig)
    return
end

include("../../problems/non_linear.jl")

function test_backward_transition(Wbound, E2, xnew, U, λ, ρ)
    W = UT.box(SVector(-Wbound, -Wbound), SVector(Wbound, Wbound))
    problem = NonLinear.problem(; U = U, W = W, noise = true, μ = ρ)
    sys = problem.system
    # Construct the linear approximation
    unew = zeros(sys.nu)
    wnew = zeros(sys.nw)
    X̄ = xnew .+ sys.ΔX
    Ū = unew .+ sys.ΔU
    W̄ = wnew .+ sys.ΔW
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
    sdp_opt = optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => true)
    maxδx = 100.0
    maxδu = 100.0

    result = ST.solve_transition_backward(
        affineSys,
        E2,
        xnew,
        unew,
        sys.Uformat,
        sys.Wformat,
        problem.transition_cost,
        L,
        sdp_opt;
        λ = λ,
        maxδx = maxδx,
        maxδu = maxδu,
    )
    E1, cont, cost = result.source, result.controller, result.cost

    # Get results
    cost_eval(x, u) = problem.transition_cost(x, u)
    ETilde = LazySets.affine_map(
        affineSys.A + affineSys.B * cont.A,
        E1,
        affineSys.B * cont.c + affineSys.c,
    )
    U_used = LazySets.affine_map(cont.A, E1, cont.c)
    # Display results
    println()
    println("Max cost : ", cost)
    println("Volume of initial ellipsoid : ", UT.get_volume(E1))
    println("Input set volume : ", UT.get_volume(U_used))

    # Display the initial set, target set and the image of the initial ellipsoid under the linear model approximation
    fig1 = plot(; aspect_ratio = :equal)
    plot!(ETilde; color = :blue)
    display(fig1)

    # Display the cost of the controller
    fig2 = plot(;
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
    col1 = parse(ColorTypes.RGB, "#ff756f")
    col2 = parse(ColorTypes.RGB, "#64bc60")
    plot!(E2; color = col1)
    plot!(ETilde; color = col2)
    display(fig2)

    #Display the feasible input set and the input set effectively used
    col3 = RGB(0.2, 0.4, 0.8)
    fig3 = plot(;
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

    plot!(sys.U; color = col3, label = "", fillalpha = 0.4, linealpha = 1.0, linewidth = 2)
    plot!(U_used; color = :green, label = "Input set used")
    display(fig3)

    # savefig(fig1, "fig1.png")
    # savefig(fig2, "fig2.png")
    # savefig(fig3, "fig3.png")
    return
end

E2 = LazySets.Ellipsoid([4.0; 4.0], inv([2.0 0.2; 0.2 0.5]))
U = UT.LazySets.IntersectionArray([
    LazySets.Ellipsoid([0.0; 0.0], [25.0 0.0; 0.0 25.0]),
    LazySets.Ellipsoid([0.0; 0.0], [20.0 0.0; 0.0 30.0]),
    UT.box(SVector(-4.0, -5.0), SVector(4.0, 5.0)),
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

#fig 3
ρ = 0.0008
Wbound = 0.15
λ = 0.01
test_backward_transition(Wbound, E2, xnew, U, λ, ρ)
