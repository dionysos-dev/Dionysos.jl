using StaticArrays, LinearAlgebra, IntervalArithmetic, Random
using JuMP, Mosek, MosekTools
using Plots, Colors
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

include("../../problems/non_linear.jl")

const FALLBACK_URL = "mosek://solve.mosek.com:30080"

function trial(E2, c, ρ, Ubound, Wbound, λ)
    U = UT.HyperRectangle(SVector(-Ubound, -Ubound), SVector(Ubound, Ubound))
    W = UT.HyperRectangle(SVector(-Wbound, -Wbound), SVector(Wbound, Wbound))
    problem = NonLinear.problem(; U = U, W = W, noise = true, μ = ρ)
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
    sdp_opt = optimizer_with_attributes(Mosek.Optimizer, MOI.Silent() => true)
    MOI.set(sdp_opt, MOI.RawOptimizerAttribute("fallback"), FALLBACK_URL)
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

    # Get results 
    if cont == nothing
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
    return (success, max_cost, init_set_volume, input_set_volume)
end

function compute_pareto_front(E2, c, ρ, Ubound, Wbound, λ_span)
    # Data vectors
    success_vector = zeros(length(λ_span))
    max_cost_vector = zeros(length(λ_span))
    init_set_volume_vector = zeros(length(λ_span))
    input_set_volume_vector = zeros(length(λ_span))

    for i in 1:length(λ_span)
        println("Problem solved : ", i, " / ", length(λ_span))
        (success, max_cost, init_set_volume, input_set_volume) =
            trial(E2, c, ρ, Ubound, Wbound, λ_span[i])
        success_vector[i] = success
        if success
            max_cost_vector[i] = max_cost
            init_set_volume_vector[i] = init_set_volume
            input_set_volume_vector[i] = input_set_volume
        else
            max_cost_vector[i] = Inf
            init_set_volume_vector[i] = 0.0
            input_set_volume_vector[i] = 0.0
        end
    end
    return max_cost_vector, init_set_volume_vector, input_set_volume_vector
end

function plot_pareto_front!(E2, c, ρ, Ubound, Wbound, λ_span, mycolorMap; lbls = false)
    max_cost_vector, init_set_volume_vector, input_set_volume_vector =
        compute_pareto_front(E2, c, ρ, Ubound, Wbound, λ_span)
    colors = [UT.get_color(mycolorMap, λ_val) for λ_val in λ_span]
    return plot!(
        -max_cost_vector,
        init_set_volume_vector;
        marker = :circle,
        line = :path,
        color = colors,
        label = "\$\\omega_{max} = $Wbound\$",
    )
end

function plot_pareto_front(E2, c, ρ, Ubound, Wbound, λ_span)
    colormap = Colors.colormap("Blues")
    mycolorMap = UT.Colormap([λ_span[1], λ_span[end]], colormap)
    colors = [UT.get_color(mycolorMap, λ_val) for λ_val in λ_span]

    fig1 = plot(;
        xtickfontsize = 10,
        ytickfontsize = 10,
        guidefontsize = 18,
        titlefontsize = 16,
        legend = false,
    )
    xlabel!("-\$\\widetilde{\\mathcal{J}}\$")
    ylabel!("\${\\rm vol}(E_1)\$")
    title!("Pareto front (\$\\lambda\$)")
    max_cost_vector, init_set_volume_vector, input_set_volume_vector =
        compute_pareto_front(E2, c, ρ, Ubound, Wbound, λ_span)
    plot!(
        -max_cost_vector,
        init_set_volume_vector;
        marker = :circle,
        color = colors,
        line = :path,
        linecolor = :black,
    )
    plot!(mycolorMap)
    return display(fig1)
end

function plot_pareto_front_Wspan(E2, c, ρ, Ubound, Wbound_span, λ_span)
    colormap = Colors.colormap("Blues")
    mycolorMap = UT.Colormap([λ_span[1], λ_span[end]], colormap)
    colors = [UT.get_color(mycolorMap, λ_val) for λ_val in λ_span]
    lColors = [:red, :green, :black, :blue, :yellow]
    fig = plot(;
        xtickfontsize = 10,
        ytickfontsize = 10,
        guidefontsize = 18,
        titlefontsize = 16,
    )
    xlabel!("-\$\\widetilde{\\mathcal{J}}\$")
    ylabel!("\${\\rm vol}(E_1)\$")
    title!("Pareto front (\$\\lambda\$)")
    for (i, Wbound) in enumerate(Wbound_span)
        max_cost_vector, init_set_volume_vector, input_set_volume_vector =
            compute_pareto_front(E2, c, ρ, Ubound, Wbound, λ_span)
        plot!(
            -max_cost_vector,
            init_set_volume_vector;
            marker = :circle,
            markersize = 4,
            color = colors,
            line = :path,
            linecolor = lColors[i],
            label = "\$\\omega_{max} = $Wbound\$",
        )
    end
    plot!(mycolorMap)
    return display(fig)
end

function plot_pareto_front_ρspan(E2, c, ρ_span, Ubound, Wbound, λ_span)
    colormap = Colors.colormap("Blues")
    mycolorMap = UT.Colormap([λ_span[1], λ_span[end]], colormap)
    colors = [UT.get_color(mycolorMap, λ_val) for λ_val in λ_span]
    lColors = [:red, :green, :black, :blue, :yellow]
    fig = plot(;
        xtickfontsize = 10,
        ytickfontsize = 10,
        guidefontsize = 18,
        titlefontsize = 16,
    )
    xlabel!("-\$\\widetilde{\\mathcal{J}}\$")
    ylabel!("\${\\rm vol}(E_1)\$")
    title!("Pareto front (\$\\lambda\$)")
    for (i, ρ) in enumerate(ρ_span)
        max_cost_vector, init_set_volume_vector, input_set_volume_vector =
            compute_pareto_front(E2, c, ρ, Ubound, Wbound, λ_span)
        plot!(
            -max_cost_vector,
            init_set_volume_vector;
            marker = :circle,
            markersize = 4,
            color = colors,
            line = :path,
            linecolor = lColors[i],
            label = "\$\\mu = $ρ\$",
        )
    end
    plot!(mycolorMap)
    return display(fig)
end

function contour_plot(
    E2,
    c,
    Ubound,
    λ,
    ρ_span,
    Wbound_span;
    cost = true,
    levels1 = 6,
    levels2 = 6,
)
    # Data vectors
    success_vector = zeros(length(ρ_span), length(Wbound_span))
    max_cost_vector = zeros(length(ρ_span), length(Wbound_span))
    init_set_volume_vector = zeros(length(ρ_span), length(Wbound_span))
    input_set_volume_vector = zeros(length(ρ_span), length(Wbound_span))

    for i in 1:length(ρ_span)
        for j in 1:length(Wbound_span)
            println(
                "Problem solved : ",
                (i - 1) * length(Wbound_span) + j,
                " / ",
                length(ρ_span) * length(Wbound_span),
            )
            (success, max_cost, init_set_volume, input_set_volume) =
                trial(E2, c, ρ_span[i], Ubound, Wbound_span[j], λ)
            success_vector[i, j] = success
            if success
                max_cost_vector[i, j] = max_cost
                init_set_volume_vector[i, j] = init_set_volume
                input_set_volume_vector[i, j] = input_set_volume
            else
                max_cost_vector[i, j] = Inf
                init_set_volume_vector[i, j] = 0.0
                input_set_volume_vector[i, j] = 0.0
            end
        end
    end

    fig1 = contour(
        ρ_span,
        Wbound_span,
        max_cost_vector';
        levels = levels1,
        color = :viridis,
        clabels = true,
        cbar = true,
        lw = 3,
        xtickfontsize = 10,
        ytickfontsize = 10,
        guidefontsize = 18,
        titlefontsize = 19,
    )
    xlabel!("\$\\mu\$")
    ylabel!("\$\\omega_{max}\$")
    title!("\$\\widetilde{\\mathcal{J}}\$")
    display(fig1)

    fig2 = contour(
        ρ_span,
        Wbound_span,
        init_set_volume_vector';
        levels = levels2,
        color = :viridis,
        clabels = true,
        cbar = true,
        lw = 3,
        xtickfontsize = 10,
        ytickfontsize = 10,
        guidefontsize = 18,
        titlefontsize = 19,
    )
    xlabel!("\$\\mu\$")
    ylabel!("\$\\omega_{max}\$")
    title!("\${\\rm vol}(E_1)\$")
    display(fig2)

    return cost ? display(fig1) : display(fig2)
end

E2 = UT.Ellipsoid([2.0 0.2; 0.2 0.5], [3.0; 3.0])
c = SVector{2, Float64}([1.0; 1.0])

#########################################################################################
#### fig1: plot pareto front with repsect to λ #####
#########################################################################################
ρ = 0.0008
Ubound = 5.0
Wbound = 0.01
λ_span = 0.0005:0.01:1.0

# fig 1.1
plot_pareto_front(E2, c, ρ, Ubound, Wbound, λ_span)

# fig 1.2
Wbound_span = [0.0, 0.1, 0.2, 0.4, 0.5]
plot_pareto_front_Wspan(E2, c, ρ, Ubound, Wbound_span, λ_span)

# fig 1.3
ρ_span = [0.0003, 0.0006, 0.001, 0.003]
plot_pareto_front_ρspan(E2, c, ρ_span, Ubound, Wbound, λ_span)

#########################################################################################
######### fig2: plot for λ=1.0 the cost as a function of noise and non-linearity #########
#########################################################################################
Ubound = 5.0
λ = 1.0
ρ_span = 0.0:0.0004:0.004
Wbound_span = 0:0.02:0.6
contour_plot(
    E2,
    c,
    Ubound,
    λ,
    ρ_span,
    Wbound_span;
    cost = true,
    levels1 = [6.8, 7.3, 7.8, 8.4, 9, 9.4],
)

#########################################################################################
#### fig3: plot for λ=0.0 the volume value as a function of noise and non-linearity #####
#########################################################################################
Ubound = 5.0
λ = 0.0
ρ_span = 0.0:0.0004:0.004
Wbound_span = 0:0.02:0.6;
contour_plot(
    E2,
    c,
    Ubound,
    λ,
    ρ_span,
    Wbound_span;
    cost = false,
    levels2 = [1, 3, 6, 12, 18, 24, 30, 36],
)
