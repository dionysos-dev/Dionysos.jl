include(joinpath(@__DIR__, "run_marche_avant.jl"))

using CSV
using DataFrames
using JLD2
using Plots

const RESULTS_DIR = joinpath(@__DIR__, "results")
const PLOTS_DIR = joinpath(RESULTS_DIR, "plots")
const SUMMARY_PATH = joinpath(RESULTS_DIR, "baseline_vs_scaled_summary.csv")
const DETAILED_PATH = joinpath(RESULTS_DIR, "baseline_vs_scaled_detailed.csv")
const BASELINE_RESULT_PATH = joinpath(RESULTS_DIR, "baseline_result.jld2")
const SCALED_RESULT_PATH = joinpath(RESULTS_DIR, "scaled_result.jld2")
const STRESS_SUMMARY_PATH = joinpath(RESULTS_DIR, "statistical_stress_test_summary.csv")
const STRESS_ROLLOUTS_PATH = joinpath(RESULTS_DIR, "statistical_rollouts.jld2")

function _initial_ellipsoid(certification)
    steps = sort(certification.steps; by = s -> s.k)
    first(steps).ellipsoid
end

function _plot_initial_ellipsoids(baseline_cert, scaled_cert)
    E_base = _initial_ellipsoid(baseline_cert)
    E_scaled = _initial_ellipsoid(scaled_cert)
    fig = plot(;
        aspect_ratio = :equal,
        legend = :topright,
        title = "Initial ellipsoid comparison (x,y)",
        size = (850, 650),
    )
    plot!(
        fig,
        project_ellipsoid_2d(E_base; dims = (1, 2));
        color = :steelblue3,
        opacity = 0.20,
        lw = 2,
        label = "baseline E0",
    )
    plot!(
        fig,
        project_ellipsoid_2d(E_scaled; dims = (1, 2));
        color = :darkorange2,
        opacity = 0.20,
        lw = 2,
        label = "scaled E0",
    )
    savefig(fig, joinpath(PLOTS_DIR, "initial_ellipsoid_comparison.png"))
    savefig(fig, joinpath(PLOTS_DIR, "initial_ellipsoid_comparison.pdf"))
end

function _plot_stress_rates()
    isfile(STRESS_SUMMARY_PATH) || return nothing
    df = CSV.read(STRESS_SUMMARY_PATH, DataFrame)
    labels = string.(df.scenario, " / ", df.method)

    success_yerr = (
        df.success_rate .- df.success_ci_low,
        df.success_ci_high .- df.success_rate,
    )
    fig = bar(
        labels,
        df.success_rate;
        yerror = success_yerr,
        ylim = (0, 1),
        ylabel = "success rate",
        title = "Monte Carlo stress test success rate",
        legend = false,
        xrotation = 35,
        size = (760, 520),
    )
    savefig(fig, joinpath(PLOTS_DIR, "stress_test_success_rates.png"))

    left_yerr = (
        df.left_chain_rate .- df.left_chain_ci_low,
        df.left_chain_ci_high .- df.left_chain_rate,
    )
    fig = bar(
        labels,
        df.left_chain_rate;
        yerror = left_yerr,
        ylim = (0, 1),
        ylabel = "left-chain rate",
        title = "Monte Carlo stress test ellipsoid-chain exits",
        legend = false,
        xrotation = 35,
        size = (760, 520),
    )
    savefig(fig, joinpath(PLOTS_DIR, "left_chain_rates.png"))
end

function _plot_sampled_rollouts()
    isfile(STRESS_ROLLOUTS_PATH) || return nothing
    data = load(STRESS_ROLLOUTS_PATH)
    baseline_rollouts = data["baseline_rollouts"]
    scaled_rollouts = data["scaled_rollouts"]
    nplot = min(40, length(baseline_rollouts), length(scaled_rollouts))

    fig = plot(;
        aspect_ratio = :equal,
        xlabel = "x",
        ylabel = "y",
        title = "Sampled closed-loop trajectories (x,y)",
        legend = :topright,
        size = (900, 680),
    )
    for i in 1:nplot
        xb = [x[1] for x in baseline_rollouts[i]]
        yb = [x[2] for x in baseline_rollouts[i]]
        xs = [x[1] for x in scaled_rollouts[i]]
        ys = [x[2] for x in scaled_rollouts[i]]
        plot!(fig, xb, yb; color = :steelblue3, alpha = 0.25, label = i == 1 ? "baseline" : "")
        plot!(fig, xs, ys; color = :darkorange2, alpha = 0.25, label = i == 1 ? "scaled" : "")
    end
    savefig(fig, joinpath(PLOTS_DIR, "sampled_trajectories_xy.png"))
end

function main()
    isfile(SUMMARY_PATH) || error("Missing summary CSV. Run 02_run_baseline_vs_scaled.jl first.")
    isfile(DETAILED_PATH) || error("Missing detailed CSV. Run 02_run_baseline_vs_scaled.jl first.")
    mkpath(PLOTS_DIR)

    detailed = CSV.read(DETAILED_PATH, DataFrame)
    fig = plot(;
        xlabel = "transition k",
        ylabel = "physical ellipsoid volume",
        yscale = :log10,
        title = "Baseline vs scaled certified volumes",
        size = (900, 620),
    )
    for method in unique(detailed.method)
        df = detailed[detailed.method .== method, :]
        sort!(df, :k)
        plot!(fig, df.k, df.volume_Ek; label = method, lw = 2, marker = :circle, ms = 3)
    end
    savefig(fig, joinpath(PLOTS_DIR, "volume_vs_transition_baseline_scaled.png"))

    if isfile(BASELINE_RESULT_PATH) && isfile(SCALED_RESULT_PATH)
        baseline_cert = load(BASELINE_RESULT_PATH, "certification")
        scaled_cert = load(SCALED_RESULT_PATH, "certification")
        _plot_initial_ellipsoids(baseline_cert, scaled_cert)
    end

    _plot_stress_rates()
    _plot_sampled_rollouts()

    println("plots: ", PLOTS_DIR)
    return PLOTS_DIR
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
