using StaticArrays
using LinearAlgebra
using JuMP
using Clarabel
using JuMP
using JLD2
using CDDLib

import HybridSystems
using LazySets
using Plots

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const OP = DI.Optim
const AB = OP.Abstraction

const PCLF = UT.PathCompleteFramework
const PCLF_FILE = joinpath(@__DIR__, "pclf_case.jld2")
const CLF_FILE = joinpath(@__DIR__, "clf_case.jld2")

function export_optimizer_jld2(opt, filename::AbstractString)
    jldopen(filename, "w") do f
        f["format_version"] = 1
        return f["optimizer"] = opt
    end
    return nothing
end

function import_optimizer_jld2(filename::AbstractString)
    return jldopen(filename, "r") do f
        v = f["format_version"]
        v == 1 || error("Unsupported optimizer file format_version=$v")
        return f["optimizer"]
    end
end

function run_and_save_case(
    problem,
    pclf,
    filename;
    verbose::Bool = true,
    atol::Float64 = 1e-4,
    level_tol::Float64 = 1e-2,
    max_slices::Int = 20,
)
    optimizer = MOI.instantiate(AB.PCLFBisimulationQuotient.OptimizerBisimulationQuotient)

    MOI.set(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient_problem"), problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("pclf"), pclf)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), verbose ? 1 : 0)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("atol"), atol)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("level_tol"), level_tol)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("max_slices"), max_slices)

    MOI.optimize!(optimizer)
    construction_time =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("construction_time_sec"))
    println("Construction time = ", construction_time)

    export_optimizer_jld2(optimizer, filename)

    return optimizer
end

function load_bisimulation(filename)
    opt = import_optimizer_jld2(filename)
    return MOI.get(opt, MOI.RawOptimizerAttribute("bisimulation_quotient"))
end

function compute_solutions()
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

    v1 = [1.0, 0.0]
    v2 = [1.0, 1.0]
    v3 = [0.0, 1.0]
    v4 = [-1.0, 1.0]
    v5 = [-1.0, 0.0]

    C1 = hcat(v1, v2)
    C2 = hcat(v2, v3)
    C3 = hcat(v3, v4)
    C4 = hcat(v4, v5)

    partitions = Dict(
        1 => [C1, C2, C3, C4],
        2 => [C1, C2, C3, C4],
        3 => [C1, C2, C3, C4],
        4 => [C1, C2, C3, C4],
    )
    lp_optimizer = JuMP.optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => 1000)
    pclf_poly =
        PCLF.compute_polyhedral_pieces_pclf(f, graph, lp_optimizer, partitions; MLF = true)
    println("Computed JSR upper bound / contraction rate = ", pclf_poly.JSRapprox)

    #opt_pclf = run_and_save_case(problem, pclf_poly, PCLF_FILE)

    # -- Computation of the induced CLF: 

    states, trans, alphabet = PCLF.build_observer_graph(graph)
    clf_poly = PCLF.build_common_lyapunov(pclf_poly)

    return opt_clf = run_and_save_case(problem, clf_poly, CLF_FILE)
end

# -- Definition of the problem:
A1 = (1.0/10.0)*[1.5519 0.4474; 7.6412 7.4716]
A2 = (1.0/10.0)*[0.4750 9.1755; 1.8955 0.1850]
f = HybridSystems.discreteswitchedsystem([A1, A2])

p = 1.7
X = HPolytope([
    HalfSpace([1.0, 0.0], p),
    HalfSpace([-1.0, 0.0], p),
    HalfSpace([0.0, 1.0], p),
    HalfSpace([0.0, -1.0], p),
])

R1 = HPolytope([
    HalfSpace([1.0, 0.0], 1.5),
    HalfSpace([-1.0, 0.0], -0.8),
    HalfSpace([0.0, 1.0], 1.5),
    HalfSpace([0.0, -1.0], -0.8),
])

observation_regions = [R1]
problem = PR.BisimulationQuotientProblem(f, X, observation_regions)

# -- The histogram plot:
#compute_solutions()

# -- The histogram plot:

bisim_pclf = load_bisimulation(PCLF_FILE)
bisim_clf = load_bisimulation(CLF_FILE)

#fig = plot(; aspect_ratio = :equal);
#plot!(bisim_pclf; what = :states, by = :states, node = (1,), show_contours = false)
#plot!(problem; opacity = 0.2)
#display(fig)

n_parts_pclf, n_faces_pclf = AB.PCLFBisimulationQuotient.cell_complexities(bisim_pclf)
n_parts_clf, n_faces_clf = AB.PCLFBisimulationQuotient.cell_complexities(bisim_clf)

p1 = histogram(
    n_parts_pclf;
    bins = 0:2:450,
    xlabel = "# polytopes",
    ylabel = "# cells",
    title = "Pieces per semilinear cell",
    yscale = :log10,
    xlims = (0, 450),
    ylims = (:auto, 1e4),
    yticks = ([1, 10, 100, 1000], ["10⁰", "10¹", "10²", "10³"]),
    legend = :topright,
    alpha = 0.45,
    label = "PCLF",
)

histogram!(p1, n_parts_clf; bins = 0:2:450, alpha = 0.45, label = "CLF", yscale = :log10)

p2 = histogram(
    n_faces_pclf;
    bins = 0:10:2000,
    xlabel = "# faces",
    ylabel = "# cells",
    title = "Faces per semilinear cell",
    yscale = :log10,
    xlims = (0, 2000),
    ylims = (:auto, 1e4),
    yticks = ([1, 10, 100, 1000], ["10⁰", "10¹", "10²", "10³"]),
    legend = :topright,
    alpha = 0.45,
    label = "PCLF",
)

histogram!(p2, n_faces_clf; bins = 0:10:2000, alpha = 0.45, label = "CLF", yscale = :log10)

fig = plot(
    p1,
    p2;
    layout = (1, 2),
    size = (1400, 500),
    margin = 8Plots.mm,
    left_margin = 8Plots.mm,
    right_margin = 8Plots.mm,
    top_margin = 8Plots.mm,
    bottom_margin = 8Plots.mm,
    titlefontsize = 14,
    guidefontsize = 12,
    tickfontsize = 10,
)

display(fig)
