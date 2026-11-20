# E4: the running example with mode 1 scaled by 1.2, which diverges on its own (paper, Section 10.4).
#
#     julia --project=. experiments/E4_unstable_plant.jl   [CACHE=0] [PROBES=n]
#
# The plant switches freely and is not stable under arbitrary switching, so no common Lyapunov
# function exists and the construction of Gol et al. has no input. The controller restricts its own
# switching to the language of the graph that forbids playing mode 1 twice in a row, on which the
# system admits a certificate. The graph respects the free constraint graph without covering it, so
# synthesis transfers to the plant and verification does not.
include(joinpath(dirname(@__DIR__), "common.jl"))
import Clarabel
import JuMP

const PROBES = parse(Int, get(ENV, "PROBES", "211"))
const SCALING = 1.2

As = [SCALING * A1, A2]
f = switched(As; controlled = true)
problem = PR.BisimulationQuotientProblem(f, X, [R1, R2, R3])

free = PCLF.common_lyapunov_graph([1, 2])
restriction = PCLF.edgeList_to_LabDigraph([(0, 0, 2), (0, 1, 1), (1, 0, 2)])
init = Dict(:clf => [0, 1])
@assert PCLF.respects_constraint(restriction, init, free)
@assert !PCLF.covers_constraint(restriction, init, free)

println("="^84, "\nE4 — the running example with mode 1 scaled by $SCALING\n", "="^84)
spectral = maximum(abs, LinearAlgebra.eigvals(As[1]))
Printf.@printf("spectral radius  ρ(%.1f A1) = %.4f\n", SCALING, spectral)

# The single-node certificate fails, as it must: its self-loop on mode 1 needs ρ(A1) < γ < 1.
single = try
    certificate(f, PCLFProducer.generate_DeBruijn_edges(2, 0))
catch err
    err isa ErrorException || rethrow()
    nothing
end
println(
    "common Lyapunov function ",
    single === nothing ? "none, as expected" : "FOUND (unexpected)",
)
single === nothing || @warn "a common Lyapunov function was found although ρ(A1) > 1"

# The rate cannot go below 0.8936, the spectral bound of the cycle 12 of the restriction, so about ten
# slices are needed whatever the template. Among the templates and solvers tried, the uniform conic
# partition of order 3 solved by Clarabel needs the fewest, 9, with the outer level containing the
# domain of interest X, as in the paper. The template of the running example gives 13.
clarabel = JuMP.optimizer_with_attributes(
    Clarabel.Optimizer,
    "max_iter" => 4000,
    "verbose" => false,
    "tol_feas" => 1e-6,
    "tol_gap_abs" => 1e-6,
    "tol_gap_rel" => 1e-6,
)
cert = certificate(
    f,
    restriction;
    outer = :box,
    cones = PCLFProducer.conic_partitions_2d(3),
    optimizer = clarabel,
)
Printf.@printf("rate on the restriction %.6f, %d slices\n", cert.γ, cert.N)
# With 12 855 cells, the inset of 10⁻³ of the other experiments leaves 5 to 8 % of the working set
# uncovered; 10⁻⁴ is the smallest inset the solver documents as safe.
built = build("E4_restricted_order3_atol4", problem, cert, cert.pclf; atol = 1e-4)
stats = cell_stats(built.quotient)
Printf.@printf("cells            %d\n", stats.cells)

synthesis = solve_optimizer(f, built)
won = Set(MOI.get(synthesis, MOI.RawOptimizerAttribute("controllable_set")))
points, area = probe_grid(PROBES)
cells_at = Dict(nd => locate(built.quotient, nd, points) for nd in restriction.verts)
domain = covered(cells_at, init[:clf])
certified = project(cells_at, won, init[:clf]; quantifier = :exists) .& domain
Printf.@printf(
    "certified        ∃-set of area %.1f on a domain of area %.1f\n",
    area * count(certified),
    area * count(domain),
)

write_numbers(
    "E4_unstable_plant",
    Pair{String, Any}[
        "EFourSpectral" => fmt(spectral; digits = 3),
        "EFourRate" => fmt(cert.γ; digits = 3),
        "EFourN" => fmt(cert.N - 1),
        "EFourCells" => fmt(stats.cells),
        "EFourDomain" => fmt(area * count(domain)),
        "EFourExists" => fmt(area * count(certified)),
        # A single build: repeating one of this size for a median costs more than it tells.
        "EFourSeconds" => fmt_seconds(built.seconds),
    ],
)
println("numbers written to results/numbers_E4_unstable_plant.tex")

# A closed loop that uses the diverging mode, so that the restriction shows: among a sample of the
# certified probe points, the one whose synthesised word plays mode 1 most often, the nearest to the
# initial point of the running example among those.
controller = PCQ.solve_concrete_problem(synthesis)
function closed_loop(x)
    memory = try
        PCQ.initial_controller_memory(synthesis, x)
    catch
        return nothing
    end
    return PCQ.simulate_closed_loop(f, controller, x, memory; N = 40)
end
# The search simulates hundreds of closed loops and costs more than the quotient; its result is cached.
start_path = joinpath(CACHE_DIR, "E4_start.jld2")
start = if USE_CACHE && isfile(start_path)
    JLD2.load(start_path, "start")
else
    candidates = points[findall(certified)][1:7:end]
    runs = [(x, closed_loop(x)) for x in candidates]
    runs = [(x, r) for (x, r) in runs if r !== nothing]
    score((x, r)) = (count(==(1), r.U), -sum(abs2, x .- X0))
    chosen = first(runs[argmax(map(score, runs))])
    USE_CACHE && JLD2.jldsave(start_path; start = chosen)
    chosen
end
closed = closed_loop(start)
println("closed loop from ", start, ": word ", join(closed.U))
@assert !occursin("11", join(closed.U))

panel = category_panel(
    PROBES,
    points,
    [d ? (c ? 2 : 1) : 0 for (d, c) in zip(domain, certified)],
    [LOST, WON];
    title = "synthesis under the restriction",
)
Plots.plot!(
    panel,
    first.(closed.X),
    last.(closed.X);
    color = :white,
    linewidth = 1.6,
    marker = :circle,
    markersize = 3.5,
    markerstrokecolor = :black,
)
figure = Plots.plot(panel; size = (520, 520), dpi = 200)
Plots.savefig(figure, joinpath(RESULTS, "fig_E4_unstable_plant.png"))
write_numbers(
    "E4_unstable_plant_trajectory",
    Pair{String, Any}[
        "EFourWord" => join(closed.U),
        "EFourStart" => Printf.@sprintf("(%.1f,%.1f)", start[1], start[2]),
    ],
)
println("figure written to results/fig_E4_unstable_plant.png")
