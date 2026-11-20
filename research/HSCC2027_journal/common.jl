# Shared setup of the experiments: the running example of the paper, the certificate search, the
# construction of a quotient, the solution of the specification on it, and the projection of the
# certified cells onto the plant. Included by every script in `experiments/`.
import LinearAlgebra
import HiGHS
import LazySets
import Printf
import Spot
import JLD2
import Statistics
import HybridSystems
import MathOptInterface as MOI
import Plots
import Dionysos

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const PCQ = DI.Optim.Abstraction.PCLFBisimulationQuotient
const PCLF = UT.PathCompleteFramework
@isdefined(PCLFProducer) ||
    include(joinpath(dirname(@__DIR__), "common", "pclf_producer.jl"))

const RESULTS = joinpath(@__DIR__, "results")
const CACHE_DIR = joinpath(@__DIR__, "cache")
const USE_CACHE = get(ENV, "CACHE", "1") == "1"
# Build times follow the conference experiments: every quotient of an experiment is built first, then
# each is rebuilt after an untimed warm-up, SAMPLES times, and the median is reported. TIMING=0 skips.
const TIMING = get(ENV, "TIMING", "1") == "1"
const SAMPLES = parse(Int, get(ENV, "SAMPLES", "3"))
mkpath(RESULTS)

# The conventions of the conference experiments (research/HSCC2027), kept so that the unconstrained
# arm reproduces their numbers: the template, the inset of every set difference, and the scaling of
# the certificate that puts every region inside the level SCALE.
const ATOL = 1e-3
const SCALE = 6.0
const RAYS_DEG = [0.0, 32.6, 65.8, 91.5, 112.6, 131.4, 152.0, 180.0]
const CONES = let a = deg2rad.(RAYS_DEG)
    [hcat([cos(a[i]), sin(a[i])], [cos(a[i + 1]), sin(a[i + 1])]) for i in 1:(length(a) - 1)]
end

# ── the running example: Gol, Ding, Lazar and Belta, Example 3.1 ────────────────────────────────────
const A1 = [-0.65 0.32; -0.42 -0.92]
const A2 = [0.65 0.32; -0.42 -0.92]

"The switched system of `As`, with the modes owned by the controller or by the environment."
function switched(As; controlled::Bool)
    owner =
        controlled ? HybridSystems.ControlledSwitching() :
        HybridSystems.AutonomousSwitching()
    return ST.with_switching(HybridSystems.discreteswitchedsystem(As), owner)
end

const XBOX = 5.9
const X = LazySets.HPolytope([
    LazySets.HalfSpace([1.0, 0.0], XBOX),
    LazySets.HalfSpace([-1.0, 0.0], XBOX),
    LazySets.HalfSpace([0.0, 1.0], XBOX),
    LazySets.HalfSpace([0.0, -1.0], XBOX),
])
const R1 = LazySets.HPolytope(
    [-0.9869 -0.1615; -0.0931 0.9957; 0.9659 0.2587; 0.0825 -0.9966],
    [6.6767, 9.2315, 2.3700, -5.9038],
)
const R2 = LazySets.HPolytope(
    [0.9993 0.0363; -0.7743 -0.6329; 0.5463 0.8376],
    [-2.1809, 6.3754, -4.8983],
)
const R3 = LazySets.HPolytope(
    [-0.9946 -0.1041; 0.5277 0.8494; 0.9999 0.0146; -0.1191 -0.9929],
    [-5.5771, 5.3510, 9.1600, 6.2406],
)
const REGIONS = [("R1", R1), ("R2", R2), ("R3", R3)]
# The solver takes co-safe formulas in which the terminal set D is observed. On the behaviours of the
# construction, which are eventually region-free, this formula is equivalent to
# G !R2 & F R1 & G(R3 -> X !R1), the form the paper uses (Remark 7.4).
const PHI = Spot.@ltl_str "((!R2 U D) & F(R1) & ((R3 -> X(!R1)) U D))"
const AP_TO_OBS = Dict(:D => -1, :R1 => 1, :R2 => 2, :R3 => 3)
const X0 = [-3.8, -7.25]

# ── certificate, quotient, specification ────────────────────────────────────────────────────────────
"""
    certificate(f, graph; with_common = false, outer = :regions, cones = CONES, optimizer = HiGHS.Optimizer)

A polyhedral certificate on `graph` with the conic template `cones`, by default that of the conference
experiments, scaled so that the level `SCALE` contains every region, with its rate `γ`, its outer
level `ΓX` and its number of slices `N`. The LP `optimizer` matters beyond its tolerances: the
template admits many optimal certificates, and each solver returns its own, with its own level sets.

The outer level is `SCALE`, which contains every region, with `outer = :regions`, the rule of the
conference experiments. With `outer = :box`, it is the smallest level whose sublevel sets contain the
domain of interest `X`, up to a factor 1.01, which is the definition of the paper (Section 6.1); a
region may then stick out of the working set, and the answer covers what lies inside.

`N` is the first level at which every sublevel set clears every region, for the certificate and,
when `with_common`, for the common Lyapunov function it induces, as the conference experiments do so
that both constructions share their slices. It counts the levels τ_0, …, τ_N of the paper, so the
paper's number of slices is `N - 1`.
"""
function certificate(
    f,
    graph;
    with_common::Bool = false,
    outer::Symbol = :regions,
    cones = CONES,
    optimizer = HiGHS.Optimizer,
)
    outer in (:regions, :box) ||
        throw(ArgumentError("outer must be :regions or :box, got :$outer"))
    nodes = collect(graph.verts)
    pclf = PCLFProducer.compute_polyhedral_pieces_pclf(
        f,
        graph,
        optimizer,
        Dict(nd => cones for nd in nodes);
        MLF = true,
    )
    isinf(pclf.JSRapprox) && error("no certificate for this graph and template")
    level(p) = maximum(
        PCQ.gamma_cover_set(p.pieces[nd], UT._as_hpolytope(R)) for nd in nodes for
        (_, R) in REGIONS
    )
    α = level(pclf) / SCALE
    for nd in nodes
        pclf.pieces[nd].G = pclf.pieces[nd].G / α
    end
    common = with_common ? PCLF.build_common_lyapunov(pclf) : nothing
    γ = pclf.JSRapprox
    ΓX = outer === :regions ? level(pclf) : PCQ.compute_tau_X(pclf, X; safety_factor = 1.01)
    regions_h = [UT._as_hpolytope(R) for (_, R) in REGIONS]
    certs = with_common ? (pclf, common) : (pclf,)
    clears(j) = all(
        c -> PCQ.all_nodes_clear_regions(c, ΓX * γ^(j - 1), regions_h; tol = 1e-2),
        certs,
    )
    N = something(findfirst(clears, 1:80), -1)
    N > 0 || error("no region-free terminal level among the first 80")
    return (; pclf, common, γ, ΓX, N, graph)
end

"""
    build(name, problem, cert, pieces; atol = ATOL)

The quotient of `problem` on `pieces`, with the levels of `cert`, or its cached copy. `atol` is the
inset of every set difference: the cells leave uncovered strips of that width between them, whose
total area grows with the number of cells.
"""
function build(name, problem, cert, pieces; atol::Float64 = ATOL)
    path = joinpath(CACHE_DIR, "$(name)_$(cert.N).jld2")
    if USE_CACHE && isfile(path)
        return JLD2.jldopen(path, "r") do file
            return (; quotient = file["quotient"], D = file["D"], seconds = file["seconds"])
        end
    end
    built = construct(problem, cert, pieces; atol)
    if USE_CACHE
        mkpath(CACHE_DIR)
        JLD2.jldopen(path, "w") do file
            file["quotient"] = built.quotient
            file["D"] = built.D
            return file["seconds"] = built.seconds
        end
    end
    return built
end

"""
    build_seconds(problem, cert, pieces; atol = ATOL)

The build time of the quotient, as the conference experiments measure it: the median of `SAMPLES`
rebuilds after an untimed warm-up, each preceded by a garbage collection. Call it once every
quotient of the experiment exists, since a timed build pays garbage collection in proportion to the
live heap. `NaN` when `TIMING=0`.
"""
function build_seconds(problem, cert, pieces; atol::Float64 = ATOL)
    TIMING || return NaN
    construct(problem, cert, pieces; atol)
    times = Float64[]
    for _ in 1:SAMPLES
        GC.gc()
        push!(times, construct(problem, cert, pieces; atol).seconds)
    end
    return Statistics.median(times)
end

"Build the quotient without the cache, timing the construction alone."
function construct(problem, cert, pieces; atol::Float64 = ATOL)
    opt = MOI.instantiate(PCQ.OptimizerBisimulationQuotient)
    MOI.set(opt, MOI.RawOptimizerAttribute("bisimulation_quotient_problem"), problem)
    MOI.set(opt, MOI.RawOptimizerAttribute("pclf"), pieces)
    MOI.set(opt, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.set(opt, MOI.RawOptimizerAttribute("atol"), atol)
    MOI.set(opt, MOI.RawOptimizerAttribute("nb_levels"), cert.N)
    MOI.set(opt, MOI.RawOptimizerAttribute("max_slices"), cert.N)
    MOI.set(opt, MOI.RawOptimizerAttribute("ΓX"), cert.ΓX)
    seconds = @elapsed MOI.optimize!(opt)
    return (;
        quotient = MOI.get(opt, MOI.RawOptimizerAttribute("bisimulation_quotient")),
        D = MOI.get(opt, MOI.RawOptimizerAttribute("D")),
        seconds,
    )
end

"""
    solve(system, built; semantics = :language)

The cells of `built` from which the specification holds: for some switching signal when `system`
owns its modes (synthesis), for every one when the environment owns them (verification). The
adversary of verification plays the edges of the certificate graph, which is the augmented system
of the paper; a mode a node does not enable is a move it cannot make.
"""
function solve(system, built; semantics::Symbol = :language)
    opt = solve_optimizer(system, built; semantics)
    return Set(MOI.get(opt, MOI.RawOptimizerAttribute("controllable_set")))
end

"The solved optimizer behind [`solve`](@ref), from which a synthesis run also yields a controller."
function solve_optimizer(system, built; semantics::Symbol = :language)
    regions = Dict{Symbol, Any}(:D => built.D, :R1 => R1, :R2 => R2, :R3 => R3)
    problem = PR.CoSafeLTLProblem(
        system,
        LazySets.Hyperrectangle(; low = X0, high = X0),
        DI.spot_stepper(PHI),
        regions,
        Dict{Symbol, Any}(ap => DI.Mapping.INNER for ap in keys(regions)),
    )
    opt = MOI.instantiate(PCQ.OptimizerCoSafeLTLOnQuotient)
    MOI.set(opt, MOI.RawOptimizerAttribute("concrete_problem"), problem)
    MOI.set(opt, MOI.RawOptimizerAttribute("bisimulation_quotient"), built.quotient)
    MOI.set(opt, MOI.RawOptimizerAttribute("ap_to_obs"), AP_TO_OBS)
    MOI.set(opt, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.set(opt, MOI.RawOptimizerAttribute("switching_semantics"), semantics)
    # The default stops once the initial state is decided and certifies a handful of cells only.
    MOI.set(opt, MOI.RawOptimizerAttribute("early_stop"), false)
    MOI.optimize!(opt)
    return opt
end

"Number of cells, polytopes and facets of a quotient."
function cell_stats(quotient)
    parts, faces = PCQ.cell_complexities(quotient)
    return (;
        cells = length(parts),
        polytopes = sum(parts),
        mean_facets = Statistics.mean(faces),
        max_facets = maximum(faces),
    )
end

# ── projection onto the plant ───────────────────────────────────────────────────────────────────────
# The quotients extend to the level SCALE of the pieces, beyond the box `X`: region R1 and the
# initial point X0 lie outside it. The probe box encloses every region and the working sets.
const PROBE_BOX = 10.5

"Probe points on a regular `n × n` grid over `[-PROBE_BOX, PROBE_BOX]²`, and the area of each."
function probe_grid(n::Int)
    ticks = range(-PROBE_BOX, PROBE_BOX; length = n)
    return [[a, b] for b in ticks for a in ticks], (step(ticks))^2
end

"The cell of `quotient` at `node` containing each point, `0` where no cell does."
function locate(quotient, node, points)
    return [something(PCQ._find_qid_in_node(quotient, node, x), 0) for x in points]
end

"""
    project(cells_at, certified, nodes; quantifier)

`PCQ.certified_at` evaluated on a whole grid at once, the point locations `cells_at[node]` being
shared between the sets: a point is certified under `:exists` when its cell at one of `nodes` is in
`certified`, and under `:forall` when its cells at all of them are.
"""
function project(cells_at, certified, nodes; quantifier::Symbol)
    holds(node) = [q != 0 && q in certified for q in cells_at[node]]
    columns = [holds(node) for node in nodes]
    combine = quantifier === :exists ? (a, b) -> a .| b : (a, b) -> a .& b
    return reduce(combine, columns)
end

"Points located in a cell at every one of `nodes`."
covered(cells_at, nodes) =
    reduce((a, b) -> a .& b, [cells_at[node] .!= 0 for node in nodes])

"Write `pairs` as LaTeX macros `\\name{value}` to `results/numbers_<name>.tex`."
function write_numbers(name, pairs)
    open(joinpath(RESULTS, "numbers_$name.tex"), "w") do io
        println(io, "% Written by experiments/$name.jl. Do not edit.")
        for (key, value) in pairs
            println(io, "\\newcommand{\\", key, "}{", value, "}")
        end
    end
    return nothing
end

fmt(x::Integer) = string(x)
fmt(x::Real; digits = 1) = Printf.@sprintf("%.*f", digits, x)
"A build time in seconds with three significant digits, whole seconds from 100 s, `--` when untimed."
function fmt_seconds(t::Real)
    isnan(t) && return "--"
    t >= 100 && return string(round(Int, t))
    return Printf.@sprintf("%.*f", max(0, 2 - floor(Int, log10(t))), t)
end

# ── figures ─────────────────────────────────────────────────────────────────────────────────────────
const WON = Plots.RGB(0.36, 0.70, 0.36)
const LOST = Plots.RGB(0.85, 0.33, 0.31)
const UNSOUND = Plots.RGB(0.98, 0.65, 0.13)
const MISSED = Plots.RGB(0.40, 0.62, 0.95)

"The outline of a planar polytope, as a closed polyline."
function outline(P)
    vs = LazySets.convex_hull(LazySets.vertices_list(P))
    return [v[1] for v in [vs; vs[1:1]]], [v[2] for v in [vs; vs[1:1]]]
end

"""
    fill_gaps(category, n)

For display only: every blank point of the `n × n` grid with at least five drawn neighbours takes
their most frequent category. The cells leave gaps of width `ATOL` between them, and a probe point
falling in one is located in no cell; areas are computed before this, on located points only.
"""
function fill_gaps(category::Vector{Int}, n::Int)
    grid = reshape(category, n, n)
    filled = copy(grid)
    for j in 2:(n - 1), i in 2:(n - 1)
        grid[i, j] == 0 || continue
        neighbours = [grid[i + a, j + b] for a in -1:1 for b in -1:1 if (a, b) != (0, 0)]
        drawn = filter(!=(0), neighbours)
        length(drawn) >= 5 || continue
        filled[i, j] = argmax(k -> count(==(k), drawn), unique(drawn))
    end
    return vec(filled)
end

"""
    category_panel(n, points, category, colours; title)

The probe points drawn by category on the `n × n` grid: category `0` is left blank and category `k`
takes `colours[k]`, isolated gaps being filled as in [`fill_gaps`](@ref). The regions are outlined
and the window is cropped to the points drawn.
"""
function category_panel(n, points, category::Vector{Int}, colours; title)
    ticks = range(-PROBE_BOX, PROBE_BOX; length = n)
    category = fill_gaps(category, n)
    Z = permutedims(reshape(Float64.(category), n, n))
    Z[Z .== 0] .= NaN
    drawn = points[category .!= 0]
    lo = minimum(x -> minimum(x), drawn) - 0.3
    hi = maximum(x -> maximum(x), drawn) + 0.3
    k = length(colours)
    p = Plots.heatmap(
        ticks,
        ticks,
        Z;
        color = Plots.cgrad(colours, k; categorical = true),
        clims = (0.5, k + 0.5),
        colorbar = false,
        aspect_ratio = :equal,
        xlims = (lo, hi),
        ylims = (lo, hi),
        framestyle = :box,
        legend = false,
        title,
        titlefontsize = 10,
        tickfontsize = 7,
    )
    for (_, R) in REGIONS
        Plots.plot!(p, outline(R)...; color = :black, linewidth = 1.4)
    end
    return p
end
