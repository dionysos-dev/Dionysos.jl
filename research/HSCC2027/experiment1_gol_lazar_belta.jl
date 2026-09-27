# EXPERIMENT 1 — the Gol-Lazar-Belta benchmark, both approaches, both quantifiers.
#
#     julia --project=. experiment1_gol_lazar_belta.jl     [SAMPLES=n] [FIGURES=0] [CACHE=0]
#
# The system, the three observation regions, the co-safe LTL formula and the initial point are those
# of Example 3.1 of Gol, Ding, Lazar and Belta, "Finite bisimulations for switched linear systems",
# IEEE TAC 59(12):3122-3134, 2014. The Lyapunov certificate is NOT theirs: theirs is a common
# polyhedral function built by hand, certifying 0.94; the path-complete framework searches instead,
# and returns a markedly tighter one on the same system.
#
# Both approaches build a quotient and then answer the same specification under both quantifiers:
# synthesis, where the modes belong to the controller, and verification, where they belong to the
# environment. The two quotients differ in cost, and the timed lines do not all favour the same arm,
# which is the point of the experiment.
#
# CACHE=0 rebuilds the quotients instead of loading them from `cache/`, and is what the reported
# build time requires.

using StaticArrays,
    LinearAlgebra, JuMP, HiGHS, LazySets, Plots, Printf, Spot, JLD2, LaTeXStrings
import HybridSystems, Statistics
import MathOptInterface as MOI

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const PCQ = DI.Optim.Abstraction.PCLFBisimulationQuotient
const PCLF = UT.PathCompleteFramework
const OPDS = DI.Optim.DiscreteSystems

gr()

const SAMPLES = parse(Int, get(ENV, "SAMPLES", "1"))
const FIGURES = get(ENV, "FIGURES", "1") == "1"
const USE_CACHE = get(ENV, "CACHE", "1") == "1"
const ATOL = 1e-3
const SCALE = 6.0
const RAYS_DEG = [0.0, 32.6, 65.8, 91.5, 112.6, 131.4, 152.0, 180.0]
const WON, LOST = :green, :red
# All three in black. These outlines sit on saturated red and green, and darkorange vanished into
# the red while navy vanished into the dark green; the regions are told apart by where they are.
const REGION_COLOURS = [:black, :black, :black]

# ── the problem ──────────────────────────────────────────────────────────────────────────────────
A1 = @SMatrix [-0.65 0.32; -0.42 -0.92]
A2 = @SMatrix [0.65 0.32; -0.42 -0.92]
# The modes are declared the CONTROLLER's. `discreteswitchedsystem` would otherwise make them
# autonomous, i.e. the environment's, which silently turns every synthesis into a verification.
f = ST.with_switching(
    HybridSystems.discreteswitchedsystem([Matrix(A1), Matrix(A2)]),
    HybridSystems.ControlledSwitching(),
)
# The same system with the modes handed to the environment: the unchanged solver then answers ∀.
f_forall = ST.with_switching(
    HybridSystems.discreteswitchedsystem([Matrix(A1), Matrix(A2)]),
    HybridSystems.AutonomousSwitching(),
)
X = LazySets.HPolytope([
    LazySets.HalfSpace([1.0, 0.0], 5.9),
    LazySets.HalfSpace([-1.0, 0.0], 5.9),
    LazySets.HalfSpace([0.0, 1.0], 5.9),
    LazySets.HalfSpace([0.0, -1.0], 5.9),
])
R1 = LazySets.HPolytope(
    [-0.9869 -0.1615; -0.0931 0.9957; 0.9659 0.2587; 0.0825 -0.9966],
    [6.6767, 9.2315, 2.3700, -5.9038],
)
R2 = LazySets.HPolytope(
    [0.9993 0.0363; -0.7743 -0.6329; 0.5463 0.8376],
    [-2.1809, 6.3754, -4.8983],
)
R3 = LazySets.HPolytope(
    [-0.9946 -0.1041; 0.5277 0.8494; 0.9999 0.0146; -0.1191 -0.9929],
    [-5.5771, 5.3510, 9.1600, 6.2406],
)
const REGIONS = [("R1", R1), ("R2", R2), ("R3", R3)]
problem = PR.BisimulationQuotientProblem(f, X, [R1, R2, R3])
φ = ltl"((!R2 U D) & F(R1) & ((R3 -> X(!R1)) U D))"
x0 = SVector(-4.0, -7.0)          # point `a` of the paper
ap_to_obs = Dict(:D => -1, :R1 => 1, :R2 => 2, :R3 => 3)

# ── the certificate ──────────────────────────────────────────────────────────────────────────────
println("="^84, "\nEXPERIMENT 1 — Gol-Lazar-Belta Example 3.1\n", "="^84)

cones = let a = deg2rad.(RAYS_DEG)
    [hcat([cos(a[i]), sin(a[i])], [cos(a[i + 1]), sin(a[i + 1])]) for i in 1:(length(a) - 1)]
end
graph = PCLF.generate_DeBruijn_edges(2, 1)
nodes = sort(collect(graph.verts); by = string)
pclf = PCLF.compute_polyhedral_pieces_pclf(
    f,
    graph,
    HiGHS.Optimizer,
    Dict(nd => cones for nd in nodes);
    MLF = true,
)
isinf(pclf.JSRapprox) && error("no certificate for this graph and template")

"The smallest level at which every node's piece contains every observation region."
region_level(p) = maximum(
    PCQ.gamma_cover_set(p.pieces[nd], UT._as_hpolytope(R)) for nd in nodes for
    (_, R) in REGIONS
)

# One common factor, so the level lands near 6. The certificate LP is free to return any scaling and
# here returns one near 0.013, where `all_nodes_clear_regions`'s absolute tolerance of 1e-2 is the
# size of the quantity it perturbs and the test answers "clear" for a terminal set that meets a
# region. Rescaling is not free, since it changes the row norms and hence the geometric thickness of
# `atol`, so it is used here and deliberately not in experiment 2.
let α = region_level(pclf) / SCALE
    for nd in keys(pclf.pieces)
        pclf.pieces[nd].G = pclf.pieces[nd].G / α
    end
end
common = PCLF.build_common_lyapunov(pclf)
γ = pclf.JSRapprox
ΓX = region_level(pclf)

# The construction's own stopping rule: the first level whose terminal set no longer meets a region.
# A terminal set that still meets one labels every point of the overlap with the terminal
# observation instead of its own.
regions_h = [UT._as_hpolytope(R) for (_, R) in REGIONS]
NB = something(
    findfirst(
        j ->
            PCQ.all_nodes_clear_regions(pclf, ΓX * γ^(j - 1), regions_h; tol = 1e-2) &&
            PCQ.all_nodes_clear_regions(common, ΓX * γ^(j - 1), regions_h; tol = 1e-2),
        1:80,
    ),
    -1,
)
NB > 0 || error("no region-free terminal level within 80 rungs")

common_parts = let S = PCLF.get_sublevel_set(common.pieces[:clf], 1.0; atol = 1e-6)
    S isa LazySets.UnionSetArray ? length(S.array) : 1
end

@printf(
    "\ngraph            primal De Bruijn order 1, %d nodes, complete: %s\n",
    length(nodes),
    PCLF.is_complete(graph, 1:2)
)
@printf(
    "template         %d placed cones, %d rows per piece\n",
    length(cones),
    size(pclf.pieces[nodes[1]].G, 1)
)
@printf("CERTIFIED RATE   %.6f\n", γ)
@printf(
    "induced common   %d disjoint polytopes, the union approach 1 works on\n",
    common_parts
)
@printf("ladder           ΓX = %.4f, %d slices, τD = %.4f\n", ΓX, NB, ΓX * γ^(NB - 1))

# The certificate itself, which every number below depends on. A polyhedral piece is the gauge
# `V_s(x) = max_i |(G_s x)_i| / w_i`, so its Γ-sublevel set is the symmetric polytope
# `{x : |G_s x| ≤ Γ w}`.
println("\nthe certificate, one polyhedral piece per node:")
for nd in nodes
    p = pclf.pieces[nd]
    @printf("  node %-6s  V(x) = max_i |(G x)_i| / w_i\n", string(nd))
    for i in 1:size(p.G, 1)
        @printf(
            "      G[%d,:] = %9.4f %9.4f      w[%d] = %8.4f\n",
            i,
            p.G[i, 1],
            p.G[i, 2],
            i,
            p.w[i]
        )
    end
end

println("\nevery node's piece must contain every region, or ∀ is not well posed:")
for nd in nodes, (nm, R) in REGIONS
    S = PCLF.get_sublevel_set(pclf.pieces[nd], ΓX)
    ok = try
        UT.is_included(UT._as_hpolytope(R), S)
    catch
        missing
    end
    @printf("  node %-6s %-3s : %s\n", string(nd), nm, ok === true ? "ok" : "NOT CONTAINED")
end

# ── the two arms ─────────────────────────────────────────────────────────────────────────────────
"Build the quotient from `cert`, on the ladder the certificate fixed above."
function build(cert)
    opt = MOI.instantiate(PCQ.OptimizerBisimulationQuotient)
    MOI.set(opt, MOI.RawOptimizerAttribute("bisimulation_quotient_problem"), problem)
    MOI.set(opt, MOI.RawOptimizerAttribute("pclf"), cert)
    MOI.set(opt, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.set(opt, MOI.RawOptimizerAttribute("atol"), ATOL)
    MOI.set(opt, MOI.RawOptimizerAttribute("nb_levels"), NB)
    MOI.set(opt, MOI.RawOptimizerAttribute("max_slices"), NB)
    MOI.set(opt, MOI.RawOptimizerAttribute("ΓX"), ΓX)
    MOI.optimize!(opt)
    return (;
        quotient = MOI.get(opt, MOI.RawOptimizerAttribute("bisimulation_quotient")),
        D = MOI.get(opt, MOI.RawOptimizerAttribute("D")),
    )
end

"Solve φ on `built`. `system` decides the game: modes owned by the controller give synthesis, modes
owned by the environment give verification."
function solve(system, built)
    regions = Dict{Symbol, Any}(:D => built.D, :R1 => R1, :R2 => R2, :R3 => R3)
    p = PR.CoSafeLTLProblem(
        system,
        LazySets.Hyperrectangle(; low = [x0[1], x0[2]], high = [x0[1], x0[2]]),
        Dionysos.spot_stepper(φ),
        regions,
        Dict{Symbol, Any}(ap => DI.Mapping.INNER for ap in keys(regions)),
    )
    opt = MOI.instantiate(PCQ.OptimizerCoSafeLTLOnQuotient)
    MOI.set(opt, MOI.RawOptimizerAttribute("concrete_problem"), p)
    MOI.set(opt, MOI.RawOptimizerAttribute("bisimulation_quotient"), built.quotient)
    MOI.set(opt, MOI.RawOptimizerAttribute("ap_to_obs"), ap_to_obs)
    MOI.set(opt, MOI.RawOptimizerAttribute("print_level"), 0)
    # Defaults to `true`, which stops the fixed point once the initial state is decided and returns
    # a certified set of a handful of cells instead of the whole one.
    MOI.set(opt, MOI.RawOptimizerAttribute("early_stop"), false)
    MOI.optimize!(opt)
    return (;
        opt,
        won = MOI.get(opt, MOI.RawOptimizerAttribute("controllable_set")),
        lost = MOI.get(opt, MOI.RawOptimizerAttribute("uncontrollable_set")),
    )
end

function timed(g, n)
    g()
    ts = Float64[]
    for _ in 1:n
        GC.gc()
        t0 = time_ns()
        g()
        push!(ts, (time_ns() - t0) / 1e9)
    end
    return Statistics.median(ts)
end

function stats(q)
    parts, faces = PCQ.cell_complexities(q)
    return (;
        cells = length(parts),
        parts,
        faces,
        sum_p = sum(parts),
        mean_p = Statistics.mean(parts),
        max_p = maximum(parts),
        sum_f = sum(faces),
        mean_f = Statistics.mean(faces),
        max_f = maximum(faces),
    )
end

function run(name, title, cert)
    println("\n", "-"^84, "\n", title, "\n", "-"^84)
    path = joinpath(@__DIR__, "cache", "exp1_$(name)_$(NB).jld2")
    built = if USE_CACHE && isfile(path)
        println("  (quotient loaded from cache; run with CACHE=0 to time the build)")
        JLD2.jldopen(path, "r") do file
            return (; quotient = file["quotient"], D = file["D"])
        end
    else
        b = build(cert)
        if USE_CACHE
            mkpath(dirname(path))
            JLD2.jldopen(path, "w") do file
                file["quotient"] = b.quotient
                return file["D"] = b.D
            end
        end
        b
    end

    parts = built.D isa LazySets.UnionSetArray ? built.D.array : [UT._as_hpolytope(built.D)]
    for (nm, R) in REGIONS, P in parts
        LazySets.isdisjoint(P, UT._as_hpolytope(R)) ||
            error("[$name] the terminal set D meets region $nm")
    end

    s = stats(built.quotient)
    syn, ver = solve(f, built), solve(f_forall, built)
    @printf(
        "  quotient         %d cells, %d polytopes, D in %d part(s), clears every region\n",
        s.cells,
        s.sum_p,
        length(parts)
    )
    @printf(
        "  certified        %d cells under ∃, %d under ∀\n",
        length(syn.won),
        length(ver.won)
    )
    return (; built, s, cert, syn, ver)
end

a1 = run(
    "approach1",
    "APPROACH 1 — determinise the PCLF, then build on the induced common",
    common,
)
a2 = run(
    "approach2",
    "APPROACH 2 — build on the PCLF directly, one partition per node",
    pclf,
)

# ── the two costs, measured side by side ─────────────────────────────────────────────────────────
# Timed here, after both quotients exist, rather than inside `run`. A timed region pays garbage
# collection proportional to the live heap, so an arm measured while the other does not yet exist is
# measured on a lighter heap: doing it inside `run` favours whichever arm runs first by a factor
# large enough to reverse a verdict. Here both arms pay the same memory.
#
# The quotient's `Σ polytopes` is the other half of the story, and it needs no timing: abstracting a
# concrete set tests it against every cell, a cell is a semilinear set, so the sweep costs one
# feasibility LP per polytope. Solving, by contrast, touches no set at all.
t_build =
    USE_CACHE ? (NaN, NaN) :
    (timed(() -> build(a1.cert), SAMPLES), timed(() -> build(a2.cert), SAMPLES))
t_syn = (timed(() -> solve(f, a1.built), SAMPLES), timed(() -> solve(f, a2.built), SAMPLES))
t_ver = (
    timed(() -> solve(f_forall, a1.built), SAMPLES),
    timed(() -> solve(f_forall, a2.built), SAMPLES),
)

println("\n", "="^84, "\nRESULT — ratio above 1 favours approach 2\n", "="^84)
@printf("%-26s %13s %13s %10s\n", "", "approach 1", "approach 2", "1 / 2")
println("-"^84)
for (nm, u, v) in (
    ("cells", float(a1.s.cells), float(a2.s.cells)),
    ("Σ polytopes", float(a1.s.sum_p), float(a2.s.sum_p)),
    ("mean polytopes / cell", a1.s.mean_p, a2.s.mean_p),
    ("MAX polytopes / cell", float(a1.s.max_p), float(a2.s.max_p)),
    ("mean facets / cell", a1.s.mean_f, a2.s.mean_f),
    ("MAX facets / cell", float(a1.s.max_f), float(a2.s.max_f)),
    ("build (s)", t_build[1], t_build[2]),
    ("solve ∃, a graph fixed point (s)", t_syn[1], t_syn[2]),
    ("solve ∀, a graph fixed point (s)", t_ver[1], t_ver[2]),
)
    if isnan(u)
        @printf("%-26s %13s %13s %10s\n", nm, "cached", "cached", "—")
    else
        @printf("%-26s %13.2f %13.2f %10.2f\n", nm, u, v, u / v)
    end
end
println("="^84)

# ── figures ──────────────────────────────────────────────────────────────────────────────────────
if FIGURES
    d = joinpath(@__DIR__, "figures")
    mkpath(d)
    out(fig, n) = (savefig(fig, joinpath(d, n)); println("wrote ", n))

    # `linealpha` is set explicitly on every overlay in this file. Plotting the quotient with
    # `show_contours = true` leaves a subplot-level line alpha of 0.5 behind, which every later
    # series inherits and which makes these outlines vanish entirely.
    function outline!(p)
        for (i, (_, R)) in enumerate(REGIONS)
            plot!(
                p,
                R;
                fillalpha = 0.0,
                linecolor = REGION_COLOURS[i],
                linealpha = 1.0,
                linewidth = 2.5,
                label = "",
            )
        end
        return p
    end
    function align!(ps)
        xs = reduce(vcat, [collect(Plots.xlims(p)) for p in ps])
        ys = reduce(vcat, [collect(Plots.ylims(p)) for p in ps])
        l = (extrema(xs), extrema(ys))
        for p in ps
            xlims!(p, l[1]...)
            ylims!(p, l[2]...)
        end
        return l
    end
    # Under equal aspect Plots derives the height from the width and the data, so a height fixed
    # independently of the data over-constrains the layout and Plots aborts.
    function rowsize(l, n; w = 430, chrome = 95)
        (xlo, xhi), (ylo, yhi) = l
        return (
            w * n,
            round(Int, w * clamp((yhi - ylo) / max(xhi - xlo, eps()), 0.35, 2.2)) + chrome,
        )
    end
    function certified(q, res, title)
        p = plot(; aspect_ratio = :equal, legend = false, title = title, titlefontsize = 9)
        for (ids, c) in ((res.lost, LOST), (res.won, WON))
            isempty(ids) && continue
            plot!(
                p,
                q;
                what = :states,
                state_ids = collect(ids),
                show_contours = false,
                user_color = c,
                fillalpha = 1.0,
            )
        end
        return outline!(p)
    end

    # (1) the induced common's quotient, its cells
    p = plot(;
        aspect_ratio = :equal,
        legend = false,
        titlefontsize = 9,
        title = "Approach 1 — the induced common's quotient, $(a1.s.cells) cells",
    )
    plot!(
        p,
        a1.built.quotient;
        what = :states,
        by = :state,
        show_contours = true,
        linewidth = 0.3,
        fillalpha = 0.9,
        merge_series = false,
    )
    outline!(p)
    l = align!([p])
    plot!(p; size = rowsize(l, 1))
    out(p, "exp1_approach1_quotient.png")

    # ── the witness runs ─────────────────────────────────────────────────────────────────────────
    # The certified sets say which points satisfy φ; they do not show what the two quantifiers mean.
    # From one point certified under both, the controller's single run and the environment's whole
    # tree of runs do, and ∀ certifies a subset of ∃, so a ∀-certified point serves both.
    #
    # Labels are read off the concrete regions rather than off a cell, so the picture is about the
    # system and not about either quotient.
    spec = Dionysos.spot_stepper(φ)
    acc = Set(OPDS.accepting_states(spec))
    A = UT.mode_matrices(f)

    function labels_at(z, D)
        ls = Symbol[]
        z ∈ D && push!(ls, :D)
        for (nm, R) in REGIONS
            z ∈ R && push!(ls, Symbol(nm))
        end
        return Tuple(ls)
    end

    "Every word the environment can play from `z`, each branch cut where the good prefix completes."
    function environment_tree(D, z; nmax)
        branches = NamedTuple{(:X, :word), Tuple{Vector{Vector{Float64}}, Vector{Int}}}[]
        function walk(path, word, qa, k)
            if qa in acc
                push!(branches, (; X = copy(path), word = copy(word)))
                return true
            end
            k == nmax && return false
            for m in eachindex(A)
                yn = A[m] * path[end]
                push!(path, yn)
                push!(word, m)
                ok = walk(path, word, OPDS.step(spec, qa, labels_at(yn, D)), k + 1)
                pop!(path)
                pop!(word)
                ok || return false
            end
            return true
        end
        z = collect(Float64, z)
        closed =
            walk([z], Int[], OPDS.step(spec, OPDS.init_state(spec), labels_at(z, D)), 0)
        depth = isempty(branches) ? 0 : maximum(length(b.word) for b in branches)
        return (; branches, closed, depth, leaves = length(branches))
    end

    "An interior point of a quotient cell; for a union, of its first convex part."
    function cell_point(S)
        V = LazySets.vertices_list(
            UT._as_hpolytope(S isa LazySets.UnionSetArray ? S.array[1] : S),
        )
        return sum(V) / length(V)
    end

    "The terminal set, the formula's target, outlined on a panel."
    function target!(p, D)
        for P in (D isa LazySets.UnionSetArray ? D.array : [D])
            plot!(
                p,
                P;
                fillalpha = 0.0,
                linecolor = :black,
                linestyle = :dash,
                linealpha = 1.0,
                linewidth = 1.5,
                label = "",
            )
        end
        return p
    end

    # The paper's point `a` is preferred, and a deeper tree over the ∀-certified cells is taken when
    # it is shallow: a tree of one or two levels shows nothing the certified set does not. The point
    # must also carry a controller on BOTH quotients, since the same point is drawn on both figures.
    const TREE_MAX = 7
    ctrl1 = PCQ.solve_concrete_problem(a1.syn.opt)
    ctrl2 = PCQ.solve_concrete_problem(a2.syn.opt)

    "The controlled run from `z`, or `nothing` if no controller starts there."
    function controlled_run(opt, ctrl, z)
        mem = try
            PCQ.initial_controller_memory(opt, z)
        catch
            return nothing
        end
        return PCQ.simulate_closed_loop(f, ctrl, z, mem; N = 2 * TREE_MAX)
    end

    "Whether the concrete run `X` completes a good prefix of φ."
    function accepts(X, D)
        qa = OPDS.step(spec, OPDS.init_state(spec), labels_at(X[1], D))
        for z in X[2:end]
            qa = OPDS.step(spec, qa, labels_at(z, D))
        end
        return qa in acc
    end

    witness, sim1, sim2 = let q = a1.built.quotient, D = a1.built.D
        cands = [(x = collect(Float64, x0), from = "the paper's point a")]
        for qid in a1.ver.won
            push!(cands, (x = cell_point(q.states[qid].set), from = "a ∀-certified cell"))
        end
        # Only the first candidate needs the membership test; the others are interior to a cell the
        # ∀ fixed point returned.
        ranked = NamedTuple[]
        for (i, c) in enumerate(cands)
            i == 1 && !any(c.x ∈ q.states[qid].set for qid in a1.ver.won) && continue
            t = environment_tree(D, c.x; nmax = TREE_MAX)
            t.closed && push!(ranked, (; c..., t))
        end
        isempty(ranked) &&
            error("no ∀-certified point whose tree closes within $TREE_MAX steps")
        sort!(
            ranked;
            by = w -> (w.from == "the paper's point a" && w.t.depth ≥ 3, w.t.depth),
            rev = true,
        )

        pick = nothing
        for w in ranked
            z = SVector{2}(w.x)
            s1 = controlled_run(a1.syn.opt, ctrl1, z)
            s1 === nothing && continue
            s2 = controlled_run(a2.syn.opt, ctrl2, z)
            s2 === nothing && continue
            pick = (w, s1, s2)
            break
        end
        pick === nothing &&
            error("no witness point carries a controller on both quotients")
        pick
    end
    xw = SVector{2}(witness.x)

    # `simulate_closed_loop` stops when the controller is no longer defined. For a satisfied good
    # prefix that is acceptance, but leaving the domain looks the same from outside, and only the
    # first is a witness.
    for (nm, s, D) in (("approach 1", sim1, a1.built.D), ("approach 2", sim2, a2.built.D))
        accepts(s.X, D) || error(
            "[$nm] the controlled run does not satisfy φ; it left the controller's domain",
        )
    end
    @printf("\nwitness point    (%.3f, %.3f), %s\n", xw[1], xw[2], witness.from)
    @printf(
        "  synthesis ∃    approach 1: %d steps, approach 2: %d steps\n",
        length(sim1.U),
        length(sim2.U)
    )
    # The same runs appear on both figures: they are runs of the system, not of a quotient.
    @printf(
        "  verification ∀ %d runs, deepest %d steps, identical on both\n",
        witness.t.leaves,
        witness.t.depth
    )

    # Only two of the environment's runs are drawn. The whole tree was a tangle in which no single
    # run could be followed, and following one is the point; the shortest and the longest bracket
    # the spread, and the title carries the count.
    const RUN_COLOURS = [:black, :blue]
    word_label(w) = isempty(w) ? "ε" : join(w, "")
    function run!(p, X, col, lab)
        xs, ys = first.(X), last.(X)
        plot!(
            p,
            xs,
            ys;
            color = col,
            linealpha = 1.0,
            linewidth = 2.2,
            marker = :circle,
            markersize = 3.5,
            markerstrokewidth = 0,
            label = lab,
        )
        return p
    end

    shown = let bs = sort(witness.t.branches; by = b -> length(b.word))
        length(bs) ≤ 2 ? bs : [bs[1], bs[end]]
    end

    # (0) the problem as the paper states it: a working pair X₁ ⊆ X₂, pairwise disjoint regions of
    # interest, and a terminal pair D₁ ⊆ D₂ with D₂ disjoint from every region. Each pair comes from
    # one level of the PCLF pieces, the intersection giving the inner set and the union the outer.
    # Neither D is invariant; what makes the pair terminal is that a run reaching D₁ stays in D₂
    # forever, and D₂ meets no region, so no region is ever observed again.
    #
    # R1 away from black, which this figure spends on D₁ and D₂.
    const PROBLEM_COLOURS = [:darkgreen, :navy, :darkorange]
    τD = ΓX * γ^(NB - 1)
    D_pieces =
        [UT._as_hpolytope(PCLF.get_sublevel_set(pclf.pieces[nd], τD)) for nd in nodes]
    D1 = reduce(LazySets.intersection, D_pieces)
    # The working pair, one level up: X₂ is the outer domain the abstraction is built on, X₁ the
    # inner set verification and synthesis are restricted to.
    X_pieces =
        [UT._as_hpolytope(PCLF.get_sublevel_set(pclf.pieces[nd], ΓX)) for nd in nodes]
    X1 = reduce(LazySets.intersection, X_pieces)

    p = plot(;
        aspect_ratio = :equal,
        legend = false,
        titlefontsize = 9,
        title = "Gol-Lazar-Belta Example 3.1 — the problem",
    )
    for P in X_pieces
        plot!(
            p,
            P;
            fillalpha = 1.0,
            fillcolor = :grey93,
            linecolor = :grey70,
            linealpha = 1.0,
            linewidth = 1.0,
            label = "",
        )
    end
    plot!(
        p,
        X1;
        fillalpha = 0.0,
        linecolor = :black,
        linealpha = 1.0,
        linewidth = 2.0,
        label = "",
    )
    # R₁ and R₂ are nudged off their centroids: the run crosses both, and a label sitting on the
    # line is unreadable. A nudge that would leave its region is dropped.
    LABEL_NUDGE = Dict("R1" => [0.9, 0.0], "R2" => [1.0, -0.9])
    for (i, (nm, R)) in enumerate(REGIONS)
        plot!(
            p,
            R;
            fillalpha = 0.30,
            fillcolor = PROBLEM_COLOURS[i],
            linealpha = 1.0,
            linecolor = PROBLEM_COLOURS[i],
            linewidth = 2.5,
            label = "",
        )
        c = cell_point(R)
        c_lab = c + get(LABEL_NUDGE, nm, [0.0, 0.0])
        c_lab ∈ UT._as_hpolytope(R) || (c_lab = c)
        annotate!(
            p,
            c_lab[1],
            c_lab[2],
            Plots.text(latexstring("R_", i), 12, PROBLEM_COLOURS[i], :center),
        )
    end
    # D₂ filled, then D₁ opaque on top: the pieces are nested and differ by 4% in area, so the pair
    # is only legible as the ring between them that this leaves.
    for P in D_pieces
        plot!(
            p,
            P;
            fillalpha = 0.40,
            fillcolor = :red,
            linecolor = :red,
            linealpha = 1.0,
            linewidth = 1.2,
            label = "",
        )
    end
    plot!(
        p,
        D1;
        fillalpha = 1.0,
        fillcolor = :grey78,
        linecolor = :black,
        linealpha = 1.0,
        linewidth = 1.5,
        label = "",
    )
    run!(p, sim2.X, :black, "")
    scatter!(
        p,
        [xw[1]],
        [xw[2]];
        color = :white,
        markerstrokecolor = :black,
        markerstrokewidth = 1.5,
        markersize = 6,
        label = "",
    )
    let c = cell_point(D1)
        annotate!(p, c[1], c[2], Plots.text(L"\mathcal{D}_1", 12, :black, :center))
    end
    # D₂ is named like X₂: just outside a vertex that lies strictly outside D₁, the most top-right
    # of them, so the label points at the ring rather than at where the two boundaries touch.
    let V = reduce(vcat, LazySets.vertices_list.(D_pieces)),
        cs_D = LazySets.constraints_list(D1)

        outside_D(w) = maximum(dot(c.a, w) - c.b for c in cs_D)
        far = filter(w -> outside_D(w) > 1e-8, V)
        cand = isempty(far) ? V : far
        v = cand[argmax([w[1] + w[2] for w in cand])]
        annotate!(
            p,
            1.10 * v[1] + 0.1,
            1.10 * v[2] + 0.4,
            Plots.text(L"\mathcal{D}_2", 12, :red, :center),
        )
    end
    # The two working sets are named on opposite sides so the labels cannot be confused: X₁ just
    # inside its dark contour at the bottom, X₂ just outside the outer boundary at the top.
    let V = LazySets.vertices_list(X1)
        v = V[argmin([w[2] for w in V])]
        annotate!(
            p,
            0.88 * v[1],
            0.88 * v[2],
            Plots.text(L"\mathcal{X}_1", 12, :black, :center),
        )
    end
    # X₂ is named at a vertex that lies strictly outside X₁, so the label sits where the two sets
    # actually differ rather than where their boundaries touch.
    let V = reduce(vcat, LazySets.vertices_list.(X_pieces)),
        cs = LazySets.constraints_list(X1)

        outside(w) = maximum(dot(c.a, w) - c.b for c in cs)
        far = filter(w -> outside(w) > 1e-8, V)
        cand = isempty(far) ? V : far
        v = cand[argmax([w[1] + w[2] for w in cand])]
        annotate!(
            p,
            1.07 * v[1] + 0.5,
            1.07 * v[2],
            Plots.text(L"\mathcal{X}_2", 12, :black, :center),
        )
    end
    # To the right of the point, not below it: below, the label fell across the X₁ contour.
    annotate!(p, xw[1] + 1.0, xw[2] - 0.5, Plots.text(L"x_0", 12, :black, :center))
    l = align!([p])
    plot!(p; size = rowsize(l, 1))
    out(p, "exp1_problem.png")

    # (2) approach 1, both games on one figure, each carrying the runs from the witness point: the
    # controller needs one of them to reach the target, the environment must be answered on all.
    ps = [
        certified(
            a1.built.quotient,
            a1.syn,
            "Approach 1 — synthesis ∃: the controlled run",
        ),
        certified(
            a1.built.quotient,
            a1.ver,
            "Approach 1 — verification ∀: $(witness.t.leaves) runs, $(length(shown)) drawn",
        ),
    ]
    run!(ps[1], sim1.X, RUN_COLOURS[1], "word " * word_label(Int.(sim1.U)))
    for (i, b) in enumerate(shown)
        run!(ps[2], b.X, RUN_COLOURS[i], "word " * word_label(b.word))
    end
    # The regions go back on top: the runs cross them, and a boundary the eye cannot follow
    # stops being a reference.
    for p in ps
        outline!(p)
        target!(p, a1.built.D)
        scatter!(
            p,
            [xw[1]],
            [xw[2]];
            color = :white,
            markerstrokecolor = :black,
            markerstrokewidth = 1.5,
            markersize = 6,
            label = "",
        )
        plot!(
            p;
            legend = :bottomright,
            legendfontsize = 7,
            background_color_legend = :white,
        )
    end
    l = align!(ps)
    out(plot(ps...; layout = (1, 2), size = rowsize(l, 2)), "exp1_approach1_spec.png")

    # (5) approach 2, the same two games from the same point. Every layer's cells are drawn, so a
    # state certified on one layer and not on the other appears in both colours; the 3-D view
    # disambiguates it.
    #
    # The environment's runs are identical to approach 1's, because they are runs of the system and
    # not of a quotient, which is the claim: the cheaper abstraction answers the same question. The
    # controlled run need not be, since the two quotients admit different strategies.
    ps = [
        certified(
            a2.built.quotient,
            a2.syn,
            "Approach 2 — synthesis ∃: the controlled run",
        ),
        certified(
            a2.built.quotient,
            a2.ver,
            "Approach 2 — verification ∀: $(witness.t.leaves) runs, $(length(shown)) drawn",
        ),
    ]
    run!(ps[1], sim2.X, RUN_COLOURS[1], "word " * word_label(Int.(sim2.U)))
    for (i, b) in enumerate(shown)
        run!(ps[2], b.X, RUN_COLOURS[i], "word " * word_label(b.word))
    end
    for p in ps
        outline!(p)
        target!(p, a2.built.D)
        scatter!(
            p,
            [xw[1]],
            [xw[2]];
            color = :white,
            markerstrokecolor = :black,
            markerstrokewidth = 1.5,
            markersize = 6,
            label = "",
        )
        plot!(
            p;
            legend = :bottomright,
            legendfontsize = 7,
            background_color_legend = :white,
        )
    end
    l = align!(ps)
    out(plot(ps...; layout = (1, 2), size = rowsize(l, 2)), "exp1_approach2_spec.png")

    # (6) facets per cell. Log counts: most cells are simple, the tail is not, and the tail is what
    # the construction pays for. `fillrange` pins both series to the same baseline, which a log axis
    # otherwise picks per series, leaving the second suspended at 10^0.
    bins = range(0, 1.02 * max(maximum(a1.s.faces), maximum(a2.s.faces)); length = 45)
    h = histogram(
        a1.s.faces;
        bins,
        yscale = :log10,
        fillrange = 0.7,
        alpha = 0.55,
        color = :firebrick,
        linecolor = :firebrick,
        label = "approach 1 (determinise first)",
        xlabel = "facets per cell",
        ylabel = "cells",
        legend = :topright,
        framestyle = :box,
        grid = false,
        ylims = (0.7, 5e4),
        size = (660, 410),
    )
    histogram!(
        h,
        a2.s.faces;
        bins,
        yscale = :log10,
        fillrange = 0.7,
        alpha = 0.55,
        color = :steelblue,
        linecolor = :steelblue,
        label = "approach 2 (on the PCLF)",
    )
    vline!(h, [a1.s.max_f]; color = :firebrick, ls = :dash, lw = 1.5, label = "")
    vline!(h, [a2.s.max_f]; color = :steelblue, ls = :dash, lw = 1.5, label = "")
    out(h, "exp1_facets_histogram.png")

    # ── 3-D. Everything Plots-based must come first: both packages export `plot`, and the import
    # is deliberate, since `using CairoMakie` makes every bare `plot` in this session ambiguous,
    # which breaks redrawing a figure from a REPL that has the script loaded.
    import CairoMakie
    mk_out(mk, n) =
        (CairoMakie.save(joinpath(d, n), mk; px_per_unit = 3); println("wrote ", n))
    node_z = Dict(nd => z for (nd, z) in zip(nodes, (0.0, 1.0)))

    function layered(title; ids = nothing, traj = nothing)
        mk = CairoMakie.Figure(; size = (900, 700))
        ax = CairoMakie.Axis3(
            mk[1, 1];
            xlabel = "x₁",
            ylabel = "x₂",
            zlabel = "memory (graph node)",
            zticks = ([0.0, 1.0], ["node $(nodes[1])", "node $(nodes[2])"]),
            azimuth = 1.2π,
            elevation = 0.16π,
            title = title,
        )
        q = a2.built.quotient
        if ids === nothing
            DI.plot_augmented_bisimulation!(
                ax,
                q;
                node_z = node_z,
                color_by = :state,
                alpha = 0.35,
                show_contours = false,
            )
        else
            DI.plot_augmented_bisimulation!(
                ax,
                q;
                state_ids = ids.lost,
                node_z = node_z,
                color_by = LOST,
                alpha = 0.30,
                show_contours = false,
            )
            DI.plot_augmented_bisimulation!(
                ax,
                q;
                state_ids = ids.won,
                node_z = node_z,
                color_by = WON,
                alpha = 0.45,
                show_contours = false,
            )
        end
        # The observation regions on every layer. The specification is about them, and without them
        # the layers read as two coloured discs. `vertices_list` does not promise a cyclic order, so
        # the outline is sorted by angle around the centroid; the small z offset keeps the line off
        # the cell meshes it would otherwise z-fight with.
        for (_, R) in REGIONS
            V = LazySets.vertices_list(UT._as_hpolytope(R))
            c = sum(V) / length(V)
            V = V[sortperm([atan(v[2] - c[2], v[1] - c[1]) for v in V])]
            push!(V, V[1])
            for z in values(node_z)
                CairoMakie.lines!(
                    ax,
                    first.(V),
                    last.(V),
                    fill(z + 0.005, length(V));
                    color = :black,
                    linewidth = 2.0,
                )
            end
        end
        traj === nothing ||
            DI.plot_augmented_trajectory!(ax, q, traj.X, traj.M; node_z = node_z)
        return mk
    end

    # (3) approach 2's quotient in 3-D, one layer per graph node
    mk_out(
        layered("Approach 2 — the quotient, one layer per graph node"),
        "exp1_approach2_3d_quotient.png",
    )

    # (4) the certified set in 3-D with the closed loop, as on the poster. It is `sim2`, the same
    # run approach 2's flat figure draws from the same witness point, so the three figures can be
    # read as one story rather than three unrelated starts.
    mk_out(
        layered(
            "Approach 2 — certified set (∃) and the closed loop";
            ids = (; won = a2.syn.won, lost = a2.syn.lost),
            traj = (; sim2.X, sim2.M),
        ),
        "exp1_approach2_3d_certified.png",
    )
end
