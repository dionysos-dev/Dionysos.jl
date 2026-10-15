# Does the abstract transition relation really contain the concrete one?
#
# An independent check of the simulation property that does not trust the algorithm that produced
# the quotient: sample interior points of an abstract cell, step them forward under the concrete
# dynamics, and verify that *some* recorded successor of that cell contains the image.
#
# Three subtleties, each of which produced a false alarm before being handled:
#
#  - cells are closed, so a point on a shared boundary lies in several of them. Testing "the cell I
#    happened to find first" against the recorded successors reports violations that are not real.
#    Test against ALL recorded successors for that mode instead.
#  - samples must be strictly interior. Points are pulled toward the cell centroid for that reason.
#  - a point whose image LEAVES the working set is owed no successor at all. Failures must be
#    classified before they are counted.
#
# Measured result: on 1 443 samples of the observer-graph quotient, every apparent violation is an
# escape from the working set, and **zero** are genuine. The Gol-Lazar-Belta quotients are clean on
# 1 902 and 5 022 samples. What still escapes this check is the converse direction -- that no
# recorded transition is spurious -- which is its mirror image and worth adding.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))

using Random
using LinearAlgebra

# How far outside a polytope does `y` fall?  Normalised constraint violation, 0 when inside.
function violation(P::LazySets.HPolytope, y)
    return maximum(
        c -> (LinearAlgebra.dot(c.a, y) - c.b) / LinearAlgebra.norm(c.a),
        LazySets.constraints_list(P),
    )
end

violation(S::UT.SemiLinearSet, y) = minimum(P -> violation(P, y), S.array)

"""
Strictly interior points of a cell: each part's centroid, plus `k` samples pulled toward it.

`pull` is how far toward the centroid a random convex combination of the vertices is dragged;
at 0 the samples sit on the boundary, which is what makes a naive version of this check report
violations that are only shared facets.
"""
function interior_points(S::UT.SemiLinearSet, k::Int, rng; pull::Float64 = 0.8)
    pts = Vector{Float64}[]
    for P in S.array
        vertices = try
            LazySets.vertices_list(P)
        catch
            continue
        end
        isempty(vertices) && continue
        centroid = sum(vertices) / length(vertices)
        push!(pts, centroid)
        for _ in 1:k
            w = rand(rng, length(vertices))
            w ./= sum(w)
            x = sum(w[i] * vertices[i] for i in eachindex(vertices))
            push!(pts, pull .* centroid .+ (1 - pull) .* x)
        end
    end
    return pts
end

"""
Audit `optimizer`'s quotient, classifying every apparent violation by its cause.

Returns `(checked, escaped, uncovered, genuine)`: images that left the working set, images inside
it but outside the covered slices, and images inside a cell of the target node that no recorded
transition reaches. Only the last is a soundness failure.
"""
function audit(optimizer; nstates::Int = 60, npts::Int = 2, seed::Int = 7)
    quotient = MOI.get(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient"))
    problem = optimizer.bisimulation_quotient_problem
    graph = optimizer.pclf.graph
    A = UT.mode_matrices(problem.system)
    X = problem.state_set

    by_node = Dict{Any, Vector{Int}}()
    for (id, state) in quotient.states
        push!(get!(by_node, state.node, Int[]), id)
    end

    rng = Random.MersenneTwister(seed)
    ids = collect(keys(quotient.states))
    Random.shuffle!(rng, ids)

    checked, escaped, uncovered, genuine = 0, 0, 0, 0
    for qid in ids[1:min(nstates, length(ids))]
        state = quotient.states[qid]
        for (src, dst, m) in graph.edges
            src == state.node || continue
            successors = [t for (mode, t) in state.next if mode == m]
            isempty(successors) && continue
            for x in interior_points(state.set, npts, rng)
                y = A[m] * x
                checked += 1
                minimum(t -> violation(quotient.states[t].set, y), successors) <= 1e-8 &&
                    continue
                if !(y in X)
                    escaped += 1
                elseif any(cid -> y in quotient.states[cid].set, get(by_node, dst, Int[]))
                    genuine += 1
                else
                    uncovered += 1
                end
            end
        end
    end
    return checked, escaped, uncovered, genuine
end

function report(name, path)
    if !isfile(path)
        println(rpad(name, 18), "cache missing -- run the experiment that writes it first")
        return nothing
    end
    checked, escaped, uncovered, genuine = audit(import_optimizer_jld2(path))
    println("\n", name)
    println("  samples checked                       : ", checked)
    println("  image leaves the working set          : ", escaped, "   (legitimate)")
    println(
        "  image inside X, outside covered slices: ",
        uncovered,
        "   (legitimate, scope)",
    )
    println(
        "  image inside X, in a target cell      : ",
        genuine,
        genuine == 0 ? "   -> SOUND on this sample" : "   *** GENUINE VIOLATION ***",
    )
    return nothing
end

println("Simulation check: does some recorded successor contain the concrete image?")
report("E2 PCLF", joinpath(dirname(dirname(@__DIR__)), "experiments", "pclf_case.jld2"))
report(
    "E8 PCLF",
    joinpath(dirname(dirname(@__DIR__)), "examples", "gol_lazar_belta_pclf.jld2"),
)
report(
    "E8 induced CLF",
    joinpath(dirname(dirname(@__DIR__)), "examples", "gol_lazar_belta_clf.jld2"),
)
