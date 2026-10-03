# N8 step 3: how does the cost of a refinement primitive actually scale with facet count?
#
# `rem:complexity` asserts that the refinement operations are "superlinear in the facet count of the
# cells involved". That assertion is the entire basis of mechanism B -- route 2 wins because its cells
# are SIMPLER, even when the totals match -- and it is currently unmeasured. This script measures it.
#
# Why it matters commercially, not just editorially. The campaign's cost story is stuck: route 2
# builds 2.18x fewer cells on N1's family, yet E2 times a tie and E8 times a 1.30x LOSS. Either the
# facet-dependent term is small (so mechanism B is real but dormant at these sizes) or it is large and
# something else eats the gain. An exponent measured on the primitives decides which, and no quotient
# has to be built to get it.
#
# Design notes, because a careless version of this measures the wrong thing:
#
#   random polytopes at a CONTROLLED facet count. Sampling normals on the circle/sphere and taking
#       the intersection of half-spaces gives exactly F facets in general position; a random V-polytope
#       would not control F, which is the independent variable.
#   overlapping, not nested or disjoint. `set_difference_decompose(P, Q)` emits one piece per
#       constraint of Q, so a Q that misses P entirely returns P and measures nothing, while a Q
#       containing P returns nothing. Both are the easy cases and both would understate the exponent.
#   the same polytopes across all timings of a given size, so the fit is not confounded by draw luck.
#   n = 2 and n = 3 -- where every benchmark in the folder lives, and where `_screen_vertices` is
#       active -- with n = 4 as the contrast, since the screen is off there.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))

using Random
using Printf
import Statistics

"""
A random bounded polytope in `n` dimensions with exactly `F` facets, in general position.

Built as an intersection of half-spaces with normals spread over the sphere and offsets at 1, so the
unit ball is contained and the result is bounded for `F >= n + 1`. Facet count is the knob being
swept, so it must be what it says: the returned polytope is cleaned and its actual facet count is
reported by the caller rather than assumed.
"""
function random_polytope(rng, n::Int, F::Int; offset::Float64 = 1.0)
    F >= 2n || error("need F >= 2n = $(2n) to stay bounded, got $F")
    cons = LazySets.HalfSpace{Float64, Vector{Float64}}[]
    # The 2n axis half-spaces first. Purely random normals can all fall in one half-space, leaving an
    # UNBOUNDED polytope -- which `clean_poly` and the emptiness LP handle on a different code path,
    # so the sweep would silently stop measuring the intended thing. Boxing first makes boundedness
    # structural rather than a property of the draw.
    for i in 1:n, s in (1.0, -1.0)
        a = zeros(n)
        a[i] = s
        push!(cons, LazySets.HalfSpace(a, offset))
    end
    # The remaining F - 2n cuts shave corners: a unit normal at distance `offset` cuts the box, whose
    # corners sit at `offset * sqrt(n)`, so these are non-redundant in general position.
    for _ in 1:(F - 2n)
        a = randn(rng, n)
        a ./= LinearAlgebra.norm(a)
        push!(cons, LazySets.HalfSpace(a, offset))
    end
    return LazySets.HPolytope(cons)
end

"""
Seconds per call of `f`, measured to nanosecond resolution and batched until the total is meaningful.

`time()` returns a `Float64` of seconds whose granularity on Windows is about a millisecond, which is
larger than several of the primitives being measured: a first version of this script reported
`clean_poly` at 0.0000 ms for F = 64 and 1.0002 ms for F = 32, i.e. pure quantization noise, and then
fitted an exponent through it. `time_ns()` has nanosecond resolution, and batching until the run
exceeds `min_total` keeps even the cheapest primitive well above the clock's granularity.
"""
function timeit(f; min_total::Float64 = 0.05, max_batches::Int = 7)
    f()
    best = Inf
    batch = 1
    for _ in 1:max_batches
        t0 = time_ns()
        for _ in 1:batch
            f()
        end
        elapsed = (time_ns() - t0) / 1e9
        best = min(best, elapsed / batch)
        elapsed >= min_total && return best
        # Grow the batch toward whatever would have taken `min_total`, capped so one slow primitive
        # cannot explode into a run that never finishes.
        batch = min(
            batch * 8,
            max(batch + 1, ceil(Int, batch * min_total / max(elapsed, 1e-9))),
        )
    end
    return best
end

"""
Least-squares exponent `p` in `t ≈ c F^p`, fitted on logs.

Reported with the r² so a fit that does not describe the data cannot be quoted as though it did.
"""
function fit_exponent(Fs, ts)
    keep = [i for i in eachindex(ts) if ts[i] > 0]
    length(keep) < 3 && return (NaN, NaN)
    x = log.(Float64.(Fs[keep]))
    y = log.(ts[keep])
    x̄, ȳ = Statistics.mean(x), Statistics.mean(y)
    p = sum((x .- x̄) .* (y .- ȳ)) / sum((x .- x̄) .^ 2)
    ŷ = ȳ .+ p .* (x .- x̄)
    ss_res = sum((y .- ŷ) .^ 2)
    ss_tot = sum((y .- ȳ) .^ 2)
    return (p, ss_tot == 0 ? NaN : 1 - ss_res / ss_tot)
end

# ---------------------------------------------------------
# Sweep 1 — set difference against facet count
# ---------------------------------------------------------

"""
Time `set_difference_decompose(P, Q)` with both operands at `F` facets.

`Q` is shifted so it genuinely overlaps `P` without containing it: the subtrahend's constraints each
spawn a piece, so overlap is what puts the primitive on its real path.
"""
function sweep_set_difference(n::Int, Fs::Vector{Int}; seed = 7)
    rng = Random.MersenneTwister(seed)
    rows = NamedTuple[]
    for F in Fs
        P = random_polytope(rng, n, F)
        Qraw = random_polytope(rng, n, F; offset = 0.8)
        shift = vcat([0.6], zeros(n - 1))
        Q = LazySets.HPolytope([
            LazySets.HalfSpace(c.a, c.b + LinearAlgebra.dot(c.a, shift)) for
            c in LazySets.constraints_list(Qraw)
        ])

        nP = length(LazySets.constraints_list(UT.clean_poly(P)))
        nQ = length(LazySets.constraints_list(UT.clean_poly(Q)))
        pieces = UT.set_difference_decompose(P, Q; atol = 1e-6)
        t = timeit(() -> UT.set_difference_decompose(P, Q; atol = 1e-6))
        push!(rows, (; F, nP, nQ, npieces = length(pieces), t))
    end
    return rows
end

# ---------------------------------------------------------
# Sweep 2 — the other primitives on the same polytopes
# ---------------------------------------------------------

function sweep_primitives(n::Int, Fs::Vector{Int}; seed = 11)
    rng = Random.MersenneTwister(seed)
    rows = NamedTuple[]
    A = Matrix{Float64}(LinearAlgebra.I, n, n)
    A[1, min(2, n)] += 0.4
    for F in Fs
        P = random_polytope(rng, n, F)
        t_clean = timeit(() -> UT.clean_poly(P))
        t_pre = timeit(() -> UT.preimage_linear(P, A))
        t_empty = timeit(() -> isempty(P))
        push!(rows, (; F, t_clean, t_pre, t_empty))
    end
    return rows
end

# ---------------------------------------------------------
# Report
# ---------------------------------------------------------

# Capped per dimension. The set difference emits one piece per constraint of the subtrahend and each
# piece is itself intersected, so cost climbs steeply in BOTH F and n: 256 facets took 15 s at n = 2,
# and the first run of this script was still grinding through n = 3 long afterwards. The exponent is
# what is wanted, and six points determine it as well as seven.
const FS_BY_DIM =
    Dict(2 => [4, 8, 16, 32, 64, 128, 256], 3 => [6, 12, 24, 48, 96], 4 => [8, 16, 32, 64])

function report(n::Int)
    println("\n", "="^78)
    println(
        "n = $n   (`_screen_vertices` is active at n <= 3, so n = 4 is the contrast case)",
    )
    println("="^78)

    println("\nset_difference_decompose(P, Q), both operands at F facets:\n")
    @printf(
        "%-8s %-8s %-8s %-10s %-14s %s\n",
        "F",
        "|P|",
        "|Q|",
        "pieces",
        "time (ms)",
        "vs F=4"
    )
    Fs = filter(>=(2n), FS_BY_DIM[n])
    rows = sweep_set_difference(n, Fs)
    base = rows[1].t
    for r in rows
        @printf(
            "%-8d %-8d %-8d %-10d %-14.3f %.1fx\n",
            r.F,
            r.nP,
            r.nQ,
            r.npieces,
            1e3 * r.t,
            r.t / base
        )
    end
    p, r2 = fit_exponent([r.F for r in rows], [r.t for r in rows])
    @printf("\n  fitted t ~ F^%.2f   (r² = %.3f)\n", p, r2)

    println("\nthe other primitives, same polytopes:\n")
    @printf("%-8s %-16s %-16s %s\n", "F", "clean (ms)", "preimage (ms)", "isempty (ms)")
    prows = sweep_primitives(n, Fs)
    for r in prows
        @printf(
            "%-8d %-16.4f %-16.4f %.4f\n",
            r.F,
            1e3 * r.t_clean,
            1e3 * r.t_pre,
            1e3 * r.t_empty
        )
    end
    for (name, key) in
        [("clean_poly", :t_clean), ("preimage_linear", :t_pre), ("isempty", :t_empty)]
        pe, re = fit_exponent([r.F for r in prows], [getfield(r, key) for r in prows])
        @printf("  %-18s t ~ F^%.2f   (r² = %.3f)\n", name, pe, re)
    end
    return nothing
end

println(
    """
N8 step 3 — the exponent behind `rem:complexity`.

Superlinear (p > 1) means a cell's facet count is worth paying to reduce, which is mechanism B's
premise. Near-linear means route 2's simpler cells buy little and only its CELL COUNT (mechanism A)
can produce a speedup. The distinction decides whether N10 is worth running at all.""",
)

for n in (2, 3, 4)
    report(n)
end
