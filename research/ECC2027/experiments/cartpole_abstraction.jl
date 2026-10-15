# The two abstractions of E2, built directly: the velocity model (F) in (θ, p) and the force model
# (T) in (θ, ẋ, θ̇), with one successor set per (cell, input) from a growth bound along the
# constant input (Definition "Abstraction" of the paper).
#
# For (F), the impulse budget |α(θ)(u − v)| ≤ I_max, bounded over the cell, decides which command
# may follow which. The synthesis is Algorithm 1 of the paper: a min–max Dijkstra search on the
# pairs (cell, command) that resolves the augmented states (cell, last command) without building
# the augmented abstraction. The controller is replayed on the exact step dynamics, which are the
# full model by the factorization.
module CartPoleAbstraction

import Dionysos
import StaticArrays: SVector, SMatrix

const PROBLEMS = joinpath(dirname(dirname(pathof(Dionysos))), "problems")
include(joinpath(PROBLEMS, "CartPendulum", "velocity_cart_pendulum.jl"))
import .VelocityCartPendulum as VC

# ------------------------------------------------------------------------------------- grid
struct Grid
    nθ::Int
    hθ::Float64
    np::Int
    hp::Float64
    p_max::Float64
end
Grid(nθ, hp, p_max) = Grid(nθ, 2π / nθ, round(Int, 2p_max / hp), hp, p_max)
ncells(G::Grid) = G.nθ * G.np
_split(G::Grid, c) = ((c - 1) % G.nθ + 1, (c - 1) ÷ G.nθ + 1)
center(G::Grid, c) =
    ((iθ, ip) = _split(G, c); SVector(-π + (iθ - 0.5) * G.hθ, -G.p_max + (ip - 0.5) * G.hp))
function bounds(G::Grid, c)
    z = center(G, c)
    return (z[1] - G.hθ / 2, z[1] + G.hθ / 2), (z[2] - G.hp / 2, z[2] + G.hp / 2)
end
function locate(G::Grid, z)
    iθ = mod(floor(Int, (z[1] + π) / G.hθ), G.nθ) + 1
    ip = floor(Int, (z[2] + G.p_max) / G.hp) + 1
    return 1 <= ip <= G.np ? iθ + G.nθ * (ip - 1) : 0
end

# Cells met by the box z ± r: an unwrapped θ range and a p range, or nothing if p leaves the domain.
function cover(G::Grid, z, r)
    θlo = floor(Int, (z[1] - r[1] + π) / G.hθ)
    θhi = floor(Int, (z[1] + r[1] + π) / G.hθ)
    plo = floor(Int, (z[2] - r[2] + G.p_max) / G.hp) + 1
    phi = floor(Int, (z[2] + r[2] + G.p_max) / G.hp) + 1
    (plo < 1 || phi > G.np) && return nothing
    return (Int32(θlo), Int32(min(θhi - θlo + 1, G.nθ)), Int32(plo), Int32(phi))
end
function cells_of(G::Grid, b)
    θlo, kθ, plo, phi = b
    return (mod(θlo + i, G.nθ) + 1 + G.nθ * (ip - 1) for ip in plo:phi for i in 0:(kθ - 1))
end
ncover(b) = Int(b[2]) * Int(b[4] - b[3] + 1)

# ------------------------------------------------------------------------------- integration
# RK4 under the constant input u.
function flow(f, z, u, τ; K = 20)
    h = τ / K
    for _ in 1:K
        k1 = f(z, u)
        k2 = f(z + h / 2 * k1, u)
        k3 = f(z + h / 2 * k2, u)
        k4 = f(z + h * k3, u)
        z += h / 6 * (k1 + 2k2 + 2k3 + k4)
    end
    return z
end

# The growth bound ṙ = J r (Reissig et al.), integrated with RK4.
function radius(J, r, τ; K = 20)
    h = τ / K
    g(r) = J * r
    for _ in 1:K
        k1 = g(r)
        k2 = g(r + h / 2 * k1)
        k3 = g(r + h / 2 * k2)
        k4 = g(r + h * k3)
        r += h / 6 * (k1 + 2k2 + 2k3 + k4)
    end
    return r
end

# ------------------------------------------------------------------------------- abstraction
"""
    Abstraction

The stored abstraction: its pairs `pairs[k] = (cell, command)`, the cover `boxes[k]` of the
successors of each, and the number `stored` of its transitions. The last command is not stored:
the synthesis adds it.
"""
struct Abstraction{GT, B}
    G::GT
    Σ::Vector{Float64}
    pairs::Vector{NTuple{2, Int32}}
    boxes::Vector{B}
    stored::Int
end

# Upper bound of the inertia α(θ) = M + m sin²θ over a cell.
inertia_bound(P, Θ) = P.M + P.m * VC.isq(VC.isin(Θ))[2]

"""
    impulse_budget(P, G, Σ, I_max) -> compatible

The impulse budget carried to the cells: `compatible(c, iv, iu)` holds when the command `Σ[iu]`
may follow `Σ[iv]` at every state of cell `c`, that is ᾱ |Σ[iu] − Σ[iv]| ≤ I_max with ᾱ the upper
bound of α(θ) over the cell.
"""
function impulse_budget(P, G::Grid, Σ, I_max)
    ᾱ = [inertia_bound(P, bounds(G, c)[1]) for c in 1:ncells(G)]
    return (c, iv, iu) -> ᾱ[c] * abs(Σ[iu] - Σ[iv]) <= I_max + 1e-9
end

function build(P, G::Grid, Σ, τ, target::BitVector)
    f = VC.dynamics(P)
    J = VC.jacobian_bound(P; p_max = G.p_max + 0.5)
    nc = ncells(G)
    r0 = SVector(G.hθ / 2, G.hp / 2)
    chunks = collect(Iterators.partition(1:nc, cld(nc, 4 * Threads.nthreads())))
    perchunk = [(NTuple{2, Int32}[], NTuple{4, Int32}[], Ref(0)) for _ in chunks]
    Threads.@threads for i in eachindex(chunks)
        mine = perchunk[i]
        for c in chunks[i]
            build_cell!(mine..., G, Σ, τ, target, c, f, J, r0)
        end
    end
    return Abstraction(
        G,
        collect(Σ),
        reduce(vcat, (t[1] for t in perchunk)),
        reduce(vcat, (t[2] for t in perchunk)),
        sum(t[3][] for t in perchunk),
    )
end

function build_cell!(pairs, boxes, stored, G, Σ, τ, target, c, f, J, r0)
    target[c] && return
    z0 = center(G, c)
    for (iu, u) in enumerate(Σ)
        z = flow(f, z0, u, τ)
        r = radius(J(SVector(abs(u))), r0, τ)
        b = cover(G, z, r)
        b === nothing && continue
        push!(pairs, (Int32(c), Int32(iu)))
        push!(boxes, b)
        stored[] += ncover(b)
    end
    return
end

# ------------------------------------------------------------------------------- synthesis
# A binary min-heap of (key, pair), enough for Knuth's algorithm.
struct Heap
    keys::Vector{Float64}
    vals::Vector{Int32}
end
Heap() = Heap(Float64[], Int32[])
Base.isempty(h::Heap) = isempty(h.keys)
function Base.push!(h::Heap, key, val)
    push!(h.keys, key)
    push!(h.vals, val)
    i = length(h.keys)
    while i > 1
        j = i ÷ 2
        h.keys[j] <= h.keys[i] && break
        h.keys[i], h.keys[j] = h.keys[j], h.keys[i]
        h.vals[i], h.vals[j] = h.vals[j], h.vals[i]
        i = j
    end
    return h
end
function Base.pop!(h::Heap)
    key, val = h.keys[1], h.vals[1]
    h.keys[1], h.vals[1] = h.keys[end], h.vals[end]
    pop!(h.keys)
    pop!(h.vals)
    n, i = length(h.keys), 1
    while true
        l, r, m = 2i, 2i + 1, i
        l <= n && h.keys[l] < h.keys[m] && (m = l)
        r <= n && h.keys[r] < h.keys[m] && (m = r)
        m == i && break
        h.keys[i], h.keys[m] = h.keys[m], h.keys[i]
        h.vals[i], h.vals[m] = h.vals[m], h.vals[i]
        i = m
    end
    return key, val
end

"""
    synthesize(A, target; compatible = nothing) -> (value, choice)

Min–max reachability (unit cost), Algorithm 1 of the paper: Knuth's generalization of Dijkstra's
algorithm run on the pairs `(cell, command)`. A pair is queued once all its successors are
resolved, at 1 plus their worst value, and popping it resolves every augmented state
`(cell, last command)` at which `compatible(cell, last command, command)` allows it. Without
`compatible` there is no memory and the states are the cells, as for (T). `value[s] = Inf` where
the target cannot be guaranteed; `choice[s]` is the index of the command to play.
"""
function synthesize(A::Abstraction, target::BitVector; compatible = nothing)
    G = A.G
    memory = compatible !== nothing
    nΣ = memory ? length(A.Σ) : 1
    nc = ncells(G)
    nstate = nc * nΣ
    state(c, iv) = memory ? c + nc * (iv - 1) : c
    # Predecessor lists in compressed rows: the pair k = (c, u) reaches the states (c′, u), so
    # each transition of the stored abstraction is listed once.
    counts = zeros(Int, nstate + 1)
    for (k, (_, iu)) in enumerate(A.pairs)
        for c′ in cells_of(G, A.boxes[k])
            counts[state(c′, iu) + 1] += 1
        end
    end
    offsets = cumsum(counts) .+ 1
    preds = Vector{Int32}(undef, offsets[end] - 1)
    fill_at = copy(offsets)
    for (k, (_, iu)) in enumerate(A.pairs)
        for c′ in cells_of(G, A.boxes[k])
            s′ = state(c′, iu)
            preds[fill_at[s′]] = k
            fill_at[s′] += 1
        end
    end
    remaining = Int32[ncover(b) for b in A.boxes]
    worst = zeros(Float64, length(A.pairs))
    value = fill(Inf, nstate)
    choice = zeros(Int32, nstate)
    heap = Heap()
    # Resolve the state s at the value θ, reached by playing the command iu, and credit the pairs
    # that lead to it: a pair whose successors are all resolved is queued at its final value.
    function resolve!(s, θ, iu)
        value[s] = θ
        choice[s] = iu
        for i in offsets[s]:(offsets[s + 1] - 1)
            k = preds[i]
            worst[k] = max(worst[k], θ)
            remaining[k] -= 1
            remaining[k] == 0 && push!(heap, 1.0 + worst[k], k)
        end
    end
    for c in 1:nc, iv in 1:nΣ
        target[c] && resolve!(state(c, iv), 0.0, 0)
    end
    while !isempty(heap)
        θ, k = pop!(heap)
        c, iu = A.pairs[k]
        for iv in 1:nΣ
            s = state(c, iv)
            value[s] < Inf && continue
            memory && !compatible(c, iv, iu) && continue
            resolve!(s, θ, iu)
        end
    end
    return value, choice
end

# ------------------------------------------------------------------------------- replay
"""
    replay(A, choice, P, τ, z0, target; nmax, K) -> NamedTuple

Run the controller on the exact step dynamics from `z0` at rest, sampling each step `K` times:
the times, angle, angular velocity, cart position, the hold force `m sin θ (g cos θ − l θ̇²)`
between samples, the impulse `α(θ)(u − v)` of each step, and whether the target was reached.
"""
function replay(A::Abstraction, choice, P, τ, z0, target; nmax = 2000, K = 20)
    G, Σ = A.G, A.Σ
    nc = ncells(G)
    f = VC.dynamics(P)
    out = (t = Float64[], θ = Float64[], θdot = Float64[], x = Float64[], F = Float64[])
    impulses = Float64[]
    z, iv, x, t = z0, findfirst(==(0.0), Σ), 0.0, 0.0
    reached = false
    for _ in 1:nmax
        c = locate(G, z)
        c == 0 && break
        if target[c]
            reached = true
            break
        end
        s = c + nc * (iv - 1)
        iu = choice[s]
        iu == 0 && break
        v, u = Σ[iv], Σ[iu]
        push!(impulses, (P.M + P.m * sin(z[1])^2) * (u - v))
        h = τ / K
        for j in 0:(K - 1)
            θdot = VC.angular_velocity(P, z[1], z[2], u)
            push!(out.t, t + j * h)
            push!(out.θ, z[1])
            push!(out.θdot, θdot)
            push!(out.x, x + u * j * h)
            push!(out.F, VC.force(P, z[1], θdot, 0.0))
            z = flow(f, z, u, h; K = 1)
        end
        x += u * τ
        t += τ
        iv = iu
    end
    return (; out..., impulses, reached)
end

# -------------------------------------------------------------------------- force model (T)
# The classical abstraction: state y = (θ, ẋ, θ̇), input the cart force, held over the step. The
# cart position is left out, as in (F). The growth bound uses a Jacobian bound computed cell by
# cell, by interval arithmetic over the tube of the step, which is tighter than a global one.
struct Grid3
    nθ::Int
    hθ::Float64
    nv::Int
    hv::Float64
    v_max::Float64
    nω::Int
    hω::Float64
    ω_max::Float64
end
Grid3(nθ, hv, v_max, hω, ω_max) = Grid3(
    nθ,
    2π / nθ,
    round(Int, 2v_max / hv),
    hv,
    v_max,
    round(Int, 2ω_max / hω),
    hω,
    ω_max,
)
ncells(G::Grid3) = G.nθ * G.nv * G.nω
function center(G::Grid3, c)
    iθ = (c - 1) % G.nθ + 1
    rest = (c - 1) ÷ G.nθ
    iv, iω = rest % G.nv + 1, rest ÷ G.nv + 1
    return SVector(
        -π + (iθ - 0.5) * G.hθ,
        -G.v_max + (iv - 0.5) * G.hv,
        -G.ω_max + (iω - 0.5) * G.hω,
    )
end
halfwidth(G::Grid3) = SVector(G.hθ / 2, G.hv / 2, G.hω / 2)
function locate(G::Grid3, y)
    iθ = mod(floor(Int, (y[1] + π) / G.hθ), G.nθ) + 1
    iv = floor(Int, (y[2] + G.v_max) / G.hv) + 1
    iω = floor(Int, (y[3] + G.ω_max) / G.hω) + 1
    return (1 <= iv <= G.nv && 1 <= iω <= G.nω) ?
           iθ + G.nθ * (iv - 1) + G.nθ * G.nv * (iω - 1) : 0
end
function cover(G::Grid3, y, r)
    θlo = floor(Int, (y[1] - r[1] + π) / G.hθ)
    θhi = floor(Int, (y[1] + r[1] + π) / G.hθ)
    vlo = floor(Int, (y[2] - r[2] + G.v_max) / G.hv) + 1
    vhi = floor(Int, (y[2] + r[2] + G.v_max) / G.hv) + 1
    ωlo = floor(Int, (y[3] - r[3] + G.ω_max) / G.hω) + 1
    ωhi = floor(Int, (y[3] + r[3] + G.ω_max) / G.hω) + 1
    (vlo < 1 || vhi > G.nv || ωlo < 1 || ωhi > G.nω) && return nothing
    return (
        Int32(θlo),
        Int32(min(θhi - θlo + 1, G.nθ)),
        Int32(vlo),
        Int32(vhi),
        Int32(ωlo),
        Int32(ωhi),
    )
end
cells_of(G::Grid3, b::NTuple{6, Int32}) = (
    mod(b[1] + i, G.nθ) + 1 + G.nθ * (iv - 1) + G.nθ * G.nv * (iω - 1) for
    iω in b[5]:b[6] for iv in b[3]:b[4] for i in 0:(b[2] - 1)
)
ncover(b::NTuple{6, Int32}) = Int(b[2]) * Int(b[4] - b[3] + 1) * Int(b[6] - b[5] + 1)

function force_dynamics(P)
    return function (y, F)
        s, c = sincos(y[1])
        a = (F + P.m * s * (P.l * y[3]^2 - P.g * c)) / (P.M + P.m * s^2)
        return SVector(y[3], a, (P.g * s - c * a) / P.l)
    end
end

_div(a, d) = extrema((a[1] / d[1], a[1] / d[2], a[2] / d[1], a[2] / d[2])) # d > 0
_hull0(a) = (min(0.0, a[1]), max(0.0, a[2]))
_absmax(a) = max(abs(a[1]), abs(a[2]))

# Interval rates of (θ, ẋ, θ̇) over a box, and the pieces the Jacobian reuses.
function _rates_T(P, Θ, Ω, F)
    S, C = VC.isin(Θ), VC.icos(Θ)
    D = VC.iadd((P.M, P.M), VC.iscale(P.m, VC.isq(S)))
    N = VC.iadd(
        (F, F),
        VC.iscale(P.m, VC.imul(S, VC.iadd(VC.iscale(P.l, VC.isq(Ω)), VC.iscale(-P.g, C)))),
    )
    A = _div(N, D)
    Θdd = VC.iscale(1 / P.l, VC.iadd(VC.iscale(P.g, S), VC.iscale(-1.0, VC.imul(C, A))))
    return A, Θdd, (S, C, D, N)
end

function _tube_T(P, Y0, F, τ)
    widen(B, k) = (
        (B[1][1] - k, B[1][2] + k),
        (B[2][1] - k, B[2][2] + k),
        (B[3][1] - 4k, B[3][2] + 4k),
    )
    B = widen(Y0, 0.25)
    found = nothing
    for _ in 1:60
        A, Θdd, _ = _rates_T(P, B[1], B[3], F)
        next = (
            VC.iadd(Y0[1], VC.iscale(τ, _hull0(B[3]))),
            VC.iadd(Y0[2], VC.iscale(τ, _hull0(A))),
            VC.iadd(Y0[3], VC.iscale(τ, _hull0(Θdd))),
        )
        if all(VC.iwithin(next[i], B[i]) for i in 1:3)
            found = next
            B = next
        elseif found === nothing
            B = widen(
                (
                    (min(B[1][1], next[1][1]), max(B[1][2], next[1][2])),
                    (min(B[2][1], next[2][1]), max(B[2][2], next[2][2])),
                    (min(B[3][1], next[3][1]), max(B[3][2], next[3][2])),
                ),
                0.25,
            )
        else
            break
        end
    end
    return found # nothing: no enclosure, the transition is not added
end

# Jacobian bound of the force model over a box: signed upper bounds on the diagonal, absolute
# values off it.
function _jacobian_T(P, Θ, Ω, F)
    _, _, (S, C, D, N) = _rates_T(P, Θ, Ω, F)
    m, l, g = P.m, P.l, P.g
    A = _div(N, D)
    cos2 = VC.icos((2Θ[1], 2Θ[2]))
    num = VC.iadd(
        VC.imul(
            VC.iscale(m, VC.iadd(VC.iscale(l, VC.imul(C, VC.isq(Ω))), VC.iscale(-g, cos2))),
            D,
        ),
        VC.iscale(-2m, VC.imul(VC.imul(S, C), N)),
    )
    dA_dθ = _div(num, VC.isq(D))
    dA_dω = _div(VC.iscale(2m * l, VC.imul(S, Ω)), D)
    dB_dθ = VC.iscale(
        1 / l,
        VC.iadd(
            VC.iadd(VC.iscale(g, C), VC.imul(S, A)),
            VC.iscale(-1.0, VC.imul(C, dA_dθ)),
        ),
    )
    dB_dω = VC.iscale(-1 / l, VC.imul(C, dA_dω))
    return SMatrix{3, 3}(
        0.0,
        _absmax(dA_dθ),
        _absmax(dB_dθ),
        0.0,
        0.0,
        0.0,
        1.0,
        _absmax(dA_dω),
        dB_dω[2],
    )
end

function build_T(P, G::Grid3, Fs, τ, target::BitVector)
    f = force_dynamics(P)
    nc = ncells(G)
    r0 = halfwidth(G)
    chunks = collect(Iterators.partition(1:nc, cld(nc, 4 * Threads.nthreads())))
    perchunk = [(NTuple{2, Int32}[], NTuple{6, Int32}[], Ref(0)) for _ in chunks]
    Threads.@threads for i in eachindex(chunks)
        mine = perchunk[i]
        for c in chunks[i]
            target[c] && continue
            y0 = center(G, c)
            Y0 = Tuple((y0[d] - r0[d], y0[d] + r0[d]) for d in 1:3)
            for (iF, F) in enumerate(Fs)
                Yt = _tube_T(P, Y0, F, τ)
                Yt === nothing && continue
                J = _jacobian_T(P, Yt[1], Yt[3], F)
                y = flow(f, y0, F, τ)
                r = radius(J, r0, τ)
                b = cover(G, y, r)
                b === nothing && continue
                push!(mine[1], (Int32(c), Int32(iF)))
                push!(mine[2], b)
                mine[3][] += ncover(b)
            end
        end
    end
    return Abstraction(
        G,
        collect(Float64, Fs),
        reduce(vcat, (t[1] for t in perchunk)),
        reduce(vcat, (t[2] for t in perchunk)),
        sum(t[3][] for t in perchunk),
    )
end

"""
    replay_T(A, choice, P, τ, y0, target; nmax, K)

Run the controller of the force model (T) on the force dynamics, the force held over each step.
"""
function replay_T(A::Abstraction, choice, P, τ, y0, target; nmax = 2000, K = 20)
    G = A.G
    f = force_dynamics(P)
    out = (t = Float64[], θ = Float64[], θdot = Float64[], x = Float64[], F = Float64[])
    y, x, t, reached = y0, 0.0, 0.0, false
    for _ in 1:nmax
        c = locate(G, y)
        c == 0 && break
        if target[c]
            reached = true
            break
        end
        iF = choice[c]
        iF == 0 && break
        F = A.Σ[iF]
        h = τ / K
        for j in 0:(K - 1)
            push!(out.t, t + j * h)
            push!(out.θ, y[1])
            push!(out.θdot, y[3])
            push!(out.x, x)
            push!(out.F, F)
            ynext = flow(f, y, F, h; K = 1)
            x += (y[2] + ynext[2]) * h / 2
            y = ynext
        end
        t += τ
    end
    return (; out..., reached)
end

end # module
