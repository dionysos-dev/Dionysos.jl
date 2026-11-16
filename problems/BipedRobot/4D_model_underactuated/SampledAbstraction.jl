# The builder below is `BipedAbstraction.build_abstraction` with a different approximation: both files
# must share one copy of that module (the provenance of a sampled approximation is registered in it),
# so it is included here only if the file that includes this one has not done it already.
isdefined(@__MODULE__, :BipedAbstraction) ||
    include(joinpath(@__DIR__, "BipedAbstraction.jl"))

"""
    SampledAbstraction

The "data-driven" way to compute the successors of a cell: instead of bounding the dynamics (a growth
bound, interval arithmetic, symbolic Jacobians), **simulate `N` points of the cell** with the RK4 map of
the system, find the cell that each image belongs to, and take the union of those cells as the
successors of the cell. Nothing has to be known about the dynamics but how to evaluate them, so it
works for any black-box system, continuous or discrete in time, in any dimension.

The set of successors is *not* guaranteed to contain the image of every point of the cell: it is an
over-approximation in the empirical sense (the more samples, the closer to the true successors). The
points are chosen by the flag `sampling`:

* `:random` (the default): `N` points independent and uniformly distributed in the cell
  ([`RandomSampler`](@ref)), a Monte Carlo sample. This is the hypothesis that a probabilistic
  certificate needs (computing one is not done here).
* `:uniform`: a regular, deterministic design ([`UniformSampler`](@ref)): the center for `N = 1`,
  the center and the `2ⁿ` corners for `N = 2ⁿ + 1` (17 in four dimensions), the lattice of three
  points per side for `N = 3ⁿ`, then finer lattices, and the points in between for the values of `N`
  that are not one of these. Nothing is random, so there is no seed and no probability to speak of.

Nothing here belongs to Dionysos itself, which already has an `approx_mode = RANDOM_SIMULATION`
doing the random version; this module differs where that one cannot serve a comparison of methods or
a computation on a cluster:

* **Reproducible.** The random points of a cell are a function of the seed and of the cell only
  (counter based, not a stateful random number generator): the same in every order of the cells,
  number of threads, chunk of a Slurm array, Julia version and machine, so the same abstraction comes
  out (up to the rounding of the dynamics, which matters only for an image that lands on the edge of
  a cell). `RANDOM_SIMULATION` draws from the global generator, so no two runs agree.
* **The same samples for every input** of a cell, so that the inputs are compared on the same points.
* **Cheap sampling.** Drawing the points of a cell costs a few nanoseconds each, against microseconds
  for `RANDOM_SIMULATION` (`LazySets.sample` allocates and rejects), which on a pendulum is more than
  simulating them: the abstraction of the simple pendulum with 16 samples is several times faster.
* **Provenance.** Cache files and Slurm chunks record `N`, the seed (if the sampling is random), the
  sampler, the sampling period and a fingerprint of the dynamics, and are refused by anything that was
  computed otherwise (see `BipedAbstraction`).

It plugs into Dionysos through the extension points that exist already, as `CertifiedRK4` does: the
`USER_DEFINED` mode of the `UniformGridAbstraction` optimizer, a point-cloud image type
([`SampledImage`](@ref)) with its own `MP.get_states_from_set_strict` method, and the quantizer of the
state grid to find the cell of a point (the one that `CENTER_SIMULATION` uses, so that one sample at
the center gives exactly that abstraction). [`build_sampled_abstraction`](@ref) is
`BipedAbstraction.build_abstraction` for it, so the cache file, the Slurm chunks and the execution
backends work the same way.

# Example

```julia
include("SampledAbstraction.jl"); using .SampledAbstraction

opt = build_sampled_abstraction(
    system;                       # continuous or discrete time
    state_grid, input_grid,
    time_step = 0.1,              # continuous time only
    nsamples = 50, seed = 1,      # 50 random points per cell
)
# ... or the center and the 16 corners of every cell (four states), without randomness:
opt = build_sampled_abstraction(
    system; state_grid, input_grid, time_step = 0.1, nsamples = 17, sampling = :uniform,
)
abstract_system = MOI.get(opt, MOI.RawOptimizerAttribute("abstract_system"))
```
"""
module SampledAbstraction

import Dionysos
const ST = Dionysos.System
const MP = Dionysos.Mapping
const UGA = Dionysos.Optim.Abstraction.UniformGridAbstraction

import LazySets
import MathematicalSystems as MS
import StaticArrays: SVector

import ..BipedAbstraction as BA

export SampledSimulation, RandomSampler, UniformSampler, SampledImage
export build_sampled_abstraction

# ------------------------------------------------------------------------------------------
# Reproducible samples: a counter-based generator (SplitMix64)
# ------------------------------------------------------------------------------------------

# The `i`-th output of the stream `key`: a pure function, so that nothing depends on the order in which
# the samples are drawn (several threads, chunks, a different `N`: sample `k` is always the same).
# SplitMix64 passes BigCrush, is the same on every machine and Julia version (`Random.Xoshiro` is not
# promised to be), and costs a few nanoseconds.
const _GOLDEN_GAMMA = 0x9e3779b97f4a7c15

@inline function _mix64(z::UInt64)
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return z ⊻ (z >> 31)
end

@inline _stream(key::UInt64, i::Integer) = _mix64(key + (UInt64(i) + 1) * _GOLDEN_GAMMA)

# A number in [0, 1), from the 53 high bits.
@inline _unit(key::UInt64, i::Integer) = (_stream(key, i) >> 11) * 0x1.0p-53

# The stream of a cell: the seed and the coordinates of its center, hashed (the fingerprint of
# `BipedAbstraction` is the one that the provenance uses) and mixed once more so that nearby cells
# get streams that are far apart.
function _cell_key(seed::Int, center::SVector)
    return _mix64(BA.CertifiedRK4.fingerprint((seed, center .+ 0.0)))
end

"""
    RandomSampler()

The sampler of `sampling = :random`: `nsamples` points, independent and uniformly distributed in the
box of a cell (a Monte Carlo sample), reproducible from the seed. A sampler is any callable

    sampler(cell, nsamples::Int, key::UInt64) -> AbstractVector of points

where `cell` is a `LazySets.AbstractHyperrectangle`, and `key` identifies the cell and the seed: a
sampler that wants to be reproducible must be a function of `(cell, nsamples, key)` alone. Its `string`
is recorded in the provenance of the transitions, so it must identify it. It is given to a
[`SampledSimulation`](@ref) as its `sampling`. (Two optional functions describe a sampler to it:
[`prepare`](@ref), to do once what depends only on the dimension and `nsamples`, and the traits
`_uses_seed` and `_contains_center`, which are `true` and `false` unless it says otherwise.)
"""
struct RandomSampler end

function (::RandomSampler)(
    cell::LazySets.AbstractHyperrectangle,
    nsamples::Int,
    key::UInt64,
)
    n = LazySets.dim(cell)
    c = SVector{n, Float64}(LazySets.center(cell))
    ρ = SVector{n, Float64}(LazySets.radius_hyperrectangle(cell))
    return [_sample_box(c, ρ, key, k) for k in 1:nsamples]
end

# Sample `k` (counting from 1): the coordinate `j` of it is the output `(k - 1) n + j` of the stream.
@inline function _sample_box(c::SVector{N, Float64}, ρ, key::UInt64, k::Int) where {N}
    return SVector{N, Float64}(
        ntuple(j -> c[j] + ρ[j] * (2 * _unit(key, (k - 1) * N + j - 1) - 1), Val(N)),
    )
end

"""
    prepare(sampler, n, nsamples) -> sampler

Called once by [`SampledSimulation`](@ref) with the state dimension `n` and the number of points: a
sampler that has something to compute that depends on nothing else (the design of a
[`UniformSampler`](@ref)) returns itself with it computed, so that nothing is computed per cell.
The default returns `sampler`.
"""
prepare(sampler, n, nsamples) = sampler

# What a sampler says about itself: whether its points depend on the seed (if not, the seed is not
# part of the provenance: another seed is the same abstraction), and whether the center of the cell
# is among them (if so, `include_center` would simulate it twice).
_uses_seed(sampler) = true
_contains_center(sampler) = false

# ------------------------------------------------------------------------------------------
# The uniform design: the center, the corners, then finer and finer lattices
# ------------------------------------------------------------------------------------------

# The points are in the cube [-1, 1]ⁿ, as the columns of a matrix, in the order in which they are
# selected: a design for `N` points is the first `N` columns of a design for more (so the abstraction
# only grows with `N`).
#
# Level 0 is the center; level 1 the corners {-1, 1}ⁿ; level `j ≥ 2` the points of the lattice with
# `2^(j-1) + 1` points per side that the lower levels do not have: the lattices {-1, 0, 1}ⁿ,
# {-1, -½, 0, ½, 1}ⁿ, ... Every level is completed before the next one starts, so that `N = 1`, `2ⁿ + 1`,
# `3ⁿ`, `5ⁿ`, ... give a regular lattice (the center and the corners for `2ⁿ + 1`). Inside a level that
# `N` ends halfway, the points are taken one at a time, the one that is the farthest from all those
# already selected (the first one in the order of the lattice on a tie): a partial level is spread over
# the cube, and not a corner of it.
const _MAX_LATTICE = 5_000_000      # points of a lattice that can be listed
const _MAX_DESIGN_WORK = 2.0e9      # coordinates compared to choose the points of a level

function _uniform_design(n::Int, N::Int)
    n >= 1 || throw(ArgumentError("the uniform design needs at least one dimension"))
    N >= 1 || throw(ArgumentError("the uniform design needs at least one point, got $N"))
    design = zeros(n, 1)
    level = 0
    previous = Float64[]                    # the coordinates of the lattice of the level before
    while size(design, 2) < N
        level += 1
        side =
            level == 1 ? [-1.0, 1.0] :
            [-1.0 + 2.0 * i / 2^(level - 1) for i in 0:(2 ^ (level - 1))]
        Float64(length(side))^n <= _MAX_LATTICE || throw(
            ArgumentError(
                "the uniform design of $N points in $n dimensions needs a lattice of " *
                "$(length(side))^$n points: use sampling = :random",
            ),
        )
        candidates = _new_lattice_points(side, n, level, previous)
        m = size(candidates, 2)
        picks = min(N - size(design, 2), m)
        m * (picks + size(design, 2)) * n <= _MAX_DESIGN_WORK || throw(
            ArgumentError(
                "the uniform design of $N points in $n dimensions is too large to build: " *
                "use sampling = :random",
            ),
        )
        design = hcat(design, _farthest_first(candidates, design, picks))
        previous = side
    end
    return design[:, 1:N]
end

# The points of the lattice `side`ⁿ that the levels before `level` do not have, in the order of the
# lattice (the first coordinate varying the fastest).
function _new_lattice_points(side, n, level, previous)
    s = length(side)
    points = Vector{Float64}[]
    digit = zeros(Int, n)
    for index in 0:(s ^ n - 1)
        rest = index
        for j in 1:n
            digit[j] = rest % s
            rest ÷= s
        end
        p = [side[d + 1] for d in digit]
        known =
            level == 1 ? false :
            level == 2 ? (all(x -> abs(x) == 1, p) || all(iszero, p)) :
            all(x -> x in previous, p)
        known || push!(points, p)
    end
    return isempty(points) ? zeros(n, 0) : reduce(hcat, points)
end

# The `picks` columns of `candidates`, each the farthest (in the Euclidean sense) from the columns of
# `selected` and from the ones chosen before it.
function _farthest_first(candidates, selected, picks)
    n, m = size(candidates)
    distance = fill(Inf, m)               # squared distance of a candidate to the selected points
    for k in axes(selected, 2), i in 1:m
        distance[i] = min(distance[i], _squared_distance(candidates, i, selected, k))
    end
    # The points of the level are mostly as far as the points of the levels below allow, which leaves
    # many ties: among the farthest, the one farthest from the points that this level already chose
    # (so the second midpoint of an edge is the opposite one), and then the first in the lattice.
    within = fill(Inf, m)
    chosen = zeros(n, picks)
    for t in 1:picks
        i = argmax(distance)
        for k in 1:m
            distance[k] == distance[i] && within[k] > within[i] && (i = k)
        end
        chosen[:, t] = candidates[:, i]
        distance[i] = -1.0                # taken
        for k in 1:m
            distance[k] >= 0 || continue
            d = _squared_distance(candidates, k, chosen, t)
            distance[k] = min(distance[k], d)
            within[k] = min(within[k], d)
        end
    end
    return chosen
end

function _squared_distance(A, i, B, k)
    d = 0.0
    for j in axes(A, 1)
        d += (A[j, i] - B[j, k])^2
    end
    return d
end

"""
    UniformSampler()

The sampler of `sampling = :uniform`: a regular design, the same in every cell, with no randomness.
The points are the center of the cell, then its `2ⁿ` corners, then the points of finer and finer
lattices (3, 5, 9, ... points per side), each level completed before the next; inside a level that
`nsamples` ends halfway, the points are chosen one at a time as the farthest from the ones already
chosen. So `nsamples = 1` is the center, `2ⁿ + 1` the center and the corners (17 for four states, 5 for
two), `3ⁿ` the lattice of three points per side, and the design of `nsamples` points is the first
`nsamples` points of that of more.

The corners lie on the boundary of the cell, which belongs to its neighbors too: the successors of a
corner are those of a state of the closed cell, which is what an over-approximation of a cell is. A
state exactly on a boundary between two cells is assigned by the quantizer of the grid (ties to even).

The design is computed once by [`prepare`](@ref), when the [`SampledSimulation`](@ref) is made; called
by hand, `UniformSampler()(cell, nsamples, key)` computes it each time. The `key` is ignored.
"""
struct UniformSampler{D}
    design::D          # `nothing`, or the points of the cube [-1, 1]ⁿ that `prepare` computed
end

UniformSampler() = UniformSampler(nothing)

Base.show(io::IO, ::UniformSampler) = print(io, "UniformSampler()")

prepare(::UniformSampler, n, nsamples) =
    UniformSampler([SVector{n, Float64}(p) for p in eachcol(_uniform_design(n, nsamples))])

_uses_seed(::UniformSampler) = false
_contains_center(::UniformSampler) = true

_design(::UniformSampler{Nothing}, n, nsamples) =
    [SVector{n, Float64}(p) for p in eachcol(_uniform_design(n, nsamples))]
_design(sampler::UniformSampler, n, nsamples) = sampler.design

function (sampler::UniformSampler)(
    cell::LazySets.AbstractHyperrectangle,
    nsamples::Int,
    key::UInt64,
)
    n = LazySets.dim(cell)
    c = SVector{n, Float64}(LazySets.center(cell))
    ρ = SVector{n, Float64}(LazySets.radius_hyperrectangle(cell))
    design = _design(sampler, n, nsamples)
    (length(design) == nsamples && eltype(design) == SVector{n, Float64}) || throw(
        ArgumentError(
            "this sampler was prepared for $(length(design)) points of dimension " *
            "$(length(first(design))), not $nsamples points of dimension $n",
        ),
    )
    return [c .+ ρ .* p for p in design]
end

# ------------------------------------------------------------------------------------------
# The image of a cell: a cloud of points
# ------------------------------------------------------------------------------------------

"""
    SampledImage(points)

The images of the sample points of a cell, as a `LazySets.LazySet` so that it is what an
over-approximation map of Dionysos may return. The cells that it covers are the cells of its points
(`MP.get_state_by_coord`, the quantizer that Dionysos uses to say which abstract state a concrete state
belongs to), not the cells that touch their bounding box. If a point is outside the mapped domain, or
is not a finite number, the lookup reports that the image leaves it (`allin = false`) and Dionysos
disables the input for that cell, as it does for every over-approximation.
"""
struct SampledImage{N} <: LazySets.LazySet{Float64}
    points::Vector{SVector{N, Float64}}
end

LazySets.dim(::SampledImage{N}) where {N} = N

# The bounding box, in case a set needs one.
LazySets.low(S::SampledImage{N}) where {N} =
    SVector{N, Float64}(ntuple(j -> minimum(p[j] for p in S.points), Val(N)))
LazySets.high(S::SampledImage{N}) where {N} =
    SVector{N, Float64}(ntuple(j -> maximum(p[j] for p in S.points), Val(N)))

# `round(Int, ·)` of the quantizer throws on a NaN or a number beyond the integers: such a successor
# is outside every domain.
function _abstract_state(mapping, y)
    all(isfinite, y) || return nothing
    try
        return MP.get_state_by_coord(mapping, y)
    catch err
        err isa InexactError || rethrow()
        return nothing
    end
end

function MP.get_states_from_set_strict(
    mapping::MP.GridMapping,
    image::SampledImage,
    ::MP.INCL_MODE,
)
    targets = Int[]
    for y in image.points
        q = _abstract_state(mapping, y)
        q === nothing && return targets, false
        push!(targets, q)
    end
    return unique!(targets), true
end

# ------------------------------------------------------------------------------------------
# The approximation
# ------------------------------------------------------------------------------------------

"""
    SampledSimulation(system::ConstrainedBlackBoxControlContinuousSystem, tstep; nsamples,
                      sampling = :random, seed = 0, include_center = false, nsystem = 5)
    SampledSimulation(system::ConstrainedBlackBoxControlDiscreteSystem; nsamples, ...)

Successors of a cell by simulation of `nsamples` of its points: the image of the cell under an input
is the [`SampledImage`](@ref) of the images of its sample points, found by the RK4 map of the system
(`nsystem` substeps of the sampling period `tstep`, as Dionysos discretizes every continuous system) or
by the map of a discrete-time system.

`sampling` chooses the points: `:random` (or `"random"`), independent points uniformly distributed in
the cell, reproducible from `seed` and the cell ([`RandomSampler`](@ref)); `:uniform`, the regular
design of [`UniformSampler`](@ref) (center, corners, then finer lattices), which has no seed; or a
sampler of your own (see [`RandomSampler`](@ref)). With `include_center = true` the center of the cell
is simulated in addition to the `nsamples` points (`nsamples = 0` and `include_center = true` is
`CENTER_SIMULATION`); a uniform design has the center already, as its first point (`nsamples = 1` is
`CENTER_SIMULATION`), and refuses it.

It is a `DiscreteTimeSystemOverApproximation` of Dionysos, **in the empirical sense only**: its images
are not guaranteed to contain those of all the points of the cell. `ST.get_system` is the discrete-time
system (the RK4 map) to use as the concrete system for simulation and synthesis, and
`ST.get_over_approximation_map` gives the function that the `USER_DEFINED` mode of the
`UniformGridAbstraction` optimizer takes; [`build_sampled_abstraction`](@ref) does both.
"""
struct SampledSimulation{N, S, F} <: ST.DiscreteTimeSystemOverApproximation
    system::S
    nsamples::Int
    seed::Int
    include_center::Bool
    sampler::F
    tstep::Union{Nothing, Float64}
    nsystem::Union{Nothing, Int}
    dynamics::String
end

function SampledSimulation(
    system::MS.ConstrainedBlackBoxControlDiscreteSystem;
    nsamples::Integer,
    sampling = :random,
    seed::Integer = 0,
    include_center::Bool = false,
    tstep::Union{Nothing, Real} = nothing,
    nsystem::Union{Nothing, Integer} = nothing,
)
    sampler = _sampler(sampling)
    nsamples >= 0 || throw(ArgumentError("nsamples must be nonnegative, got $nsamples"))
    if _contains_center(sampler)
        include_center && throw(
            ArgumentError(
                "$sampler has the center among its points: `include_center` would simulate it twice",
            ),
        )
        nsamples >= 1 ||
            throw(ArgumentError("$sampler needs at least one point (the center), got 0"))
    else
        nsamples + include_center >= 1 || throw(
            ArgumentError("no point to simulate: nsamples = 0 and include_center = false"),
        )
    end
    sampler = prepare(sampler, MS.statedim(system), Int(nsamples))
    # The state dimension is a type parameter, so that the points of a cell are `SVector`s of a known
    # size: as a runtime integer it made the map ten times slower than everything it computes.
    return SampledSimulation{MS.statedim(system), typeof(system), typeof(sampler)}(
        system,
        Int(nsamples),
        Int(seed),
        include_center,
        sampler,
        tstep === nothing ? nothing : Float64(tstep),
        nsystem === nothing ? nothing : Int(nsystem),
        _dynamics_fingerprint(MS.mapping(system), MS.stateset(system), MS.inputset(system)),
    )
end

function SampledSimulation(
    system::MS.ConstrainedBlackBoxControlContinuousSystem,
    tstep::Real;
    nsystem::Integer = ST.DEFAULT_NUM_SUBSTEPS,
    kwargs...,
)
    isfinite(tstep) && tstep > 0 ||
        throw(ArgumentError("the sampling period must be positive and finite, got $tstep"))
    nsystem > 0 || throw(ArgumentError("nsystem must be positive, got $nsystem"))
    discrete =
        ST.discretize_continuous_system(system, Float64(tstep); num_substeps = Int(nsystem))
    return SampledSimulation(discrete; tstep, nsystem, kwargs...)
end

# The flag `sampling` is a name or a sampler.
function _sampler(sampling::Union{Symbol, AbstractString})
    name = Symbol(sampling)
    name === :random && return RandomSampler()
    name === :uniform && return UniformSampler()
    return throw(
        ArgumentError("sampling must be :random, :uniform or a sampler, got :$name"),
    )
end
_sampler(sampler) = sampler

ST.get_system(approx::SampledSimulation) = approx.system

function ST.get_over_approximation_map(approx::SampledSimulation{N}) where {N}
    F = ST.get_system_map(approx)
    return (cell, u) -> begin
        center = SVector{N, Float64}(LazySets.center(cell))
        xs = approx.sampler(cell, approx.nsamples, _cell_key(approx.seed, center))
        points = Vector{SVector{N, Float64}}(undef, length(xs) + approx.include_center)
        k = 0
        if approx.include_center
            points[k += 1] = SVector{N, Float64}(F(center, u))
        end
        for x in xs
            points[k += 1] = SVector{N, Float64}(F(x, u))
        end
        return SampledImage(points)
    end
end

# ------------------------------------------------------------------------------------------
# Provenance of the stored transitions
# ------------------------------------------------------------------------------------------

"""
    KERNEL_VERSION

Version of the sampling kernel. It is part of [`get_provenance`](@ref): bump it whenever a change of
the kernel (the generator, the way a cell is keyed) can change the transitions it returns.
"""
const KERNEL_VERSION = 1

# What identifies the dynamics of a black-box system: its values at a few states and inputs, always the
# same ones (a stream, not a random generator), rounded to ten digits so that a last-bit difference
# between two machines does not tell two runs apart. A system that is written differently but computes
# the same function has the same fingerprint, and one whose parameters differ does not.
function _dynamics_fingerprint(F, X, U; nprobes::Int = 16)
    n, m = LazySets.dim(X), LazySets.dim(U)
    fp = BA.CertifiedRK4.fingerprint
    xbox, ubox = _probe_box(X, n), _probe_box(U, m)
    h = fp((n, m, nprobes))
    for k in 1:nprobes
        x = _probe_point(xbox, k, 0x1)
        u = _probe_point(ubox, k, 0x2)
        values = try
            collect(Float64, F(x, u))
        catch
            Float64[]                 # a domain error is part of the behavior too
        end
        h = fp(map(v -> isfinite(v) ? round(v; sigdigits = 10) + 0.0 : v, values), h)
    end
    return BA.CertifiedRK4.fingerprint_string(h)
end

# The box where the dynamics are probed: that of the set if it is bounded, [-1, 1] otherwise.
function _probe_box(S, n)
    try
        if LazySets.isbounded(S)
            box = LazySets.box_approximation(S)
            return (
                SVector{n, Float64}(LazySets.low(box)),
                SVector{n, Float64}(LazySets.high(box)),
            )
        end
    catch
    end
    return (-ones(SVector{n, Float64}), ones(SVector{n, Float64}))
end

function _probe_point(
    (lo, hi)::Tuple{SVector{N, Float64}, SVector{N, Float64}},
    k,
    tag,
) where {N}
    key = _mix64(UInt64(k) * _GOLDEN_GAMMA + UInt64(tag))
    return SVector{N, Float64}(
        ntuple(j -> lo[j] + (hi[j] - lo[j]) * _unit(key, j - 1), Val(N)),
    )
end

"""
    get_provenance(approx::SampledSimulation) -> NamedTuple

Plain data that identifies **how** the successors of `approx` are computed: the kernel version, the
number of samples, the seed (`nothing` for a sampler whose points do not depend on it, as the uniform
design), whether the center is simulated, the sampler, the sampling period and RK4 substeps, and a
fingerprint of the dynamics. `BipedAbstraction` stores it with cached transitions and chunks, and
refuses transitions that do not carry the same one.
"""
get_provenance(approx::SampledSimulation) = (;
    approximation = "SampledSimulation",
    kernel_version = KERNEL_VERSION,
    nsamples = approx.nsamples,
    seed = _uses_seed(approx.sampler) ? approx.seed : nothing,
    include_center = approx.include_center,
    sampler = string(approx.sampler),
    tstep = approx.tstep,
    nsystem = approx.nsystem,
    dynamics = approx.dynamics,
)

BA.provenance(approx::SampledSimulation) = get_provenance(approx)

# ------------------------------------------------------------------------------------------
# Building the transition system with Dionysos
# ------------------------------------------------------------------------------------------

"""
    build_sampled_abstraction(system; state_grid, input_grid, nsamples, sampling = :random,
                              seed = 0, include_center = false, time_step = nothing,
                              nsystem = 5, attributes...) -> optimizer
    build_sampled_abstraction(approx::SampledSimulation; state_grid, input_grid, attributes...)

Transition system (symbolic model) of `system`, with the successors of a cell computed by simulating
`nsamples` of its points ([`SampledSimulation`](@ref)). It is
[`BipedAbstraction.build_abstraction`](@ref) for this approximation: the same grids, the same
attributes (`print_level`, `execution_backend`, `state_input_filter`, ...), the same optimizer
returned with the abstraction cached, and the same `abstraction_file`, `abstraction_info` and
`merge_chunks` for a cache file and the chunks of a Slurm array, with the provenance of the
approximation stored and checked (a cache file or chunk computed with another `nsamples`, `sampling`,
`seed` (if the sampling is random), sampling period or dynamics is refused). `sampling` is `:random`
or `:uniform`, see [`SampledSimulation`](@ref).

`system` is a continuous-time system, with its sampling period `time_step`, or a discrete-time one.
The concrete system for simulation and synthesis is `ST.get_system` of the approximation, the discrete
system of the RK4 map: the optimizer is given that one, not `system`, and its `discrete_time_system`
is that map. `overapproximation_map`, `approx_mode` and `approx` are set by this function and cannot be
overridden.

```julia
opt = build_sampled_abstraction(
    system; state_grid, input_grid, time_step = 0.1, nsamples = 50, seed = 1,
)
abstract_system = MOI.get(opt, MOI.RawOptimizerAttribute("abstract_system"))
```
"""
function build_sampled_abstraction(
    system::Union{
        MS.ConstrainedBlackBoxControlContinuousSystem,
        MS.ConstrainedBlackBoxControlDiscreteSystem,
    };
    nsamples,
    sampling = :random,
    seed = 0,
    include_center = false,
    time_step = nothing,
    nsystem = ST.DEFAULT_NUM_SUBSTEPS,
    attributes...,
)
    approx = if system isa MS.ConstrainedBlackBoxControlContinuousSystem
        time_step === nothing && throw(
            ArgumentError("a continuous-time system needs its sampling period `time_step`"),
        )
        SampledSimulation(
            system,
            time_step;
            nsamples,
            sampling,
            seed,
            include_center,
            nsystem,
        )
    else
        time_step === nothing || throw(
            ArgumentError(
                "`time_step` is for continuous-time systems; this one is discrete",
            ),
        )
        SampledSimulation(system; nsamples, sampling, seed, include_center)
    end
    return build_sampled_abstraction(approx; attributes...)
end

function build_sampled_abstraction(approx::SampledSimulation; state_grid, attributes...)
    for name in (:overapproximation_map, :approx_mode, :approx, :time_step, :nsystem)
        haskey(attributes, name) &&
            throw(ArgumentError("cannot override $name in a sampled abstraction"))
    end
    return BA.build_abstraction(
        ST.get_system(approx);
        state_grid,
        approx_mode = UGA.USER_DEFINED,
        overapproximation_map = ST.get_over_approximation_map(approx),
        approx,
        attributes...,
    )
end

end # module
