"""
    LocalGrowthAbstraction

Discrete-time transition systems from **cell-wise (local) growth bounds**.

The classical growth bound of Dionysos (`ST.ContinuousTimeGrowthBound`) uses one Jacobian bound
`L(u)` valid on the whole state space `X`: every cell of radius `r` is mapped to
`Φ(c, u) ± ρ(r, u)` with `ρ̇ = L(u) ρ`, `ρ(0) = r`. Whenever the dynamics are only strongly
non-linear in a small part of `X`, `L` is dictated by the worst point and the image of every other
cell is needlessly large, which creates spurious transitions.

This module computes, for each pair `(cell R, input u)`, a **local** bound:

1. start from the cell `R = c ± r` and the input `u`;
2. compute a conservative image `Φ(c, u) ± ρ_glob` with the *global* growth bound, whose
   Jacobian bound is obtained with interval arithmetic over all of `X`;
3. take a larger box `B ⊇ R ∪ image` and bound the Jacobian over `B` only, again with interval
   arithmetic: the true (local) growth bound `L_B`;
4. recompute the image with it: `Φ(c, u) ± ρ_loc`, `ρ̇_loc = L_B ρ_loc`.

The growth-bound theorem needs `L_B` to bound the Jacobian on a convex
set containing the trajectories of **all** points of `R` during the **whole** sampling period, not
just at its end points. So `B` is verified, with an interval Picard–Lindelöf check, to be an
a-priori enclosure of that tube (`R + [0, τ]·f(B, u) ⊆ interior(B)` with outward rounding).
If the first refinement cannot verify a box, it falls back to the global radius. The global
Jacobian bound is only valid on the outer box of `X`. The computed radius is capped by the
global radius.

The tube check applies to the exact ODE flow. The endpoint image still uses ordinary RK4 for
the centre and radius, the global radius cap, and optional radius deflation. It does not
certify successor coverage for the exact flow or the RK4 map.
Use [`CertifiedRK4GrowthBound`](@ref) (file `CertifiedRK4.jl`, next to this one) to build a
certified abstraction of its own Float64 RK4 map.

The result, [`LocalGrowthBound`](@ref), is a `Dionysos.System.DiscreteTimeSystemOverApproximation`,
so it plugs into Dionysos' abstraction machinery through the `USER_DEFINED` approximation mode of
the `UniformGridAbstraction` optimizer (every execution backend works). [`build_abstraction`](@ref)
wraps that. Nothing in Dionysos itself had to change for it.

# Example

```julia
using Symbolics, IntervalArithmetic          # loads Dionysos' interval extension too
include("LocalGrowthAbstraction.jl"); using .LocalGrowthAbstraction

@variables x[1:2] u[1:1]
f = [x[2], -9.81 * sin(x[1]) - 0.1 * x[2] + u[1]]          # pendulum
X = LazySets.Hyperrectangle(; low = [-π, -8.0], high = [π, 8.0])
U = LazySets.Hyperrectangle(; low = [-2.0], high = [2.0])

lgb = LocalGrowthBound(f, collect(x), collect(u), X, U, 0.1)
opt = build_abstraction(
    lgb;
    state_grid = MP.GridFree(SVector(0.0, 0.0), SVector(0.1, 0.2)),
    input_grid = MP.GridFree(SVector(0.0), SVector(1.0)),
)
abstract_system = MOI.get(opt, MOI.RawOptimizerAttribute("abstract_system"))
```
"""
module LocalGrowthAbstraction

import Dionysos
const UT = Dionysos.Utils
const ST = Dionysos.System
const PR = Dionysos.Problem
const MP = Dionysos.Mapping
const SY = Dionysos.Symbolic
const OP = Dionysos.Optim
const AB = OP.Abstraction

import MathematicalSystems as MS
import LazySets
import MathOptInterface as MOI
import JLD2
import Symbolics
import IntervalArithmetic as IA
using StaticArrays

# The certified successors of the Float64 RK4 map (a submodule of this file's folder, not of Dionysos).
include(joinpath(@__DIR__, "CertifiedRK4.jl"))
using .CertifiedRK4: CertifiedRK4GrowthBound
import .CertifiedRK4: overapproximation_map

export SymbolicDynamics, LocalGrowthBound, LocalGrowthBoundStats, CertifiedRK4GrowthBound
export overapproximation_map
export jacobian_bound, flow_enclosure, build_abstraction, stats, reset_stats!
export save_abstraction, load_abstraction!, merge_chunks!, merge_chunks_to_file

# ------------------------------------------------------------------------------------------
# Symbolic dynamics, compiled for floats and for intervals
# ------------------------------------------------------------------------------------------

"""
    SymbolicDynamics{N, M}

Continuous-time dynamics `ẋ = f(x, u)` (`N` states, `M` inputs) given symbolically, compiled once
into three callables:

* `f(x, u) -> SVector{N}`: plain floating-point evaluation;
* `f_interval(xbox, ubox) -> Vector`: interval evaluation (an enclosure of `f` over the boxes);
* `jacobian_interval(xbox, ubox) -> Matrix`: interval enclosure of `∂f/∂x` over the boxes.

The Jacobian is differentiated symbolically and evaluated with interval arithmetic, exactly as
the Dionysos Symbolics extension does for `compute_jacobian_bound`.

    SymbolicDynamics(fsym, x, u)
    SymbolicDynamics(system::MS.ConstrainedBlackBoxControlContinuousSystem)

`fsym` is a vector of `Symbolics` expressions in the variable vectors `x` and `u`. The second form
traces the black-box `f` of a Dionysos system with symbolic variables (it must be written with
generic operations). Interval evaluation assumes those operations have sound interval
implementations. Numeric constants denote their supplied values, including the stored binary
values of floating-point constants.
"""
struct SymbolicDynamics{N, M, F, FI, JI}
    f::F
    f_interval::FI
    jacobian_interval::JI
end

function _compile(expr, x, u)
    # `nanmath = false`: the default NaNMath-backed `sin`/`cos`/... do not accept intervals.
    # `cse = true`: shared sub-expressions are computed once; without it the large, deeply
    # nested expressions of a multibody normal form take minutes to compile.
    compiled = Symbolics.build_function(
        expr,
        x,
        u;
        expression = Val(false),
        nanmath = false,
        cse = true,
    )
    return compiled isa Tuple ? first(compiled) : compiled
end

function SymbolicDynamics(fsym::AbstractVector, x::AbstractVector, u::AbstractVector)
    N, M = length(x), length(u)
    length(fsym) == N ||
        throw(DimensionMismatch("f has $(length(fsym)) components but there are $N states"))
    fvec = Symbolics.Num[Symbolics.Num(e) for e in fsym]
    xs = Symbolics.Num[Symbolics.Num(e) for e in x]
    us = Symbolics.Num[Symbolics.Num(e) for e in u]
    # One compiled function serves both float and interval arguments (it is generic code).
    f_interval = _compile(fvec, xs, us)
    f = (xv, uv) -> SVector{N, Float64}(f_interval(collect(xv), collect(uv)))
    jacobian_interval = _compile(Symbolics.jacobian(fvec, xs), xs, us)
    return SymbolicDynamics{N, M, typeof(f), typeof(f_interval), typeof(jacobian_interval)}(
        f,
        f_interval,
        jacobian_interval,
    )
end

function SymbolicDynamics(system::MS.ConstrainedBlackBoxControlContinuousSystem)
    nx, nu = MS.statedim(system), MS.inputdim(system)
    xs = [Symbolics.variable(:x, i) for i in 1:nx]
    us = [Symbolics.variable(:u, i) for i in 1:nu]
    fsym = try
        collect(MS.mapping(system)(xs, us))
    catch err
        error(
            "SymbolicDynamics: could not trace the dynamics symbolically " *
            "($(sprint(showerror, err))). Build the dynamics from symbolic expressions instead.",
        )
    end
    return SymbolicDynamics(fsym, xs, us)
end

# Intervals ---------------------------------------------------------------------------------

_to_interval(v::IA.Interval{Float64}) = v
_to_interval(v::Real) = IA.interval(Float64, v)

# Plain `Vector`s (not `SizedVector`s, which a comprehension over an `SVector`'s indices gives):
# the compiled functions build their output with the container type of their first argument.
_box(lo, hi) =
    IA.Interval{Float64}[IA.interval(Float64, lo[i], hi[i]) for i in 1:length(lo)]
_point_box(u) = IA.Interval{Float64}[IA.interval(Float64, u[i]) for i in 1:length(u)]
_cell_box(c, r) = IA.Interval{Float64}[
    IA.interval(Float64, c[i], r[i]; format = :midpoint) for i in 1:length(c)
]

# A finite interval can still represent a partial domain, e.g. sqrt([-1, 1]).
# Symbolics may mix interval arguments with plain literals. Their stored values are the
# model coefficients, so IntervalArithmetic's NG flag for literal promotion is expected.
_valid_interval_enclosure(v) = IA.decoration(v) >= IA.dac && IA.isbounded(v)

# Form a fresh box from validated bounds. IntervalArithmetic's hull operation can
# lower the decoration even when both operands are finite and continuous.
_interval_hull(x, y) =
    IA.interval(Float64, min(IA.inf(x), IA.inf(y)), max(IA.sup(x), IA.sup(y)))

# Metzler majorant of an interval Jacobian: the sup of the diagonal (signed) and the magnitude of
# the off-diagonal entries, which is what the growth-bound ODE `ρ̇ = L ρ` requires.
function _metzler(J, ::Val{N}) where {N}
    M = Matrix{Float64}(undef, N, N)
    for j in 1:N, i in 1:N
        v = _to_interval(J[i, j])
        M[i, j] = if !_valid_interval_enclosure(v)
            Inf
        else
            i == j ? IA.sup(v) : IA.mag(v)
        end
    end
    return SMatrix{N, N, Float64}(M)
end

"""
    jacobian_bound(dynamics::SymbolicDynamics, xlow, xhigh, u) -> SMatrix

Growth-bound matrix `L` (Metzler: signed supremum on the diagonal, magnitude off the diagonal)
of the Jacobian `∂f/∂x`, bounded with interval arithmetic for all states in the box
`[xlow, xhigh]` and the input `u`. An entry is `Inf` if its interval evaluation fails
to give a finite, defined and continuous enclosure.
"""
function jacobian_bound(
    dynamics::SymbolicDynamics{N},
    xlow::AbstractVector,
    xhigh::AbstractVector,
    u::AbstractVector,
) where {N}
    return _metzler(dynamics.jacobian_interval(_box(xlow, xhigh), _point_box(u)), Val(N))
end

# ------------------------------------------------------------------------------------------
# A-priori enclosure of the flow tube
# ------------------------------------------------------------------------------------------

"""
    flow_enclosure(dynamics, cell_low, cell_high, guess_low, guess_high, u, tstep;
                   max_iterations = 10, inflation = 0.1, shrink = false)
    flow_enclosure(dynamics, cell::LazySets.AbstractHyperrectangle, guess_low, guess_high,
                   u, tstep; kwargs...)

Return the lower and upper bounds of a box containing the exact ODE trajectories of **all**
initial states in the cell throughout `[0, tstep]` with constant input `u`. Return `nothing`
if the test cannot verify an enclosure. The hyperrectangle overload rounds the source
endpoints outward from its stored centre and half-widths.

The Picard inclusion test computes `E = R + [0, tstep]·[f](B, u)` with outward rounding and
requires `E` to lie strictly inside `B`. A trajectory reaching the boundary would still lie
inside `E` by the integral equation, so it cannot leave `B`. This assumes a locally Lipschitz
vector field on a neighborhood of the accepted box and a sound interval evaluation of it.

The initial candidate is the hull of the source and the guess. On failure the hull with `E`
is inflated by `inflation` times its width, with at least one representable step of slack on
each side. With `shrink = true`, return `E` from the successful test as a tighter enclosure.
"""
function flow_enclosure(
    dynamics::SymbolicDynamics,
    cell_low,
    cell_high,
    guess_low,
    guess_high,
    u,
    tstep;
    kwargs...,
)
    length(cell_low) == length(cell_high) ||
        throw(DimensionMismatch("cell bounds must have the same length"))
    source = _box(cell_low, cell_high)
    return _flow_enclosure(dynamics, source, guess_low, guess_high, u, tstep; kwargs...)
end

function flow_enclosure(
    dynamics::SymbolicDynamics,
    cell::LazySets.AbstractHyperrectangle,
    guess_low,
    guess_high,
    u,
    tstep;
    kwargs...,
)
    source = _cell_box(LazySets.center(cell), LazySets.radius_hyperrectangle(cell))
    return _flow_enclosure(dynamics, source, guess_low, guess_high, u, tstep; kwargs...)
end

function _flow_enclosure(
    dynamics::SymbolicDynamics{N, M},
    source,
    guess_low,
    guess_high,
    u,
    tstep;
    max_iterations::Int = 10,
    inflation::Real = 0.1,
    shrink::Bool = false,
) where {N, M}
    length(source) == length(guess_low) == length(guess_high) == N ||
        throw(DimensionMismatch("cell and guess must have $N state coordinates"))
    length(u) == M || throw(DimensionMismatch("input must have $M coordinates"))
    isfinite(tstep) && tstep >= 0 ||
        throw(ArgumentError("tstep must be finite and nonnegative"))
    isfinite(inflation) && inflation >= 0 ||
        throw(ArgumentError("inflation must be finite and nonnegative"))
    # Detect the library's non-rigorous rounding mode without relying on its configuration
    # internals. The upper endpoint of this sum must round up from halfway between floats.
    one_interval = IA.interval(Float64, 1.0)
    half_ulp = IA.interval(Float64, eps(1.0) / 2)
    IA.sup(one_interval + half_ulp) >= nextfloat(1.0) ||
        throw(ArgumentError("flow_enclosure requires outward interval rounding"))

    guess = _box(guess_low, guess_high)
    time = IA.interval(Float64, 0, tstep)
    ubox = _point_box(u)
    all(_valid_interval_enclosure, source) &&
    all(_valid_interval_enclosure, guess) &&
    all(_valid_interval_enclosure, ubox) &&
    _valid_interval_enclosure(time) || return nothing
    B = _interval_hull.(source, guess)
    for _ in 1:max_iterations
        F = _to_interval.(dynamics.f_interval(B, ubox))
        all(_valid_interval_enclosure, F) || return nothing
        E = source .+ Ref(time) .* F
        all(_valid_interval_enclosure, E) || return nothing
        if all(i -> IA.isinterior(E[i], B[i]), eachindex(B))
            enclosure = shrink ? E : B
            return IA.inf.(enclosure), IA.sup.(enclosure)
        end
        for i in eachindex(B)
            lo = min(IA.inf(B[i]), IA.inf(E[i]))
            hi = max(IA.sup(B[i]), IA.sup(E[i]))
            # Inflation only proposes a new box. The next inclusion test validates it.
            padding = iszero(inflation) ? 0.0 : inflation * (hi - lo)
            inflated_lo = min(prevfloat(lo), lo - padding)
            inflated_hi = max(nextfloat(hi), hi + padding)
            isfinite(inflated_lo) && isfinite(inflated_hi) || return nothing
            B[i] = IA.interval(Float64, inflated_lo, inflated_hi)
        end
    end
    return nothing
end

# ------------------------------------------------------------------------------------------
# The over-approximation
# ------------------------------------------------------------------------------------------

"""
    LocalGrowthBoundStats

Counters kept by a [`LocalGrowthBound`](@ref) (thread safe): `calls` (cell/input pairs processed),
`enclosure_failures` (no verified flow enclosure found, global bound used), `unbounded` (the local
Jacobian bound was not finite, global bound used) and `improved` (pairs whose local radius is
strictly smaller than the global one in at least one dimension).
"""
struct LocalGrowthBoundStats
    calls::Threads.Atomic{Int}
    enclosure_failures::Threads.Atomic{Int}
    unbounded::Threads.Atomic{Int}
    improved::Threads.Atomic{Int}
end

LocalGrowthBoundStats() = LocalGrowthBoundStats(
    Threads.Atomic{Int}(0),
    Threads.Atomic{Int}(0),
    Threads.Atomic{Int}(0),
    Threads.Atomic{Int}(0),
)

"""
    LocalGrowthBound <: Dionysos.System.DiscreteTimeSystemOverApproximation

Over-approximation of the sampled dynamics (sampling period `tstep`) whose growth bound is
computed locally for every cell, see the [module documentation](@ref LocalGrowthAbstraction).

    LocalGrowthBound(system, dynamics::SymbolicDynamics, tstep; kwargs...)
    LocalGrowthBound(system::MS.ConstrainedBlackBoxControlContinuousSystem, tstep; kwargs...)
    LocalGrowthBound(fsym, x, u, X, U, tstep; kwargs...)

`system` is a Dionysos continuous-time system (state set `X` a box, input set `U`): its `f` is what
the cell centre is simulated with (RK4, `nsystem` substeps, via `ST.discretize`). The three
constructors differ in where the symbolic form used for the interval bounds comes from: given
explicitly, traced from `system`, or `fsym` itself (then `system` is built from it).

Keywords: `refinements = 1` (how many times steps 3-4 are repeated, each time starting from the
previous, smaller image), `ngrowthbound = nsystem = ST.DEFAULT_NUM_SUBSTEPS`, and the
[`flow_enclosure`](@ref) options `max_enclosure_iterations = 10`, `enclosure_inflation = 0.1`,
`shrink_enclosure = false`.

`boundary_tolerance = 1e-9` is the fraction of the cell radius removed from every image to suppress
face contacts. It can also remove true successors near a face. Set it to `0` to disable the
deflation and keep the closed-cell lookup of the computed box.

Fields of interest: `global_growth` (the Dionysos `DiscreteTimeGrowthBound` used for step 2 and as
fallback) and `stats` (see [`LocalGrowthBoundStats`](@ref)).
"""
struct LocalGrowthBound{S, G, D} <: ST.DiscreteTimeSystemOverApproximation
    system::S
    global_growth::G
    dynamics::D
    tstep::Float64
    refinements::Int
    ngrowthbound::Int
    max_enclosure_iterations::Int
    enclosure_inflation::Float64
    shrink_enclosure::Bool
    boundary_tolerance::Float64
    stats::LocalGrowthBoundStats
end

function LocalGrowthBound(
    system::MS.ConstrainedBlackBoxControlContinuousSystem,
    dynamics::SymbolicDynamics{N},
    tstep::Real;
    refinements::Int = 1,
    ngrowthbound::Int = ST.DEFAULT_NUM_SUBSTEPS,
    nsystem::Int = ST.DEFAULT_NUM_SUBSTEPS,
    max_enclosure_iterations::Int = 10,
    enclosure_inflation::Real = 0.1,
    shrink_enclosure::Bool = false,
    boundary_tolerance::Real = 1e-9,
) where {N}
    MS.statedim(system) == N ||
        throw(DimensionMismatch("system and dynamics differ in dimension"))
    refinements >= 0 || throw(ArgumentError("`refinements` must be nonnegative"))
    boundary_tolerance >= 0 ||
        throw(ArgumentError("`boundary_tolerance` must be nonnegative"))
    # Step 2's ingredient: the *global* growth bound of Dionysos, with the Jacobian bound taken
    # with interval arithmetic over the whole state set and the given input.
    X = UT._outer_box(MS.stateset(system))
    xlow, xhigh = LazySets.low(X), LazySets.high(X)
    global_L = u -> jacobian_bound(dynamics, xlow, xhigh, u)
    isfinite(sum(global_L(LazySets.center(MS.inputset(system))))) || error(
        "the global Jacobian bound over the state set is unbounded; the dynamics are not " *
        "Lipschitz on `X` (a singularity inside the box?)",
    )
    continuous_growth =
        ST.ContinuousTimeGrowthBound(system; jacobian_bound = global_L, ngrowthbound)
    global_growth = ST.discretize(continuous_growth, Float64(tstep); num_substeps = nsystem)
    return LocalGrowthBound(
        ST.get_system(global_growth),
        global_growth,
        dynamics,
        Float64(tstep),
        refinements,
        ngrowthbound,
        max_enclosure_iterations,
        Float64(enclosure_inflation),
        shrink_enclosure,
        Float64(boundary_tolerance),
        LocalGrowthBoundStats(),
    )
end

LocalGrowthBound(
    system::MS.ConstrainedBlackBoxControlContinuousSystem,
    tstep::Real;
    kwargs...,
) = LocalGrowthBound(system, SymbolicDynamics(system), tstep; kwargs...)

function LocalGrowthBound(
    fsym::AbstractVector,
    x::AbstractVector,
    u::AbstractVector,
    X,
    U,
    tstep::Real;
    kwargs...,
)
    dynamics = SymbolicDynamics(fsym, x, u)
    system = MS.ConstrainedBlackBoxControlContinuousSystem(
        dynamics.f,
        length(x),
        length(u),
        X,
        U,
    )
    return LocalGrowthBound(system, dynamics, tstep; kwargs...)
end

ST.get_system(lgb::LocalGrowthBound) = lgb.system

"""
    stats(lgb) -> NamedTuple

Snapshot of the counters of `lgb` (see [`LocalGrowthBoundStats`](@ref)).
"""
function stats(lgb::LocalGrowthBound)
    s = lgb.stats
    return (
        calls = s.calls[],
        enclosure_failures = s.enclosure_failures[],
        unbounded = s.unbounded[],
        improved = s.improved[],
    )
end

"""
    reset_stats!(lgb)

Zero the counters of `lgb`.
"""
function reset_stats!(lgb::LocalGrowthBound)
    for c in (
        lgb.stats.calls,
        lgb.stats.enclosure_failures,
        lgb.stats.unbounded,
        lgb.stats.improved,
    )
        c[] = 0
    end
    return lgb
end

# The per-input hoisting of Dionysos: the *global* radius depends on the (uniform) cell radius and
# the input only, so it is integrated once per input and reused for every cell.
ST.input_cache(lgb::LocalGrowthBound, r, u) = ST.input_cache(lgb.global_growth, r, u)

function ST.reach_set(lgb::LocalGrowthBound, cell, u, global_radius)
    c = LazySets.center(cell)
    r0 = LazySets.radius_hyperrectangle(cell)
    Fc = ST.get_system_map(lgb.global_growth)(c, u)      # centre simulation
    Threads.atomic_add!(lgb.stats.calls, 1)

    ρ = global_radius                                    # step 2: Φ(c, u) ± ρ_glob
    for _ in 1:(lgb.refinements)
        # Step 3: a larger box around the cell and the current conservative image, verified to
        # enclose the whole flow tube, and the Jacobian bounded over it.
        B = flow_enclosure(
            lgb.dynamics,
            cell,
            Fc .- ρ,
            Fc .+ ρ,
            u,
            lgb.tstep;
            max_iterations = lgb.max_enclosure_iterations,
            inflation = lgb.enclosure_inflation,
            shrink = lgb.shrink_enclosure,
        )
        if B === nothing
            Threads.atomic_add!(lgb.stats.enclosure_failures, 1)
            break
        end
        L = jacobian_bound(lgb.dynamics, B[1], B[2], u)
        if !all(isfinite, L)
            Threads.atomic_add!(lgb.stats.unbounded, 1)
            break
        end
        # Step 4: recompute the radius from the local Jacobian and keep the smaller estimate.
        ρ_local = ST.runge_kutta4((r, _) -> L * r, r0, u, lgb.tstep, lgb.ngrowthbound)
        ρ = min.(ρ, ρ_local)
    end
    any(ρ .< global_radius) && Threads.atomic_add!(lgb.stats.improved, 1)
    # Suppress face contacts, with the loss of coverage described in the constructor docs.
    ρ = max.(ρ .- lgb.boundary_tolerance .* r0, 0.0)
    return LazySets.Hyperrectangle(Fc, ρ)
end

ST.get_over_approximation_map(lgb::LocalGrowthBound) =
    (elem, u) -> ST.reach_set(
        lgb,
        elem,
        u,
        ST.input_cache(lgb, LazySets.radius_hyperrectangle(elem), u),
    )

"""
    overapproximation_map(approx::ST.DiscreteTimeSystemOverApproximation, state_grid) -> Function

The function `(cell, u) -> reach set` that the `USER_DEFINED` mode of the `UniformGridAbstraction`
optimizer takes as its `overapproximation_map`. It reproduces what the kernel does for the
approximations that Dionysos knows: `ST.input_cache` is evaluated **once per input** (for the uniform
cell radius `h/2` of `state_grid`), not once per cell, and `ST.reach_set` is then called with it. The
table of caches is shared by the threads of the build, behind a lock (a lookup is negligible next to
the reach set of one cell).

[`CertifiedRK4GrowthBound`](@ref) has its own method, which needs the grid to find the cells that
its image can quantize into.
"""
function overapproximation_map(approx::ST.DiscreteTimeSystemOverApproximation, state_grid)
    r = MP.get_h(state_grid) / 2.0
    caches = Dict{Any, Any}()
    lock = ReentrantLock()
    return (elem, u) -> begin
        cache = @lock lock get!(() -> ST.input_cache(approx, r, u), caches, u)
        return ST.reach_set(approx, elem, u, cache)
    end
end

# ------------------------------------------------------------------------------------------
# Building the transition system with Dionysos
# ------------------------------------------------------------------------------------------

"""
    build_abstraction(approx::ST.DiscreteTimeSystemOverApproximation; state_grid, input_grid,
                      print_level = 0, abstraction_file = nothing, abstraction_info = nothing,
                      merge_chunks = nothing, attributes...) -> optimizer

Discrete-time transition system (symbolic model) of the system of `approx`, computed by the
Dionysos `UniformGridAbstraction` optimizer in `USER_DEFINED` mode, with the function that
[`overapproximation_map`](@ref) builds from `approx` and `state_grid`. Its `ST.input_cache` is
computed once per input per build. For a [`LocalGrowthBound`](@ref) this caches the global starting
radius across cells. Pass `lgb.global_growth` to build the global baseline.

`state_grid` / `input_grid` are Dionysos grids (`MP.GridFree(x0, h)`); further keyword arguments
are forwarded as optimizer attributes (`execution_backend`, `use_implicit_mapping`, ...).
Returns the optimizer, which caches the abstraction so that a control problem can be set on it
afterwards:

```julia
abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
```

With `abstraction_file`, the transitions are **cached on disk**: if the file exists they are
loaded instead of computed ([`load_abstraction!`](@ref), which checks `abstraction_info`), otherwise
they are computed and written there ([`save_abstraction`](@ref)). This is how an abstraction computed
on a server is used on another machine.

With a `SY.SlurmArrayBackend` as `execution_backend`, the call computes **one chunk** of the
transitions of a Slurm array and writes it to the backend's directory (nothing is cached in the
optimizer). With `merge_chunks = (; outdir, nchunks)` the transitions are not computed but read from
those chunks: into the optimizer ([`merge_chunks!`](@ref)), or, if there is an `abstraction_file`,
straight into that file ([`merge_chunks_to_file`](@ref)), which needs far less memory and leaves the
optimizer **without** an abstraction (run again with the file to load it).

**Certified approximations.** A `CertifiedRK4GrowthBound` can use cache files and chunks. They store
its provenance (`CertifiedRK4.get_provenance`: the kernel version, the sampling period, the substeps
and a fingerprint of the dynamics, plus a fingerprint of the cells and inputs): in the cache file,
and for a chunk in a small file next to it (`provenance_part_<id>_of_<N>.jld2`, with the size and a
fingerprint of the chunk file, written once the chunk is complete and removed before any chunk is
computed). Loading or merging **refuses** anything that does not carry the same provenance: a file or
chunk computed by an uncertified approximation, by another version of the kernel, for other
dynamics, grids, state box or carving, or a chunk modified after its record was written. The
`abstraction_info` you pass should still identify what the provenance cannot see (the
specification, `state_input_filter`, ...).
"""
function build_abstraction(
    approx::ST.DiscreteTimeSystemOverApproximation;
    state_grid,
    input_grid,
    print_level::Int = 0,
    abstraction_file = nothing,
    abstraction_info = nothing,
    merge_chunks = nothing,
    attributes...,
)
    if approx isa CertifiedRK4GrowthBound
        # Cache files and chunks are allowed: they record how their transitions were computed
        # (`provenance`), and are refused unless it is this construction.
        for name in (:overapproximation_map, :approx_mode, :out_of_domain_handler)
            haskey(attributes, name) &&
                throw(ArgumentError("cannot override $name in a certified RK4 abstraction"))
        end
    end
    backend = get(attributes, :execution_backend, nothing)
    chunk_job = backend isa SY.SlurmArrayBackend && merge_chunks === nothing
    # A chunk and its provenance record go together: the record of an older chunk must not outlive
    # it, whoever overwrites the chunk file (it is written again when this one is complete).
    chunk_job && _forget_chunk_provenance(backend)
    system = ST.get_system(approx)
    problem = PR.AlternatingSimulationProblem(system, system.X)

    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("approx_mode"),
        AB.UniformGridAbstraction.USER_DEFINED,
    )
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("overapproximation_map"),
        overapproximation_map(approx, state_grid),
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), print_level)
    for (name, value) in attributes
        MOI.set(optimizer, MOI.RawOptimizerAttribute(String(name)), value)
    end
    if abstraction_file !== nothing && isfile(abstraction_file)
        print_level >= 1 &&
            @info "Loading the abstraction from $abstraction_file (not computed)"
        load_abstraction!(
            optimizer,
            abstraction_file;
            info = abstraction_info,
            print_level,
            approx,
        )
    else
        abstraction_file === nothing ||
            @warn "$abstraction_file does not exist (relative paths are relative to $(pwd())): the abstraction is built (computed or merged from chunks), then written there"
        if merge_chunks === nothing
            MOI.optimize!(optimizer)
            chunk_job && _record_chunk_provenance(backend, optimizer, approx)
            abstraction_file === nothing || save_abstraction(
                abstraction_file,
                optimizer;
                info = abstraction_info,
                approx,
            )
        elseif abstraction_file === nothing
            merge_chunks!(
                optimizer,
                merge_chunks.outdir;
                nchunks = merge_chunks.nchunks,
                print_level = max(print_level, 1),
                approx,
            )
        else
            # Straight from the chunks to the file: no automaton, hence a small fraction of the memory.
            merge_chunks_to_file(
                optimizer,
                merge_chunks.outdir,
                abstraction_file;
                nchunks = merge_chunks.nchunks,
                info = abstraction_info,
                print_level = max(print_level, 1),
                approx,
            )
        end
    end
    return optimizer
end

# ------------------------------------------------------------------------------------------
# Saving and loading the transitions (to compute them on a server)
# ------------------------------------------------------------------------------------------

# A few concrete states spread over the state ids: if the same ids are not the same cells (another
# grid or state box), a loaded abstraction would be silently wrong, so these are stored and checked.
function _id_probe(abstract_system)
    n = SY.get_n_state(abstract_system)
    ids = unique([1, (n + 1) ÷ 2, n])
    return [collect(SY.get_concrete_state(abstract_system, id)) for id in ids]
end

_peak_gb() = round(Sys.maxrss() / 2^30; digits = 1)

# ------------------------------------------------------------------------------------------
# Provenance of the stored transitions
# ------------------------------------------------------------------------------------------

# Transitions that are written to a file (a cache, or the chunk of a Slurm array) and read back later
# are only as trustworthy as the record of how they were computed. An approximation that keeps such a
# record (`provenance`) has it stored with the transitions and checked when they are read, so that
# they cannot be mixed with transitions computed another way. Dionysos itself stores no such record:
# the cache file is this module's, and a chunk gets a small file of its own next to it.

"""
    provenance(approx) -> Union{Nothing, NamedTuple}

Plain data that identifies **how** the successors of `approx` are computed: for a
[`CertifiedRK4GrowthBound`](@ref), `CertifiedRK4.get_provenance`. `nothing`, the default, means that
`approx` keeps no record and that stored transitions are not checked (the local growth bound).
"""
provenance(approx) = nothing
provenance(approx::CertifiedRK4GrowthBound) = CertifiedRK4.get_provenance(approx)

function _ids_fingerprint(abstract_system)
    fp = CertifiedRK4.fingerprint
    h = fp((
        SY.get_n_state(abstract_system),
        SY.get_n_allowed_state(abstract_system),
        SY.get_n_input(abstract_system),
    ))
    for q in SY.enum_allowed_states(abstract_system)
        h = fp(q, h)
        h = fp(SY.get_concrete_state(abstract_system, q), h)
    end
    for u in SY.enum_inputs(abstract_system)
        h = fp(u, h)
        h = fp(SY.get_concrete_input(abstract_system, u), h)
    end
    return CertifiedRK4.fingerprint_string(h)
end

"""
    abstraction_provenance(abstract_system, approx) -> Union{Nothing, NamedTuple}

What a cache file or a chunk must record so that the transitions computed for `(abstract_system,
approx)` are only read back by the same construction: [`provenance`](@ref)`(approx)`, and a
fingerprint of the cells (the retained states, with their ids) and of the inputs behind the abstract
ids. The latter catches another grid, state box or carving that happens to have the same number of
states. `nothing` if `approx` keeps no record. `abstract_system` is the whole model, not a partition.
"""
function abstraction_provenance(abstract_system, approx)
    approximation = provenance(approx)
    approximation === nothing && return nothing
    return (; approximation, ids = _ids_fingerprint(abstract_system))
end

# Error unless the provenance `found` in `what` (a cache file or a chunk) is the `expected` one.
function check_provenance(found, expected, what)
    found == expected && return nothing
    reason = if found === nothing
        "it has no provenance record: it was computed by an approximation that keeps none (not " *
        "certified), or before provenance was recorded"
    elseif expected === nothing
        "it records the provenance of a certified construction, and this abstraction's " *
        "approximation has none"
    else
        "it was computed by another construction (other dynamics, sampling period, kernel " *
        "version, grid, state box or carving)"
    end
    return error(
        "$what cannot be used for this abstraction: $reason.\n" *
        "stored:    $found\nrequested: $expected",
    )
end

_provenance_file(outdir, id, nchunks) =
    joinpath(outdir, "provenance_part_$(id)_of_$(nchunks).jld2")

function _forget_chunk_provenance(backend::SY.SlurmArrayBackend)
    backend.chunk_id === nothing && return nothing
    file = _provenance_file(backend.outdir, backend.chunk_id, backend.nchunks)
    isfile(file) && rm(file)
    return nothing
end

# The size and fingerprint of a chunk file, to bind its record to its content.
_chunk_digest(file) =
    (; bytes = filesize(file), fingerprint = CertifiedRK4.fingerprint_file(file))

# Once the chunk is complete: its provenance, and the digest of the file it describes.
function _record_chunk_provenance(backend::SY.SlurmArrayBackend, optimizer, approx)
    expected = abstraction_provenance(
        MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system")),
        approx,
    )
    expected === nothing && return nothing
    chunk = SY.slurm_chunk_file(backend.outdir, backend.chunk_id, backend.nchunks)
    file = _provenance_file(backend.outdir, backend.chunk_id, backend.nchunks)
    JLD2.jldsave(file * ".tmp"; provenance = expected, digest = _chunk_digest(chunk))
    mv(file * ".tmp", file; force = true)
    return nothing
end

# Chunk `id` must have been computed by the construction `expected` (`nothing`: one without a
# record, which a chunk with a record is not), and must not have changed since.
function _check_chunk_provenance(outdir, id, nchunks, expected)
    what = "chunk $id of $nchunks"
    file = _provenance_file(outdir, id, nchunks)
    record = isfile(file) ? JLD2.load(file) : nothing
    check_provenance(record === nothing ? nothing : record["provenance"], expected, what)
    expected === nothing && return nothing
    chunk = SY.slurm_chunk_file(outdir, id, nchunks)
    isfile(chunk) || return nothing           # reading it will say that it is missing
    record["digest"] == _chunk_digest(chunk) || error(
        "$what was modified after its provenance was recorded ($file): it is not the chunk " *
        "that the certified construction wrote",
    )
    return nothing
end

# File format: `n_state`, `n_input`, `info`, `probe`, then the transitions in `nparts` parts (one per
# chunk of a Slurm array, one in all for `save_abstraction`) as the integer vectors `target_k`,
# `source_k`, `symbol_k`, and `n_transitions`. The files of the first version (one part, keys `target`,
# `source`, `symbol`, no `nparts`) are still read. `parts` is an iterator, so that the parts of a big
# abstraction are never all in memory. `provenance` (`abstraction_provenance`) is stored under its
# own key when the approximation has one, and is absent otherwise.
function _write_abstraction(file, abstract_system, info, parts; provenance = nothing)
    mkpath(dirname(abspath(file)))
    tmp = file * ".tmp"
    nparts, n = 0, 0
    JLD2.jldopen(tmp, "w") do f
        f["n_state"] = SY.get_n_state(abstract_system)
        f["n_input"] = SY.get_n_input(abstract_system)
        f["info"] = info
        f["probe"] = _id_probe(abstract_system)
        provenance === nothing || (f["provenance"] = provenance)
        for (target, source, symbol) in parts
            nparts += 1
            f["target_$nparts"], f["source_$nparts"], f["symbol_$nparts"] =
                target, source, symbol
            n += length(target)
        end
        f["nparts"] = nparts
        f["n_transitions"] = n
        return
    end
    mv(tmp, file; force = true)   # a job killed while writing leaves no half-written file
    return n
end

"""
    save_abstraction(file, optimizer; info = nothing, approx = nothing)

Write the transitions of the abstraction cached in `optimizer` (see [`build_abstraction`](@ref)) to
`file` (JLD2, three integer vectors `target`, `source`, `symbol` per part: 12 bytes per transition). The
compiled symbolic dynamics cannot be saved, so only the transitions are: they are all the control
synthesis needs. `info` is any plain data that identifies the problem (grids, state box, `tstep`, ...):
it is stored, and [`load_abstraction!`](@ref) refuses a file whose `info` differs.

`approx` is the approximation the transitions were computed with. If it keeps a record of how it
computes them ([`provenance`](@ref), as a `CertifiedRK4GrowthBound` does), that record is stored too,
with a fingerprint of the cells, and a file without it is refused when loaded for such an
approximation (and the other way round).
"""
function save_abstraction(file::AbstractString, optimizer; info = nothing, approx = nothing)
    abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    n_state = SY.get_n_state(abstract_system)
    n_input = SY.get_n_input(abstract_system)
    n = SY.get_n_transitions(abstract_system)
    T = max(n_state, n_input) < typemax(Int32) ? Int32 : Int
    target, source, symbol = Vector{T}(undef, n), Vector{T}(undef, n), Vector{T}(undef, n)
    k = 0
    for (t, s, u) in SY.enum_transitions(abstract_system)
        k += 1
        target[k], source[k], symbol[k] = t, s, u
    end
    k == n || error("the abstraction has $k transitions, expected $n")
    _write_abstraction(
        file,
        abstract_system,
        info,
        ((target, source, symbol),);
        provenance = abstraction_provenance(abstract_system, approx),
    )
    return file
end

"""
    merge_chunks!(optimizer, outdir; nchunks, print_level = 0, approx = nothing) -> abstract_system

Fill the abstraction of `optimizer` (built with the same problem, grids and attributes as the
chunks, but not computed) with the transitions of the `nchunks` files that the jobs of a Slurm array
wrote in `outdir` (each one computed by [`build_abstraction`](@ref) with
`execution_backend = SY.SlurmArrayBackend(nchunks, chunk_id, outdir, :contiguous, true)`), and cache
it in the optimizer. Errors if a chunk is missing or was computed for another number of chunks,
states or inputs, or (if `approx` keeps a record of how it computes the transitions,
`provenance`) by another construction than `approx`: dynamics, sampling period, kernel,
grids, state box or carving. All the transitions end up in an automaton: for a big abstraction use
[`merge_chunks_to_file`](@ref) instead.
"""
function merge_chunks!(
    optimizer,
    outdir::AbstractString;
    nchunks::Int,
    print_level::Int = 0,
    approx = nothing,
)
    abstract_system =
        AB.UniformGridAbstraction.build_empty_abstraction!(optimizer.abstraction_solver)
    expected = abstraction_provenance(abstract_system, approx)
    for id in 1:nchunks
        _check_chunk_provenance(outdir, id, nchunks, expected)
    end
    SY.merge_transition_chunks!(abstract_system, outdir; nchunks, print_level)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("abstract_system"), abstract_system)
    return abstract_system
end

"""
    merge_chunks_to_file(optimizer, outdir, file; nchunks, info = nothing, print_level = 0,
                         approx = nothing) -> n

Write the transitions of the `nchunks` chunks in `outdir` (see [`merge_chunks!`](@ref)) to `file`
([`save_abstraction`](@ref) format) **without building the abstraction in memory**: the chunks are read
and written one after the other, so the memory is that of one chunk (about 24 bytes per transition of
it), instead of the hundreds of bytes per transition of the automaton. Returns the number of
transitions. Errors if a chunk is missing or was computed for another number of chunks, states or
inputs, or by another construction than `approx` (see [`merge_chunks!`](@ref)); the provenance of
`approx`, if it has one, is stored in the file. The optimizer is left without an abstraction:
[`load_abstraction!`](@ref) fills it from `file`.
With `print_level ≥ 1` every chunk is logged with the peak memory of the process so far.
"""
function merge_chunks_to_file(
    optimizer,
    outdir::AbstractString,
    file::AbstractString;
    nchunks::Int,
    info = nothing,
    print_level::Int = 0,
    approx = nothing,
)
    abstract_system =
        AB.UniformGridAbstraction.build_empty_abstraction!(optimizer.abstraction_solver)
    expected = abstraction_provenance(abstract_system, approx)
    n_state, n_input = SY.get_n_state(abstract_system), SY.get_n_input(abstract_system)
    T = max(n_state, n_input) < typemax(Int32) ? Int32 : Int
    seen = Ref(0)
    function part(id)
        _check_chunk_provenance(outdir, id, nchunks, expected)
        data = try
            SY.read_transition_chunk(outdir, id, nchunks)
        catch e
            error(
                "chunk $id of $nchunks cannot be read ($(SY.slurm_chunk_file(outdir, id, nchunks))): " *
                "missing, or its task was killed while writing it: " *
                sprint(showerror, e),
            )
        end
        data.chunk_id == id ||
            error("chunk $id: its file says it is chunk $(data.chunk_id)")
        data.nchunks == nchunks ||
            error("chunk $id was computed for $(data.nchunks) chunks, not $nchunks")
        (data.nstates, data.ninputs) == (n_state, n_input) || error(
            "chunk $id has $(data.nstates) states and $(data.ninputs) inputs, the problem has " *
            "$n_state and $n_input: it was computed for another grid or state box",
        )
        n = length(data.transitions)
        target, source, symbol =
            Vector{T}(undef, n), Vector{T}(undef, n), Vector{T}(undef, n)
        for (k, (t, s, u)) in enumerate(data.transitions)
            target[k], source[k], symbol[k] = t, s, u
        end
        seen[] += n
        print_level >= 1 && @info(
            "Merging chunk $id / $nchunks",
            transitions = n,
            total = seen[],
            peak_memory_GB = _peak_gb(),
        )
        return (target, source, symbol)
    end
    # The optimizer must not keep the empty abstraction: it would look computed.
    optimizer.abstraction_solver.abstract_system = nothing
    n = _write_abstraction(
        file,
        abstract_system,
        info,
        (part(id) for id in 1:nchunks);
        provenance = expected,
    )
    print_level >= 1 && @info(
        "Merged $nchunks chunks to $file",
        transitions = n,
        size_GB = round(filesize(file) / 2^30; digits = 2),
        peak_memory_GB = _peak_gb(),
    )
    return n
end

"""
    load_abstraction!(optimizer, file; info = nothing, print_level = 0, approx = nothing) -> abstract_system

Fill the abstraction of `optimizer` (built with the same problem, grids and attributes as the one that
was saved, but not computed) with the transitions of `file` written by [`save_abstraction`](@ref) or
[`merge_chunks_to_file`](@ref), and cache it in the optimizer so that a control problem can be set
and solved on it. Errors if `info`, the number of states and inputs or the cells behind the state ids
differ from the saved ones, or if the provenance stored in the file is not that of `approx` (see
[`save_abstraction`](@ref)): in particular a file computed without a record of how, by an
uncertified approximation, is refused for a certified one. With `print_level ≥ 1` the loading is
logged part by part, with the peak memory of the process.
"""
function load_abstraction!(
    optimizer,
    file::AbstractString;
    info = nothing,
    print_level::Int = 0,
    approx = nothing,
)
    return JLD2.jldopen(file, "r") do f
        f["info"] == info || error(
            "$file was computed for another problem.\nsaved:     $(f["info"])\nrequested: $info",
        )
        abstract_system =
            AB.UniformGridAbstraction.build_empty_abstraction!(optimizer.abstraction_solver)
        check_provenance(
            haskey(f, "provenance") ? f["provenance"] : nothing,
            abstraction_provenance(abstract_system, approx),
            file,
        )
        (SY.get_n_state(abstract_system), SY.get_n_input(abstract_system)) ==
        (f["n_state"], f["n_input"]) || error(
            "$file has $(f["n_state"]) states and $(f["n_input"]) inputs, the problem has " *
            "$(SY.get_n_state(abstract_system)) and $(SY.get_n_input(abstract_system))",
        )
        _id_probe(abstract_system) == f["probe"] ||
            error("$file numbers the cells differently: another state grid or state box")
        nparts = haskey(f, "nparts") ? f["nparts"] : 0       # 0: the first version of the format
        total = haskey(f, "n_transitions") ? f["n_transitions"] : nothing
        print_level >= 1 &&
            @info("Loading $file", parts = max(nparts, 1), transitions = total)
        loaded = 0
        # In blocks: a single vector of all the `(target, source, symbol)` tuples (24 bytes each) would
        # be a needless peak for a big abstraction.
        block = 1 << 22
        transitions = SY.TransitionKey[]
        for k in 1:max(nparts, 1)
            suffix = nparts == 0 ? "" : "_$k"
            target, source, symbol =
                f["target" * suffix], f["source" * suffix], f["symbol" * suffix]
            for first in 1:block:length(target)
                last = min(first + block - 1, length(target))
                resize!(transitions, last - first + 1)
                for j in first:last
                    transitions[j - first + 1] =
                        (Int(target[j]), Int(source[j]), Int(symbol[j]))
                end
                SY.add_transitions!(abstract_system, transitions)
            end
            loaded += length(target)
            print_level >= 1 && @info(
                "Loaded part $k / $(max(nparts, 1))",
                transitions = loaded,
                of = total,
                peak_memory_GB = _peak_gb(),
            )
        end
        MOI.set(optimizer, MOI.RawOptimizerAttribute("abstract_system"), abstract_system)
        return abstract_system
    end
end

end # module
