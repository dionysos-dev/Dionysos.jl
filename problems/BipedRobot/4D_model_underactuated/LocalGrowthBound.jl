"""
    LocalGrowthBound

What `Dionysos.System.CellwiseBound` needs to be used with dynamics given symbolically: a growth bound
**per cell** for the abstraction of a system, instead of one for the whole state space.

The classical growth bound of Dionysos (`ST.ContinuousTimeGrowthBound`) uses one Jacobian bound
`L(u)` valid on the whole state space `X`: every cell of radius `r` is mapped to
`Φ(c, u) ± ρ(r, u)` with `ρ̇ = L(u) ρ`, `ρ(0) = r`. Whenever the dynamics are only strongly
non-linear in a small part of `X`, `L` is dictated by the worst point and the image of every other
cell is needlessly large, which creates spurious transitions. A `CellwiseBound` replaces `L(u)` by
a bound over an *enclosure of the step from the cell*. Dionysos asks the caller for two functions;
this module supplies them for symbolic dynamics ([`cellwise_bound`](@ref)):

1. `enclosure(cell, u, τ)`: a box `B` that contains the trajectories of **all** points of the cell
   during the **whole** sampling period, not just at its end points. It is verified with an
   interval Picard–Lindelöf check (`R + [0, τ]·f(B, u) ⊆ interior(B)` with outward rounding,
   [`picard_enclosure`](@ref)), and the cell gets no transition under `u` if none is found;
2. `bound(B, u)`: the Metzler matrix `L_B` of the Jacobian `∂f/∂x` over `B`, bounded with interval
   arithmetic ([`jacobian_bound_over_box`](@ref)): signed supremum on the diagonal, magnitude off it.

Both need `f` and `∂f/∂x` evaluated on intervals, which is what [`SymbolicDynamics`](@ref) compiles.

Dionysos then simulates the centre with RK4, integrates `ρ̇ = L_B ρ` with RK4 and finds the cells
that the box `Φ(c, u) ± ρ` touches. The tube check applies to the exact ODE flow, whereas the
centre and the radius are integrated with RK4: this does not certify successor coverage for the
exact flow or for the RK4 map (for the RK4 map, see `CertifiedRK4GrowthBound` in `BipedAbstraction`).

The module builds no abstraction: give `cellwise_bound(dynamics)` as the `jacobian_bound` of
`ST.ContinuousTimeGrowthBound`, or of the `GROWTH` mode of the `UniformGridAbstraction` optimizer.
Nothing in Dionysos itself had to change for it.

# Example

```julia
using Symbolics, IntervalArithmetic
include("LocalGrowthBound.jl"); using .LocalGrowthBound

@variables x[1:2] u[1:1]
f = [x[2], -9.81 * sin(x[1]) - 0.5 * x[2] + u[1]]          # pendulum
X = LazySets.Hyperrectangle(; low = [-3.0, -6.0], high = [3.0, 6.0])
U = LazySets.Hyperrectangle(; low = [-2.0], high = [2.0])

dynamics = SymbolicDynamics(f, collect(x), collect(u))
system = MS.ConstrainedBlackBoxControlContinuousSystem(dynamics.f, 2, 1, X, U)
growth = ST.discretize(
    ST.ContinuousTimeGrowthBound(system; jacobian_bound = cellwise_bound(dynamics)),
    0.1,
)
cell = LazySets.Hyperrectangle(SVector(1.5, 0.0), SVector(0.1, 0.2))
ST.get_over_approximation_map(growth)(cell, SVector(0.0))   # the image of the cell
```
"""
module LocalGrowthBound

import Dionysos
const ST = Dionysos.System

import MathematicalSystems as MS
import LazySets
import Symbolics
import IntervalArithmetic as IA
using StaticArrays

export SymbolicDynamics, FlowEnclosure
export jacobian_bound_over_box, picard_enclosure, cellwise_bound

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
    jacobian_bound_over_box(dynamics::SymbolicDynamics, xlow, xhigh, u) -> SMatrix

Growth-bound matrix `L` (Metzler: signed supremum on the diagonal, magnitude off the diagonal)
of the Jacobian `∂f/∂x`, bounded with interval arithmetic for all states in the box
`[xlow, xhigh]` and the input `u`. An entry is `Inf` if its interval evaluation fails
to give a finite, defined and continuous enclosure.
"""
function jacobian_bound_over_box(
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
    picard_enclosure(dynamics, cell_low, cell_high, guess_low, guess_high, u, tstep;
                   max_iterations = 10, inflation = 0.1, shrink = false)
    picard_enclosure(dynamics, cell::LazySets.AbstractHyperrectangle, guess_low, guess_high,
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
function picard_enclosure(
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
    return _picard_enclosure(dynamics, source, guess_low, guess_high, u, tstep; kwargs...)
end

function picard_enclosure(
    dynamics::SymbolicDynamics,
    cell::LazySets.AbstractHyperrectangle,
    guess_low,
    guess_high,
    u,
    tstep;
    kwargs...,
)
    source = _cell_box(LazySets.center(cell), LazySets.radius_hyperrectangle(cell))
    return _picard_enclosure(dynamics, source, guess_low, guess_high, u, tstep; kwargs...)
end

function _picard_enclosure(
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
        throw(ArgumentError("picard_enclosure requires outward interval rounding"))

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
# The cell-wise bound of Dionysos, supplied for symbolic dynamics
# ------------------------------------------------------------------------------------------

"""
    FlowEnclosure

What [`cellwise_bound`](@ref) returns as the enclosure of the step from a cell: the box `[low, high]`
(an a-priori enclosure of the exact trajectories from the cell, see [`picard_enclosure`](@ref)) and the
bound `L` of the Jacobian over it ([`jacobian_bound_over_box`](@ref)), which is what
`Dionysos.System.CellwiseBound` asks of its `bound` function.
"""
struct FlowEnclosure{V, M}
    low::V
    high::V
    L::M
end

"""
    cellwise_bound(dynamics::SymbolicDynamics; max_iterations = 10, inflation = 0.1,
                   shrink = false) -> Dionysos.System.CellwiseBound

The two functions of a `Dionysos.System.CellwiseBound`, for the symbolic dynamics `dynamics`:

* `enclosure(cell, u, tstep)`: a box that contains the exact trajectories from `cell` under the
  constant input `u` over `[0, tstep]` ([`picard_enclosure`](@ref), starting from the cell itself; the
  keywords are its options), with the Jacobian bound over it, as a [`FlowEnclosure`](@ref). It is
  `nothing` if no box could be verified or if the Jacobian is not bounded over it: Dionysos then
  disables the input for that cell;
* `bound(E, u)`: the Metzler matrix of the Jacobian over `E` (signed supremum on the diagonal,
  magnitude off it), which Dionysos integrates, `ṙ = L r`, to obtain the radius of the image.

The growth bound is then taken over the enclosure of the step from each cell, not over the whole
state set. Everything else (the RK4 simulation of the centre, the integration of the radius, the
lookup of the successor cells) is done by Dionysos. Use it as the `jacobian_bound` of
`ST.ContinuousTimeGrowthBound`, or of the `GROWTH` mode of the `UniformGridAbstraction` optimizer.
"""
function cellwise_bound(
    dynamics::SymbolicDynamics;
    max_iterations::Int = 10,
    inflation::Real = 0.1,
    shrink::Bool = false,
)
    function enclosure(cell, u, tstep)
        c, r = LazySets.center(cell), LazySets.radius_hyperrectangle(cell)
        box = picard_enclosure(
            dynamics,
            cell,
            c .- r,
            c .+ r,
            u,
            tstep;
            max_iterations,
            inflation,
            shrink,
        )
        box === nothing && return nothing
        L = jacobian_bound_over_box(dynamics, box[1], box[2], u)
        all(isfinite, L) || return nothing
        return FlowEnclosure(box[1], box[2], L)
    end
    return ST.CellwiseBound(enclosure, (E, u) -> E.L)
end

end # module
