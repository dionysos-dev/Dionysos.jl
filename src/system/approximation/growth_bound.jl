# --------------------------------------------------
#  GROWTH OVERAPPROXIMATION IMPLEMENTATION
# --------------------------------------------------

"""
    JacobianBoundPrecision

How much of the state and input space a single evaluation of a derived Jacobian bound has to
cover at once. Every level is rigorous; they differ in tightness and in how often the bound is
recomputed. See [`compute_jacobian_bound`](@ref).

- `GLOBAL_BOUND` — one matrix for all of `X` and `U`.
- `INPUT_BOUND` — over all of `X`, re-derived per input.
- `REGIONWISE_BOUND` — over one region of `X` at a time, per (region, input).
"""
@enum JacobianBoundPrecision GLOBAL_BOUND INPUT_BOUND REGIONWISE_BOUND

"""
    RegionwiseBound(bound, region_of, nregions)

A Jacobian bound that is **piecewise constant over the state space**: `bound(i, u)` is the
matrix valid throughout region `i`, and `region_of(x)` says which region a point falls in.

A bound covering the whole state set must be as large as its worst point. On a nonlinear
system that is usually far larger than any one region needs — `|-(g/l)cos x₁|` is `g/l` over a
full turn but near zero on most slices of it — so splitting the state space buys tightness,
and tightness is what decides whether synthesis succeeds at a given grid size.

Crucially the regions are *few and fixed*, not one per cell: the radius can still be
integrated once per (region, input) and hoisted out of the cell loop by
[`input_cache`](@ref), leaving only an array lookup on the per-cell hot path.
"""
struct RegionwiseBound{F, R}
    bound::F        # (region::Int, u) -> SMatrix
    region_of::R    # x -> region::Int
    nregions::Int
end

"""
    RegionwiseGrowth(radius, region_of, nregions)

The growth-bound map produced by a [`RegionwiseBound`](@ref): `radius(i, r, u[, tstep])` is the
inflated radius for a cell of radius `r` sitting in region `i`.
"""
struct RegionwiseGrowth{F, R}
    radius::F
    region_of::R
    nregions::Int
end

"""
    DiscreteTimeGrowthBound <: DiscreteTimeSystemOverApproximation

A discrete-time overapproximation based on **growth bounds**.

Given a system and a `growthbound_map`, this approximation inflates the center trajectory by a radius that depends on the current state set's size and the input.

# Fields
- `system`: A `ConstrainedBlackBoxControlDiscreteSystem` from `MathematicalSystems.jl`.
- `growthbound_map`: A function  
    `f(radius::SVector, u::SVector) -> SVector`  
    that computes how uncertainty in state evolves under the system.

# Overapproximation Map
Returns a function of the form:
    `f(rect::LazySets.AbstractHyperrectangle, u::SVector) -> LazySets.Hyperrectangle`
This function simulates the image of the center and inflates it using the computed growth bound.
"""
struct DiscreteTimeGrowthBound{S <: MS.ConstrainedBlackBoxControlDiscreteSystem, F} <:
       DiscreteTimeSystemOverApproximation
    system::S
    growthbound_map::F
end

get_system(approx::DiscreteTimeGrowthBound) = approx.system

# The cell radius is uniform across the grid, so a regionwise bound is still hoistable: the
# radius is integrated once per region here, not once per cell, and the per-cell work drops to
# a region lookup and an array index.
input_cache(approx::DiscreteTimeGrowthBound, r, u) =
    _input_cache(approx.growthbound_map, r, u)

_input_cache(g, r, u) = g(r, u)
_input_cache(g::RegionwiseGrowth, r, u) = [g.radius(i, r, u) for i in 1:(g.nregions)]

# For a plain map the cache *is* the radius; for a regionwise one it is the table to index.
_radius_from_cache(g, cache, x) = cache
_radius_from_cache(g::RegionwiseGrowth, cache, x) = cache[g.region_of(x)]

function reach_set(approx::DiscreteTimeGrowthBound, elem, u, Fr)
    x = LazySets.center(elem)
    Fx = get_system_map(approx)(x, u)
    return LazySets.Hyperrectangle(Fx, _radius_from_cache(approx.growthbound_map, Fr, x))
end

function get_over_approximation_map(approx::DiscreteTimeGrowthBound)
    return (elem, u) -> begin
        r = LazySets.radius_hyperrectangle(elem)
        return reach_set(approx, elem, u, input_cache(approx, r, u))
    end
end

"""
    ContinuousTimeGrowthBound <: ContinuousTimeSystemOverApproximation

A continuous-time overapproximation based on **growth bounds** for reachable set propagation.

It estimates how uncertainty evolves through time using a `growthbound_map` which depends on the radius, input, and time step.

# Fields
- `system`: A `ConstrainedBlackBoxControlContinuousSystem` from `MathematicalSystems.jl`.
- `growthbound_map`: A function  
    `f(radius::SVector, u::SVector, tstep::Real) -> SVector`  
    that estimates how uncertainty grows over a time step.

# Overapproximation Map
Returns a function of the form:
    `f(rect::LazySets.AbstractHyperrectangle, u::SVector, tstep::Real) -> LazySets.Hyperrectangle`
This function simulates the image of the center and inflates it using the computed growth bound.
"""
struct ContinuousTimeGrowthBound{S <: MS.ConstrainedBlackBoxControlContinuousSystem, F} <:
       ContinuousTimeSystemOverApproximation
    system::S
    growthbound_map::F
end

get_system(approx::ContinuousTimeGrowthBound) = approx.system
function get_over_approximation_map(approx::ContinuousTimeGrowthBound)
    return (rect, u, tstep) -> begin
        x = LazySets.center(rect)
        r = LazySets.radius_hyperrectangle(rect)
        Fx = get_system_map(approx)(x, u, tstep)
        g = approx.growthbound_map
        Fr =
            g isa RegionwiseGrowth ? g.radius(g.region_of(x), r, u, tstep) : g(r, u, tstep)
        return LazySets.Hyperrectangle(Fx, Fr)
    end
end
function discretize(
    approx::ContinuousTimeGrowthBound,
    tstep::Float64;
    num_substeps::Int = DEFAULT_NUM_SUBSTEPS,
)
    discretized_system =
        discretize_continuous_system(get_system(approx), tstep; num_substeps = num_substeps)
    g = approx.growthbound_map
    discretized_growthbound_map = if g isa RegionwiseGrowth
        RegionwiseGrowth((i, r, u) -> g.radius(i, r, u, tstep), g.region_of, g.nregions)
    else
        (r, u) -> g(r, u, tstep)
    end
    return DiscreteTimeGrowthBound(discretized_system, discretized_growthbound_map)
end

"""
    ContinuousTimeGrowthBound(system; jacobian_bound = nothing, ngrowthbound = DEFAULT_NUM_SUBSTEPS)

Build a continuous-time growth-bound overapproximation from a Jacobian bound: the
radius dynamics `ṙ = jacobian_bound(u) * r` are integrated with `ngrowthbound` RK4
substeps per time step. When `jacobian_bound` is not provided, it is derived from
the system via `compute_jacobian_bound` (requires an extension providing it).
"""
function ContinuousTimeGrowthBound(
    system::MS.ConstrainedBlackBoxControlContinuousSystem;
    jacobian_bound = nothing,
    ngrowthbound::Int = DEFAULT_NUM_SUBSTEPS,
)
    if jacobian_bound === nothing
        jacobian_bound = compute_jacobian_bound(system)
    end
    if jacobian_bound isa RegionwiseBound
        # `L` differs from region to region, so the radius ODE is integrated per region — a
        # handful of integrations, not one per cell.
        growthbound_map = RegionwiseGrowth(
            (i, r, u, tstep) -> runge_kutta4(
                (rr, uu) -> jacobian_bound.bound(i, uu) * rr,
                r,
                u,
                tstep,
                ngrowthbound,
            ),
            jacobian_bound.region_of,
            jacobian_bound.nregions,
        )
        return ContinuousTimeGrowthBound(system, growthbound_map)
    end
    modified_jacobian_bound = (r, u) -> jacobian_bound(u) * r
    growthbound_map =
        (r, u, tstep) -> runge_kutta4(modified_jacobian_bound, r, u, tstep, ngrowthbound)
    return ContinuousTimeGrowthBound(system, growthbound_map)
end

# Deprecated: use `ContinuousTimeGrowthBound(system; jacobian_bound, ngrowthbound)`.
function ContinuousTimeGrowthBound_from_jacobian_bound(
    system::MS.ConstrainedBlackBoxControlContinuousSystem,
    jacobian_bound;
    ngrowthbound::Int = DEFAULT_NUM_SUBSTEPS,
)
    return ContinuousTimeGrowthBound(
        system;
        jacobian_bound = jacobian_bound,
        ngrowthbound = ngrowthbound,
    )
end

"""
    compute_jacobian_bound(system; precision = INPUT_BOUND, nsplit = 4)

Derive a Jacobian bound from the dynamics instead of writing one by hand: the extension traces
`f` symbolically and bounds each entry with interval arithmetic, which is a proof rather than a
sample. Requires `Symbolics` to be loaded.

A hand-written bound is still worth having — it can exploit structure the interval arithmetic
cannot see, and costs nothing at runtime — but it is also the single input whose being wrong
invalidates every guarantee downstream, silently. This is the alternative.

`precision` (a [`JacobianBoundPrecision`](@ref)) trades tightness against work:

| | state ranged over | re-derived per | returns |
| :--- | :--- | :--- | :--- |
| `GLOBAL_BOUND` | all of `X`, all of `U` | never | `u -> SMatrix` |
| `INPUT_BOUND` | all of `X` | input | `u -> SMatrix` |
| `REGIONWISE_BOUND` | one of `nsplit^n` sub-boxes of `X` | (region, input) | [`RegionwiseBound`](@ref) |

`GLOBAL_BOUND` is enough when the Jacobian barely varies. `INPUT_BOUND` — the default — is free
in the abstraction's hot loop, because the radius map is hoisted out of the cell loop by
[`input_cache`](@ref). `REGIONWISE_BOUND` stays hoisted too: the regions are few and fixed, so
the radius is integrated once per (region, input) and the per-cell work is a lookup. Reach for
it when a bound taken over the whole state space is too conservative for synthesis to succeed —
on a pendulum it cut spurious transitions by 5–10% at `nsplit ≥ 8`.

Note `nsplit^n` grows with the state dimension `n`, so lower `nsplit` above 2–3 states.
"""
function compute_jacobian_bound(system; kwargs...)
    return error(
        "Automatic Jacobian-bound computation requires Symbolics.jl " *
        "(load it with `using Symbolics`), or provide the bound explicitly via " *
        "`ContinuousTimeGrowthBound(system; jacobian_bound = jacobian_bound)`.",
    )
end

# ------------------------------------------------------------
# Perturbed systems
# ------------------------------------------------------------
# Reissig, Weber & Rungger, "Feedback Refinement Relations for the Synthesis of Symbolic
# Controllers" (IEEE TAC 2017) define the growth bound for the perturbed system
# `ẋ ∈ f(x, u) + [-w̄, w̄]` from the start: the nominal trajectory is simulated at the centre of
# the disturbance set and the radius dynamics gain a constant, `ṙ = L(u)·r + z`. The disturbance
# is folded into the successor set here, at construction — it never becomes part of the input
# alphabet, and it costs no extra reach-set computations: the same one nominal and one radius
# integration per (cell, input) as the unperturbed kernel.
#
# The folding produces a *nominal* system plus a disturbed radius map, so everything downstream —
# `discretize`, the abstraction loop, the per-input hoisting — sees the types it always saw, and
# the unperturbed path is untouched.

function _fold_noise(system, noise_bound)
    W = MS.noiseset(system)
    box = LazySets.box_approximation(W)
    w_c = LazySets.center(box)

    z = if noise_bound !== nothing
        noise_bound
    else
        # The additive reading `ẋ = f(x, u) + w` — the RWR'17 setting — where the disturbance's
        # effect on the radius is exactly the disturbance set's half-width. A non-additive
        # `f(x, u, w)` needs `zᵢ ≥ sup |fᵢ(x, u, w) − fᵢ(x, u, w_c)|`, which cannot be read off
        # a black box; it is derived symbolically when the extension is loaded
        # (`compute_noise_bound`), and must be supplied otherwise.
        if MS.noisedim(system) == MS.statedim(system)
            LazySets.radius_hyperrectangle(box)
        else
            compute_noise_bound(system)
        end
    end
    return w_c, z
end

"""
    ContinuousTimeGrowthBound(system::MS.NoisyConstrainedBlackBoxControlContinuousSystem; jacobian_bound = nothing, noise_bound = nothing, ngrowthbound = DEFAULT_NUM_SUBSTEPS)

The perturbed growth bound of Reissig–Weber–Rungger: the nominal trajectory at the centre `w_c`
of the disturbance set, and the radius dynamics `ṙ = jacobian_bound(u)·r + z`.

`z` bounds the disturbance's effect per state dimension. When `noise_bound` is not given it is
the half-width of the disturbance set's box approximation, which is exact for an additive
disturbance `ẋ = f(x, u) + w` and requires `noisedim == statedim`; a non-additive disturbance
must supply it. A `jacobian_bound` given explicitly must hold over `w ∈ W`, not only at `w_c` —
for an additive disturbance the two coincide, which is why deriving it from the nominal system
is sound there.

Robust abstraction at the cost of nominal abstraction: `w` is never enumerated, and the folded
kernel runs the same two integrations per (cell, input) as the unperturbed one.
"""
function ContinuousTimeGrowthBound(
    system::MS.NoisyConstrainedBlackBoxControlContinuousSystem;
    jacobian_bound = nothing,
    noise_bound = nothing,
    ngrowthbound::Int = DEFAULT_NUM_SUBSTEPS,
)
    w_c, z = _fold_noise(system, noise_bound)

    f = system.f
    nominal = MS.ConstrainedBlackBoxControlContinuousSystem(
        (x, u) -> f(x, u, w_c),
        MS.statedim(system),
        MS.inputdim(system),
        MS.stateset(system),
        MS.inputset(system),
    )

    if jacobian_bound === nothing
        # Derived from the *noisy* system, so the extension bounds ∂f/∂x over w ∈ W as well —
        # for a non-additive disturbance the Jacobian sees w, and a bound taken at w_c alone
        # would under-approximate it.
        jacobian_bound = compute_jacobian_bound(system)
    end
    jacobian_bound isa RegionwiseBound && error(
        "A regionwise Jacobian bound with a disturbance is not implemented; use an input-only " *
        "bound, or fold the disturbance into a hand-written growthbound_map.",
    )

    disturbed_radius_dynamics = (r, u) -> jacobian_bound(u) * r .+ z
    growthbound_map =
        (r, u, tstep) -> runge_kutta4(disturbed_radius_dynamics, r, u, tstep, ngrowthbound)
    return ContinuousTimeGrowthBound(nominal, growthbound_map)
end

function ContinuousTimeGrowthBound(
    ::MS.NoisyConstrainedBlackBoxControlContinuousSystem,
    growthbound_map,
)
    return error(
        "A hand-written growthbound_map owns the whole over-approximation, disturbance " *
        "included — build the nominal system yourself and make the map account for every " *
        "w ∈ W, or drop the map and let the keyword constructor fold the disturbance.",
    )
end

"""
    DiscreteTimeGrowthBound(system::MS.NoisyConstrainedBlackBoxControlDiscreteSystem, growthbound_map; noise_bound = nothing)

The discrete-time twin of the perturbed growth bound: for `x⁺ ∈ f(x, u) + [-w̄, w̄]`, the nominal
map at the centre of the disturbance set and the inflated radius `growthbound_map(r, u) .+ z`.
"""
function DiscreteTimeGrowthBound(
    system::MS.NoisyConstrainedBlackBoxControlDiscreteSystem,
    growthbound_map;
    noise_bound = nothing,
)
    w_c, z = _fold_noise(system, noise_bound)

    f = system.f
    nominal = MS.ConstrainedBlackBoxControlDiscreteSystem(
        (x, u) -> f(x, u, w_c),
        MS.statedim(system),
        MS.inputdim(system),
        MS.stateset(system),
        MS.inputset(system),
    )

    return DiscreteTimeGrowthBound(nominal, (r, u) -> growthbound_map(r, u) .+ z)
end

"""
    compute_noise_bound(system; precision = INPUT_BOUND)

Derive the per-dimension disturbance bound `zᵢ ≥ sup_{x, u, w ∈ W} |fᵢ(x, u, w) − fᵢ(x, u, w_c)|`
of a noisy system symbolically, instead of writing it by hand: the extension traces the
difference, simplifies it — for an additive disturbance the dynamics cancel and the bound is
exact — and bounds what remains with interval arithmetic, which is a proof rather than a guess.
Requires Symbolics.jl (`using Symbolics`); dependency-heavy expressions may bound loosely, and
an explicit `noise_bound` always overrides.
"""
function compute_noise_bound(system; kwargs...)
    return error(
        "Automatic disturbance-bound computation requires Symbolics.jl " *
        "(load it with `using Symbolics`), or provide the bound explicitly via " *
        "`noise_bound = z` with zᵢ ≥ sup over x, u, w ∈ W of |fᵢ(x, u, w) − fᵢ(x, u, w_c)|.",
    )
end
