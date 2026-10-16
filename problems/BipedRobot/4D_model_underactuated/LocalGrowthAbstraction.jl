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

Step 3 is made rigorous: the growth-bound theorem needs `L_B` to bound the Jacobian on a convex
set containing the trajectories of **all** points of `R` during the **whole** sampling period, not
just at its end points. So `B` is verified, with an interval Picard–Lindelöf check, to be an
a-priori enclosure of that tube (`R + [0, τ]·f(B, u) ⊆ B`), and is enlarged until it is. If no such
box is found the cell falls back to the global bound, which is always valid. Both images are sound,
so their intersection `Φ(c, u) ± min(ρ_glob, ρ_loc)` is too, and the local image is never worse
than the global one.

The result, [`LocalGrowthBound`](@ref), is a `Dionysos.System.DiscreteTimeSystemOverApproximation`,
so it plugs straight into Dionysos' abstraction machinery (`SY.compute_abstract_system_from_concrete_system!`,
the `UniformGridAbstraction` optimizer, every execution backend). [`build_abstraction`](@ref) wraps
that.

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

export SymbolicDynamics, LocalGrowthBound, LocalGrowthBoundStats
export jacobian_bound, flow_enclosure, build_abstraction, stats, reset_stats!
export save_abstraction, load_abstraction!, merge_chunks!

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
generic operations).
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

_to_interval(v::IA.Interval) = v
_to_interval(v::Real) = IA.interval(float(v))

# Plain `Vector`s (not `SizedVector`s, which a comprehension over an `SVector`'s indices gives):
# the compiled functions build their output with the container type of their first argument.
_box(lo, hi) =
    IA.Interval{Float64}[IA.interval(Float64(lo[i]), Float64(hi[i])) for i in 1:length(lo)]
_point_box(u) = IA.Interval{Float64}[IA.interval(Float64(u[i])) for i in 1:length(u)]

# Metzler majorant of an interval Jacobian: the sup of the diagonal (signed) and the magnitude of
# the off-diagonal entries, which is what the growth-bound ODE `ρ̇ = L ρ` requires.
function _metzler(J, ::Val{N}) where {N}
    M = Matrix{Float64}(undef, N, N)
    for j in 1:N, i in 1:N
        v = _to_interval(J[i, j])
        M[i, j] = if !IA.isbounded(v)
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
`[xlow, xhigh]` and the input `u`. May contain `Inf` if the interval evaluation is unbounded.
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

A box `B` that provably contains the trajectory `x(t)`, `t ∈ [0, tstep]`, of **every** initial
state in the cell `[cell_low, cell_high]` under the constant input `u` (or `nothing` if none was
found).

The proof is the interval Picard–Lindelöf test: if `E = R + [0, τ]·f(B, u) ⊆ B` then no
trajectory leaves `B` before `τ`. The search starts from `guess` (the user's "larger cell":
the hull of the cell and its conservative image); when the test fails the box is replaced by the
hull of itself and `E`, inflated by `inflation` times its width, and the test is repeated up to
`max_iterations` times. With `shrink = true` the (tighter, still valid) `E` is returned instead of
the verified guess.

The endpoint arithmetic is plain floating point, not outward rounded.
"""
function flow_enclosure(
    dynamics::SymbolicDynamics,
    cell_low,
    cell_high,
    guess_low,
    guess_high,
    u,
    tstep;
    max_iterations::Int = 10,
    inflation::Real = 0.1,
    shrink::Bool = false,
)
    lo = min.(cell_low, guess_low)
    hi = max.(cell_high, guess_high)
    ubox = _point_box(u)
    for _ in 1:max_iterations
        F = dynamics.f_interval(_box(lo, hi), ubox)
        flo = [IA.inf(_to_interval(v)) for v in F]
        fhi = [IA.sup(_to_interval(v)) for v in F]
        (all(isfinite, flo) && all(isfinite, fhi)) || return nothing
        Elo = cell_low .+ tstep .* min.(0.0, flo)
        Ehi = cell_high .+ tstep .* max.(0.0, fhi)
        if all(Elo .>= lo) && all(Ehi .<= hi)
            return shrink ? (Elo, Ehi) : (lo, hi)
        end
        nlo, nhi = min.(lo, Elo), max.(hi, Ehi)
        w = nhi .- nlo
        lo, hi = nlo .- inflation .* w, nhi .+ inflation .* w
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

`boundary_tolerance = 1e-9` is the fraction of the cell radius by which the radius of every image
is deflated, so that an image that merely *touches* a neighbouring cell along a face does not
count as overlapping it. This matters whenever some coordinate evolves exactly on the grid, e.g.
`q̇ = u` with `τ u` a multiple of the cell size: its image is exactly one cell, and with closed
cells (the `OUTER` lookup of Dionysos has `tol = 0`) it would otherwise be reported as three. It is
the half-open-cell convention of grid abstractions; the tolerance is orders of magnitude above
floating-point noise and orders of magnitude below any actual growth, but it does mean that an
image protruding into a neighbour by less than `boundary_tolerance · r` is ignored. Use `0` for
the strict closed-cell behaviour.

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
    cell_low, cell_high = LazySets.low(cell), LazySets.high(cell)
    Fc = ST.get_system_map(lgb.global_growth)(c, u)      # centre simulation
    Threads.atomic_add!(lgb.stats.calls, 1)

    ρ = global_radius                                    # step 2: Φ(c, u) ± ρ_glob
    for _ in 1:(lgb.refinements)
        # Step 3: a larger box around the cell and the current conservative image, verified to
        # enclose the whole flow tube, and the Jacobian bounded over it.
        B = flow_enclosure(
            lgb.dynamics,
            cell_low,
            cell_high,
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
        # Step 4: the image with the local growth bound. Both radii are valid, take the smaller.
        ρ_local = ST.runge_kutta4((r, _) -> L * r, r0, u, lgb.tstep, lgb.ngrowthbound)
        ρ = min.(ρ, ρ_local)
    end
    any(ρ .< global_radius) && Threads.atomic_add!(lgb.stats.improved, 1)
    # Half-open cells: an image touching a neighbour along a face does not overlap it.
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

# ------------------------------------------------------------------------------------------
# Building the transition system with Dionysos
# ------------------------------------------------------------------------------------------

"""
    build_abstraction(approx::ST.DiscreteTimeSystemOverApproximation; state_grid, input_grid,
                      print_level = 0, abstraction_file = nothing, abstraction_info = nothing,
                      merge_chunks = nothing, attributes...) -> optimizer

Discrete-time transition system (symbolic model) of the system of `approx`, computed by the
Dionysos `UniformGridAbstraction` optimizer in `USER_DEFINED` mode with the reach sets of
`approx`: any `DiscreteTimeSystemOverApproximation` works, in particular a [`LocalGrowthBound`](@ref)
and its `global_growth` (the baseline to compare with).

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
those chunks ([`merge_chunks!`](@ref)).
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
        ST.get_over_approximation_map(approx),
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), print_level)
    for (name, value) in attributes
        MOI.set(optimizer, MOI.RawOptimizerAttribute(String(name)), value)
    end
    if abstraction_file !== nothing && isfile(abstraction_file)
        print_level >= 1 &&
            @info "Loading the abstraction from $abstraction_file (not computed)"
        load_abstraction!(optimizer, abstraction_file; info = abstraction_info)
    else
        abstraction_file === nothing ||
            @warn "$abstraction_file does not exist (relative paths are relative to $(pwd())): the abstraction is built (computed or merged from chunks), then written there"
        if merge_chunks === nothing
            MOI.optimize!(optimizer)
        else
            merge_chunks!(
                optimizer,
                merge_chunks.outdir;
                nchunks = merge_chunks.nchunks,
                print_level = max(print_level, 1),
            )
        end
        abstraction_file === nothing ||
            save_abstraction(abstraction_file, optimizer; info = abstraction_info)
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

"""
    save_abstraction(file, optimizer; info = nothing)

Write the transitions of the abstraction cached in `optimizer` (see [`build_abstraction`](@ref)) to
`file` (JLD2, three integer vectors `target`, `source`, `symbol`: 12 bytes per transition). The
compiled symbolic dynamics cannot be saved, so only the transitions are: they are all the control
synthesis needs. `info` is any plain data that identifies the problem (grids, state box, `tstep`, ...):
it is stored, and [`load_abstraction!`](@ref) refuses a file whose `info` differs.
"""
function save_abstraction(file::AbstractString, optimizer; info = nothing)
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
    mkpath(dirname(abspath(file)))
    tmp = file * ".tmp"
    JLD2.jldsave(
        tmp;
        n_state,
        n_input,
        info,
        probe = _id_probe(abstract_system),
        target,
        source,
        symbol,
    )
    mv(tmp, file; force = true)   # a job killed while writing leaves no half-written file
    return file
end

"""
    merge_chunks!(optimizer, outdir; nchunks, print_level = 0) -> abstract_system

Fill the abstraction of `optimizer` (built with the same problem, grids and attributes as the
chunks, but not computed) with the transitions of the `nchunks` files that the jobs of a Slurm array
wrote in `outdir` (each one computed by [`build_abstraction`](@ref) with
`execution_backend = SY.SlurmArrayBackend(nchunks, chunk_id, outdir, :contiguous, true)`), and cache
it in the optimizer. Errors if a chunk is missing or was computed for another number of chunks,
states or inputs.
"""
function merge_chunks!(
    optimizer,
    outdir::AbstractString;
    nchunks::Int,
    print_level::Int = 0,
)
    abstract_system =
        AB.UniformGridAbstraction.build_empty_abstraction!(optimizer.abstraction_solver)
    SY.merge_transition_chunks!(abstract_system, outdir; nchunks, print_level)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("abstract_system"), abstract_system)
    return abstract_system
end

"""
    load_abstraction!(optimizer, file; info = nothing) -> abstract_system

Fill the abstraction of `optimizer` (built with the same problem, grids and attributes as the one that
was saved, but not computed) with the transitions of `file` written by [`save_abstraction`](@ref),
and cache it in the optimizer so that a control problem can be set and solved on it. Errors if
`info`, the number of states and inputs or the cells behind the state ids differ from the saved ones.
"""
function load_abstraction!(optimizer, file::AbstractString; info = nothing)
    data = JLD2.load(file)
    data["info"] == info || error(
        "$file was computed for another problem.\nsaved:     $(data["info"])\nrequested: $info",
    )
    abstract_system =
        AB.UniformGridAbstraction.build_empty_abstraction!(optimizer.abstraction_solver)
    (SY.get_n_state(abstract_system), SY.get_n_input(abstract_system)) ==
    (data["n_state"], data["n_input"]) || error(
        "$file has $(data["n_state"]) states and $(data["n_input"]) inputs, the problem has " *
        "$(SY.get_n_state(abstract_system)) and $(SY.get_n_input(abstract_system))",
    )
    _id_probe(abstract_system) == data["probe"] ||
        error("$file numbers the cells differently: another state grid or state box")
    target, source, symbol = data["target"], data["source"], data["symbol"]
    transitions = Vector{SY.TransitionKey}(undef, length(target))
    for k in eachindex(transitions)
        transitions[k] = (Int(target[k]), Int(source[k]), Int(symbol[k]))
    end
    SY.add_transitions!(abstract_system, transitions)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("abstract_system"), abstract_system)
    return abstract_system
end

end # module
