"""
    LocalGrowthAbstraction

Tools to build the abstraction of a system with a growth bound **per cell**, on top of what Dionysos
provides (`Dionysos.System.CellwiseBound`, `ContinuousTimeGrowthBound` and the `GROWTH` mode of the
`UniformGridAbstraction` optimizer), and to compute and store it on a server.

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
   [`flow_enclosure`](@ref)), and the cell gets no transition under `u` if none is found;
2. `bound(B, u)`: the Metzler matrix `L_B` of the Jacobian `∂f/∂x` over `B`, bounded with interval
   arithmetic ([`jacobian_bound`](@ref)): signed supremum on the diagonal, magnitude off it.

Dionysos then simulates the centre with RK4, integrates `ρ̇ = L_B ρ` with RK4 and finds the cells
that the box `Φ(c, u) ± ρ` touches. The tube check applies to the exact ODE flow, whereas the
centre and the radius are integrated with RK4: this does not certify successor coverage for the
exact flow or for the RK4 map. Use [`CertifiedRK4GrowthBound`](@ref) (file `CertifiedRK4.jl`, next
to this one) to build a certified abstraction of the Float64 RK4 map.

[`build_abstraction`](@ref) runs the `UniformGridAbstraction` optimizer for such a system, and
stores or merges the transitions to compute them on a server ([`save_abstraction`](@ref),
[`load_abstraction!`](@ref), [`merge_chunks!`](@ref), [`merge_chunks_to_file`](@ref)). Nothing in
Dionysos itself had to change for it.

# Example

```julia
using Symbolics, IntervalArithmetic
include("LocalGrowthAbstraction.jl"); using .LocalGrowthAbstraction
const UGA = Dionysos.Optim.Abstraction.UniformGridAbstraction

@variables x[1:2] u[1:1]
f = [x[2], -9.81 * sin(x[1]) - 0.1 * x[2] + u[1]]          # pendulum
X = LazySets.Hyperrectangle(; low = [-π, -8.0], high = [π, 8.0])
U = LazySets.Hyperrectangle(; low = [-2.0], high = [2.0])

dynamics = SymbolicDynamics(f, collect(x), collect(u))
system = MS.ConstrainedBlackBoxControlContinuousSystem(dynamics.f, 2, 1, X, U)
opt = build_abstraction(
    system;
    state_grid = MP.GridFree(SVector(0.0, 0.0), SVector(0.1, 0.2)),
    input_grid = MP.GridFree(SVector(0.0), SVector(1.0)),
    approx_mode = UGA.GROWTH,
    jacobian_bound = cellwise_bound(dynamics),
    time_step = 0.1,
)
abstract_system = MOI.get(opt, MOI.RawOptimizerAttribute("abstract_system"))
```
"""
module LocalGrowthAbstraction

import Dionysos
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

export SymbolicDynamics, StepEnclosure, CertifiedRK4GrowthBound
export jacobian_bound, flow_enclosure, cellwise_bound, build_abstraction
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
# The cell-wise bound of Dionysos, supplied for symbolic dynamics
# ------------------------------------------------------------------------------------------

"""
    StepEnclosure

What [`cellwise_bound`](@ref) returns as the enclosure of the step from a cell: the box `[low, high]`
(an a-priori enclosure of the exact trajectories from the cell, see [`flow_enclosure`](@ref)) and the
bound `L` of the Jacobian over it ([`jacobian_bound`](@ref)), which is what
`Dionysos.System.CellwiseBound` asks of its `bound` function.
"""
struct StepEnclosure{V, M}
    low::V
    high::V
    L::M
end

"""
    cellwise_bound(dynamics::SymbolicDynamics; max_iterations = 10, inflation = 0.1,
                   shrink = false) -> Dionysos.System.CellwiseBound

The two functions of a `Dionysos.System.CellwiseBound`, for the symbolic dynamics `dynamics`:

* `enclosure(cell, u, tstep)`: a box that contains the exact trajectories from `cell` under the
  constant input `u` over `[0, tstep]` ([`flow_enclosure`](@ref), starting from the cell itself; the
  keywords are its options), with the Jacobian bound over it, as a [`StepEnclosure`](@ref). It is
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
        box = flow_enclosure(
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
        L = jacobian_bound(dynamics, box[1], box[2], u)
        all(isfinite, L) || return nothing
        return StepEnclosure(box[1], box[2], L)
    end
    return ST.CellwiseBound(enclosure, (E, u) -> E.L)
end

# ------------------------------------------------------------------------------------------
# Building the transition system with Dionysos
# ------------------------------------------------------------------------------------------

"""
    build_abstraction(system; state_grid, input_grid, print_level = 0, abstraction_file = nothing,
                      abstraction_info = nothing, merge_chunks = nothing, approx = nothing,
                      attributes...) -> optimizer
    build_abstraction(approx::CertifiedRK4GrowthBound; state_grid, input_grid, kwargs...)

Transition system (symbolic model) of `system`, a Dionysos continuous-time or discrete-time system,
computed by the `UniformGridAbstraction` optimizer. `state_grid` / `input_grid` are Dionysos grids
(`MP.GridFree(x0, h)`); the further keyword arguments are forwarded as optimizer attributes, which
choose how the successors are computed: for a growth bound per cell,

```julia
build_abstraction(
    system;                                   # continuous time
    state_grid, input_grid,
    approx_mode = UGA.GROWTH,                 # UGA = Dionysos.Optim.Abstraction.UniformGridAbstraction
    jacobian_bound = cellwise_bound(SymbolicDynamics(system)),
    time_step = 0.05,
)
```

The second method is for a [`CertifiedRK4GrowthBound`](@ref): it builds the abstraction of its
system with the function of [`overapproximation_map`](@ref), which is also what registers `approx`
for the provenance below (the attributes it sets cannot be overridden).

The optimizer is returned with the abstraction cached, so that a control problem can be set on it
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
    system::Union{
        MS.ConstrainedBlackBoxControlContinuousSystem,
        MS.ConstrainedBlackBoxControlDiscreteSystem,
    };
    state_grid,
    input_grid,
    print_level::Int = 0,
    abstraction_file = nothing,
    abstraction_info = nothing,
    merge_chunks = nothing,
    approx = nothing,
    attributes...,
)
    backend = get(attributes, :execution_backend, nothing)
    chunk_job = backend isa SY.SlurmArrayBackend && merge_chunks === nothing
    # A chunk and its provenance record go together: the record of an older chunk must not outlive
    # it, whoever overwrites the chunk file (it is written again when this one is complete).
    chunk_job && _forget_chunk_provenance(backend)
    problem = PR.AlternatingSimulationProblem(system, system.X)

    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
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

function build_abstraction(approx::CertifiedRK4GrowthBound; state_grid, attributes...)
    for name in (:overapproximation_map, :approx_mode, :out_of_domain_handler, :approx)
        haskey(attributes, name) &&
            throw(ArgumentError("cannot override $name in a certified RK4 abstraction"))
    end
    return build_abstraction(
        ST.get_system(approx);
        state_grid,
        approx_mode = AB.UniformGridAbstraction.USER_DEFINED,
        overapproximation_map = overapproximation_map(approx, state_grid),
        approx,
        attributes...,
    )
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
`approx` keeps no record and that stored transitions are not checked (any growth bound of Dionysos).
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
