# ============================================================
# Optimal control with bounded input variation.
#
# Consecutive inputs along the closed loop must satisfy d(u⁻, u) ≤ Δ — an input
# slew-rate constraint. For velocity-controlled systems this is an acceleration
# limit, which is what makes the synthesized profile trackable by the low-level
# motor controllers.
#
# This is the classical turn-restricted shortest path (Caldwell 1961): the
# dynamic program lives on the *line graph* — nodes are (state, input) pairs,
# not states — exactly like the Δu state augmentation of MPC. Collapsing the
# value table onto states alone loses paths whose upstream compatibility needs
# a locally more expensive input, so the value function here is keyed on pairs
# and the controller keeps the previous input as memory.
# ============================================================

import DataStructures: PriorityQueue, dequeue_pair!

"""
    AbstractInputVariation

An input slew-rate constraint, consulted only as `_compatible(constraint, u⁻, u)`.

Two forms: [`BoundedInputVariation`](@ref), written by the user over concrete
inputs and holding a distance *function*, and [`TabulatedInputVariation`](@ref),
the plain-data form over abstract input symbols that a synthesized controller can
actually carry through serialization.
"""
abstract type AbstractInputVariation end

"""
    BoundedInputVariation(input_distance, max_variation; target_input = nothing, initial_input = nothing)

Input slew-rate constraint for [`compute_bounded_input_variation_controller`](@ref):
consecutive inputs must satisfy `input_distance(u⁻, u) ≤ max_variation`.
`target_input` additionally constrains the *last* input before entering the
target (e.g. the rest input, so velocities ramp down), `initial_input` the
*first* one; `nothing` leaves them free.

At the discrete level inputs are abstract symbols (`Int`); the
UniformGridAbstraction front-end accepts the same struct expressed on
*concrete* inputs and lifts it to a [`TabulatedInputVariation`](@ref) with
`SY.get_concrete_input` / `SY.get_abstract_input`.
"""
struct BoundedInputVariation{D, TU, IU} <: AbstractInputVariation
    input_distance::D
    max_variation::Float64
    target_input::TU
    initial_input::IU
end

function BoundedInputVariation(
    input_distance,
    max_variation::Real;
    target_input = nothing,
    initial_input = nothing,
)
    return BoundedInputVariation(
        input_distance,
        Float64(max_variation),
        target_input,
        initial_input,
    )
end

_compatible(constraint::BoundedInputVariation, u_prev::Int, u::Int) =
    constraint.input_distance(u_prev, u) <= constraint.max_variation

# `nothing` boundary input = unconstrained.
_compatible(constraint::BoundedInputVariation, ::Nothing, u::Int) = true

# Why this type exists, at length (kept out of the docstring, which renders into
# an API page already at Documenter's size limit):
#
# A synthesized controller keeps a reference to its constraint, so whatever the
# constraint holds ends up inside every serialized controller. A
# `BoundedInputVariation` holds `input_distance` -- a function, and after lifting
# an anonymous closure over the abstract system. JLD2 cannot reconstruct a closure
# type in a fresh session, which is where a deployed controller is loaded: it
# substitutes a non-callable placeholder, and the controller then throws from
# `is_defined` while `output_control` keeps answering. That is a controller still
# emitting commands having silently lost its domain check.
#
# The constraint is only ever consulted as `_compatible(constraint, u⁻, u)` on
# input symbols, so tabulating it loses nothing. 625² bits = 49 KB for the biped.
"""
    TabulatedInputVariation(compatible; target_input = nothing, initial_input = nothing)

Plain-data form of [`BoundedInputVariation`](@ref) over *abstract* input symbols:
`compatible[u⁻, u]` says whether playing `u` after `u⁻` is allowed. Unlike a
distance function, it survives serialization, so a controller holding it can be
saved and reloaded.

Built by `lift_bounded_input_variation`; users write a
[`BoundedInputVariation`](@ref) on concrete inputs and never construct this.
"""
struct TabulatedInputVariation{TU, IU} <: AbstractInputVariation
    compatible::BitMatrix
    target_input::TU
    initial_input::IU
end

function TabulatedInputVariation(
    compatible::AbstractMatrix{Bool};
    target_input = nothing,
    initial_input = nothing,
)
    size(compatible, 1) == size(compatible, 2) ||
        error("the compatibility table must be square, got $(size(compatible))")
    return TabulatedInputVariation(BitMatrix(compatible), target_input, initial_input)
end

_compatible(constraint::TabulatedInputVariation, u_prev::Int, u::Int) =
    constraint.compatible[u_prev, u]

_compatible(constraint::TabulatedInputVariation, ::Nothing, u::Int) = true

# Per-state candidate inputs, cheapest first, flat (CSR): state `q` owns
# `inputs[offsets[q] : offsets[q+1]-1]`. Internal to this solver -- the
# controller captures it, so it is also what a saved controller carries.
#
# A `Dict{Int, Vector{Tuple{Int, Float64}}}` costs a heap object per state plus
# the dict, and keeps a cost per candidate that nothing reads once the lists are
# sorted. Both go: 743 MB -> a flat pair of arrays on the 4-D biped.
# `inputs` is Int32: abstract input symbols index the model's input mapping, so
# they are far below 2^31, and this array is the dominant allocation. `select`
# widens the one symbol it returns, so nothing downstream sees the narrow type --
# unlike `ControlTable`, which hands out its storage wholesale and so has to stay
# `Int` for `get_coord_by_state`.
struct _RankedCandidates
    offsets::Vector{Int32}
    inputs::Vector{Int32}
end

function (C::_RankedCandidates)(q::Int)
    # Out of range reads as "no candidate", matching a state the table was never
    # sized for rather than throwing.
    1 <= q < length(C.offsets) || return view(C.inputs, 1:0)
    return @inbounds view(C.inputs, C.offsets[q]:(C.offsets[q + 1] - 1))
end

"""
    compute_bounded_input_variation_controller(
        autom::SY.AbstractAutomatonList,
        target_set,
        constraint::AbstractInputVariation;
        initial_set = SY.enum_states(autom),
        safe_set = nothing,
        cost_function = nothing,
    )

Optimal reach(-avoid) controller under the input slew-rate `constraint`
(consecutive inputs at most `max_variation` apart in `input_distance`).

Dijkstra on the pair graph: the value of `(q, u)` is the optimal cost-to-target
from `q` when playing `u` now, subject to every later consecutive pair — and
the final input, when `target_input` is set — being compatible. The returned
controller is *dynamic* — its memory is the previously played input — and at `q`
with memory `u⁻` it plays the compatible input of least value.

The controller holds on to `constraint`, so it is serializable exactly when the
constraint is — pass a [`TabulatedInputVariation`](@ref), not a distance function.

Only **deterministic** automata are supported for now — the exact-lattice
abstractions this constraint is designed for are deterministic by construction.

Returns `(controller, controllable_set, uncontrollable_set, value_fun_tab)`
with `value_fun_tab[q] = min_u value(q, u)` (unconstrained start at `q`).
"""
function compute_bounded_input_variation_controller(
    autom::SY.AbstractAutomatonList,
    target_set,
    constraint::AbstractInputVariation;
    initial_set = SY.enum_states(autom),
    safe_set = nothing,
    cost_function = nothing,
)
    SY.is_deterministic(autom) || error(
        "compute_bounded_input_variation_controller requires a deterministic automaton; " *
        "for nondeterministic abstractions the worst-case pair-graph fixed point is not " *
        "implemented yet.",
    )

    nstates = SY.get_n_state(autom)
    cost = cost_function === nothing ? ((q, u) -> 1.0) : cost_function
    safe_bits = _safe_bits(safe_set, nstates)

    value = Dict{Tuple{Int, Int}, Float64}()
    pq = PriorityQueue{Tuple{Int, Int}, Float64}()

    # Boundary condition on the last edge of every path: it must be compatible
    # with `target_input` (when set) and start from a safe state.
    for q_target in _safe_targets(target_set, safe_bits)
        for (q, u) in SY.pre(autom, q_target)
            safe_bits[q] || continue
            _compatible(constraint, constraint.target_input, u) || continue
            c = cost(q, u)
            if c < get(value, (q, u), Inf)
                value[(q, u)] = c
                pq[(q, u)] = c
            end
        end
    end

    target_bits = _bitset_from_states(target_set, nstates)
    init_states = collect(initial_set)
    init_bits = _bitset_from_states(init_states, nstates)
    init_covered = falses(nstates)
    num_init_uncovered = count(q -> !target_bits[q], init_states)

    settled = Set{Tuple{Int, Int}}()
    while !isempty(pq) && num_init_uncovered > 0
        (pair, v) = dequeue_pair!(pq)
        pair in settled && continue
        push!(settled, pair)
        q′, u′ = pair

        if init_bits[q′] &&
           !init_covered[q′] &&
           _compatible(constraint, constraint.initial_input, u′)
            init_covered[q′] = true
            num_init_uncovered -= 1
        end

        # Backward relaxation: play `u` at `q`, land (deterministically) on q′,
        # then continue with u′ — allowed only when the pair is compatible.
        for (q, u) in SY.pre(autom, q′)
            safe_bits[q] || continue
            target_bits[q] && continue # paths restart at the target boundary
            _compatible(constraint, u, u′) || continue
            candidate = cost(q, u) + v
            if candidate < get(value, (q, u), Inf)
                value[(q, u)] = candidate
                pq[(q, u)] = candidate
            end
        end
    end

    # Candidate inputs per state, cheapest first — the controller's lookup table,
    # and the whole size of a saved controller. Stored flat (CSR) rather than as a
    # `Dict{Int, Vector{Tuple{Int, Float64}}}`: that costs a heap-allocated inner
    # vector per state plus the dict itself, 743 MB of a 745 MB controller on the
    # 4-D biped. The cost only orders each list — `select` below discards it — so
    # it is dropped once the ordering is baked in, halving what remains.
    counts = zeros(Int32, nstates)
    for ((q, _), _) in value
        counts[q] += 1
    end
    offsets = Vector{Int32}(undef, nstates + 1)
    let k = Int32(1)
        for q in 1:nstates
            offsets[q] = k
            k += counts[q]
        end
        offsets[nstates + 1] = k
    end
    total = Int(offsets[nstates + 1]) - 1
    cand_inputs = Vector{Int32}(undef, total)
    cand_values = Vector{Float64}(undef, total)   # temporary: ordering only
    let fill_at = copy(offsets)
        for ((q, u), v) in value
            i = fill_at[q]
            cand_inputs[i] = u
            cand_values[i] = v
            fill_at[q] = i + Int32(1)
        end
    end
    for q in 1:nstates
        lo, hi = Int(offsets[q]), Int(offsets[q + 1]) - 1
        hi > lo || continue
        perm = sortperm(view(cand_values, lo:hi))
        cand_inputs[lo:hi] = cand_inputs[lo:hi][perm]
    end
    candidates = _RankedCandidates(offsets, cand_inputs)

    # Deterministic selection shared by output and memory update: cheapest
    # compatible input at `q` given the previously played input (0 = none yet).
    select = function (q::Int, u_prev::Int)
        target_bits[q] && return nothing # reached: the controller is done
        for u32 in candidates(q)
            u = Int(u32)   # widen here so callers only ever see `Int` symbols
            if u_prev == 0
                _compatible(constraint, constraint.initial_input, u) && return u
            else
                _compatible(constraint, u_prev, u) && return u
            end
        end
        return nothing
    end

    controller = ST.DiscreteDynamicController(
        0, # memory: previously played input, 0 before the first step
        ST.PredicateDomain(((mem, q),) -> select(q, mem) !== nothing),
        (mem, q) -> select(q, mem),
        (mem, q) -> select(q, mem),
        false,
    )

    value_fun_tab = fill(Inf, nstates)
    for ((q, _), v) in value
        value_fun_tab[q] = min(value_fun_tab[q], v)
    end
    for q in _safe_targets(target_set, safe_bits)
        value_fun_tab[q] = 0.0
    end

    # A state is controllable when the controller can actually start there:
    # some finite-value input compatible with `initial_input` (or the state is
    # already in the target).
    controllable_set = Set{Int}()
    for q in 1:nstates
        if target_bits[q] && safe_bits[q]
            push!(controllable_set, q)
            # `candidates(q)` is empty for a state with no finite-value input, so
            # `any` already answers false — no membership test needed.
        elseif any(
            u -> _compatible(constraint, constraint.initial_input, Int(u)),
            candidates(q),
        )
            push!(controllable_set, q)
        end
    end
    uncontrollable_set = setdiff(Set(SY.enum_states(autom)), controllable_set)

    return controller, controllable_set, uncontrollable_set, value_fun_tab
end
