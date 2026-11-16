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
#
# Non-deterministic automata are handled by the min-max version (Knuth 1977): a
# pair's value is its cost plus the worst, over its successors, of the best
# compatible continuation there.
# ============================================================

import DataStructures: BinaryMinHeap

"""
    AbstractInputVariation

A constraint on consecutive inputs, consulted as `_compatible(constraint, q, u⁻, u)`:
may `u` be played at state `q` right after `u⁻`?

Three forms: [`BoundedInputVariation`](@ref), written by the user over concrete
inputs and holding a distance *function*; [`TabulatedInputVariation`](@ref), the
plain-data form over abstract input symbols that a synthesized controller can carry
through serialization; and [`StateInputVariation`](@ref), whose check also depends
on the state, e.g. an actuator torque limit.
"""
abstract type AbstractInputVariation end

# A constraint that ignores the state answers from the inputs alone, and a
# `nothing` boundary input (no initial or final input prescribed) is unconstrained.
_compatible(c::AbstractInputVariation, ::Int, u_prev::Int, u::Int) =
    _compatible(c, u_prev, u)
_compatible(::AbstractInputVariation, ::Int, ::Nothing, ::Int) = true
_compatible(::AbstractInputVariation, ::Int, ::Int, ::Nothing) = true
_compatible(::AbstractInputVariation, ::Int, ::Nothing, ::Nothing) = true

"""
    BoundedInputVariation(input_distance, max_variation; target_input = nothing, initial_input = nothing)

Input slew-rate constraint for [`compute_bounded_input_variation_controller`](@ref):
consecutive inputs must satisfy `input_distance(u⁻, u) ≤ max_variation`.
`target_input` additionally constrains where a run may stop: it ends at a target
only if `target_input` may follow its last input (e.g. the rest input, so
velocities ramp down), and otherwise passes through. `initial_input` constrains
the *first* input; `nothing` leaves them free.

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

"""
    StateInputVariation(compatible; target_input = nothing, initial_input = nothing)

Constraint on consecutive inputs that depends on the state: `compatible(q, u⁻, u)`
says whether `u` may be played at abstract state `q` right after `u⁻`, all three
being abstract symbols. It must hold for every concrete state of `q`'s cell; the
canonical use is an actuator torque limit during the switch from `u⁻` to `u`.

`target_input` and `initial_input` play the same role as in
[`BoundedInputVariation`](@ref), and are checked through `compatible` at the state
where the motion ends or starts.

Unlike [`TabulatedInputVariation`](@ref), it holds a function, so a controller
carrying it cannot be reloaded in a fresh session.
"""
struct StateInputVariation{F, TU, IU} <: AbstractInputVariation
    compatible::F
    target_input::TU
    initial_input::IU
end

StateInputVariation(compatible; target_input = nothing, initial_input = nothing) =
    StateInputVariation(compatible, target_input, initial_input)

_compatible(c::StateInputVariation, q::Int, u_prev::Int, u::Int) =
    c.compatible(q, u_prev, u)

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

# The symbols leaving (or entering) each state, deduplicated and sorted, flat (CSR):
# state `k` owns `symbols[offsets[k] : offsets[k+1]-1]`, and `counts` holds how many
# transitions share each (state, symbol). Grouped by source, an entry is a pair and
# its count is its number of successors; grouped by target, an entry is an incoming
# symbol. `keyof` picks the grouping state from a `(target, source, symbol)` triple.
function _grouped_symbols(autom, nstates::Int, keyof, keep)
    deg = zeros(Int, nstates)
    for t in SY.enum_transitions(autom)
        keep(t) && (deg[keyof(t)] += 1)
    end
    offs = Vector{Int}(undef, nstates + 1)
    acc = 1
    for k in 1:nstates
        offs[k] = acc
        acc += deg[k]
    end
    offs[nstates + 1] = acc
    raw = Vector{Int32}(undef, acc - 1)
    fill_at = copy(offs)
    for t in SY.enum_transitions(autom)
        keep(t) || continue
        k = keyof(t)
        raw[fill_at[k]] = t[3]
        fill_at[k] += 1
    end

    offsets = Vector{Int}(undef, nstates + 1)
    symbols = Int32[]
    counts = Int32[]
    sizehint!(symbols, length(raw))
    sizehint!(counts, length(raw))
    for k in 1:nstates
        offsets[k] = length(symbols) + 1
        for s in sort!(view(raw, offs[k]:(offs[k + 1] - 1)))
            if length(symbols) >= offsets[k] && symbols[end] == s
                counts[end] += Int32(1)
            else
                push!(symbols, s)
                push!(counts, Int32(1))
            end
        end
    end
    offsets[nstates + 1] = length(symbols) + 1
    return offsets, symbols, counts
end

# Position of symbol `u` in state `k`'s slice, or 0 when `k` has no such entry.
function _entry(offsets, symbols, k::Int, u::Int)
    lo, hi = offsets[k], offsets[k + 1] - 1
    i = searchsortedfirst(view(symbols, lo:hi), u)
    return i <= hi - lo + 1 && symbols[lo + i - 1] == u ? lo + i - 1 : 0
end

# Mark `q′` resolved at value `v` for every incoming symbol that `accepts`, then
# credit the predecessor pairs reaching `q′` with one of those symbols: a pair whose
# successors are now all resolved is queued, once, at its final worst-case value.
function _resolve!(heap, s, autom, cost, q′::Int, v::Float64, accepts)
    empty!(s.newly_list)
    for e in s.in_offsets[q′]:(s.in_offsets[q′ + 1] - 1)
        s.resolved[e] && continue
        u = Int(s.in_inputs[e])
        accepts(u) || continue
        s.resolved[e] = true
        s.newly[u] = true
        push!(s.newly_list, u)
    end
    isempty(s.newly_list) && return
    for (q, u) in SY.pre(autom, q′)
        s.newly[u] || continue
        p = _entry(s.pair_offsets, s.pair_inputs, q, u)
        p == 0 && continue # a pair from an unsafe state carries no value
        s.worst[p] = max(s.worst[p], v)
        s.remaining[p] -= Int32(1)
        s.remaining[p] == 0 && push!(heap, (cost(q, u) + s.worst[p], p))
    end
    for u in s.newly_list
        s.newly[u] = false
    end
    return
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

Optimal worst-case reach(-avoid) controller under a constraint on consecutive inputs.

The search runs on the pairs `(q, u)` of the automaton. The value of a pair is the
optimal worst-case cost from `q` when playing `u` now: `V(q, u) = cost(q, u) + max
over successors q′ of W(q′, u)`, where `W(q′, u)` is `0` if `q′` is a target at which
`target_input` may follow `u`, and otherwise the least `V(q′, u′)` over the `u′`
compatible with `u` at `q′`. A run stops at the first target it reaches with such a
final input; a target it reaches too fast to stop is crossed like any safe state.

For a deterministic automaton this is Dijkstra on the line graph (Caldwell 1961);
otherwise it is its min-max generalization (Knuth 1977). The returned controller is
*dynamic*: its memory is the previously played input, and at `q` with memory `u⁻`
it plays the compatible input of least value.

The controller holds on to `constraint`, so it is serializable exactly when the
constraint is: pass a [`TabulatedInputVariation`](@ref), not a distance function.

Returns `(controller, controllable_set, uncontrollable_set, value_fun_tab)`
with `value_fun_tab[q] = min_u V(q, u)` (unconstrained start at `q`).
"""
function compute_bounded_input_variation_controller(
    autom::SY.AbstractAutomatonList,
    target_set,
    constraint::AbstractInputVariation;
    initial_set = SY.enum_states(autom),
    safe_set = nothing,
    cost_function = nothing,
)
    nstates = SY.get_n_state(autom)
    cost = cost_function === nothing ? ((q, u) -> 1.0) : cost_function
    safe_bits = _safe_bits(safe_set, nstates)
    target_bits = _bitset_from_states(target_set, nstates)

    # A run stops at a safe target only if `target_input` may follow its last input
    # (`initial_input` before the first step); elsewhere it continues.
    prev_or_initial(u_prev::Int) = u_prev == 0 ? constraint.initial_input : u_prev
    stops(q::Int, u_prev::Int) =
        target_bits[q] &&
        safe_bits[q] &&
        _compatible(constraint, q, prev_or_initial(u_prev), constraint.target_input)

    # Pairs from every safe state carry a value, targets included: a run that enters a
    # target too fast to stop leaves it again.
    pair_offsets, pair_inputs, remaining =
        _grouped_symbols(autom, nstates, t -> t[2], t -> safe_bits[t[2]])
    in_offsets, in_inputs, _ =
        _grouped_symbols(autom, nstates, t -> t[1], t -> safe_bits[t[2]])
    npairs = length(pair_inputs)
    pair_source = Vector{Int32}(undef, npairs)
    for q in 1:nstates, p in pair_offsets[q]:(pair_offsets[q + 1] - 1)
        pair_source[p] = q
    end

    search = (;
        pair_offsets,
        pair_inputs,
        remaining,
        in_offsets,
        in_inputs,
        worst = zeros(npairs),
        resolved = falses(length(in_inputs)),
        newly = falses(SY.get_n_input(autom)),
        newly_list = Int[],
    )
    value = fill(Inf, npairs)
    heap = BinaryMinHeap{Tuple{Float64, Int}}()

    # Last step of every run: entering a safe target with a final input compatible
    # with `target_input`. The other arrivals at a target are resolved later, by the
    # pairs leaving it.
    for q′ in _safe_targets(target_set, safe_bits)
        _resolve!(
            heap,
            search,
            autom,
            cost,
            q′,
            0.0,
            u -> _compatible(constraint, q′, u, constraint.target_input),
        )
    end

    init_states = collect(initial_set)
    init_bits = _bitset_from_states(init_states, nstates)
    # An initial state is covered once a run can start there: it stops at once, or a
    # settled pair there is compatible with `initial_input`.
    init_covered = falses(nstates)
    for q in init_states
        stops(q, 0) && (init_covered[q] = true)
    end
    num_init_uncovered = count(q -> !init_covered[q], unique(init_states))

    while !isempty(heap) && num_init_uncovered > 0
        v, p = pop!(heap)
        value[p] = v
        q′ = Int(pair_source[p])
        u′ = Int(pair_inputs[p])

        if init_bits[q′] &&
           !init_covered[q′] &&
           _compatible(constraint, q′, constraint.initial_input, u′)
            init_covered[q′] = true
            num_init_uncovered -= 1
        end

        # `u′` at `q′` is now settled: every input `u` arriving at `q′` from which the
        # switch to `u′` is allowed has found its best continuation.
        _resolve!(heap, search, autom, cost, q′, v, u -> _compatible(constraint, q′, u, u′))
    end

    # Candidate inputs per state, cheapest first: the controller's lookup table, and
    # the whole size of a saved controller. The values only order each list.
    cand_offsets = Vector{Int32}(undef, nstates + 1)
    cand_inputs = Int32[]
    buffer = Int[]
    for q in 1:nstates
        cand_offsets[q] = length(cand_inputs) + 1
        empty!(buffer)
        for p in pair_offsets[q]:(pair_offsets[q + 1] - 1)
            isfinite(value[p]) && push!(buffer, p)
        end
        sort!(buffer; by = p -> value[p])
        for p in buffer
            push!(cand_inputs, pair_inputs[p])
        end
    end
    cand_offsets[nstates + 1] = length(cand_inputs) + 1
    candidates = _RankedCandidates(cand_offsets, cand_inputs)

    # Shared by output and memory update: the cheapest input compatible with the
    # previously played one (0 = none yet, so `initial_input` applies).
    select = function (q::Int, u_prev::Int)
        stops(q, u_prev) && return nothing # reached: the controller is done
        for u32 in candidates(q)
            u = Int(u32)
            if u_prev == 0
                _compatible(constraint, q, constraint.initial_input, u) && return u
            else
                _compatible(constraint, q, u_prev, u) && return u
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
    for p in 1:npairs
        q = Int(pair_source[p])
        value_fun_tab[q] = min(value_fun_tab[q], value[p])
    end
    for q in _safe_targets(target_set, safe_bits)
        value_fun_tab[q] = 0.0
    end

    # A state is controllable when a run can start there: it stops at once, or some
    # finite-value input there is compatible with `initial_input`.
    controllable_set = Set{Int}()
    for q in 1:nstates
        if stops(q, 0)
            push!(controllable_set, q)
        elseif any(
            u -> _compatible(constraint, q, constraint.initial_input, Int(u)),
            candidates(q),
        )
            push!(controllable_set, q)
        end
    end
    uncontrollable_set = setdiff(Set(SY.enum_states(autom)), controllable_set)

    return controller, controllable_set, uncontrollable_set, value_fun_tab
end
