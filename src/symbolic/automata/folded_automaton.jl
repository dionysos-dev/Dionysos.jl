"""
    FoldedAutomaton{A <: AbstractAutomatonList} <: AbstractAutomatonList

Read-only view of an automaton with part of its input alphabet handed to the environment.

`fold` maps each inner symbol to the folded symbol the controller sees; every inner symbol of a
folded group is the environment's refinement of that single controller choice. The view therefore
presents `Post(q, u) = ⋃_{s : fold(s) = u} Post_inner(q, s)`, which is what turns the solvers'
existing controlled predecessor — every successor of the pair must be acceptable, some symbol must
survive — into the robust one `∃u ∀w` without touching them: the `∀w` is the union, the `∃u` is
the folded alphabet.

`FoldedAutomaton(inner)` folds the whole alphabet into one symbol — the pure-verification view,
where nothing is left to choose and the winning set is the answer.

The view is a faithful re-labelling, deliberately without deduplication: two inner transitions
that collapse onto the same `(target, source, symbol)` triple are enumerated twice, by
`enum_transitions` and by `pre` alike. The fixed-point counters initialise and decrement over the
same enumeration, so they stay consistent; deduplicating in one place but not the other would
desynchronise them.

Mutation goes through the underlying automaton, never through the view: `add_transition!`,
`HybridSystems.add_state!` and `Base.empty!` throw.

The inner symbol is the environment's move, so it is the witness a counterexample needs:
`inner_symbols(view, u)` lists the moves folded into `u`, and the transitions enumerated from
`inner(view)` name the move behind every folded edge.
"""
struct FoldedAutomaton{A <: AbstractAutomatonList} <: AbstractAutomatonList
    inner::A
    fold::Vector{Int}
    unfold::Vector{Vector{Int}}

    function FoldedAutomaton(inner::A, fold::Vector{Int}) where {A <: AbstractAutomatonList}
        length(fold) == get_n_input(inner) || throw(
            ArgumentError(
                "`fold` maps every inner symbol: expected length $(get_n_input(inner)), " *
                "got $(length(fold)).",
            ),
        )
        nfolded = isempty(fold) ? 0 : maximum(fold)
        unfold = [Int[] for _ in 1:nfolded]
        for (s, u) in enumerate(fold)
            1 <= u ||
                throw(ArgumentError("folded symbols are 1-based; got $u for symbol $s."))
            push!(unfold[u], s)
        end
        for u in 1:nfolded
            isempty(unfold[u]) && throw(
                ArgumentError(
                    "folded symbol $u has no inner symbol: an empty controller choice would " *
                    "be vacuously winning, which is exactly the unsoundness folding exists " *
                    "to avoid.",
                ),
            )
        end
        return new{A}(inner, fold, unfold)
    end
end

FoldedAutomaton(inner::AbstractAutomatonList) =
    FoldedAutomaton(inner, ones(Int, get_n_input(inner)))

"The automaton the view re-labels."
inner(a::FoldedAutomaton) = a.inner

"The environment's moves folded into controller symbol `u` — the witness set for counterexamples."
inner_symbols(a::FoldedAutomaton, u::Int) = a.unfold[u]

get_n_state(a::FoldedAutomaton) = get_n_state(a.inner)
get_n_input(a::FoldedAutomaton) = length(a.unfold)

enum_transitions(a::FoldedAutomaton) =
    ((q′, q, a.fold[s]) for (q′, q, s) in enum_transitions(a.inner))

pre(a::FoldedAutomaton, target::Int) =
    ((source, a.fold[s]) for (source, s) in pre(a.inner, target))

function compute_post!(targetlist, a::FoldedAutomaton, source::Int, symbol::Int)
    for s in a.unfold[symbol]
        compute_post!(targetlist, a.inner, source, s)
    end
    return targetlist
end

function post(a::FoldedAutomaton, source::Int, symbol::Int)
    targets = Int[]
    return compute_post!(targets, a, source, symbol)
end

const _FOLDED_READ_ONLY = "FoldedAutomaton is a read-only view; mutate the inner automaton and \
                           rebuild the view."

add_transition!(::FoldedAutomaton, ::Int, ::Int, ::Int) = error(_FOLDED_READ_ONLY)
add_transitions!(::FoldedAutomaton, translist) = error(_FOLDED_READ_ONLY)
HybridSystems.add_state!(::FoldedAutomaton) = error(_FOLDED_READ_ONLY)
Base.empty!(::FoldedAutomaton) = error(_FOLDED_READ_ONLY)

"""
    complete_with_sink(autom; groups = collect(1:get_n_input(autom))) -> (completed, sink, ncompleted)

The automaton extended so that every `(state, symbol group)` has a successor, by routing the
groups with none to a fresh absorbing `sink` state.

An empty group in an abstraction is a *missing behaviour*, and missing behaviour is precisely
what a fold must not tolerate: under `∀` a state whose environment move was dropped is vacuously
easier, so a quotient with empty pairs would report verified states it has no right to. Routing
those pairs to a sink that satisfies nothing is the pessimistic repair — the unmodelled move is
assumed to lose — which restores soundness at the price of conservatism, and `ncompleted` says
how much was repaired so the caller can report it.

`groups[s]` names the environment move symbol `s` realises. The default puts every symbol in its
own group. The distinction matters when the alphabet is finer than the environment's choices —
edge symbols `(mode, announcement)` over a lifted quotient, say: a layer that never carried some
announcement is a *structural non-edge*, not a dropped behaviour, and completing it per symbol
would poison the fold; only a **mode** with no successor at all is genuinely missing.

The sink is state `get_n_state(autom) + 1` of `completed`; it carries every symbol as a self-loop
and corresponds to no state of the original automaton, so translations back must skip it.
"""
function complete_with_sink(
    autom::AbstractAutomatonList;
    groups::Vector{Int} = collect(1:get_n_input(autom)),
)
    n = get_n_state(autom)
    m = get_n_input(autom)
    length(groups) == m ||
        throw(ArgumentError("`groups` maps every symbol: expected length $m."))
    sink = n + 1

    completed = IndexedAutomatonList(sink, m)
    for (q′, q, u) in enum_transitions(autom)
        add_transition!(completed, q, q′, u)
    end

    group_ids = unique(groups)
    ncompleted = 0
    for q in 1:n, g in group_ids
        symbols = findall(==(g), groups)
        if all(u -> isempty(post(autom, q, u)), symbols)
            add_transition!(completed, q, sink, first(symbols))
            ncompleted += 1
        end
    end
    for u in 1:m
        add_transition!(completed, sink, sink, u)
    end

    return completed, sink, ncompleted
end
