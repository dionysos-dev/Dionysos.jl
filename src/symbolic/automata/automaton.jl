abstract type AbstractAutomatonList <: HybridSystems.AbstractAutomaton end

# === Required Interface ===
get_n_state(autom::AbstractAutomatonList) =
    error("implement `get_n_state` for $(typeof(autom))")
get_n_input(autom::AbstractAutomatonList) =
    error("implement `get_n_input` for $(typeof(autom))")

# Transition-tuple convention (see `TransitionKey` in metadata.jl):
#   - stored / enumerated as `(target, source, symbol)`  ← target first
#   - added as `add_transition!(autom, source, target, symbol)`  ← source first
# The ordering flip is deliberate but easy to get wrong; use the
# `transition_{target,source,symbol}` accessors on enumerated tuples.
enum_transitions(autom::AbstractAutomatonList) =
    error("implement `enum_transitions` for $(typeof(autom))")
add_transition!(autom::AbstractAutomatonList, source::Int, target::Int, symbol::Int) =
    error("implement `add_transition!` for $(typeof(autom))")
pre(autom::AbstractAutomatonList, target::Int) =
    error("implement `pre` for $(typeof(autom))")
post(autom::AbstractAutomatonList, source::Int, symbol::Int) =
    error("implement `post` for $(typeof(autom))")
Base.empty!(autom::AbstractAutomatonList) =
    error("implement `Base.empty!` for $(typeof(autom))")
# Concrete subtypes implement `HybridSystems.add_state!(autom)` (returns the new state id);
# a separate stub here would shadow that method and silently return `nothing`.

finalize!(autom::AbstractAutomatonList) = autom

# === Common Default Implementations ===
enum_states(autom::AbstractAutomatonList) = 1:get_n_state(autom)
enum_inputs(autom::AbstractAutomatonList) = 1:get_n_input(autom)

function HybridSystems.ntransitions(autom::AbstractAutomatonList)
    return length(enum_transitions(autom))
end

function add_transitions!(autom::AbstractAutomatonList, translist)
    for (q′, q, u) in translist
        add_transition!(autom, q, q′, u)
    end
end

"Append the successors of `(source, symbol)` to `targetlist` in place."
function compute_post!(targetlist, autom::AbstractAutomatonList, source::Int, symbol::Int)
    return append!(targetlist, post(autom, source, symbol))
end

function is_deterministic(autom::AbstractAutomatonList)
    seen = Dict{Tuple{Int, Int}, Int}()
    for (q′, q, u) in enum_transitions(autom)
        key = (q, u)
        seen[key] = get(seen, key, 0) + 1
        if seen[key] > 1
            return false
        end
    end
    return true
end

function nondeterminism_counts(autom::AbstractAutomatonList)
    count = Dict{Tuple{Int, Int}, Int}()
    for (q′, q, u) in enum_transitions(autom)
        count[(q, u)] = get(count, (q, u), 0) + 1
    end
    return collect(values(count))
end

function count_self_loops(autom::AbstractAutomatonList)
    count = 0
    for (q′, q, u) in enum_transitions(autom)
        if q′ == q
            count += 1
        end
    end
    return count
end
