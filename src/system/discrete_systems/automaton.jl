abstract type AbstractAutomatonList{N, M} <: HybridSystems.AbstractAutomaton end

# === Required Interface ===
function get_n_state(autom::AbstractAutomatonList{N, M}) where {N, M} end
function get_n_input(autom::AbstractAutomatonList{N, M}) where {N, M} end

# The transition enumeration should return a list of tuples (target, source, symbol)
function enum_transitions(autom::AbstractAutomatonList{N, M}) where {N, M} end
function add_transition!(
    autom::AbstractAutomatonList{N, M},
    source::Int,
    target::Int,
    symbol::Int,
) where {N, M} end
function pre(autom::AbstractAutomatonList{N, M}, target::Int) where {N, M} end
function post(autom::AbstractAutomatonList{N, M}, source::Int, symbol::Int) where {N, M} end
function Base.empty!(autom::AbstractAutomatonList{N, M}) where {N, M} end
function add_state!(autom::AbstractAutomatonList{N, M}) where {N, M} end

# === Common Default Implementations ===
enum_states(autom::AbstractAutomatonList) = 1:get_n_state(autom)
enum_inputs(autom::AbstractAutomatonList) = 1:get_n_input(autom)

function HybridSystems.ntransitions(autom::AbstractAutomatonList{N, M}) where {N, M}
    return length(enum_transitions(autom))
end

function add_transitions!(autom::AbstractAutomatonList{N, M}, translist) where {N, M}
    for (q′, q, u) in translist
        add_transition!(autom, q, q′, u)
    end
end

function is_deterministic(autom::AbstractAutomatonList{N, M}) where {N, M}
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

function nondeterminism_counts(autom::AbstractAutomatonList{N, M}) where {N, M}
    count = Dict{Tuple{Int, Int}, Int}()
    for (q′, q, u) in enum_transitions(autom)
        count[(q, u)] = get(count, (q, u), 0) + 1
    end
    return collect(values(count))
end

function count_self_loops(autom::AbstractAutomatonList{N, M}) where {N, M}
    count = 0
    for (q′, q, u) in enum_transitions(autom)
        if q′ == q
            count += 1
        end
    end
    return count
end
