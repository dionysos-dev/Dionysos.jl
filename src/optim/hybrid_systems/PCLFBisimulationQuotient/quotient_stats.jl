
# ============================================================
# Basic statistics on the bisimulation quotient
# ============================================================

num_nodes(T::PCBisimulationQuotient) = length(T.slices)
function num_slices(T::PCBisimulationQuotient)
    isempty(T.slices) && return 0
    return length(first(values(T.slices)))
end

num_states(T::PCBisimulationQuotient) = length(T.states)
num_transitions(T::PCBisimulationQuotient) = sum(length(q.next) for q in values(T.states))

function states_by_obs(T::PCBisimulationQuotient)
    d = Dict{Int, Int}()
    for q in values(T.states)
        d[q.obs] = get(d, q.obs, 0) + 1
    end
    return sort(collect(d); by = first) |> Dict
end

function states_by_slice(T::PCBisimulationQuotient)
    d = Dict{Int, Int}()
    for q in values(T.states)
        d[q.slice] = get(d, q.slice, 0) + 1
    end
    return sort(collect(d); by = first) |> Dict
end

function states_by_node(T::PCBisimulationQuotient)
    d = Dict{typeof(first(keys(T.part_ids))), Int}()
    for q in values(T.states)
        d[q.node] = get(d, q.node, 0) + 1
    end
    return d
end

function transitions_by_mode(T::PCBisimulationQuotient)
    d = Dict{Int, Int}()
    for q in values(T.states)
        for (m, _) in q.next
            d[m] = get(d, m, 0) + 1
        end
    end
    return sort(collect(d); by = first) |> Dict
end

function outgoing_degree_by_state(T::PCBisimulationQuotient)
    d = Dict{Int, Int}()
    for (qid, q) in T.states
        d[qid] = length(q.next)
    end
    return d
end

function outgoing_degree_stats(T::PCBisimulationQuotient)
    degs = [length(q.next) for q in values(T.states)]
    isempty(degs) && return Dict(:min => 0, :max => 0, :mean => 0.0, :median => 0.0)
    return Dict(
        :min => minimum(degs),
        :max => maximum(degs),
        :mean => sum(degs) / length(degs),
        :median => sort(degs)[cld(length(degs), 2)],
    )
end

"""
    state_ids_in_node(quotient, node; state_ids)

The ids among `state_ids` whose state belongs to `node`.

`state_ids` is meant to be narrowed by the caller -- a controllable set, say -- so it may name
states the quotient no longer holds; those are skipped rather than treated as an error. The
by-slice and by-observation variants that used to sit beside this one had no callers and were
removed; `git` has them if a use appears.
"""
function state_ids_in_node(T::PCBisimulationQuotient, node; state_ids = keys(T.states))
    return [qid for qid in state_ids if haskey(T.states, qid) && T.states[qid].node == node]
end

function transitions_from_mode(T::PCBisimulationQuotient, mode::Int)
    out = Vector{Tuple{Int, Int}}()
    for (qid, q) in T.states
        for (m, dst) in q.next
            if m == mode
                push!(out, (qid, dst))
            end
        end
    end
    return out
end

function reachable_target_states(T::PCBisimulationQuotient)
    targets = Set{Int}()
    for q in values(T.states)
        for (_, dst) in q.next
            push!(targets, dst)
        end
    end
    return targets
end

function deadend_states(T::PCBisimulationQuotient)
    return [q.id for q in values(T.states) if isempty(q.next)]
end

function self_loop_count(T::PCBisimulationQuotient)
    c = 0
    for (qid, q) in T.states
        for (_, dst) in q.next
            c += (dst == qid)
        end
    end
    return c
end

function num_parts(S::UT.SemiLinearSet)
    return length(S.array)
end

function num_faces(P::Poly)
    Q = UT.clean_poly(copy(P))
    return length(LazySets.constraints_list(Q))
end

function num_faces(S::UT.SemiLinearSet)
    total = 0
    for P in S.array
        try
            total += num_faces(P)
        catch
            continue
        end
    end
    return total
end

function cell_complexities(T::PCBisimulationQuotient)
    n_parts = Int[]
    n_faces = Int[]

    for q in values(T.states)
        S = q.set
        push!(n_parts, num_parts(S))
        push!(n_faces, num_faces(S))
    end

    return n_parts, n_faces
end

function bisimulation_stats(T::PCBisimulationQuotient)
    return Dict(
        :num_states => num_states(T),
        :num_transitions => num_transitions(T),
        :states_by_obs => states_by_obs(T),
        :states_by_slice => states_by_slice(T),
        :states_by_node => states_by_node(T),
        :transitions_by_mode => transitions_by_mode(T),
        :outgoing_degree_stats => outgoing_degree_stats(T),
        :deadend_states => deadend_states(T),
        :num_deadend_states => length(deadend_states(T)),
        :self_loop_count => self_loop_count(T),
    )
end

function print_bisimulation_stats(T::PCBisimulationQuotient)
    println("Bisimulation quotient statistics")
    println("--------------------------------")
    println("Number of nodes        : ", num_nodes(T))
    println("Number of slices       : ", num_slices(T))
    println("Number of states       : ", num_states(T))
    println("Number of transitions  : ", num_transitions(T))

    println("States by observation  : ", sort(collect(states_by_obs(T)); by = first))
    println("States by slice        : ", sort(collect(states_by_slice(T)); by = first))
    println("States by node         : ", sort(collect(states_by_node(T)); by = first))
    println("Transitions by mode    : ", sort(collect(transitions_by_mode(T)); by = first))

    println("Outgoing degree stats  : ", outgoing_degree_stats(T))
    println("Deadend states         : ", length(deadend_states(T)))
    println("Self-loops             : ", self_loop_count(T))
    return nothing
end
