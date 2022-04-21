mutable struct ProbaAutomaton{S<:AbstractSet{Tuple{Int,Int,Int,Float64}}} <: HybridSystems.AbstractAutomaton
    nstates::Int #not used
    nsymbols::Int
    transitions::S
end

function print_automaton(autom::ProbaAutomaton)
    println()
    println("nstates: ",autom.nstates)
    println("nsymbols: ",autom.nsymbols)
    println("n_transitions: ",ntransitions(autom))
    println(autom.transitions.data)
    println()
end
function ProbaAutomaton{S}(nstates, nsymbols) where {S}
    transitions = S()
    return ProbaAutomaton(nstates, nsymbols, transitions)
end
NewProbaAutomaton(nstates, nsymbols) = ProbaAutomaton{SortedTupleSet{4,Tuple{Int,Int,Int,Float64}}}(nstates, nsymbols)

function HybridSystems.ntransitions(autom::ProbaAutomaton)
    return length(autom.transitions)
end

# Assumes not add twice same transition !!!
function HybridSystems.add_transition!(autom::ProbaAutomaton, source, symbol, target, proba)
    push_new!(autom.transitions, (target, source, symbol, proba))
end

function comparison(t1::Tuple{Int,Int,Int,Float64},t2::Tuple{Int,Int})
    return t1[2:3] == t2
end
function delete_transition_post!(autom::ProbaAutomaton, source, symbol)
    delete!(autom.transitions, (source, symbol), comparison)
end

# translist is an iterable of Tuple{Int,Int,Int}
append_new!(s::Set, translist) = union!(s, translist)
function add_transitions!(autom::ProbaAutomaton, translist)
    append_new!(autom.transitions, translist)
end
Base.empty!(autom::ProbaAutomaton) = empty!(autom.transitions)

function compute_post!(targetlist, autom::ProbaAutomaton, source, symbol)
    fix_and_eliminate_tail!(targetlist, autom.transitions, (source, symbol))
end

function pre(autom::ProbaAutomaton, target)
    return fix_and_eliminate_first(autom.transitions, target)
end

function post(autom::ProbaAutomaton, source)
    translist = []
    for el in autom.transitions.data
        if el[2] == source
            push!(translist,(el[1],el[3],el[4]))
        end
    end
    return translist
end

function post(autom::ProbaAutomaton, source, symbol)
    targetlist = []
    fix_and_eliminate_tail!(targetlist, autom.transitions, (source, symbol))
    return targetlist
end

function HybridSystems.add_state!(autom::ProbaAutomaton)
    autom.nstates += 1
end

function get_transitions()

end
