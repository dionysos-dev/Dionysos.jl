mutable struct ListGridController{S} <: Controller
    pairs::S
end

struct ListGridControl{ST,LT}
    state::ST
    label::LT
end

const ListGridControl_seed = hash("ListGridControl")
Base.hash(pair::ListGridControl, h::UInt) =
    hash(pair.state, hash(pair.label, h + ListGridControl_seed))

_controltype(::Type{<:ListGridManager}, ::Any, ::Any) = ListGridControl{ListGridSymbol,ListGridSymbol}

AddController!(mng::ListGridManager) = ListGridController(Set{controltype(mng)}())

# ---
Base.empty!(mng::ListGridManager, contr::Controller) = (empty!(contr.pairs); contr)
get_npairs(mng::ListGridManager, contr) = length(contr.pairs)
# Return an iterable of all transitions
enum_pairs(mng::ListGridManager, contr) = contr.pairs
# Return of iterable of enabled labels
function enum_enabled_labels(mng::ListGridManager, contr, state)
    LT = labeltype(mng)
    label_iter = Set{LT}()
    for pair in contr.pairs
        if pair.state === state
            push!(label_iter, pair.label)
        end
    end
    return label_iter
end
# ---

function add_controls!(mng::ListGridManager, contr::Controller,
        stateset::StateSet, labelset::LabelSet)
    pairs = contr.pairs
    for state in stateset.elems
        for label in labelset.elems
            push!(pairs, ListGridControl(state, label))
        end
    end
end

function add_symbols!(mng::ListGridManager, labelset::LabelSet,
        contr::Controller, stateset::StateSet)
    selems = stateset.elems
    lelems = labelset.elems
    for pair in contr.pairs
        if pair.state ∈ selems
            push!(lelems, pair.label)
        end
    end
end

## *

struct ReachSpec{IS,TS} <: Specification
    initset::IS
    targetset::TS
end

ReachSpec(mng::Manager, initset, targetset) = ReachSpec(initset, targetset)

# Return an iterable of all transitions whose target = `target`
function _incoming_transitions(transitions, target)
    idxlist = searchsorted(transitions, target, lt = _isless)
    return view(transitions, idxlist)
end

# Compute a table with for each pair of integers, the number of transitions with
# source and label corresponding to the associated integers
function _set_nposts(transitions, nstates, nlabels)
    nposts = zeros(UInt32, nstates, nlabels)
    for trans in transitions
        nposts[trans.source.val, trans.label.val] += 1
    end
    return nposts
end

# In the following UNDRV means "undrivable", and DRV means "drivable"

function _initialize_reach(mng::ListGridManager, autom, initset, targetset)
    ensure_sorted!(autom)
    nstates = _nsymbols(mng.XMap)
    nlabels = _nsymbols(mng.UMap)
    npostsUNDRV = _set_nposts(autom.transitions, nstates, nlabels)
    allDRV = Set{ListGridSymbol}(targetset.elems)
    currentDRV = collect(targetset.elems)
    nextDRV = ListGridSymbol[]
    ninitsDRV = 0
    telems = targetset.elems
    for state in initset.elems
        if state ∉ telems
            ninitsDRV += 1
        end
    end
    return npostsUNDRV, ninitsDRV, allDRV, currentDRV, nextDRV
end

function _add_controls_reach!(mng::ListGridManager, contr, autom, initset,
        npostsUNDRV, ninitsDRV, allDRV, currentDRV, nextDRV)
    pairs = contr.pairs
    transitions = autom.transitions
    ielems = initset.elems
    @inline decr_nposts(trans) = npostsUNDRV[trans.source.val, trans.label.val] -= 1
    while !isempty(currentDRV) && ninitsDRV > 0
        empty!(nextDRV)
        for target in currentDRV
            for trans in _incoming_transitions(transitions, target)
                source = trans.source
                source ∈ allDRV && continue
                # source was not drivable yet
                iszero(decr_nposts(trans)) || continue
                push!(allDRV, source)
                push!(nextDRV, source)
                push!(pairs, ListGridControl(source, trans.label))
                if source ∈ ielems
                    ninitsDRV -= 1
                end
            end
        end
        currentDRV, nextDRV = nextDRV, currentDRV
    end
    return iszero(ninitsDRV)
end

function add_controls!(mng::ListGridManager, contr::Controller,
        autom::Automaton, spec::ReachSpec)
    println("add_controls_reach! started")
    # TODO: try to infer whether npostsUNDRV is sparse or not,
    # and if sparse, use a Dict instead
    # 2020-09-09: I tried with a Dict and it was much slower for pathplanning.
    data = _initialize_reach(mng, autom, spec.initset, spec.targetset)
    init_covered = _add_controls_reach!(mng, contr, autom, spec.initset, data...)
    if init_covered
        println("\nadd_controls_reach! terminated with success")
    else
        println("\nadd_controls_reach! terminated without covering init set")
    end
    return init_covered
end
