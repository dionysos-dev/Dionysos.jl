mutable struct ListGridAutomaton{VT} <: Automaton
    transitions::VT
    issorted::Bool
end

struct ListGridTransition{ST,LT}
    source::ST
    label::LT
    target::ST
end

# Base.isless(symb1::ListGridSymbol, symb2::ListGridSymbol) =
#     symb1.val < symb2.val
# Base.isless(trans1::ListGridTransition, trans2::ListGridTransition) =
#     trans1.target < trans2.target ? true :
#     trans1.target > trans2.target ? false :
#     trans1.source < trans2.source ? true :
#     trans1.source > trans2.source ? false :
#     trans1.label < trans2.label

_isless(symb1::ListGridSymbol, symb2::ListGridSymbol) = symb1.val < symb2.val
_isless(trans1::ListGridTransition, trans2::ListGridTransition) = _isless(trans1.target, trans2.target)
_isless(target::ListGridSymbol, trans::ListGridTransition) = _isless(target, trans.target)
_isless(trans::ListGridTransition, target::ListGridSymbol,) = _isless(trans.target, target)

_transitiontype(::Type{<:ListGridManager}, ::Any, ::Any) = ListGridTransition{ListGridSymbol,ListGridSymbol}

AddAutomaton!(mng::ListGridManager) = ListGridAutomaton(transitiontype(mng)[], true)

# ---
Base.empty!(mng::ListGridManager, autom::Automaton) =
    (empty!(autom.transitions); autom.issorted = true; autom)
get_ntransitions(mng::ListGridManager, autom) = length(autom.transitions)
# Return an iterable of all transitions
enum_transitions(mng::ListGridManager, autom) = autom.transitions
# ---

function ensure_sorted!(autom::ListGridAutomaton)
    if !autom.issorted
        sort!(autom.transitions, lt = _isless)
        autom.issorted = true
    end
end

add_transition!(mng::ListGridManager, autom::Automaton, trans) =
    (push!(autom.transitions, trans); autom.issorted = false)

function add_transitions!(mng::ListGridManager, autom::Automaton,
        sourceset::StateSet, labelset::LabelSet, targetset::StateSet)
    transitions = autom.transitions
    for source in sourceset.elems
        for label in labelset.elems
            for target in targetset.elems
                push!(transitions, ListGridTransition(source, label, target))
            end
        end
    end
    autom.issorted = false
end

function add_symbols!(mng::ListGridManager, targetset::StateSet,
        autom::Automaton, sourceset::StateSet, labelset::LabelSet)
    selems = sourceset.elems
    lelems = labelset.elems
    telems = targetset.elems
    for trans in autom.transitions
        if trans.source ∈ selems && trans.label ∈ lelems
            push!(telems, trans.target)
        end
    end
end

function _pos_lims_autom(grid, rect)
    lbI = _ceil_trunc.((rect.lb - grid.orig)./grid.h .- 0.5)
    ubI = _floor_trunc.((rect.ub - grid.orig)./grid.h .+ 0.5)
    return HyperRectangle(lbI, ubI)
end

function add_transitions!(mng::ListGridManager, autom::Automaton,
        sourceset::StateSet, labelset::LabelSet, targetset::StateSet,
        contsys::ControlSystemGrowth)
    @assert xspacetype(mng) === xspacetype(contsys)
    was_empty = isempty(autom.transitions)
    println("add_transitions! started")
    tstep = contsys.tstep
    ntrans = 0
    XDisc = mng.XDisc
    UDisc = mng.UDisc
    XMap = mng.XMap
    UMap = mng.UMap
    #
    get_xcoord_center = cell2coord(XDisc)
    get_xcell_by_pos = pos2cell(XDisc)
    get_state_try_by_cell = cell2symbol_try(XMap)
    get_xcell_by_state = symbol2cell(XMap)
    get_ucoord_center = cell2coord(UDisc)
    get_ucell_by_label = symbol2cell(UMap)
    #
    telems = targetset.elems
    transitions = autom.transitions
    # Vector to store transitions
    translist = ListGridTransition{ListGridSymbol,ListGridSymbol}[]
    r = XDisc.grid.h/2 + contsys.measnoise

    # Updates every 1 seconds
    # Commented because it changes the number of allocations
    # @showprogress 1 "Computing symbolic control system: " (
    for label in labelset.elems
        ucell = get_ucell_by_label(label)
        u = get_ucoord_center(ucell)
        Fr = contsys.growthbound_map(r, u, contsys.tstep) + contsys.measnoise
        for source in sourceset.elems
            empty!(translist)
            xcell = get_xcell_by_state(source)
            x = get_xcoord_center(xcell)
            Fx = contsys.sys_map(x, u, tstep)
            rect = HyperRectangle(Fx - Fr, Fx + Fr)
            rectI = _pos_lims_autom(XDisc.grid, rect)
            ypost_iter = Iterators.product(_ranges(rectI)...)
            allin = true
            for ypost in ypost_iter
                ycell = get_xcell_by_pos(ypost)
                target, is_valid = get_state_try_by_cell(ycell)
                if !(is_valid && target ∈ telems)
                    allin = false
                    break
                end
                push!(translist, ListGridTransition(source, label, target))
            end
            if allin
                append!(transitions, translist)
                ntrans += length(translist)
            end
        end
    end
    # )
    autom.issorted = false
    if !was_empty
        unique!(transitions)
    end
    println("add_transitions! ended with success: $(ntrans) transitions added")
end
