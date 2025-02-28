"""
    LazySymbolicModel{N, M, S1, S2, A} <: GridBasedSymbolicModel{N, M}

A symbolic model using lazy abstraction where the automaton is computed only for a subset of the state space.
"""
mutable struct LazySymbolicModelList{
    N,
    M,
    S1 <: DO.GridDomainType{N},
    S2 <: DO.CustomList{M},
    A,
} <: GridBasedSymbolicModel{N, M}
    Xdom::S1
    Udom::S2
    autom::A
    xpos2int::Dict{NTuple{N, Int}, Int}
    xint2pos::Vector{NTuple{N, Int}}
    ucoord2int::Any
    uint2coord::Any
end

function LazySymbolicModelList(
    Xdom::DO.GeneralDomainList{N, DO.RectangularObstacles{NTuple{N, T1}}},
    Udom::DO.DomainType{N, T2},
) where {N, T1, T2}
    customDomainList = DO.convert_to_custom_domain(Udom)
    nu = DO.get_ncells(customDomainList)
    uint2coord = [coord for coord in DO.enum_elems(customDomainList)]
    ucoord2int =
        Dict((coord, i) for (i, coord) in enumerate(DO.enum_elems(customDomainList)))
    return symmodel = LazySymbolicModelList(
        Xdom,
        customDomainList,
        AutomatonList{Set{NTuple{3, Int}}}(0, nu),
        Dict{NTuple{N, T1}, Int}(),
        NTuple{N, T1}[],
        ucoord2int,
        uint2coord,
    )
end

get_n_state(symmodel::LazySymbolicModelList) = length(symmodel.xint2pos)
get_n_input(symmodel::LazySymbolicModelList) = length(symmodel.uint2coord)
enum_states(symmodel::LazySymbolicModelList) = 1:get_n_state(symmodel)
enum_inputs(symmodel::LazySymbolicModelList) = 1:get_n_input(symmodel)
get_state_domain(symmodel::LazySymbolicModelList) = symmodel.Xdom
get_input_domain(symmodel::LazySymbolicModelList) = symmodel.Udom

get_xpos_by_state(symmodel::LazySymbolicModelList, state) = symmodel.xint2pos[state]
function get_state_by_xpos(symmodel::LazySymbolicModelList, pos)
    id = get(symmodel.xpos2int, pos, nothing)
    created = false
    if id === nothing
        if pos in get_state_domain(symmodel)
            created = true
            push!(symmodel.xint2pos, pos)
            id = length(symmodel.xint2pos)
            symmodel.xpos2int[pos] = id
            i = HybridSystems.add_state!(symmodel.autom)
            @assert i == id
        else
            error("$pos is not in state domain $(symmodel.Xdom)")
        end
    end
    return id::Int
end
is_xpos(symmodel::LazySymbolicModelList, xpos) = haskey(symmodel.xpos2int, xpos)

get_concrete_input(symmodel::LazySymbolicModelList, input) = symmodel.uint2coord[input]
get_abstract_input(symmodel::LazySymbolicModelList, u) = symmodel.ucoord2int[u]
add_transitions!(symmodel::LazySymbolicModelList, translist) =
    add_transitions!(symmodel.autom, translist)

@recipe function f(
    symmodel::LazySymbolicModelList;
    dims = [1, 2],
    arrowsB = false,
    cost = false,
    lyap_fun = [],
)
    dom = symmodel.Xdom
    grid = DO.get_grid(dom)
    if cost
        LyapMax = max(filter(isfinite, getfield.([lyap_fun...], :second))...)
        colormap = Colors.colormap("Blues")
        mycolorMap = UT.Colormap([0.0, LyapMax], colormap)
        cost_ordered =
            reverse(sort(hcat([(lyap, state) for (state, lyap) in lyap_fun]...); dims = 2))
        for (lyap, state) in cost_ordered
            pos = get_xpos_by_state(symmodel, state)
            @series begin
                lyap ≠ Inf ? color := UT.get_color(mycolorMap, lyap) : color := :yellow
                return grid, pos
            end
        end
        @series begin
            mycolorMap
        end
    else
        dict = Dict{NTuple{2, Int}, Any}()
        for s in enum_states(symmodel)
            pos = get_xpos_by_state(symmodel, s)
            if !haskey(dict, pos[dims])
                dict[pos[dims]] = true
                @series begin
                    legend := false
                    return grid, pos
                end
            end
        end
    end
    # Display the arrows
    if arrowsB
        for (target, source, symbol) in symmodel.autom.transitions
            if source == target
                @series begin
                    p1 = DO.get_coord_by_pos(grid, get_xpos_by_state(symmodel, source))
                    return UT.DrawPoint(p1)
                end
            else
                @series begin
                    p1 = DO.get_coord_by_pos(grid, get_xpos_by_state(symmodel, source))
                    p2 = DO.get_coord_by_pos(grid, get_xpos_by_state(symmodel, target))
                    return UT.DrawArrow(p1, p2)
                end
            end
        end
    end
end
