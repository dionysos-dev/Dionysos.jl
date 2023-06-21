mutable struct LazySymbolicModel{N, M, S1 <: DO.DomainType{N}, S2 <: DO.DomainType{M}, A} <:
               SymbolicModel{N, M}
    Xdom::S1
    Udom::S2
    autom::A
    xpos2int::Dict{NTuple{N, Int}, Int}
    xint2pos::Vector{NTuple{N, Int}}
    upos2int::Any
    uint2pos::Any
end

function LazySymbolicModel(
    Xdom::DO.GeneralDomainList{N, DO.RectangularObstacles{NTuple{N, T1}}},
    Udom::DO.DomainType{N, T2},
) where {N, T1, T2}
    nu = DO.get_ncells(Udom)
    uint2pos = [pos for pos in DO.enum_pos(Udom)]
    upos2int = Dict((pos, i) for (i, pos) in enumerate(DO.enum_pos(Udom)))
    return symmodel = LazySymbolicModel(
        Xdom,
        Udom,
        AutomatonList{Set{NTuple{3, Int}}}(0, nu),
        Dict{NTuple{N, T1}, Int}(),
        NTuple{N, T1}[],
        upos2int,
        uint2pos,
    )
end

function get_xpos_by_state(symmodel::LazySymbolicModel, state)
    return symmodel.xint2pos[state]
end

function get_state_by_xpos(symmodel::LazySymbolicModel, pos)
    id = get(symmodel.xpos2int, pos, nothing)
    created = false
    if id === nothing
        if pos in symmodel.Xdom
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
    return id::Int, created
end

function get_state_by_coord(symmodel::LazySymbolicModel, coord)
    pos = DO.get_pos_by_coord(symmodel.Xdom, coord)
    state, created = get_state_by_xpos(symmodel, pos)
    return state
end

function get_symbol(symmodel, subset, incl_mode::DO.INCL_MODE)
    Xdom = symmodel.Xdom
    grid = Xdom.grid
    posL = DO.get_subset_pos(Xdom, subset, incl_mode)
    symbolsList = [get_state_by_xpos(symmodel, pos)[1] for pos in posL]
    return symbolsList
end

function get_symbols(symmodel, subsetList, incl_mode::DO.INCL_MODE)
    symbols = Int[]
    for subset in subsetList
        append!(symbols, get_symbol(symmodel, subset, incl_mode))
    end
    return symbols
end

function get_ncells(symmodel::LazySymbolicModel)
    return symmodel.autom.nstates
end

function enum_cells(symmodel::LazySymbolicModel)
    return 1:get_ncells(symmodel)
end

@recipe function f(symmodel::LazySymbolicModel; dims = [1, 2], arrowsB = false, cost = false, lyap_fun = [])
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
        for s in enum_cells(symmodel)
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

