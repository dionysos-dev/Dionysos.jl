abstract type SymbolicModel end

struct SymbolicModelList{NX, NU, S<:Automaton} <: SymbolicModel
    Xgrid::GridSpaceList{NX}
    Ugrid::GridSpaceList{NU}
    autom::S
    xpos2int::Dict{NTuple{NX, Int}, Int}
    xint2pos::Vector{NTuple{NX, Int}}
    upos2int::Dict{NTuple{NU, Int}, Int}
    uint2pos::Vector{NTuple{NU, Int}}
end

function NewSymbolicModelListList(
    Xgrid::GridSpaceList{NX}, Ugrid::GridSpaceList{NU}) where {NX, NU}
    #---------------------------------------------------------------------------
    nx = get_ncells(Xgrid)
    nu = get_ncells(Ugrid)
    xint2pos = [pos for pos in enum_pos(Xgrid)]
    xpos2int = Dict((pos, i) for (i, pos) in enumerate(enum_pos(Xgrid)))
    uint2pos = [pos for pos in enum_pos(Ugrid)]
    upos2int = Dict((pos, i) for (i, pos) in enumerate(enum_pos(Ugrid)))
    autom = NewAutomatonList(nx, nu)
    return SymbolicModelList(
        Xgrid, Ugrid, autom, xpos2int, xint2pos, upos2int, uint2pos)
end

function get_xpos_by_state(symmodel::SymbolicModelList, state)
    return symmodel.xint2pos[state]
end

function get_state_by_xpos(symmodel::SymbolicModelList, xpos)
    return symmodel.xpos2int[xpos]
end

function get_upos_by_symbol(symmodel::SymbolicModelList, symbol)
    return symmodel.uint2pos[symbol]
end

function get_symbol_by_upos(symmodel::SymbolicModelList, upos)
    return symmodel.upos2int[upos]
end

# Assumes that automaton is "empty"
function compute_symmodel_from_controlsystem!(symmodel, contsys)
    println("compute_symmodel_from_controlsystem! started")
    Xgrid = symmodel.Xgrid
    Ugrid = symmodel.Ugrid
    tstep = contsys.tstep
    ntrans = 0

    # Updates every 1 seconds
    @showprogress 1 "Computing symbolic control system: " (
    for upos in enum_pos(Ugrid)
        symbol = get_symbol_by_upos(symmodel, upos)
        u = get_coord_by_pos(Ugrid, upos)
        r = Xgrid.h./2 .+ contsys.measnoise
        r = contsys.bound_map(r, u, contsys.tstep)
        r = r .+ contsys.measnoise
        for xpos in enum_pos(Xgrid)
            source = get_state_by_xpos(symmodel, xpos)
            x = get_coord_by_pos(Xgrid, xpos)
            x = contsys.sys_map(x, u, tstep)
            rectI = get_pos_lims_outer(Xgrid, HyperRectangle(x .- r, x .+ r))
            ypos_iter = Iterators.product(_ranges(rectI)...)
            any(x -> !(x ∈ Xgrid), ypos_iter) && continue
            for ypos in ypos_iter
                target = get_state_by_xpos(symmodel, ypos)
                add_transition!(symmodel.autom, source, symbol, target)
            end
            ntrans += length(ypos_iter)
        end
    end)
    println("compute_symmodel_from_controlsystem! terminated with success: ",
        "$(ntrans) transitions created")
end
