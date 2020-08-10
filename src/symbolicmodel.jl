abstract type SymbolicModel{N,M} end

struct SymbolicModelList{N,M,S1<:GridSpace{N},S2<:GridSpace{M},A<:Automaton} <: SymbolicModel{N,M}
    Xgrid::S1
    Ugrid::S2
    autom::A
    xpos2int::Dict{NTuple{N,Int},Int}
    xint2pos::Vector{NTuple{N,Int}}
    upos2int::Dict{NTuple{M,Int},Int}
    uint2pos::Vector{NTuple{M,Int}}
end

# ListList refers to List for SymbolicModel, and List for automaton
function NewSymbolicModelListList(Xgrid, Ugrid)
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
# Compare to OLD implementation (see below), we do not make a first check before:
# we go through the list only once; this requires to store the transitions in a
# vector (trans_list). This approach uses a bit more allocations than the OLD one
# (29 vs 24/26) on pathplanning-simple/hard but is faster in both cases.
function compute_symmodel_from_controlsystem!(symmodel, contsys::ControlSystemGrowth)
    println("compute_symmodel_from_controlsystem! started")
    Xgrid = symmodel.Xgrid
    Ugrid = symmodel.Ugrid
    tstep = contsys.tstep
    r = Xgrid.h/2.0 + contsys.measnoise
    ntrans = 0
    # Vector to store transitions
    trans_list = Tuple{Int,Int,Int}[]

    # Updates every 1 seconds
    # Commented because it changes the number of allocations
    # @showprogress 1 "Computing symbolic control system: " (
    for upos in enum_pos(Ugrid)
        symbol = get_symbol_by_upos(symmodel, upos)
        u = get_coord_by_pos(Ugrid, upos)
        Fr = contsys.growthbound_map(r, u, contsys.tstep) + contsys.measnoise
        for xpos in enum_pos(Xgrid)
            empty!(trans_list)
            source = get_state_by_xpos(symmodel, xpos)
            x = get_coord_by_pos(Xgrid, xpos)
            Fx = contsys.sys_map(x, u, tstep)
            rectI = get_pos_lims_outer(Xgrid, HyperRectangle(Fx - Fr, Fx + Fr))
            ypos_iter = Iterators.product(_ranges(rectI)...)
            allin = true
            for ypos in ypos_iter
                if !(ypos in Xgrid)
                    allin = false
                    break
                end
                target = get_state_by_xpos(symmodel, ypos)
                push!(trans_list, (source, symbol, target))
            end
            if allin
                for trans in trans_list
                    add_transition!(symmodel.autom, trans...)
                end
                ntrans += length(trans_list)
            end
        end
    end
    # )
    println("compute_symmodel_from_controlsystem! terminated with success: ",
        "$(ntrans) transitions created")
end

# Assumes that automaton is "empty"
function compute_symmodel_from_controlsystem_OLD!(symmodel, contsys::ControlSystemGrowth)
    println("compute_symmodel_from_controlsystem! started")
    Xgrid = symmodel.Xgrid
    Ugrid = symmodel.Ugrid
    tstep = contsys.tstep
    r = Xgrid.h/2.0 + contsys.measnoise
    ntrans = 0
    # Define the function out of loop. This allowed to recudes the allocations
    # from 1.6M (on pathplanning-simple) to 24!
    not_in_Xgrid = x -> !(x ∈ Xgrid)

    # Updates every 1 seconds
    # @showprogress 1 "Computing symbolic control system: " (
    for upos in enum_pos(Ugrid)
        symbol = get_symbol_by_upos(symmodel, upos)
        u = get_coord_by_pos(Ugrid, upos)
        Fr = contsys.growthbound_map(r, u, contsys.tstep) + contsys.measnoise
        for xpos in enum_pos(Xgrid)
            source = get_state_by_xpos(symmodel, xpos)
            x = get_coord_by_pos(Xgrid, xpos)
            Fx = contsys.sys_map(x, u, tstep)
            rectI = get_pos_lims_outer(Xgrid, HyperRectangle(Fx - Fr, Fx + Fr))
            ypos_iter = Iterators.product(_ranges(rectI)...)
            any(not_in_Xgrid, ypos_iter) && continue
            for ypos in ypos_iter
                target = get_state_by_xpos(symmodel, ypos)
                add_transition!(symmodel.autom, source, symbol, target)
            end
            ntrans += length(ypos_iter)
        end
    end
    # )
    println("compute_symmodel_from_controlsystem! terminated with success: ",
        "$(ntrans) transitions created")
end

# TODO: check where to place contsys.measnoise (for pathplanning, it is equal to zero)
# So not critical for the moment...
function compute_symmodel_from_controlsystem!(
        symmodel::SymbolicModel{N}, contsys::ControlSystemLinearized) where N
    println("compute_symmodel_from_controlsystem! started")
    Xgrid = symmodel.Xgrid
    Ugrid = symmodel.Ugrid
    tstep = contsys.tstep
    r = Xgrid.h/2.0 + contsys.measnoise
    _H_ = SMatrix{N,N}(I).*r
    _ONE_ = ones(SVector{N})
    e = norm(r, Inf)
    ntrans = 0
    trans_list = Tuple{Int,Int,Int}[]

    # Updates every 1 seconds
    # Commented because it changes the number of allocations
    # @showprogress 1 "Computing symbolic control system: " (
    for upos in enum_pos(Ugrid)
        symbol = get_symbol_by_upos(symmodel, upos)
        u = get_coord_by_pos(Ugrid, upos)
        Fe = contsys.error_map(e, u, contsys.tstep)
        Fr = r .+ Fe
        for xpos in enum_pos(Xgrid)
            empty!(trans_list)
            source = get_state_by_xpos(symmodel, xpos)
            x = get_coord_by_pos(Xgrid, xpos)
            Fx, DFx = contsys.linsys_map(x, _H_, u, tstep)
            A = inv(DFx)
            b = abs.(A)*Fr .+ 1.0
            HP = CenteredPolyhedron(A, b)
            # TODO: can we improve abs.(DFx)*_ONE_?
            rad = contsys.measnoise + abs.(DFx)*_ONE_ .+ Fe
            rectI = get_pos_lims_outer(Xgrid, HyperRectangle(Fx - rad, Fx + rad))
            ypos_iter = Iterators.product(_ranges(rectI)...)
            allin = true
            for ypos in ypos_iter
                y = get_coord_by_pos(Xgrid, ypos) - Fx
                !(y in HP) && continue
                if !(ypos in Xgrid)
                    allin = false
                    break
                end
                target = get_state_by_xpos(symmodel, ypos)
                push!(trans_list, (source, symbol, target))
            end
            if allin
                for trans in trans_list
                    add_transition!(symmodel.autom, trans...)
                end
                ntrans += length(trans_list)
            end
        end
    end
    # )
    println("compute_symmodel_from_controlsystem! terminated with success: ",
        "$(ntrans) transitions created")
end
