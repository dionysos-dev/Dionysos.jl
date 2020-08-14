abstract type SymbolicModel{N,M} end

struct SymbolicModelList{N,M,S1<:Domain{N},S2<:Domain{M},A<:Automaton} <: SymbolicModel{N,M}
    Xdom::S1
    Udom::S2
    autom::A
    xpos2int::Dict{NTuple{N,Int},Int}
    xint2pos::Vector{NTuple{N,Int}}
    upos2int::Dict{NTuple{M,Int},Int}
    uint2pos::Vector{NTuple{M,Int}}
end

# ListList refers to List for SymbolicModel, and List for automaton
function NewSymbolicModelListList(Xdom, Udom)
    nx = get_ncells(Xdom)
    nu = get_ncells(Udom)
    xint2pos = [pos for pos in enum_pos(Xdom)]
    xpos2int = Dict((pos, i) for (i, pos) in enumerate(enum_pos(Xdom)))
    uint2pos = [pos for pos in enum_pos(Udom)]
    upos2int = Dict((pos, i) for (i, pos) in enumerate(enum_pos(Udom)))
    autom = NewAutomatonList(nx, nu)
    return SymbolicModelList(
        Xdom, Udom, autom, xpos2int, xint2pos, upos2int, uint2pos)
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
# vector (translist). This approach uses a bit more allocations than the OLD one
# (29 vs 24/26) on pathplanning-simple/hard but is faster in both cases.
function compute_symmodel_from_controlsystem!(symmodel::SymbolicModel{N},
        contsys::ControlSystemGrowth{N}) where N
    println("compute_symmodel_from_controlsystem! started")
    Xdom = symmodel.Xdom
    Udom = symmodel.Udom
    tstep = contsys.tstep
    r = Xdom.grid.h/2.0 + contsys.measnoise
    ntrans = 0
    # Vector to store transitions
    translist = Tuple{Int,Int,Int}[]

    # Updates every 1 seconds
    # Commented because it changes the number of allocations
    # @showprogress 1 "Computing symbolic control system: " (
    for upos in enum_pos(Udom)
        symbol = get_symbol_by_upos(symmodel, upos)
        u = get_coord_by_pos(Udom.grid, upos)
        Fr = contsys.growthbound_map(r, u, contsys.tstep) + contsys.measnoise
        for xpos in enum_pos(Xdom)
            empty!(translist)
            source = get_state_by_xpos(symmodel, xpos)
            x = get_coord_by_pos(Xdom.grid, xpos)
            Fx = contsys.sys_map(x, u, tstep)
            rectI = get_pos_lims_outer(Xdom.grid, HyperRectangle(Fx - Fr, Fx + Fr))
            ypos_iter = Iterators.product(_ranges(rectI)...)
            allin = true
            for ypos in ypos_iter
                if !(ypos in Xdom)
                    allin = false
                    break
                end
                target = get_state_by_xpos(symmodel, ypos)
                push!(translist, (target, source, symbol))
            end
            if allin
                add_transitions!(symmodel.autom, translist)
                ntrans += length(translist)
            end
        end
    end
    # )
    println("compute_symmodel_from_controlsystem! terminated with success: ",
        "$(ntrans) transitions created")
end

# Assumes that automaton is "empty"
function compute_symmodel_from_controlsystem_OLD!(symmodel::SymbolicModel{N},
        contsys::ControlSystemGrowth{N}) where N
    println("compute_symmodel_from_controlsystem! started")
    Xdom = symmodel.Xdom
    Udom = symmodel.Udom
    tstep = contsys.tstep
    r = Xdom.grid.h/2.0 + contsys.measnoise
    ntrans = 0
    # Define the function out of the loop. This allowed to recudes the allocations
    # from 1.6M (on pathplanning-simple) to 24!
    not_in_Xdom = x -> !(x ∈ Xdom)

    # Updates every 1 seconds
    # @showprogress 1 "Computing symbolic control system: " (
    for upos in enum_pos(Udom)
        symbol = get_symbol_by_upos(symmodel, upos)
        u = get_coord_by_pos(Udom.grid, upos)
        Fr = contsys.growthbound_map(r, u, contsys.tstep) + contsys.measnoise
        for xpos in enum_pos(Xdom)
            source = get_state_by_xpos(symmodel, xpos)
            x = get_coord_by_pos(Xdom.grid, xpos)
            Fx = contsys.sys_map(x, u, tstep)
            rectI = get_pos_lims_outer(Xdom.grid, HyperRectangle(Fx - Fr, Fx + Fr))
            ypos_iter = Iterators.product(_ranges(rectI)...)
            any(not_in_Xdom, ypos_iter) && continue
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
function compute_symmodel_from_controlsystem!(symmodel::SymbolicModel{N},
        contsys::ControlSystemLinearized{N}) where N
    println("compute_symmodel_from_controlsystem! started")
    Xdom = symmodel.Xdom
    Udom = symmodel.Udom
    tstep = contsys.tstep
    r = Xdom.grid.h/2.0 + contsys.measnoise
    _H_ = SMatrix{N,N}(I).*r
    _ONE_ = ones(SVector{N})
    e = norm(r, Inf)
    ntrans = 0
    translist = Tuple{Int,Int,Int}[]

    # Updates every 1 seconds
    # Commented because it changes the number of allocations
    # @showprogress 1 "Computing symbolic control system: " (
    for upos in enum_pos(Udom)
        symbol = get_symbol_by_upos(symmodel, upos)
        u = get_coord_by_pos(Udom.grid, upos)
        Fe = contsys.error_map(e, u, contsys.tstep)
        Fr = r .+ Fe
        for xpos in enum_pos(Xdom)
            empty!(translist)
            source = get_state_by_xpos(symmodel, xpos)
            x = get_coord_by_pos(Xdom.grid, xpos)
            Fx, DFx = contsys.linsys_map(x, _H_, u, tstep)
            A = inv(DFx)
            b = abs.(A)*Fr .+ 1.0
            HP = CenteredPolyhedron(A, b)
            # TODO: can we improve abs.(DFx)*_ONE_?
            rad = contsys.measnoise + abs.(DFx)*_ONE_ .+ Fe
            rectI = get_pos_lims_outer(Xdom.grid, HyperRectangle(Fx - rad, Fx + rad))
            ypos_iter = Iterators.product(_ranges(rectI)...)
            allin = true
            for ypos in ypos_iter
                y = get_coord_by_pos(Xdom.grid, ypos) - Fx
                !(y in HP) && continue
                if !(ypos in Xdom)
                    allin = false
                    break
                end
                target = get_state_by_xpos(symmodel, ypos)
                push!(translist, (target, source, symbol))
            end
            if allin
                add_transitions!(symmodel.autom, translist)
                ntrans += length(translist)
            end
        end
    end
    # )
    println("compute_symmodel_from_controlsystem! terminated with success: ",
        "$(ntrans) transitions created")
end
