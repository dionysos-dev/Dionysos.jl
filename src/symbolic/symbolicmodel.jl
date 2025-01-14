using Plots, Colors
using ProgressMeter

"""
    SymbolicModel{N, M}

is the abtract type which defines a symbolic model.
"""
abstract type SymbolicModel{N, M} end

"""
    SymbolicModelList{N, M, S1 <: DO.DomainType{N}, S2 <: DO.DomainType{M}, A} <: SymbolicModel{N, M}
    
is one implementation of the `SymbolicModel` type for classical abstraction-based methods, i.e. when the whole domain is partitioned/covered.
"""
mutable struct SymbolicModelList{N, M, S1 <: DO.DomainType{N}, S2 <: DO.DomainType{M}, A} <:
               SymbolicModel{N, M}
    Xdom::S1
    Udom::S2
    autom::A
    xpos2int::Dict{NTuple{N, Int}, Int}
    xint2pos::Vector{NTuple{N, Int}}
    upos2int::Dict{NTuple{M, Int}, Int}
    uint2pos::Vector{NTuple{M, Int}}
end

# ListList refers to List for SymbolicModel, and List for automaton
function NewSymbolicModelListList(
    Xdom,
    Udom,
    ::Type{S} = UT.SortedTupleSet{4, NTuple{4, Int}},
) where {S}
    nx = DO.get_ncells(Xdom)
    nu = DO.get_ncells(Udom)
    xint2pos = [pos for pos in DO.enum_pos(Xdom)]
    xpos2int = Dict((pos, i) for (i, pos) in enumerate(DO.enum_pos(Xdom)))
    uint2pos = [pos for pos in DO.enum_pos(Udom)]
    upos2int = Dict((pos, i) for (i, pos) in enumerate(DO.enum_pos(Udom)))
    autom = AutomatonList{S}(nx, nu)
    return SymbolicModelList(Xdom, Udom, autom, xpos2int, xint2pos, upos2int, uint2pos)
end

function with_automaton(symmodel::SymbolicModelList, autom)
    return SymbolicModelList(
        symmodel.Xdom,
        symmodel.Udom,
        autom,
        symmodel.xpos2int,
        symmodel.xint2pos,
        symmodel.upos2int,
        symmodel.uint2pos,
    )
end

function is_in(symmodel::SymbolicModelList, xpos)
    return haskey(symmodel.xpos2int, xpos)
end

function get_xpos_by_state(symmodel::SymbolicModelList, state)
    return symmodel.xint2pos[state]
end

function get_state_by_xpos(symmodel::SymbolicModelList, xpos)
    return symmodel.xpos2int[xpos]
end

function get_state_by_coord(symmodel::SymbolicModelList, x)
    xpos = DO.get_pos_by_coord(symmodel.Xdom, x)
    return get_state_by_xpos(symmodel, xpos)
end

function get_all_states_by_xpos(symmodel::SymbolicModelList, l_xpos)
    return [symmodel.xpos2int[xpos] for xpos in l_xpos]
end

function get_upos_by_symbol(symmodel::SymbolicModel, symbol)
    return symmodel.uint2pos[symbol]
end

function get_symbol_by_upos(symmodel::SymbolicModel, upos)
    return symmodel.upos2int[upos]
end

function enum_cells(symmodel::SymbolicModelList)
    return 1:length(symmodel.xint2pos)
end

function get_domain_from_symbols(symmodel::SymbolicModelList, symbols)
    newDomain = DO.DomainList(symmodel.Xdom.grid)
    for symbol in symbols
        xpos = get_xpos_by_state(symmodel, symbol)
        DO.add_pos!(newDomain, xpos)
    end
    return newDomain
end

function compute_symmodel_from_data!(
    symmodel::SymbolicModel{N},
    contsys::ST.ControlSystemGrowth{N};
    n_samples = 50,
    ε = 0.0
) where {N}
    println("compute_symmodel_from_data! started")
    Xdom = symmodel.Xdom
    dim = length(Xdom.grid.orig)
    Udom = symmodel.Udom
    tstep = contsys.tstep
    γ = 1.1                 # growth factor of timestep
    tstep_max = 0.5
    count = 0

    write = true
    transdict = Dict{Tuple{Int, Int, Int, Int}, Float64}() # {(target, source, symbol): prob}
    
    enum_u = DO.enum_pos(Udom)
    for (i, upos) in enumerate(enum_u) # for each input
        println("$i / $(length(enum_u))")
        symbol = get_symbol_by_upos(symmodel, upos)
        u = DO.get_coord_by_pos(Udom.grid, upos)
        enum_x = DO.enum_pos(Xdom)
        for (j, xpos) in enumerate(DO.enum_pos(Xdom)) # for each cell
            # if j % 10000 == 0
            #     println("  $j / $(length(enum_x))")
            # end

            # while ici, adapter time step pour chaque tuple cellule, input
            self_loops = true
            tstep_cur = tstep
            power = 0
            source = get_state_by_xpos(symmodel, xpos) # cellule source
            rec = DO.get_rec(Xdom.grid, xpos)
            # uniform sampling in the cell
            x_sampled = [SVector(rec.lb .+ (rec.ub .- rec.lb) .* rand(dim)) for _ in 1:n_samples]

            while self_loops && tstep_cur <= tstep_max
                self_loops = false
                Fx_sampled = [contsys.sys_map(x, u, tstep_cur) for x ∈ x_sampled] # x_k+1
                pos_sampled = [DO.get_pos_by_coord(Xdom.grid, Fx) for Fx ∈ Fx_sampled]
                pos_contained_sampled = [ypos ∈ Xdom for ypos in pos_sampled]
                if !all(pos_contained_sampled)
                    continue
                end
                target_sampled = [get_state_by_xpos(symmodel, pos) for pos ∈ pos_sampled]
                
                for target in target_sampled # enumère les cellules images
                    if target == source
                        self_loops = true
                        # tstep *= 1.1
                        power += 1
                        tstep_cur = tstep * γ^power
                        if tstep_cur > tstep_max
                            count += 1
                        end
                        break
                    end
                end

                # écrire les transitions finales
                if !self_loops || tstep_cur > tstep_max
                    if tstep_cur > tstep_max
                        power = 0
                        # si on arrive pas à éliminer les self-loops, on reprend le plus petit time step
                        # surement moyen de faire + efficace
                        Fx_sampled = [contsys.sys_map(x, u, tstep) for x ∈ x_sampled] # x_k+1
                        pos_sampled = [DO.get_pos_by_coord(Xdom.grid, Fx) for Fx ∈ Fx_sampled]
                        pos_contained_sampled = [ypos ∈ Xdom for ypos in pos_sampled]
                        if !all(pos_contained_sampled)
                            continue
                        end
                        target_sampled = [get_state_by_xpos(symmodel, pos) for pos ∈ pos_sampled]
                    end
                    for target in target_sampled # enumère les cellules images
                        if (target, source, symbol, power) ∈ keys(transdict)
                            transdict[target, source, symbol, power] += 1
                        else
                            transdict[target, source, symbol, power] = 1
                        end
                    end
                end
            end 
        end
    end
    println("Number of times max time step reached: $count")
    #translist = Tuple{Int, Int, Int}[]
    translist = Tuple{Int, Int, Int, Int}[]
    for (t, s, sym, p) ∈ keys(transdict)
        if transdict[(t, s, sym, p)] / n_samples >= ε
            #push!(translist, (t, s, sym))
            push!(translist, (t, s, sym, p))
        end
    end
    add_transitions!(symmodel.autom, translist)
    return println(
        "compute_symmodel_from_data! terminated with success: ",
        "$(length(translist)) transitions created",
    )
end

# Assumes that automaton is "empty"
# Compare to OLD implementation (see below), we do not make a first check before:
# we go through the list only once; this requires to store the transitions in a
# vector (translist). This approach uses a bit more allocations than the OLD one
# (29 vs 24/26) on pathplanning-simple/hard but is faster in both cases.
function compute_symmodel_from_controlsystem!(
    symmodel::SymbolicModel{N},
    contsys::ST.ControlSystemGrowth{N},
) where {N}
    println("compute_symmodel_from_controlsystem! started")
    Xdom = symmodel.Xdom
    Udom = symmodel.Udom
    tstep = contsys.tstep
    r = Xdom.grid.h / 2.0 + contsys.measnoise
    ntrans = 0
    # Vector to store transitions
    translist = Tuple{Int, Int, Int}[]

    # Updates every 1 seconds
    # Commented because it changes the number of allocations
    # @showprogress 1 "Computing symbolic control system: " (
    for upos in DO.enum_pos(Udom)
        symbol = get_symbol_by_upos(symmodel, upos)
        u = DO.get_coord_by_pos(Udom.grid, upos)
        Fr = contsys.growthbound_map(r, u, contsys.tstep) + contsys.measnoise
        for xpos in DO.enum_pos(Xdom)
            empty!(translist)
            source = get_state_by_xpos(symmodel, xpos)
            x = DO.get_coord_by_pos(Xdom.grid, xpos)
            Fx = contsys.sys_map(x, u, tstep)
            rectI = DO.get_pos_lims_outer(Xdom.grid, UT.HyperRectangle(Fx - Fr, Fx + Fr))
            ypos_iter = Iterators.product(DO._ranges(rectI)...)
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
    return println(
        "compute_symmodel_from_controlsystem! terminated with success: ",
        "$(ntrans) transitions created",
    )
end

function compute_deterministic_symmodel_from_controlsystem!(
    symmodel::SymbolicModel{N},
    contsys;
    tstep = contsys.tstep,
) where {N}
    println("compute_deterministic_symmodel_from_controlsystem! started")
    Xdom = symmodel.Xdom
    Udom = symmodel.Udom
    ntrans = 0

    for upos in DO.enum_pos(Udom)
        symbol = get_symbol_by_upos(symmodel, upos)
        u = DO.get_coord_by_pos(Udom.grid, upos)
        for xpos in DO.enum_pos(Xdom)
            source = get_state_by_xpos(symmodel, xpos)
            x = DO.get_coord_by_pos(Xdom.grid, xpos)
            Fx = contsys.sys_map(x, u, tstep)
            ypos = DO.get_pos_by_coord(Xdom.grid, Fx)
            if ypos in Xdom
                target = get_state_by_xpos(symmodel, ypos)
                HybridSystems.add_transition!(symmodel.autom, source, target, symbol)
                ntrans += 1
            end
        end
    end
    return println(
        "compute_symmodel_from_controlsystem! terminated with success: ",
        "$(ntrans) transitions created",
    )
end

# TODO: check where to place contsys.measnoise (for pathplanning, it is equal to zero)
# So not critical for the moment...
function compute_symmodel_from_controlsystem!(
    symmodel::SymbolicModel{N},
    contsys::ST.ControlSystemLinearized{N},
) where {N}
    println("compute_symmodel_from_controlsystem! started")
    Xdom = symmodel.Xdom
    Udom = symmodel.Udom
    tstep = contsys.tstep
    r = Xdom.grid.h / 2.0 + contsys.measnoise
    _H_ = SMatrix{N, N}(I) .* r
    _ONE_ = ones(SVector{N})
    e = norm(r, Inf)
    ntrans = 0
    translist = Tuple{Int, Int, Int}[]

    # Updates every 1 seconds
    # Commented because it changes the number of allocations
    # @showprogress 1 "Computing symbolic control system: " (
    for upos in DO.enum_pos(Udom)
        symbol = get_symbol_by_upos(symmodel, upos)
        u = DO.get_coord_by_pos(Udom.grid, upos)
        Fe = contsys.error_map(e, u, contsys.tstep)
        Fr = r .+ Fe
        for xpos in DO.enum_pos(Xdom)
            empty!(translist)
            source = get_state_by_xpos(symmodel, xpos)
            x = DO.get_coord_by_pos(Xdom.grid, xpos)
            Fx, DFx = contsys.linsys_map(x, _H_, u, tstep)
            A = inv(DFx)
            b = abs.(A) * Fr .+ 1.0
            HP = UT.CenteredPolyhedron(A, b)
            # TODO: can we improve abs.(DFx)*_ONE_?
            rad = contsys.measnoise + abs.(DFx) * _ONE_ .+ Fe
            rectI = DO.get_pos_lims_outer(Xdom.grid, UT.HyperRectangle(Fx - rad, Fx + rad))
            ypos_iter = Iterators.product(DO._ranges(rectI)...)
            allin = true
            for ypos in ypos_iter
                y = DO.get_coord_by_pos(Xdom.grid, ypos) - Fx
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
    return println(
        "compute_symmodel_from_controlsystem! terminated with success: ",
        "$(ntrans) transitions created",
    )
end

@recipe function f(symmodel::SymbolicModel; arrowsB = false, cost = false, lyap_fun = [])
    # Display the cells
    state_grid = symmodel.Xdom.grid
    if cost
        LyapMax = max(filter(isfinite, getfield.([lyap_fun...], :second))...)
        colormap = Colors.colormap("Blues")
        mycolorMap = UT.Colormap([0.0, LyapMax], colormap)
        cost_ordered =
            reverse(sort(hcat([(lyap, state) for (state, lyap) in lyap_fun]...); dims = 2))
        for (lyap, state) in cost_ordered
            pos = get_xpos_by_state(symmodel, state)
            elli = DO.get_elem_by_pos(state_grid, pos)
            @series begin
                lyap ≠ Inf ? color := UT.get_color(mycolorMap, lyap) : color := :yellow
                return elli
            end
        end
        @series begin
            mycolorMap
        end
    else
        @series begin
            symmodel.Xdom
        end
    end
    # Display the arrows
    if arrowsB
        for t in symmodel.autom.transitions.data
            if t[1] == t[2]
                @series begin
                    color = RGB(
                        abs(0.6 * sin(t[1])),
                        abs(0.6 * sin(t[1] + 2π / 3)),
                        abs(0.6 * sin(t[1] - 2π / 3)),
                    )
                    p1 = DO.get_coord_by_pos(state_grid, get_xpos_by_state(symmodel, t[2]))
                    return UT.DrawPoint(p1)
                end
            else
                @series begin
                    color = RGB(
                        abs(0.6 * sin(t[1])),
                        abs(0.6 * sin(t[1] + 2π / 3)),
                        abs(0.6 * sin(t[1] - 2π / 3)),
                    )
                    p1 = DO.get_coord_by_pos(state_grid, get_xpos_by_state(symmodel, t[2]))
                    p2 = DO.get_coord_by_pos(state_grid, get_xpos_by_state(symmodel, t[1]))
                    return UT.DrawArrow(p1, p2)
                end
            end
        end
    end
end