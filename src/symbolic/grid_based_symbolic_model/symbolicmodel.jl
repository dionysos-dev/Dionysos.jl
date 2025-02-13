using Plots, Colors

"""
    Abstract Type: SymbolicModel{N, M}

Defines a generic symbolic model interface, where:
- `N` is the state space dimension.
- `M` is the input space dimension.
"""
abstract type SymbolicModel{N, M} end

function get_n_state(symmodel::SymbolicModel) end
function get_n_input(symmodel::SymbolicModel) end
function enum_states(symmodel::SymbolicModel) end
function enum_inputs(symmodel::SymbolicModel) end
function get_state_domain(symmodel::SymbolicModel) end
function get_input_domain(symmodel::SymbolicModel) end

function concretize_state(symmodel::SymbolicModel, state) end
function concretize_input(symmodel::SymbolicModel, input) end
function abstract_state(symmodel::SymbolicModel, x) end
function abstract_input(symmodel::SymbolicModel, u) end

"""
    GridBasedSymbolicModel{N, M} <: SymbolicModel{N, M}

An intermediate abstract type for symbolic models that rely on a grid-based discretization.
- `N`: Dimension of the state space.
- `M`: Dimension of the input space.
"""
abstract type GridBasedSymbolicModel{N, M} <: SymbolicModel{N, M} end

function get_xpos_by_state(symmodel::GridBasedSymbolicModel, state) end
function get_state_by_xpos(symmodel::GridBasedSymbolicModel, xpos) end
function get_upos_by_symbol(symmodel::GridBasedSymbolicModel, symbol) end
function get_symbol_by_upos(symmodel::GridBasedSymbolicModel, upos) end
function is_xpos(symmodel::GridBasedSymbolicModel, xpos) end

function concretize_state(symmodel::GridBasedSymbolicModel, state)
    Xdom = get_state_domain(symmodel)
    xpos = get_xpos_by_state(symmodel, state)
    return DO.get_coord_by_pos(Xdom, xpos)
end

function concretize_input(symmodel::GridBasedSymbolicModel, input)
    Udom = get_input_domain(symmodel)
    upos = get_upos_by_symbol(symmodel, input)
    return DO.get_coord_by_pos(Udom, upos)
end

function abstract_state(symmodel::GridBasedSymbolicModel, x)
    xpos = DO.get_pos_by_coord(get_state_domain(symmodel), x)
    return get_state_by_xpos(symmodel, xpos)
end

function abstract_input(symmodel::GridBasedSymbolicModel, u)
    upos = DO.get_pos_by_coord(get_input_domain(symmodel), u)
    return get_symbol_by_upos(symmodel, upos)
end

get_state_grid(symmodel::GridBasedSymbolicModel) = DO.get_grid(get_state_domain(symmodel))

function get_state_by_coord(symmodel::GridBasedSymbolicModel, x)
    xpos = DO.get_pos_by_coord(get_state_domain(symmodel), x)
    return get_state_by_xpos(symmodel, xpos)
end

function get_all_states_by_xpos(symmodel::GridBasedSymbolicModel, l_xpos)
    return [get_state_by_xpos(symmodel, xpos) for xpos in l_xpos]
end

function get_domain_from_states(symmodel::GridBasedSymbolicModel, states)
    newDomain = DO.DomainList(get_state_grid(symmodel))
    for state in states
        DO.add_pos!(newDomain, get_xpos_by_state(symmodel, state))
    end
    return newDomain
end

function get_states_from_set(
    symmodel::GridBasedSymbolicModel,
    set::UT.HyperRectangle,
    position_in_domain,
)
    grid = get_state_grid(symmodel)
    domain_list = DO.DomainList(grid)
    DO.add_subset!(domain_list, get_state_domain(symmodel), set, position_in_domain)
    return [get_state_by_xpos(symmodel, pos) for pos in DO.enum_pos(domain_list)]
end

"""
    SymbolicModelList{N, M, S1, S2, A} <: GridBasedSymbolicModel{N, M}

A classical symbolic model where the entire domain is partitioned into grid cells.
"""
mutable struct SymbolicModelList{
    N,
    M,
    S1 <: DO.GridDomainType{N},
    S2 <: DO.GridDomainType{M},
    A,
} <: GridBasedSymbolicModel{N, M}
    Xdom::S1
    Udom::S2
    autom::A
    xpos2int::Dict{NTuple{N, Int}, Int}
    xint2pos::Vector{NTuple{N, Int}}
    upos2int::Dict{NTuple{M, Int}, Int}
    uint2pos::Vector{NTuple{M, Int}}
end

function NewSymbolicModelListList(
    Xdom,
    Udom,
    ::Type{S} = UT.SortedTupleSet{3, NTuple{3, Int}},
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

get_n_state(symmodel::SymbolicModelList) = length(symmodel.xint2pos)
get_n_input(symmodel::SymbolicModelList) = length(symmodel.uint2pos)
enum_states(symmodel::SymbolicModelList) = 1:get_n_state(symmodel)
enum_inputs(symmodel::SymbolicModelList) = 1:get_n_input(symmodel)
get_state_domain(symmodel::SymbolicModelList) = symmodel.Xdom
get_input_domain(symmodel::SymbolicModelList) = symmodel.Udom

get_xpos_by_state(symmodel::SymbolicModelList, state) = symmodel.xint2pos[state]
get_state_by_xpos(symmodel::SymbolicModelList, xpos) = symmodel.xpos2int[xpos]
get_upos_by_symbol(symmodel::SymbolicModelList, symbol) = symmodel.uint2pos[symbol]
get_symbol_by_upos(symmodel::SymbolicModelList, upos) = symmodel.upos2int[upos]
is_xpos(symmodel::SymbolicModelList, xpos) = haskey(symmodel.xpos2int, xpos)

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
    Xdom = get_state_domain(symmodel)
    Udom = get_input_domain(symmodel)
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
