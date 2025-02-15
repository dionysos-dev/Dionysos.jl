"""
    SymbolicModelList{N, M, S1, S2, A} <: GridBasedSymbolicModel{N, M}

A classical symbolic model where the entire domain is partitioned into grid cells.
"""
mutable struct SymbolicModelList{
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

function NewSymbolicModelListList(
    Xdom,
    Udom,
    ::Type{S} = UT.SortedTupleSet{3, NTuple{3, Int}},
) where {S}
    nx = DO.get_ncells(Xdom)
    xint2pos = [pos for pos in DO.enum_pos(Xdom)]
    xpos2int = Dict((pos, i) for (i, pos) in enumerate(DO.enum_pos(Xdom)))

    customDomainList = DO.convert_to_custom_domain(Udom)
    nu = DO.get_ncells(customDomainList)
    uint2coord = [coord for coord in DO.enum_elems(customDomainList)]
    ucoord2int =
        Dict((coord, i) for (i, coord) in enumerate(DO.enum_elems(customDomainList)))
    autom = AutomatonList{S}(nx, nu)

    return SymbolicModelList(
        Xdom,
        customDomainList,
        autom,
        xpos2int,
        xint2pos,
        ucoord2int,
        uint2coord,
    )
end

function with_automaton(symmodel::SymbolicModelList, autom)
    return SymbolicModelList(
        symmodel.Xdom,
        symmodel.Udom,
        autom,
        symmodel.xpos2int,
        symmodel.xint2pos,
        symmodel.ucoord2int,
        symmodel.uint2coord,
    )
end

get_n_state(symmodel::SymbolicModelList) = length(symmodel.xint2pos)
get_n_input(symmodel::SymbolicModelList) = length(symmodel.uint2coord)
enum_states(symmodel::SymbolicModelList) = 1:get_n_state(symmodel)
enum_inputs(symmodel::SymbolicModelList) = 1:get_n_input(symmodel)
get_state_domain(symmodel::SymbolicModelList) = symmodel.Xdom
get_input_domain(symmodel::SymbolicModelList) = symmodel.Udom

get_xpos_by_state(symmodel::SymbolicModelList, state) = symmodel.xint2pos[state]
get_state_by_xpos(symmodel::SymbolicModelList, xpos) = symmodel.xpos2int[xpos]
is_xpos(symmodel::SymbolicModelList, xpos) = haskey(symmodel.xpos2int, xpos)

get_concrete_input(symmodel::SymbolicModelList, input) = symmodel.uint2coord[input]
get_abstract_input(symmodel::SymbolicModelList, u) = symmodel.ucoord2int[u]

add_transitions!(symmodel::SymbolicModelList, translist) =
    add_transitions!(symmodel.autom, translist)

# for continuous time system with growth-bound function
function compute_reachable_set(
    contsys::ST.ControlSystemGrowth{N},
    rect::UT.HyperRectangle,
    u,
) where {N}
    x = UT.get_center(rect)
    r = UT.get_r(rect)
    Fx = contsys.sys_map(x, u, contsys.tstep)
    Fr = contsys.growthbound_map(r, u, contsys.tstep) + contsys.measnoise
    return UT.HyperRectangle(Fx - Fr, Fx + Fr)
end

# for continuous time system with deltaGUAS
# function compute_reachable_set(contsys::ST.ControlSystemGrowth{N}, rect::UT.HyperRectangle, u) where {N}
#     x = UT.get_center(rect)
#     Fx = contsys.sys_map(x, u, contsys.tstep)
#     return UT.HyperRectangle(Fx, Fx)
# end

function compute_abstract_transitions!(
    symmodel::GridBasedSymbolicModel,
    reachable_set::UT.HyperRectangle,
    abstract_state,
    abstract_input,
    translist,
)
    Xdom = get_state_domain(symmodel)
    rectI = DO.get_pos_lims_outer(DO.get_grid(Xdom), reachable_set)
    ypos_iter = Iterators.product(DO._ranges(rectI)...)
    allin = true
    for ypos in ypos_iter
        if !(ypos in Xdom)
            allin = false
            break
        end
        target = get_state_by_xpos(symmodel, ypos)
        push!(translist, (target, abstract_state, abstract_input))
    end
    return allin
end

# Compute the transition for abstract_system to be in feedback refinement relation with concrete_system
# -> General case
function compute_symmodel_from_controlsystem!(
    abstract_system::GridBasedSymbolicModel{N},
    concrete_system; #::ST.ControlSystemGrowth{N},
    compute_reachable_set = compute_reachable_set,
) where {N}
    println("compute_symmodel_from_controlsystem! started")
    ntrans = 0
    translist = Tuple{Int, Int, Int}[]
    for abstract_input in enum_inputs(abstract_system)
        concrete_input = get_concrete_input(abstract_system, abstract_input)
        for abstract_state in enum_states(abstract_system)
            concrete_elem = get_concrete_elem(abstract_system, abstract_state)
            reachable_set =
                compute_reachable_set(concrete_system, concrete_elem, concrete_input)
            empty!(translist)
            allin = compute_abstract_transitions!(
                abstract_system,
                reachable_set,
                abstract_state,
                abstract_input,
                translist,
            )
            if allin
                add_transitions!(abstract_system, translist)
                ntrans += length(translist)
            end
        end
    end
    return println(
        "compute_symmodel_from_controlsystem! terminated with success: ",
        "$(ntrans) transitions created",
    )
end

function get_over_approximation(concrete_system::ST.ControlSystemGrowth, concrete_input, r)
    return concrete_system.growthbound_map(r, concrete_input, concrete_system.tstep) +
           concrete_system.measnoise
end

# Compute the transition for abstract_system to be in feedback refinement relation with concrete_system
# -> If the overapproximation is independant of the state x, but only of the input u, and radius of the grid
# simulate the center and then add the overapproximation
# we should add a faster variable to decide to apply this more efficient implementation
function compute_symmodel_from_controlsystem!(
    abstract_system::GridBasedSymbolicModel{N},
    concrete_system::ST.ControlSystemGrowth{N};
    get_over_approximation = get_over_approximation, # independant of the state
) where {N}
    println("compute_symmodel_from_controlsystem! started")
    Xdom = get_state_domain(abstract_system)
    r = DO.get_h(DO.get_grid(Xdom)) / 2.0 + concrete_system.measnoise
    ntrans = 0
    translist = Tuple{Int, Int, Int}[]

    for abstract_input in enum_inputs(abstract_system)
        concrete_input = get_concrete_input(abstract_system, abstract_input)
        Fr = get_over_approximation(concrete_system, concrete_input, r)
        for abstract_state in enum_states(abstract_system)
            empty!(translist)
            x = get_concrete_state(abstract_system, abstract_state)
            Fx = concrete_system.sys_map(x, concrete_input, concrete_system.tstep) # should be adapted
            reachable_set = UT.HyperRectangle(Fx - Fr, Fx + Fr)
            rectI = DO.get_pos_lims_outer(Xdom.grid, reachable_set)
            ypos_iter = Iterators.product(DO._ranges(rectI)...)
            allin = true
            for ypos in ypos_iter
                if !(ypos in Xdom)
                    allin = false
                    break
                end
                target = get_state_by_xpos(abstract_system, ypos)
                push!(translist, (target, abstract_state, abstract_input))
            end
            if allin
                add_transitions!(abstract_system, translist)
                ntrans += length(translist)
            end
        end
    end
    return println(
        "compute_symmodel_from_controlsystem! terminated with success: ",
        "$(ntrans) transitions created",
    )
end

function compute_deterministic_symmodel_from_controlsystem!(
    abstract_system::GridBasedSymbolicModel{N},
    concrete_system;
    tstep = concrete_system.tstep,
) where {N}
    println("compute_deterministic_symmodel_from_controlsystem! started")
    Xdom = get_state_domain(abstract_system)
    Udom = get_input_domain(abstract_system)
    ntrans = 0

    for abstract_input in enum_inputs(abstract_system)
        concrete_input = get_concrete_input(abstract_system, abstract_input)
        for xpos in DO.enum_pos(Xdom)
            source = get_state_by_xpos(abstract_system, xpos)
            x = DO.get_coord_by_pos(Xdom.grid, xpos)
            Fx = concrete_system.sys_map(x, concrete_input, tstep)
            ypos = DO.get_pos_by_coord(Xdom.grid, Fx)
            if ypos in Xdom
                target = get_state_by_xpos(abstract_system, ypos)
                HybridSystems.add_transition!(
                    abstract_system.autom,
                    source,
                    target,
                    abstract_input,
                )
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
    abstract_system::GridBasedSymbolicModel{N},
    concrete_system::ST.ControlSystemLinearized{N},
) where {N}
    println("compute_symmodel_from_controlsystem! started")
    Xdom = get_state_domain(abstract_system)
    Udom = get_input_domain(abstract_system)
    tstep = concrete_system.tstep
    r = Xdom.grid.h / 2.0 + concrete_system.measnoise
    _H_ = SMatrix{N, N}(I) .* r
    _ONE_ = ones(SVector{N})
    e = norm(r, Inf)
    ntrans = 0
    translist = Tuple{Int, Int, Int}[]

    for abstract_input in enum_inputs(abstract_system)
        concrete_input = get_concrete_input(abstract_system, abstract_input)
        Fe = concrete_system.error_map(e, concrete_input, concrete_system.tstep)
        Fr = r .+ Fe
        for xpos in DO.enum_pos(Xdom)
            empty!(translist)
            source = get_state_by_xpos(abstract_system, xpos)
            x = DO.get_coord_by_pos(Xdom.grid, xpos)
            Fx, DFx = concrete_system.linsys_map(x, _H_, concrete_input, tstep)
            A = inv(DFx)
            b = abs.(A) * Fr .+ 1.0
            HP = UT.CenteredPolyhedron(A, b)
            # TODO: can we improve abs.(DFx)*_ONE_?
            rad = concrete_system.measnoise + abs.(DFx) * _ONE_ .+ Fe
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
                target = get_state_by_xpos(abstract_system, ypos)
                push!(translist, (target, source, abstract_input))
            end
            if allin
                add_transitions!(abstract_system.autom, translist)
                ntrans += length(translist)
            end
        end
    end
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
