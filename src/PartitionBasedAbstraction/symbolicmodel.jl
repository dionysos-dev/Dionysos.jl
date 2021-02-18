using LazySets: center
import .Abstraction.compute_symmodel_from_controlsystem!

#true si P intersect un des sets de L
function intersect_union(P,L)
    return any(L) do l
        !is_intersection_empty(P,l)
    end
end

# return symbols of cells of the abstraction fully contained in P or intersecting P
function get_symbols(P::AbstractPolytope{T},symmodel::AB.SymbolicModelList{N},incl_mode::AB.INCL_MODE) where {N,T}
    grid = symmodel.Xdom.grid
    B = overapproximate(P) #bounding box Hyperrectangle
    c = B.center
    r = B.radius

    rect = AB.HyperRectangle(SVector{N,T}(c.-r),SVector{N,T}(c.+r))
    rectI = AB.get_pos_lims(grid, rect, AB.OUTER)
    list = Int[]
    for pos in Iterators.product(AB._ranges(rectI)...)
        if pos ∈ symmodel.Xdom
            center = Array(AB.get_coord_by_pos(grid, pos))
            cell = Hyperrectangle(center,grid.h./2)
            if incl_mode == AB.INNER
                if issubset(cell,P)
                    push!(list,AB.get_state_by_xpos(symmodel, pos))
                end
            else
                if !is_intersection_empty(cell,P)
                    push!(list,AB.get_state_by_xpos(symmodel, pos))
                end
            end
        end
    end
    return list
end

function mesh(P::AbstractPolytope{T}, Δ::SVector{N,T},obstacles::Vector) where {N,T}
    B = overapproximate(P) #bounding box Hyperrectangle
    c = B.center
    r = B.radius

    grid = AB.GridFree(SVector{N,T}(zeros(T,N)),Δ) #arbitrary choice:SVector{N,T}(c) decided to set the center in the "center" of P
    domain = AB.DomainList(grid)
    rect = AB.HyperRectangle(SVector{N,T}(c.-r),SVector{N,T}(c.+r))
    rectI = AB.get_pos_lims(grid, rect, AB.OUTER)
    for pos in Iterators.product(AB._ranges(rectI)...)
        center = Array(AB.get_coord_by_pos(grid, pos))
        cell = Hyperrectangle(center,grid.h./2)
        if center ∈ P && !intersect_union(cell, obstacles)
            AB.add_pos!(domain, pos)
        end
    end
    return domain
end

function AB.compute_symmodel_from_controlsystem!(symmodel::AB.SymbolicModel{N},contsys::ControlSystemLipschitz{N}) where N
    println("compute_symmodel_from_controlsystem! started")
    Xdom = symmodel.Xdom
    Udom = symmodel.Udom
    tstep = contsys.tstep
    r = Xdom.grid.h/2.0 + contsys.measnoise
    Fr = r#r .+ r*exp(contsys.L*tstep)# in this case
    ntrans = 0
    translist = Tuple{Int,Int,Int}[]

    for upos in AB.enum_pos(Udom)
        symbol = AB.get_symbol_by_upos(symmodel, upos)
        u = AB.get_coord_by_pos(Udom.grid, upos)
        for xpos in AB.enum_pos(Xdom)
            AB.empty!(translist)
            source = AB.get_state_by_xpos(symmodel, xpos)
            x = AB.get_coord_by_pos(Xdom.grid, xpos)
            Fx = contsys.sys_map(x, u, tstep)
            rectI = AB.get_pos_lims_outer(Xdom.grid, AB.HyperRectangle(Fx .- Fr, Fx .+ Fr))
            ypos_iter = Iterators.product(AB._ranges(rectI)...)
            allin = true
            for ypos in ypos_iter
                if !(ypos in Xdom)
                    allin = false
                    break
                end
                target = AB.get_state_by_xpos(symmodel, ypos)
                push!(translist, (target, source, symbol))
            end
            if allin
                AB.add_transitions!(symmodel.autom, translist)
                ntrans += AB.length(translist)
            end
        end
    end
    println("compute_symmodel_from_controlsystem! terminated with success: ","$(ntrans) transitions created")
end

function build_abstraction(_X_::AbstractPolytope{T},hx::SVector,_U_::AbstractPolytope{T},hu::SVector,contsys::ControlSystem,obstacles::Vector{<:AbstractPolytope{T}}) where T
    Xfull = mesh(_X_,hx,obstacles)
    Ufull = mesh(_U_,hu,[])
    #plot_mesh(Xfull)
    #plot_mesh(Ufull)
    symmodel = AB.NewSymbolicModelListList(Xfull, Ufull)
    compute_symmodel_from_controlsystem!(symmodel, contsys)
    return symmodel
end



######################### to delete later ##############################
function get_coord_from_state(symmodel::AB.SymbolicModelList{N},state) where N
    return AB.get_coord_by_pos(symmodel.Xdom.grid, AB.get_xpos_by_state(symmodel, state))
end
function get_coord_from_symbol(symmodel::AB.SymbolicModelList{N},symbol) where N
    return AB.get_coord_by_pos(symmodel.Udom.grid, AB.get_upos_by_symbol(symmodel, symbol))
end

function print_mesh(domain::AB.DomainList)
    fig = plot(aspect_ratio = 1,legend = false)
    points = [AB.get_coord_by_pos(domain.grid, pos) for pos in domain.elems]
    plot!(Singleton.(points), color=:yellow)
    display(fig)
end

function print_symbolset(symmodel::AB.SymbolicModelList{N},initList::Vector{Int},targetList::Vector{Int}) where N
    for state in initList
        x = get_coord_from_state(symmodel,state)
        plot!(Singleton(x), color = :green)
    end
    for state in targetList
        x = get_coord_from_state(symmodel,state)
        plot!(Singleton(x), color = :red)
    end
end
