using StaticArrays, Plots

mutable struct NestedDomain{N}
    domains::Vector{DO.GeneralDomainList}
    active::Vector{Dict{NTuple{N,Int},Any}}
    levels::Int
end

function NestedDomain(dom::DO.GeneralDomainList{N}) where {N}
    dict = Dict{NTuple{N,Int},Any}()
    for pos in DO.enum_pos(dom)
        dict[pos] = true
    end
    return NestedDomain(DO.GeneralDomainList[dom], [dict], 1)
end

function is_active(Ndomain::NestedDomain, pos, l)
    if l > Ndomain.levels
        return false
    end
    return Ndomain.active[l][pos]
end

function get_depth(Ndomain::NestedDomain, x)
    for l in 1:Ndomain.levels
        pos = DO.get_pos_by_coord(Ndomain, l, x)
        if is_pos(Ndomain,pos,l)
            if is_active(Ndomain, pos, l)
                return l
            end
        else
            return 0
        end
    end
    return 0
end

function add_dom!(Ndomain::NestedDomain, dom::DO.GeneralDomainList) where {N,T}
    push!(Ndomain.domains,dom)
    Ndomain.levels += 1
end


function DO.get_pos_by_coord(Ndomain::NestedDomain, l, x)
    return  DO.get_pos_by_coord(Ndomain.domains[l], x)
end

function DO.get_coord_by_pos(Ndomain::NestedDomain, l, xpos)
    return  DO.get_coord_by_pos(Ndomain.domains[l].grid,xpos)
end

function get_grid(Ndomain::NestedDomain, l)
    return DO.get_grid(Ndomain.domains[l])
end

function DO.get_grid(Ndomain::NestedDomain, x)
    return DO.get_grid(Ndomain, get_depth(Ndomain, x))
end

#add a domain fitting inside the previous one
function add_sub_dom!(Ndomain::NestedDomain{N}) where {N}
    l = Ndomain.levels
    dom = Ndomain.domains[l]
    hx = DO.get_h(dom.grid)/2.0
    subdom = DO.GeneralDomainList(hx;elems = dom.elemsCoord,periodic=dom.periodic,periods=dom.periods,T0=dom.T0,fit=dom.fit)
    push!(Ndomain.domains, subdom)
    push!(Ndomain.active, Dict{NTuple{N,Int},Any}())
    Ndomain.levels += 1
end

function get_subpos(pos)
    lbI = Tuple([p*2 for p in pos])
    ubI = Tuple([p*2+1 for p in pos])
    rectI = UT.HyperRectangle(lbI,ubI)
    return rectI
end

function cut_pos!(Ndomain, pos, l)
    if l == Ndomain.levels
        add_sub_dom!(Ndomain)
    end
    dict = Ndomain.active[l+1]
    Ndomain.active[l][pos] = false
    subpos = Iterators.product(DO._ranges(get_subpos(pos))...)
    for spos in subpos
        dict[spos] = true
    end
    return subpos
end

function Base.empty!(Ndomain::NestedDomain)
    for dom in Ndomain.domains
        empty!(dom)
    end
end

function is_pos(Ndomain::NestedDomain, pos, l)
    return in(pos,Ndomain.domains[l])
end


function DO.get_pos_by_coord(Ndomain::NestedDomain, x)
    l = get_depth(Ndomain, x)
    return  (DO.get_pos_by_coord(Ndomain, l, x), l)
end

function Base.isempty(Ndomain::NestedDomain)
    for dom in Ndomain.domains
        if isempty(dom)
            return true
        end
    end
    return false
end

function DO.get_ncells(Ndomain::NestedDomain)
    n = 0
    for dict in Ndomain.active
        for (pos,v) in dict
            if dict[pos]
                n += 1
            end
        end
    end
    return n
end


function  DO.enum_pos(Ndomain::NestedDomain)
    L = []
    for (l,dom) in enumerate(Ndomain.domains)
         push!(L,(DO.enum_pos(dom)))
    end
    return L
end

function Plots.plot!(Ndomain::NestedDomain;dims=[1,2])
    for l in 1:Ndomain.levels
        for (pos, v) in Ndomain.active[l]
            grid = get_grid(Ndomain, l)
            if v == true
                DO.plot_elem!(grid, pos)
            end
        end
    end
end


#end
