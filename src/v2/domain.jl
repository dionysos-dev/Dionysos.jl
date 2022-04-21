#module MultiDomainList

# using ..Abstraction
# const AB = Abstraction

using StaticArrays, Plots

mutable struct MultiDomain# <: AB.Domain{N,T}
    domains::Vector{D.GeneralDomainList}
    levels::Int
end


function MultiDomain(dom::D.GeneralDomainList) where {N,T}
    return MultiDomain([dom],1)
end

function add_dom!(Mdomain::MultiDomain,dom::D.GeneralDomainList) where {N,T}
    push!(Mdomain.domains,dom)
    Mdomain.levels += 1
end

function AB.get_pos_by_coord(Mdomain::MultiDomain,l,x)
    return  AB.get_pos_by_coord(Mdomain.domains[l],x)
end

function AB.get_coord_by_pos(Mdomain::MultiDomain,l,xpos)
    return  AB.get_coord_by_pos(Mdomain.domains[l].grid,xpos)
end

#add a domain fitting inside the previous one
function add_sub_dom!(Mdomain::MultiDomain) where {N,T}
    l = Mdomain.levels
    dom = Mdomain.domains[l]
    hx = AB.get_h(dom.grid)/2.0
    subdom = D.GeneralDomainList(hx,dom.elemsCoord;periodic=dom.periodic,periods=dom.periods,T0=dom.T0,fit=dom.fit)
    push!(Mdomain.domains,subdom)
    Mdomain.levels += 1
end

function AB.add_pos!(Mdomain::MultiDomain, pos, l)
    AB.add_pos!(Mdomain.domains[l], pos)
end

function AB.add_coord!(Mdomain::MultiDomain, x, l)
    D.add_coord!(Mdomain.domains[l],x)
end

function AB.add_set!(Mdomain::MultiDomain, rect::AB.HyperRectangle, incl_mode::AB.INCL_MODE, l)
    D.add_set!(Mdomain.domains[l], rect, incl_mode)
end

#assuming that the rectangle is already in the domain for periodic dimensions
function AB.get_subset_pos(Mdomain::MultiDomain,rect::AB.HyperRectangle,incl_mode::AB.INCL_MODE,l) where {N}
    return [(pos,l) for pos in AB.get_subset_pos(Mdomain.domains[l],rect,incl_mode)]
end

#assuming that the rectangle is already in the domain for periodic dimensions
function get_subset_pos2(Mdomain::MultiDomain,rect::AB.HyperRectangle,incl_mode::AB.INCL_MODE,l) where {N}
    return [(pos,l) for pos in D.get_subset_pos2(Mdomain.domains[l],rect,incl_mode)]
end

############################################################################################
##
function set_in_period_coord(Mdomain::MultiDomain,x,l)
    return set_in_period_coord(Mdomain.domains[l],x)
end
function set_in_period_pos(Mdomain::MultiDomain,pos)
    return set_in_period_pos(Mdomain.domains[l],pos)
end
##


function AB.add_subset!(Mdomain1::MultiDomain, Mdomain2::MultiDomain, rect::AB.HyperRectangle, incl_mode::AB.INCL_MODE,l)
    AB.add_subset!(Mdomain1.domains[l], Mdomain2.domains[l], rect, incl_mode)
end

function  AB.remove_pos!(Mdomain::MultiDomain, pos, l)
    AB.remove_pos!(Mdomain.domains[l], pos)
end

function  AB.remove_coord!(Mdomain::MultiDomain, x, l)
    AB.remove_coord!(Mdomain.domains[l], x)
end

function  AB.remove_set!(Mdomain::MultiDomain, rect::AB.HyperRectangle, incl_mode::AB.INCL_MODE, l)
    AB.remove_set!(Mdomain.domains[l], rect, incl_mode)
end

function Base.union!(Mdomain1::MultiDomain, Mdomain2::MultiDomain)
    # if Mdomains
    # for dom in Mdomain1
    # union!(domain1.elems, domain2.elems)
end
#
# function Base.setdiff!(domain1::GeneralDomainList, domain2::GeneralDomainList)
#     setdiff!(domain1.elems, domain2.elems)
# end
############################################################################################
function Base.empty!(Mdomain::MultiDomain)
    for dom in Mdomain.domains
        empty!(dom)
    end
end

function is_pos(Mdomain::MultiDomain,pos,l)
    return in(pos,Mdomain.domains[l])
end

function Base.isempty(Mdomain::MultiDomain)
    for dom in Mdomain.domains
        if isempty(dom)
            return true
        end
    end
    return false
end

# function Base.issubset(domain1::MultiDomain, domain2::MultiDomain)
#     return issubset(domain1.elems, domain2.elems)
# end

function  AB.get_ncells(Mdomain::MultiDomain)
    n = 0
    for dom in Mdomain.domains
        n += AB.get_ncells(dom)
    end
    return n
end


function  AB.enum_pos(Mdomain::MultiDomain)
    L = []
    for (l,dom) in enumerate(Mdomain.domains)
         push!(L,(AB.enum_pos(dom)))
    end
    return L
end

function Plots.plot!(Mdomain::MultiDomain;dims=[1,2],l=0) where {N,T}
    if l == 0
        for dom in Mdomain.domains
            D.plot!(dom)
        end
    else
        D.plot!(Mdomain.domains[l])
    end
end
#end
