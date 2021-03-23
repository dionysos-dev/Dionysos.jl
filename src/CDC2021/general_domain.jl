
module DomainList

using ..Abstraction
const AB = Abstraction

using StaticArrays


## add the periodicity in the domain (add into Domain.jl)
struct GeneralDomainList{N,T,S<:AB.Grid{N,T}} <: AB.Domain{N,T}
    grid::S
    elems::Set{NTuple{N,Int}}
    periodic::Vector{Int} # components which are periodic
    periods::Vector{Float64}  # periods
    T0::Vector{Float64}
    nx::Vector{Int} # number of cell in the periodic directions
    lims::Union{Nothing,AB.HyperRectangle} # lower and upper bound on the non periodic dimensions
                            # can be used to be more efficient to compute the cells which belong to a
                            # given hypperrectangle and the domain if the rectangle is outside the domain
end

# hx is the desired step size, but for periodic direction, it can be changed to fit exactly to the period
# (for periodic dims, hx has to fit exactly in the period)
# The periods are [T0[i], T0[i] + periods[i]]
# for periodic dimensions, I set the origin in T0[dim],it makes it easy to manage for pos thanks to nx
function GeneralDomainList(hx;periodic=Int[],periods=Float64[],T0=zeros(length(periodic)),lims=nothing)
    N = length(hx)
    x0 = zeros(N)
    nx = zeros(Int, length(periodic))
    hx = collect(hx)
    for (i,dim) in enumerate(periodic)
        nx[i] = round(periods[i]/hx[dim])
        hx[dim] = periods[i]./nx[i]
        x0[dim] = T0[i] + hx[dim]/2.0
    end

    grid = AB.GridFree(SVector{N,Float64}(x0), SVector{N,Float64}(hx))
    return GeneralDomainList(grid, Set{NTuple{N,Int}}(),periodic,periods,T0,nx,lims)
end
# it corrects the grid to be valid (with respect periodicity)
function GeneralDomainList(grid::AB.GridFree;periodic=Int[],periods=Float64[],T0=zeros(length(periodic)),lims=nothing)
    nx = zeros(Int, length(periodic))
    x0 = collect(grid.orig)
    hx = collect(grid.h)
    N = length(hx)
    for (i,dim) in enumerate(periodic)
        nx[i] = round(periods[i]/hx[dim])
        hx[dim] = periods[i]./nx[i]
        x0[dim] = T0[i] + hx[dim]/2.0
    end

    grid = AB.GridFree(SVector{N,Float64}(x0), SVector{N,Float64}(hx))
    return GeneralDomainList(grid, Set{NTuple{N,Int}}(),periodic,periods,T0,nx,lims)
end


function GeneralDomainList(domain::GeneralDomainList)
    domain2 = GeneralDomainList(domain.hx;periodic=domain.periodic,periods=domain.periods)
    return domain2
end
##
function set_in_period_coord(domain::AB.DomainList,x)
    return x
end
function set_in_period_pos(domain::AB.DomainList,pos)
    return pos
end
##

_coord_tuple(domain, i, j, t::Tuple{}) = tuple()
function _coord_tuple(domain, i, j, t::Tuple)
    el, rest = first(t), Base.tail(t)
    if i <= length(domain.periodic) && domain.periodic[i] == j
        el = domain.T0[i] + mod(el - domain.T0[i], domain.periods[i])
        i += 1
    end
    return tuple(el, _coord_tuple(domain, i, j + 1, rest)...)
end

function set_in_period_coord(domain::GeneralDomainList,x::SVector)
    #=
    for (i,dim) in enumerate(domain.periodic)  #doesnt work because of SVector
        x[dim] = mod(x[dim],domain.periods[i])
    end
    =#
    if isempty(domain.periodic)
        return x
    else
        return SVector(_coord_tuple(domain, 1, 1, x.data))
    end
end

_pos_tuple(domain, i, j, t::Tuple{}) = tuple()
function _pos_tuple(domain, i, j, t::Tuple)
    el, rest = first(t), Base.tail(t)
    if i <= length(domain.periodic) && domain.periodic[i] == j
        el = mod(el, domain.nx[i])
        i += 1
    end
    return tuple(el, _pos_tuple(domain, i, j + 1, rest)...)
end

function set_in_period_pos(domain::GeneralDomainList,pos)
    #=for (i,dim) in enumerate(domain.periodic)
        pos[dim] = mod(pos[dim],domain.nx[i])
    end=#
    if isempty(domain.periodic)
        return pos
    else
        return _pos_tuple(domain, 1, 1, pos)
    end
end

function AB.add_pos!(domain::GeneralDomainList, pos)
    pos = set_in_period_pos(domain,pos)
    push!(domain.elems, pos)
end

function AB.add_coord!(domain::GeneralDomainList, x)
    AB.add_pos!(domain, AB.get_pos_by_coord(domain.grid, x))
end

function AB.add_set!(domain::GeneralDomainList, rect::AB.HyperRectangle, incl_mode::AB.INCL_MODE)
    rectI = AB.get_pos_lims(domain.grid, rect, incl_mode)
    for pos in Iterators.product(AB._ranges(rectI)...)
        pos = set_in_period_pos(domain,pos)
        AB.add_pos!(domain, pos)
    end
end

#assuming that the rectangle is already in the domain for periodic dimensions
function AB.get_subset_pos(domain::GeneralDomainList,rect::AB.HyperRectangle,incl_mode::AB.INCL_MODE)
    lims = domain.lims
    if lims != nothing
        if any(rect.lb .> lims.ub) ||  any(rect.ub .< lims.lb)
            return []
        else
            rect = AB.HyperRectangle(max.(rect.lb,lims.lb),min.(rect.ub,lims.ub))
        end
    end
    posL = get_subset_pos2(domain,rect,incl_mode)
    return posL
end

function get_subset_pos2(domain::GeneralDomainList,rect::AB.HyperRectangle,incl_mode::AB.INCL_MODE)
    rectI = AB.get_pos_lims(domain.grid, rect, incl_mode)
    pos_iter = Iterators.product(AB._ranges(rectI)...)
    posL = []
    for pos in pos_iter
        pos = set_in_period_pos(domain,pos)
        if pos ∈ domain
            push!(posL, pos)
        end
    end
    return posL
end


function AB.add_subset!(domain1::GeneralDomainList, domain2::GeneralDomainList, rect::AB.HyperRectangle, incl_mode::AB.INCL_MODE)
    rectI = AB.get_pos_lims(domain1.grid, rect, incl_mode)
    pos_iter = Iterators.product(AB._ranges(rectI)...)
    if length(pos_iter) < AB.get_ncells(domain2)
        for pos in pos_iter
            pos = set_in_period_pos(domain1,pos)
            if pos ∈ domain2
                AB.add_pos!(domain1, pos)
            end
        end
    else
        for pos in AB.enum_pos(domain2)
            pos = set_in_period_pos(domain2,pos)
            if pos ∈ rectI
                AB.add_pos!(domain1, pos)
            end
        end
    end
end

function  AB.remove_pos!(domain::GeneralDomainList, pos)
    pos = set_in_period_pos(domain,pos)
    delete!(domain.elems, pos)
end

function  AB.remove_coord!(domain::GeneralDomainList, x)
    x = set_in_period_coord(domain,x)
    remove_pos!(domain, AB.get_pos_by_coord(domain.grid, x))
end

function  AB.remove_set!(domain::GeneralDomainList, rect::AB.HyperRectangle, incl_mode::AB.INCL_MODE)
    rectI = AB.get_pos_lims(domain.grid, rect, incl_mode)
    pos_iter = Iterators.product(AB._ranges(rectI)...)
    if length(pos_iter) < AB.get_ncells(domain)
        for pos in pos_iter
            pos = set_in_period_pos(domain,pos)
            AB.remove_pos!(domain, pos)
        end
    else
        for pos in AB.enum_pos(domain)
            pos = set_in_period_pos(domain,pos)
            if pos ∈ rectI
                AB.remove_pos!(domain, pos)
            end
        end
    end
end

function Base.union!(domain1::GeneralDomainList, domain2::GeneralDomainList)
    union!(domain1.elems, domain2.elems)
end

function Base.setdiff!(domain1::GeneralDomainList, domain2::GeneralDomainList)
    setdiff!(domain1.elems, domain2.elems)
end

function Base.empty!(domain::GeneralDomainList)
    empty!(domain.elems)
end

function Base.in(pos, domain::GeneralDomainList)
    pos = set_in_period_pos(domain,pos)
    return in(pos, domain.elems)
end

function Base.isempty(domain::GeneralDomainList)
    return isempty(domain.elems)
end

function Base.issubset(domain1::GeneralDomainList, domain2::GeneralDomainList)
    return issubset(domain1.elems, domain2.elems)
end

function  AB.get_ncells(domain::GeneralDomainList)
    return length(domain.elems)
end

function  AB.get_somepos(domain::GeneralDomainList)
    return first(domain.elems)
end

function  AB.enum_pos(domain::GeneralDomainList)
    return domain.elems
end

# make fit the grid exactly to the rectangle
function build_grid_in_rec(X,hx)
    N = length(hx)
    x0 = zeros(N)
    hx = collect(hx)
    for i=1:N
        L = X.ub[i]-X.lb[i]
        n = round(L/hx[i])
        hx[i] = L/n
        x0[i] = X.lb[i] + hx[i]/2.0
    end
    return AB.GridFree(SVector{N,Float64}(x0), SVector{N,Float64}(hx))
end


function one_direction(lb,ub,T,T0)
    if ub-lb>=T
        return [[T0,T0+T]]
    else
        lb = T0 + mod(lb-T0,T)
        ub = T0 + mod(ub-T0,T)
        if lb<=ub
            return [[lb,ub]]
        else
            return [[T0,ub],[lb,T0+T]]
        end
    end
end
function recursive(L,rec,lb,ub,periodic,periods,T0,i)
    if i>length(periodic)
        N = length(lb)
        push!(L,AB.HyperRectangle(SVector{N}(lb), SVector{N}(ub)))
        return
    end
    dim = periodic[i]
    intervals = one_direction(rec.lb[dim],rec.ub[dim],periods[i],T0[i])
    for interval in intervals
        l = copy(lb)
        u = copy(ub)
        l[dim] = interval[1]
        u[dim] = interval[2]
        recursive(L,rec,l,u,periodic,periods,T0,i+1)
    end
end
function set_rec_in_period(periodic,periods,T0,rec)
    L = []
    lb = collect(rec.lb)
    ub = collect(rec.ub)
    recursive(L,rec,lb,ub,periodic,periods,T0,1)
    return L
end

end
