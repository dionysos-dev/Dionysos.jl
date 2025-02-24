"""
    RectangularObstacles{VT} <: AbstractSet{VT}

Struct for a rectangular domain with rectangular obstacles.
"""
struct RectangularObstacles{VT} <: AbstractSet{VT}
    X::UT.HyperRectangle{VT}
    O::Union{Any, Vector{UT.HyperRectangle{VT}}}
end

function Base.in(pos, dom::RectangularObstacles)
    return !mapreduce(Base.Fix1(in, pos), |, dom.O; init = !in(pos, dom.X))
end

function _fit_grid(elems::RectangularObstacles, grid, nx = nothing, fit = false)
    if fit == true
        N = length(nx)
        lbI = ntuple(i -> 0, Val(N))
        ubI = ntuple(i -> nx[i] - 1, Val(N))
        return RectangularObstacles(
            UT.HyperRectangle(lbI, ubI),
            [get_pos_lims_outer(grid, Oi) for Oi in elems.O],
        )
    else
        return RectangularObstacles(
            get_pos_lims_outer(grid, elems.X),
            [get_pos_lims_outer(grid, Oi) for Oi in elems.O],
        )
    end
end
_fit_grid(elems::Set, grid, nx, fit) = elems

## add the periodicity in the domain (add into Domain.jl)
"""
    GeneralDomainList{N,E<:AbstractSet{NTuple{N,Int}},T,S<:Grid{N,T},F} <: DomainType{N,T}

Struct for a rectangular domain with rectangular obstacles.
"""
struct GeneralDomainList{N, E <: AbstractSet{NTuple{N, Int}}, T, S <: Grid{N, T}, F} <:
       GridDomainType{N, T}
    grid::S
    elems::E
    periodic::Vector{Int} # components which are periodic
    periods::Vector{Float64}  # periods
    T0::Vector{Float64}
    nx::Vector{Int} # number of cell in all direction (if fit) //the periodic directions
    lims::Any#::Union{Nothing,UT.HyperRectangle} # lower and upper bound on the non periodic dimensions
    # can be used to be more efficient to compute the cells which belong to a
    # given hypperrectangle and the domain if the rectangle is outside the domain
    elemsCoord::F
    fit::Bool
end

# hx is the desired step size, but for periodic direction, it can be changed to fit exactly to the period
# (for periodic dims, hx has to fit exactly in the period)
# The periods are [T0[i], T0[i] + periods[i]]
# for periodic dimensions, I set the origin in T0[dim],it makes it easy to manage for pos thanks to nx
function GeneralDomainList(
    hx;
    elems = Set{NTuple{length(hx), Int}}(),
    periodic = Int[],
    periods = Float64[],
    T0 = zeros(length(periodic)),
    lims = nothing,
    fit = false,
    f = nothing,
    fi = nothing,
    A = nothing,
)
    N = length(hx)
    x0 = zeros(N)
    nx = zeros(Int, N)

    hx = collect(hx)
    if fit == true
        X = elems.X
        for i in 1:N
            if !(i in periodic)
                l = X.ub[i] - X.lb[i]
                nx[i] = round(l / hx[i])
                hx[i] = l ./ nx[i]
                x0[i] = X.lb[i] + hx[i] / 2.0
            end
        end
    end

    for (i, dim) in enumerate(periodic)
        nx[dim] = round(periods[i] / hx[dim])
        hx[dim] = periods[i] ./ nx[dim]
        x0[dim] = T0[i] + hx[dim] / 2.0
    end
    grid = GridFree(SVector{N, Float64}(x0), SVector{N, Float64}(hx))
    if f !== nothing
        grid = DeformedGrid(grid, f, fi, A)
    end
    return GeneralDomainList(
        grid,
        _fit_grid(elems, grid, nx, fit),
        periodic,
        periods,
        T0,
        nx,
        lims,
        elems,
        fit,
    )
end
#it corrects the grid to be valid (with respect periodicity)
function GeneralDomainList(
    grid::Grid{N},
    elems = Set{NTuple{N, Int}}();
    periodic = Int[],
    periods = Float64[],
    T0 = zeros(length(periodic)),
    lims = nothing,
) where {N}
    nx = zeros(Int, length(periodic))
    for (i, dim) in enumerate(periodic)
        nx = zeros(Int, length(periodic))
        x0 = collect(get_origin(grid))
        hx = collect(get_h(grid))
        nx[i] = round(periods[i] / hx[dim])
        hx[dim] = periods[i] ./ nx[i]
        x0[dim] = T0[i] + hx[dim] / 2.0
        grid = GridFree(SVector{N, Float64}(x0), SVector{N, Float64}(hx))
    end
    return GeneralDomainList(
        grid,
        _fit_grid(elems, grid),
        periodic,
        periods,
        T0,
        nx,
        lims,
        elems,
        false,
    )
end

function GeneralDomainList(domain::GeneralDomainList)
    domain2 =
        GeneralDomainList(domain.hx; periodic = domain.periodic, periods = domain.periods)
    return domain2
end

function get_dim(domain::GeneralDomainList)
    return get_dim(domain.grid)
end

function get_grid(domain::GeneralDomainList)
    return domain.grid
end

#
function set_in_period_coord(domain::DomainList, x)
    return x
end
function set_in_period_pos(domain::DomainList, pos)
    return pos
end
#

_coord_tuple(domain, i, j, t::Tuple{}) = tuple()
function _coord_tuple(domain, i, j, t::Tuple)
    el, rest = first(t), Base.tail(t)
    if i <= length(domain.periodic) && domain.periodic[i] == j
        el = domain.T0[i] + mod(el - domain.T0[i], domain.periods[i])
        i += 1
    end
    return tuple(el, _coord_tuple(domain, i, j + 1, rest)...)
end

function set_in_period_coord(domain::GeneralDomainList, x::SVector)
    if isempty(domain.periodic)
        return x
    else
        return SVector(_coord_tuple(domain, 1, 1, x.data))
    end
end

function _pos_tuple(domain, pos)
    pos = collect(pos)
    for i in domain.periodic
        pos[i] = mod(pos[i], domain.nx[i])
    end
    return Tuple(pos)
end

function set_in_period_pos(domain::GeneralDomainList, pos)
    if isempty(domain.periodic)
        return pos
    else
        return _pos_tuple(domain, pos)#_pos_tuple(domain, 1, 1, pos)
    end
end

#assuming that the rectangle is already in the domain for periodic dimensions
function get_subset_pos(
    domain::GeneralDomainList{N},
    rect::UT.HyperRectangle,
    incl_mode::INCL_MODE,
) where {N}
    lims = domain.lims
    if lims !== nothing
        if any(rect.lb .> lims.ub) || any(rect.ub .< lims.lb)
            return NTuple{N, Int}[]
        else
            rect = UT.HyperRectangle(max.(rect.lb, lims.lb), min.(rect.ub, lims.ub))
        end
    end
    posL = get_subset_pos2(domain, rect, incl_mode)
    return posL
end

function get_subset_pos2(
    domain::GeneralDomainList{N},
    rect::UT.HyperRectangle,
    incl_mode::INCL_MODE,
) where {N}
    rectI = get_pos_lims(domain.grid, rect, incl_mode)
    pos_iter = Iterators.product(_ranges(rectI)...)
    return NTuple{N, Int}[
        set_in_period_pos(domain, pos) for pos in pos_iter if pos ∈ domain
    ]
    return posL
end

function Base.in(pos, domain::GeneralDomainList)
    pos = set_in_period_pos(domain, pos)
    return in(pos, domain.elems)
end

function Base.issubset(domain1::GeneralDomainList, domain2::GeneralDomainList)
    return issubset(domain1.elems, domain2.elems)
end

function get_pos(domain::GeneralDomainList, elems::Set)
    return elems
end

function get_pos(domain::GeneralDomainList, elems::RectangularObstacles)
    posL = get_subset_pos(domain, domain.elemsCoord.X, INNER)
    L = []
    for pos in posL
        if pos in domain.elems
            push!(L, pos)
        end
    end
    return L
end

function enum_pos(domain::GeneralDomainList)
    posL = get_pos(domain, domain.elems)
    return posL
end

function get_ncells(domain::GeneralDomainList)
    return length(enum_pos(domain))
end

###############################################################################
#fonctionne seulement si elems quand on crée GeneralDomainList est un set classique

function add_pos!(domain::GeneralDomainList, pos)
    pos = set_in_period_pos(domain, pos)
    return push!(domain.elems, pos)
end

function add_coord!(domain::GeneralDomainList, x)
    return add_pos!(domain, get_pos_by_coord(domain.grid, x))
end

function get_pos_by_coord(domain::GeneralDomainList, x)
    pos = get_pos_by_coord(domain.grid, x)
    return set_in_period_pos(domain, pos)
end

function add_set!(domain::GeneralDomainList, rect::UT.HyperRectangle, incl_mode::INCL_MODE)
    rectI = get_pos_lims(domain.grid, rect, incl_mode)
    for pos in Iterators.product(_ranges(rectI)...)
        pos = set_in_period_pos(domain, pos)
        add_pos!(domain, pos)
    end
end

function add_subset!(
    domain1::GeneralDomainList,
    domain2::GeneralDomainList,
    rect::UT.HyperRectangle,
    incl_mode::INCL_MODE,
)
    rectI = get_pos_lims(domain1.grid, rect, incl_mode)
    pos_iter = Iterators.product(_ranges(rectI)...)
    if length(pos_iter) < get_ncells(domain2)
        for pos in pos_iter
            pos = set_in_period_pos(domain1, pos)
            if pos ∈ domain2
                add_pos!(domain1, pos)
            end
        end
    else
        for pos in enum_pos(domain2)
            pos = set_in_period_pos(domain2, pos)
            if pos ∈ rectI
                add_pos!(domain1, pos)
            end
        end
    end
end

function remove_pos!(domain::GeneralDomainList, pos)
    pos = set_in_period_pos(domain, pos)
    return delete!(domain.elems, pos)
end

function remove_coord!(domain::GeneralDomainList, x)
    x = set_in_period_coord(domain, x)
    return remove_pos!(domain, get_pos_by_coord(domain.grid, x))
end

function remove_set!(
    domain::GeneralDomainList,
    rect::UT.HyperRectangle,
    incl_mode::INCL_MODE,
)
    rectI = get_pos_lims(domain.grid, rect, incl_mode)
    pos_iter = Iterators.product(_ranges(rectI)...)
    if length(pos_iter) < get_ncells(domain)
        for pos in pos_iter
            pos = set_in_period_pos(domain, pos)
            remove_pos!(domain, pos)
        end
    else
        for pos in enum_pos(domain)
            pos = set_in_period_pos(domain, pos)
            if pos ∈ rectI
                remove_pos!(domain, pos)
            end
        end
    end
end

function Base.union!(domain1::GeneralDomainList, domain2::GeneralDomainList)
    return union!(domain1.elems, domain2.elems)
end

function Base.setdiff!(domain1::GeneralDomainList, domain2::GeneralDomainList)
    return setdiff!(domain1.elems, domain2.elems)
end

function Base.empty!(domain::GeneralDomainList)
    return empty!(domain.elems)
end

function Base.isempty(domain::GeneralDomainList)
    return isempty(domain.elems)
end

###############################################################################

# make fit the grid exactly to the rectangle
function build_grid_in_rec(X, hx)
    N = length(hx)
    x0 = zeros(N)
    hx = collect(hx)
    for i in 1:N
        L = X.ub[i] - X.lb[i]
        n = ceil(L / hx[i])
        hx[i] = L / n
        x0[i] = X.lb[i] + hx[i] / 2.0
    end
    return GridFree(SVector{N, Float64}(x0), SVector{N, Float64}(hx))
end

function one_direction(lb, ub, T, T0)
    if ub - lb >= T
        return [(T0, T0 + T)]
    else
        lb = T0 + mod(lb - T0, T)
        ub = T0 + mod(ub - T0, T)
        if lb <= ub
            return [(lb, ub)]
        else
            return [(T0, ub), (lb, T0 + T)]
        end
    end
end
function recursive(L, rec, lb, ub, periodic, periods, T0, i)
    N = length(lb)
    if i > length(periodic)
        push!(L, UT.HyperRectangle(SVector(lb), SVector(ub)))
        return
    end
    dim = periodic[i]
    intervals = one_direction(rec.lb[dim], rec.ub[dim], periods[i], T0[i])
    for interval in intervals
        l = ntuple(i -> i == dim ? interval[1] : lb[i], Val(N))
        u = ntuple(i -> i == dim ? interval[2] : ub[i], Val(N))
        recursive(L, rec, l, u, periodic, periods, T0, i + 1)
    end
end
function set_rec_in_period(periodic, periods, T0, rec::UT.HyperRectangle)
    L = typeof(rec)[]
    recursive(L, rec, rec.lb, rec.ub, periodic, periods, T0, 1)
    return L
end

# ################### symbolic model
# function _SymbolicModel(Xdom::GeneralDomainList{N,RectangularObstacles{NTuple{N,T}}}, Udom::Domain{M}) where {N,M,T}
#     nu = get_ncells(Udom)
#     uint2pos = [pos for pos in enum_pos(Udom)]
#     upos2int = Dict((pos, i) for (i, pos) in enumerate(enum_pos(Udom)))
#     symmodel = SymbolicModelList(
#         Xdom,
#         Udom,
#         AutomatonList{Set{NTuple{3,T}}}(0, nu),
#         Dict{NTuple{N,T},Int}(),
#         NTuple{N,T}[],
#         upos2int,
#         uint2pos,
#     )
# end
#
#
# function get_state_by_xpos(
#     symmodel::SymbolicModelList{N,M,<:GeneralDomainList{N,RectangularObstacles{NTuple{N,T}}}},
#     pos,
# ) where {N,M,T}
#     #pos = set_in_period_pos(domain,pos)
#     id = get(symmodel.xpos2int, pos, nothing)
#     created = false
#     if id === nothing
#         if pos in symmodel.Xdom
#             created = true
#             push!(symmodel.xint2pos, pos)
#             id = length(symmodel.xint2pos)
#             symmodel.xpos2int[pos] = id
#             i = HybridSystems.add_state!(symmodel.autom)
#             @assert i == id
#         else
#             error("$pos is not in state domain $(symmodel.Xdom)")
#         end
#     end
#     return id::Int,created
# end
