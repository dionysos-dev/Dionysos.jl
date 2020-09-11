mutable struct BDDGridXDomain <: XDomain
    Rroot::Ref{Ptr{Node}}
    ncells::Int
end
mutable struct BDDGridUDomain <: UDomain
    Rroot::Ref{Ptr{Node}}
    ncells::Int
end

function _add_domain!(BDD::BDDManager, ::Type{DT}) where DT
    root = Ref(CUDD.Cudd_ReadLogicZero(BDD))
    _Ref(root)
    Rroot = Ref(root)
    push!(BDD.Rroots, Rroot)
    return DT(Rroot, 0)
end

AddXDomain!(mng::BDDGridManager) = _add_domain!(mng.BDD, BDDGridXDomain)
AddUDomain!(mng::BDDGridManager) = _add_domain!(mng.BDD, BDDGridUDomain)

# ---
Base.empty!(mng::BDDGridManager, domain::Domain) = (empty!(domain.cells); domain)
get_ncells(mng::BDDGridManager, domain) = length(domain.cells)
# Return an iterable of all cells
enum_cells(mng::BDDGridManager, domain) = domain.cells
get_some_cell(mng::BDDGridManager, domain) =
    isempty(domain.cells) ? nothing : first(domain.cells)
# ---

add_cell!(mng::BDDGridManager, domain::Domain, cell) = push!(domain.cells, cell)

function add_cells!(mng::BDDGridManager, domain::Domain,
        rect::HyperRectangle, incl_mode::INCL_MODE)
    disc = get_disc(mng, _XUType(domain))
    get_cell = pos2cell(disc)
    cells = domain.cells
    rectI = coord2pos_set(disc)(rect, incl_mode)
    for post in Iterators.product(_ranges(rectI)...)
        push!(cells, get_cell(post))
    end
end

# Return true if size of rectI <= min(MaxVal, typemax(Int))
function _compare_size(rectI::HyperRectangle, MaxVal)
    prodsize = 1
    for i in eachindex(rectI)
        (n, s) = _diff_trunc(rectI.ub[i], rectI.lb[i])
        s || return false
        n === 0 && return true
        MaxVal ÷ prodsize < n && return false
        prodsize *= n
    end
    return true
end

function add_cells!(mng::BDDGridManager, domain::D, sup_domain::D,
        rect::HyperRectangle, incl_mode::INCL_MODE) where {D<:Domain}
    disc = get_disc(mng, _XUType(domain))
    get_cell = pos2cell(disc)
    get_pos = cell2pos(disc)
    cells = domain.cells
    sup_cells = sup_domain.cells
    rectI = coord2pos_set(disc)(rect, incl_mode)
    if _compare_size(rectI, length(sup_cells))
        post_iter = Iterators.product(_ranges(rectI)...)
        for post in post_iter
            cell = get_cell(post)
            if cell ∈ sup_cells
                push!(cells, cell)
            end
        end
    else
        for cell in sup_cells
            if get_pos(cell) ∈ rectI
                push!(cells, cell)
            end
        end
    end
end

function remove_cells!(mng::BDDGridManager, domain::Domain,
        rect::HyperRectangle, incl_mode::INCL_MODE)
    disc = get_disc(mng, _XUType(domain))
    get_cell = pos2cell(disc)
    get_pos = cell2pos(disc)
    cells = domain.cells
    rectI = coord2pos_set(disc)(rect, incl_mode)
    if _compare_size(rectI, length(cells))
        post_iter = Iterators.product(_ranges(rectI)...)
        for post in post_iter
            delete!(cells, get_cell(post))
        end
    else
        for cell in cells
            if get_pos(cell) ∈ rectI
                delete!(cells, cell)
            end
        end
    end
end
