using Polyhedra

abstract type Grid{N,T} end

"""
    GridFree{N,T} <: Grid{N,T}

Uniform grid on unbounded space, centered at `orig` and with steps set by the vector `h`.
"""
struct GridFree{N,T} <: Grid{N,T}
    orig::SVector{N,T}
    h::SVector{N,T}
end


"""
    GridFree{N,T} <: Grid{N,T}

Uniform grid on rectagular space `rect`, centered at `orig` and with steps set by the vector `h`.
"""
struct GridRectangular{N,T} <: Grid{N,T}
    orig::SVector{N,T}
    h::SVector{N,T}
    rect::Any
end


"""
    GridFree{N,T} <: Grid{N,T}

Uniform grid on rectagular space `rect`, centered at `orig` and with steps set by the vector `h`.
Cells are (possibly overlapping) ellipsoids defined at each grid point `c` as `(x-c)'P(x-c) ≤ 1`
"""
struct GridEllipsoidalRectangular{N,T} <: Grid{N,T}
    orig::SVector{N,T}
    h::SVector{N,T}
    P::SMatrix{N,N}
    rect::Any
end


function get_pos_by_coord(grid::Grid{N}, x) where N
    return ntuple(i -> round(Int, (x[i] - grid.orig[i])/grid.h[i]), Val(N))
end

function get_all_pos_by_coord(grid::GridEllipsoidalRectangular{N}, x) where N
    center = get_pos_by_coord(grid,x)
    all_pos = typeof(center)[]
    for dpos in Iterators.product(eachrow(repeat([-1 0 1],N))...)
        coord = get_coord_by_pos(grid, dpos.+center)
        if (x-coord)'grid.P*(x-coord) ≤ 1
            push!(all_pos, (dpos.+center))
        end
    end
    return all_pos
end


function get_coord_by_pos(grid, pos)
    return grid.orig + pos.*grid.h
end

function get_pos_lims_inner(grid::Grid{N}, rect; tol=1e-6) where N
    lbI = ntuple(i -> ceil(Int, (rect.lb[i] - tol - grid.orig[i])/grid.h[i] + 0.5), Val(N))
    ubI = ntuple(i -> floor(Int, (rect.ub[i] + tol - grid.orig[i])/grid.h[i] - 0.5), Val(N))
    return UT.HyperRectangle(lbI, ubI)
end

function get_pos_lims_outer(grid::Grid{N}, rect; tol=0.0) where N
    lbI = ntuple(i -> ceil(Int, (rect.lb[i] + tol - grid.orig[i])/grid.h[i] - 0.5), Val(N))
    ubI = ntuple(i -> floor(Int, (rect.ub[i] - tol - grid.orig[i])/grid.h[i] + 0.5), Val(N))
    return UT.HyperRectangle(lbI, ubI)
end

function get_pos_lims(grid, rect, incl_mode::INCL_MODE)
    if incl_mode == INNER
        return get_pos_lims_inner(grid, rect)
    else
        return get_pos_lims_outer(grid, rect)
    end
end

function _ranges(rect::UT.HyperRectangle{NTuple{N,T}}) where {N,T}
    return ntuple(i -> UnitRange(rect.lb[i], rect.ub[i]), Val(N))
end

function rectangle(c,r)
    Shape(c[1].-r[1] .+ [0,2*r[1],2*r[1],0], c[2].-r[2] .+ [0,0,2*r[2],2*r[2]])
end

function get_rec(grid::GridFree, pos)
    x = get_coord_by_pos(grid, pos)
    r = grid.h/2.0
    return UT.HyperRectangle(x-r, x+r)
end

function get_dim(grid::GridFree)
    return length(grid.orig)
end

function get_h(grid::GridFree)
    return grid.h
end

function get_origin(grid::GridFree)
    return grid.h
end

function get_volume(grid::GridFree)
    r = get_h(grid)/2.0
    return UT.volume(UT.HyperRectangle(-r,r))
end

function sample_elem(grid::GridFree, xpos, N::Int)
    x = get_coord_by_pos(grid, xpos)
    r = grid.h/2
    rec = UT.HyperRectangle(x .- r, x .+ r)
    return UT.sample_from_rec(rec,N)
end

function plot_elem!(grid::GridFree, pos; dims=[1,2], opacity=.9, color=:yellow)
    center = get_coord_by_pos(grid, pos)
    h = grid.h[dims]
    plot!(rectangle(center[dims],h./2), opacity=opacity,color=color)
end

######################## deformed grid ########################

# f is an inversible function
struct DeformedGrid{N,T} <: Grid{N,T}
    grid::GridFree{N,T}
    f::Function
    fi::Function
    A
end

function DeformedGrid(grid::GridFree{N,T},f::Function,fi::Function;A=nothing) where {N,T}
    return DeformedGrid(grid,f,fi,A)
end

function get_pos_by_coord(Dgrid::DeformedGrid{N}, x) where N
    return get_pos_by_coord(Dgrid.grid, Dgrid.fi(x))
end

function get_coord_by_pos(Dgrid::DeformedGrid{N}, pos) where N
    return Dgrid.f(get_coord_by_pos(Dgrid.grid, pos))
end

function get_pos_lims_inner(Dgrid::DeformedGrid{N}, rect; tol=1e-6) where N
    return get_pos_lims_inner(Dgrid.grid, rect; tol=tol)
end

function get_pos_lims_outer(Dgrid::DeformedGrid{N}, rect; tol=1e-6) where N
    return get_pos_lims_outer(Dgrid.grid, rect; tol=tol)
end

function get_dim(Dgrid::DeformedGrid)
    return get_dim(Dgrid.grid)
end

function get_h(Dgrid::DeformedGrid)
    return get_h(Dgrid.grid)
end

function get_origin(Dgrid::DeformedGrid)
    return get_origin(Dgrid.grid)
end

# only for linear transformation of the grid
function get_volume(Dgrid::DeformedGrid)
    if Dgrid.A != nothing
        return abs(det(Dgrid.A))*get_volume(Dgrid.grid)
    else
        println("volume is state-dependant for nonlinear transformation")
        return get_volume(Dgrid.grid)
    end
end

function sample_elem(Dgrid::DeformedGrid, xpos, N::Int)
    points = sample_elem(Dgrid.grid, xpos, N)
    return [Dgrid.f(x) for x in points]
end

function plot_deformed_rectangle!(rec,f;dims=[1,2],opacity=0.9,color=:yellow,N=2)
     lb = rec.lb[dims]
     ub = rec.ub[dims]
     vertices = [SVector(lb[1],lb[2]), SVector(lb[1],ub[2]), SVector(ub[1],lb[2]), SVector(ub[1],ub[2])]
     points = SVector[]
     for x in LinRange(lb[1],ub[1],N)
         push!(points,f(SVector(x,lb[2])))
     end
     for x in LinRange(lb[2],ub[2],N)
         push!(points,f(SVector(ub[1],x)))
     end
     for x in LinRange(ub[1],lb[1],N)
         push!(points,f(SVector(x,ub[2])))
     end
     for x in LinRange(ub[2],lb[2],N)
         push!(points,f(SVector(lb[1],x)))
     end
     unique!(points)
     x = [point[1] for point in points]
     y = [point[2] for point in points]
     plot!(Shape(x,y),opacity=opacity,color=color)
end

function plot_elem!(Dgrid::DeformedGrid, pos; dims=[1,2], opacity=1.0, color=:yellow, N=8)
    rec = get_rec(Dgrid.grid,pos)
    plot_deformed_rectangle!(rec,Dgrid.f;dims=dims,opacity=opacity,color=color, N=N)
end
