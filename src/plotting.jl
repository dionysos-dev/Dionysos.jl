# All plotting functions
module Plot

import ..Abstraction
AB = Abstraction

using LinearAlgebra
using StaticArrays
using LazySets
using Polyhedra
using CDDLib
using PyPlot

FC(c, a) =  matplotlib.colors.colorConverter.to_rgba(c, alpha = a)
const _mult = (SVector(-1,-1), SVector(-1,1), SVector(1,1), SVector(1,-1))

function verts_rect(c, h)
    return ntuple(i -> c + h.*_mult[i], 4)
end

function project(x, vars)
    return SVector(x[vars[1]], x[vars[2]])
end

# Cells
function subset!(ax, vars, subset::AB.SubSet{N,T};
        fc = "red", fa = 0.5, ec = "black", ea = 1.0, ew = 1.5) where {N,T}
    gridspace = subset.gridspace
    @assert length(vars) == 2 && N >= 2
    fca = FC(fc, fa)
    eca = FC(ec, ea)
    h = project(gridspace.h, vars)

    vertslist = NTuple{4,SVector{2,T}}[]

    for pos in unique(x -> x[vars], AB.enum_pos(subset))
        c = project(AB.get_coord_by_pos(gridspace, pos), vars)
        push!(vertslist, verts_rect(c, h/2.0))
    end

    polylist = matplotlib.collections.PolyCollection(vertslist)
    polylist.set_facecolor(fca)
    polylist.set_edgecolor(eca)
    polylist.set_linewidth(ew)
    ax.add_collection(polylist)
end

# Sets
function set!(ax, vars, rect::AB.HyperRectangle;
        fc = "green", fa = 0.5, ec = "black", ea = 1.0, ew = 1.5)
    @assert length(vars) == 2 && length(rect.lb) == length(rect.ub) >= 2
    c = (rect.lb + rect.ub)/2.0
    h = (rect.ub - rect.lb)/2.0
    polylist = matplotlib.collections.PolyCollection([verts_rect(c[vars], h[vars])])
    fca = FC(fc, fa)
    eca = FC(ec, ea)
    polylist.set_facecolor(fca)
    polylist.set_edgecolor(eca)
    polylist.set_linewidth(ew)
    ax.add_collection(polylist)
end

# Trajectory open loop
function trajectory_open_loop!(ax, vars, contsys::AB.ControlSystem{N,T}, x0, u, nstep;
        lc = "red", lw = 1.5, mc = "black", ms = 5.0, nsub = 5) where {N,T}
    @assert length(vars) == 2 && N >= 2
    tstep = contsys.tstep/nsub
    Nstep = nstep*nsub + 1
    X1list = Vector{T}(undef, Nstep)
    X2list = Vector{T}(undef, Nstep)
    X1list[1] = x0[vars[1]]
    X2list[1] = x0[vars[2]]
    x = x0

    for i = 2:Nstep
        x = contsys.sys_map(x, u, tstep)
        X1list[i] = x[vars[1]]
        X2list[i] = x[vars[2]]
    end

    idx = 1:nsub:Nstep
    ax.plot(X1list, X2list, lw = lw, c = lc)
    ax.plot(X1list[idx], X2list[idx],
        ls = "None", marker = ".", ms = ms, mfc = mc, mec = mc)
end

# Images
function cell_image!(ax, vars, Xsub, Usub, contsys::AB.ControlSystem{N,T};
        nsub = fill(5, N),
        fc = "blue", fa = 0.5, ec = "darkblue", ea = 1.0, ew = 1.5) where {N,T}
    @assert length(vars) == 2 && N >= 2
    fca = FC(fc, fa)
    eca = FC(ec, ea)
    vertslist = Vector{SVector{2,T}}[]
    ns = nsub .- 1

    for xpos in AB.enum_pos(Xsub), upos in AB.enum_pos(Usub)
        x = AB.get_coord_by_pos(Xsub.gridspace, xpos)
        u = AB.get_coord_by_pos(Usub.gridspace, upos)
        subpos_axes = ((0:ns[i])./ns[i] .- 0.5 for i = 1:length(ns))
        subpos_iter = Iterators.product(subpos_axes...)
        x_iter = (x + subpos.*Xsub.gridspace.h for subpos in subpos_iter)
        Fx_iter = (contsys.sys_map(x, u, contsys.tstep) for x in x_iter)
        push!(vertslist, convex_hull([project(Fx, vars) for Fx in Fx_iter][:]))
    end

    polylist = matplotlib.collections.PolyCollection(vertslist)
    polylist.set_facecolor(fca)
    polylist.set_edgecolor(eca)
    polylist.set_linewidth(ew)
    ax.add_collection(polylist)
end

# Outer-approximation
function cell_approx!(ax, vars, Xsub, Usub, contsys::AB.ControlSystemGrowth{N,T};
        fc = "yellow", fa = 0.5, ec = "gold", ea = 1.0, ew = 0.5) where {N,T}
    @assert length(vars) == 2 && N >= 2
    fca = FC(fc, fa)
    eca = FC(ec, ea)
    vertslist = NTuple{4,SVector{2,T}}[]
    r = Xsub.gridspace.h/2.0 + contsys.measnoise

    for xpos in AB.enum_pos(Xsub), upos in AB.enum_pos(Usub)
        x = AB.get_coord_by_pos(Xsub.gridspace, xpos)
        u = AB.get_coord_by_pos(Usub.gridspace, upos)
        Fx = contsys.sys_map(x, u, contsys.tstep)
        Fr = contsys.growthbound_map(r, u, contsys.tstep)
        Fr = Fr + contsys.measnoise
        push!(vertslist, verts_rect(project(Fx, vars), project(Fr, vars)))
    end

    polylist = matplotlib.collections.PolyCollection(vertslist)
    polylist.set_facecolor(fca)
    polylist.set_edgecolor(eca)
    polylist.set_linewidth(ew)
    ax.add_collection(polylist)
end

function cell_approx!(ax, vars, Xsub, Usub, contsys::AB.ControlSystemLinearized{N,T};
        fc = "yellow", fa = 0.5, ec = "gold", ea = 1.0, ew = 0.5) where {N,T}
    @assert length(vars) == 2 && N >= 2
    fca = FC(fc, fa)
    eca = FC(ec, ea)
    vertslist = Vector{SVector{2,T}}[]
    _H_ = SMatrix{N,N}(I).*(Xsub.gridspace.h/2.0)
    e = norm(Xsub.gridspace.h/2.0 + contsys.measnoise, Inf)

    for xpos in AB.enum_pos(Xsub), upos in AB.enum_pos(Usub)
        x = AB.get_coord_by_pos(Xsub.gridspace, xpos)
        u = AB.get_coord_by_pos(Usub.gridspace, upos)
        Fx, DFx = contsys.linsys_map(x, _H_, u, contsys.tstep)
        A = inv(DFx)
        Fr = contsys.measnoise .+ contsys.error_map(e, u, contsys.tstep)
        b = abs.(A)*Fr .+ 1.0
        HP = HPolytope([A; -A], vcat(b, b))
        verts1 = vertices_list(HP, backend=CDDLib.Library())
        verts2 = [project(Fx + v, vars) for v in verts1]
        push!(vertslist, convex_hull(verts2))
    end

    polylist = matplotlib.collections.PolyCollection(vertslist)
    polylist.set_facecolor(fca)
    polylist.set_edgecolor(eca)
    polylist.set_linewidth(ew)
    ax.add_collection(polylist)
end

# Trajectory closed loop
function trajectory_closed_loop!(ax, vars, contsys, symmodel, contr, x0, nstep;
        lc = "red", lw = 1.5, mc = "black", ms = 5.0, nsub = 5, randchoose = false)
    Xgrid = symmodel.Xgrid
    Usub = AB.NewSubSet(symmodel.Ugrid)
    Ysub = AB.NewSubSet(symmodel.Xgrid)
    for i = 1:nstep
        xpos = AB.get_pos_by_coord(Xgrid, x0)
        if !(xpos ∈ Xgrid)
            @warn("Trajectory out of domain")
            return
        end
        source = AB.get_state_by_xpos(symmodel, xpos)
        symbollist = Int[]
        AB.compute_enabled_symbols!(symbollist, contr, source)
        if isempty(symbollist)
            @warn("Uncontrollable state")
            return
        end
        if randchoose
            symbol = rand(symbollist)
        else
            symbol = symbollist[1]
        end
        upos = AB.get_upos_by_symbol(symmodel, symbol)
        u = AB.get_coord_by_pos(symmodel.Ugrid, upos)
        Plot.trajectory_open_loop!(ax, vars, contsys, x0, u, 1,
            lc = lc, lw = lw, mc = mc, ms = ms, nsub = nsub)
        x0 = contsys.sys_map(x0, u, contsys.tstep)
    end
end
end  # Plot
