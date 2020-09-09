# All plotting functions
module Plot

import ..Abstraction
ABS = Abstraction

using LinearAlgebra
using StaticArrays
using LazySets
using Polyhedra
using CDDLib
using PyPlot

#---
FC(c, a) =  matplotlib.colors.colorConverter.to_rgba(c, alpha = a)
const _mult = (SVector(-1,-1), SVector(-1,1), SVector(1,1), SVector(1,-1))

function verts_rect(c, h)
    return ntuple(i -> c + h.*_mult[i], 4)
end
#---

#--- Cells
function _plot_cells!(ax, vars::SVector{2,Int},
        disc::ABS.ListGridDisc, cell_iter, fca, eca, ew)
    h = disc.grid.h[vars]
    get_coord = ABS.cell2coord(disc)
    get_pos = ABS.cell2pos(disc)
    _, T = ABS._spacetype(typeof(disc))

    vertslist = NTuple{4,SVector{2,T}}[]
    pos2_recorder = Set{SVector{2,Int}}()

    for cell in cell_iter
        pos2 = get_pos(cell)[vars]
        pos2 ∈ pos2_recorder && continue
        push!(pos2_recorder, pos2)
        c = get_coord(cell)[vars]
        push!(vertslist, verts_rect(c, h/2.0))
    end

    polylist = matplotlib.collections.PolyCollection(vertslist)
    polylist.set_facecolor(fca)
    polylist.set_edgecolor(eca)
    polylist.set_linewidth(ew)
    ax.add_collection(polylist)
end

function cells!(ax, vars, mng::ABS.Manager, domain::ABS.Domain;
        fc = "red", fa = 0.5, ec = "black", ea = 1.0, ew = 1.5)
    fca = FC(fc, fa)
    eca = FC(ec, ea)
    _vars = SVector{2,Int}(vars)
    disc = ABS.get_disc(mng, ABS._XUType(domain))
    cell_iter = ABS.enum_cells(mng, domain)
    _plot_cells!(ax, _vars, disc, cell_iter, fca, eca, ew)
end

function cells!(ax, vars, mng::ABS.Manager, symbset::ABS.SymbolSet;
        fc = "red", fa = 0.5, ec = "black", ea = 1.0, ew = 1.5)
    fca = FC(fc, fa)
    eca = FC(ec, ea)
    _vars = SVector{2,Int}(vars)
    disc = ABS.get_disc(mng, ABS._XUType(symbset))
    get_cell = ABS.symbol2cell(ABS.get_map(mng, ABS._XUType(symbset)))
    cell_iter = Base.Generator(get_cell, ABS.enum_symbols(mng, symbset))
    _plot_cells!(ax, _vars, disc, cell_iter, fca, eca, ew)
end

#--- Images
function _plot_image!(ax, vars::SVector{2,Int},
        mng::ABS.GridManager, xcell_iter, ucell_iter, contsys,
        nsub::SVector{N,Int}, fca, eca, ew) where N
    XDisc = ABS.get_disc(mng, ABS.XType)
    UDisc = ABS.get_disc(mng, ABS.UType)
    h = XDisc.grid.h[vars]
    get_xcoord_center = ABS.cell2coord(XDisc)
    get_ucoord_center = ABS.cell2coord(UDisc)
    _, T = ABS._spacetype(typeof(XDisc))

    vertslist = Vector{SVector{2,T}}[]
    rectI = ABS.HyperRectangle(zeros(SVector{N,Int}), nsub)
    sub_iter = Iterators.product(ABS._ranges(rectI)...)
    size = 1.0 ./ (nsub .- 1.0)
    sub2pos(sub) = (sub.*size .- 0.5)

    for xcell in xcell_iter, ucell in ucell_iter
        x = get_xcoord_center(xcell)
        u = get_ucoord_center(ucell)
        sub2Fcoord(sub) = contsys.sys_map(x + sub2pos(sub).*h, u, contsys.tstep)
        Fx_iter = Base.Generator(sub2Fcoord, sub_iter)
        push!(vertslist, convex_hull([Fx[vars] for Fx in Fx_iter][:]))
    end

    polylist = matplotlib.collections.PolyCollection(vertslist)
    polylist.set_facecolor(fca)
    polylist.set_edgecolor(eca)
    polylist.set_linewidth(ew)
    ax.add_collection(polylist)
end

function images!(ax, vars, mng::ABS.Manager, XDom::ABS.XDomain, UDom::ABS.UDomain,
        contsys::ABS.ControlSystem{N,T};
        nsub = 5*ones(SVector{N,Int}),
        fc = "blue", fa = 0.5, ec = "darkblue", ea = 1.0, ew = 1.5) where {N,T}
    fca = FC(fc, fa)
    eca = FC(ec, ea)
    xcell_iter = ABS.enum_cells(mng, XDom)
    ucell_iter = ABS.enum_cells(mng, UDom)
    _vars = SVector{2,Int}(vars)
    _nsub = SVector{N,Int}(nsub)
    _plot_image!(ax, _vars, mng, xcell_iter, ucell_iter, contsys, _nsub, fca, eca, ew)
end

function images!(ax, vars, mng::ABS.Manager, stateset::ABS.StateSet, labelset::ABS.LabelSet,
        contsys::ABS.ControlSystem{N,T};
        nsub = 5*ones(SVector{N,Int}),
        fc = "blue", fa = 0.5, ec = "darkblue", ea = 1.0, ew = 1.5) where {N,T}
    fca = FC(fc, fa)
    eca = FC(ec, ea)
    _vars = SVector{2,Int}(vars)
    get_xcell = ABS.symbol2cell(ABS.get_map(mng, ABS.XType))
    get_ucell = ABS.symbol2cell(ABS.get_map(mng, ABS.UType))
    xcell_iter = Base.Generator(get_xcell, ABS.enum_symbols(mng, stateset))
    ucell_iter = Base.Generator(get_ucell, ABS.enum_symbols(mng, labelset))
    _nsub = SVector{N,Int}(nsub)
    _plot_image!(ax, _vars, mng, xcell_iter, ucell_iter, contsys, _nsub, fca, eca, ew)
end

# Trajectory open loop
function _plot_trajectory!(ax, vars::SVector{2,Int}, contsys, x0, u, nstep;
        nsub, lc, lw, mc, ms)
    tstep = contsys.tstep/nsub
    Nstep = nstep*nsub + 1
    X1list = Vector{eltype(x0)}(undef, Nstep)
    X2list = Vector{eltype(x0)}(undef, Nstep)
    X1list[1] = x0[vars[1]]
    X2list[1] = x0[vars[2]]
    x = x0

    for i in 2:Nstep
        x = contsys.sys_map(x, u, tstep)
        X1list[i] = x[vars[1]]
        X2list[i] = x[vars[2]]
    end

    idx = 1:nsub:Nstep
    ax.plot(X1list, X2list, lw = lw, c = lc)
    ax.plot(X1list[idx], X2list[idx],
        ls = "None", marker = ".", ms = ms, mfc = mc, mec = mc)
end

function trajectory!(ax, vars,
        contsys::ABS.ControlSystem{N,T}, x0::SVector{N,T}, u, nstep;
        nsub = 5, lc = "red", lw = 1.5, mc = "black", ms = 5.0) where {N,T}
    _vars = SVector{2,Int}(vars)
    _plot_trajectory!(ax, _vars, contsys, x0, u, nstep; nsub, lc, lw, mc, ms)
end

function _plot_controlled_trajectory!(ax, vars::SVector{2,Int},
        mng, contr, contsys, XDisc, XMap, UDisc, UMap, x0, nstep;
        nsub, lc, lw, mc, ms)
    get_xcell_by_coord = ABS.coord2cell(XDisc)
    get_state_by_xcell_try = ABS.cell2symbol_try(XMap)
    get_ucell_by_label = ABS.symbol2cell(UMap)
    get_ucoord_center = ABS.cell2coord(UDisc)
    for i in 1:nstep
        xcell = get_xcell_by_coord(x0)
        state, is_valid = get_state_by_xcell_try(xcell)
        if !is_valid
            @warn("Trajectory out of symbolic domain")
            return
        end
        label_iter = ABS.enum_enabled_labels(mng, contr, xcell)
        if isempty(label_iter)
            @warn("Uncontrollable state")
            return
        end
        label = first(label_iter)
        ucell = get_ucell_by_label(label)
        u = get_ucoord_center(ucell)
        Plot._plot_trajectory!(ax, vars, contsys, x0, u, 1,
            nsub = nsub, lc = lc, lw = lw, mc = mc, ms = ms)
        x0 = contsys.sys_map(x0, u, contsys.tstep)
    end
end

function trajectory!(ax, vars,
        mng::ABS.Manager, contr::ABS.Controller, contsys::ABS.ControlSystem{N,T},
        x0::SVector{N,T}, nstep;
        nsub = 5, lc = "red", lw = 1.5, mc = "black", ms = 5.0) where {N,T}
    _vars = SVector{2,Int}(vars)
    XDisc = ABS.get_disc(mng, ABS.XType)
    UDisc = ABS.get_disc(mng, ABS.UType)
    XMap = ABS.get_map(mng, ABS.XType)
    UMap = ABS.get_map(mng, ABS.UType)
    get_xcell_by_coord = ABS.coord2cell(XDisc)
    get_state_by_xcell_try = ABS.cell2symbol_try(XMap)
    get_ucell_by_label = ABS.symbol2cell(UMap)
    get_ucoord_center = ABS.cell2coord(UDisc)
    for i in 1:nstep
        xcell = get_xcell_by_coord(x0)
        state, is_valid = get_state_by_xcell_try(xcell)
        if !is_valid
            @warn("Trajectory out of symbolic domain")
            return
        end
        label_iter = ABS.enum_enabled_labels(mng, contr, state)
        if isempty(label_iter)
            @warn("Uncontrollable state")
            return
        end
        label = first(label_iter)
        ucell = get_ucell_by_label(label)
        u = get_ucoord_center(ucell)
        Plot._plot_trajectory!(ax, _vars, contsys, x0, u, 1,
            nsub = nsub, lc = lc, lw = lw, mc = mc, ms = ms)
        x0 = contsys.sys_map(x0, u, contsys.tstep)
    end
end

#=
#--- Sets
function set!(ax, vars, rect::ABS.HyperRectangle;
        fc = "green", fa = 0.5, ec = "black", ea = 1.0, ew = 1.5)
    @assert length(vars) == 2 && length(rect.lb) == length(rect.ub) >= 2
    c = (rect.lb + rect.ub)/2.0
    h = (rect.ub - rect.lb)/2.0
    polylist = matplotlib.collections.PolyCollection((verts_rect(c[vars], h[vars]),))
    fca = FC(fc, fa)
    eca = FC(ec, ea)
    polylist.set_facecolor(fca)
    polylist.set_edgecolor(eca)
    polylist.set_linewidth(ew)
    ax.add_collection(polylist)
end



# Outer-approximation
function cell_approx!(ax, vars, Xdom, Udom, contsys::ABS.ControlSystemGrowth{N,T};
        fc = "yellow", fa = 0.5, ec = "gold", ea = 1.0, ew = 0.5) where {N,T}
    @assert length(vars) == 2 && N >= 2
    fca = FC(fc, fa)
    eca = FC(ec, ea)
    vertslist = NTuple{4,SVector{2,T}}[]
    r = Xdom.grid.h/2.0 + contsys.measnoise

    for xpos in ABS.enum_pos(Xdom), upos in ABS.enum_pos(Udom)
        x = ABS.get_coord_by_pos(Xdom.grid, xpos)
        u = ABS.get_coord_by_pos(Udom.grid, upos)
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

function cell_approx!(ax, vars, Xdom, Udom, contsys::ABS.ControlSystemLinearized{N,T};
        fc = "yellow", fa = 0.5, ec = "gold", ea = 1.0, ew = 0.5) where {N,T}
    @assert length(vars) == 2 && N >= 2
    fca = FC(fc, fa)
    eca = FC(ec, ea)
    vertslist = Vector{SVector{2,T}}[]
    _H_ = SMatrix{N,N}(I).*(Xdom.grid.h/2.0)
    e = norm(Xdom.grid.h/2.0 + contsys.measnoise, Inf)

    for xpos in ABS.enum_pos(Xdom), upos in ABS.enum_pos(Udom)
        x = ABS.get_coord_by_pos(Xdom.grid, xpos)
        u = ABS.get_coord_by_pos(Udom.grid, upos)
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
    for i in 1:nstep
        xpos = ABS.get_pos_by_coord(symmodel.Xdom.grid, x0)
        if !(xpos ∈ symmodel.Xdom)
            @warn("Trajectory out of domain")
            return
        end
        source = ABS.get_state_by_xpos(symmodel, xpos)
        symbollist = Int[]
        ABS.compute_enabled_symbols!(symbollist, contr, source)
        if isempty(symbollist)
            @warn("Uncontrollable state")
            return
        end
        if randchoose
            symbol = rand(symbollist)
        else
            symbol = symbollist[1]
        end
        upos = ABS.get_upos_by_symbol(symmodel, symbol)
        u = ABS.get_coord_by_pos(symmodel.Udom.grid, upos)
        Plot.trajectory_open_loop!(ax, vars, contsys, x0, u, 1,
            lc = lc, lw = lw, mc = mc, ms = ms, nsub = nsub)
        x0 = contsys.sys_map(x0, u, contsys.tstep)
    end
end
=#

end  # Plot
