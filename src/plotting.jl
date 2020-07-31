# All plotting functions
module Plot

import ..Abstraction
AB = Abstraction

using PyPlot
using PyCall

const spatial = pyimport_conda("scipy.spatial", "scipy")
FC(c, a) =  matplotlib.colors.colorConverter.to_rgba(c, alpha = a)

function verts_rect(c, h)
    return (c .+ (-h[1], -h[2]), c .+ (-h[1], h[2]), c .+ (h[1], h[2]), c .+ (h[1], -h[2]))
end

## =============================================================================
# Cells
function subset!(ax, vars, sub_set::AB.SubSet{N};
        fc = "red", fa = 0.5, ec = "black", ea = 1.0, ew = 1.5) where N
    #---------------------------------------------------------------------------
    grid_space = sub_set.grid_space
    @assert length(vars) == 2 && N >= 2
    fca = FC(fc, fa)
    eca = FC(ec, ea)

    pos_iter = (AB.get_pos_by_ref(grid_space, ref) for ref in AB.enumerate_subset_ref(sub_set))
    verts_vec = NTuple{4, Tuple{Float64, Float64}}[]

    for pos in unique(x -> x[vars], pos_iter)
        c = AB.get_coords_by_pos(grid_space, pos)[vars]
        push!(verts_vec, verts_rect(c, grid_space.h[vars]./2))
    end

    poly_list = matplotlib.collections.PolyCollection(verts_vec)
    poly_list.set_facecolor(fca)
    poly_list.set_edgecolor(eca)
    poly_list.set_linewidth(ew)
    ax.add_collection(poly_list)
end

## =============================================================================
# Sets
function set!(ax, vars, rect::AB.HyperRectangle;
        fc = "green", fa = 0.5, ec = "black", ea = 1.0, ew = 1.5)
    #---------------------------------------------------------------------------
    @assert length(vars) == 2 && length(rect.lb) == length(rect.ub) >= 2
    c = (rect.lb .+ rect.ub)./2
    h = (rect.ub .- rect.lb)./2
    poly_list = matplotlib.collections.PolyCollection([verts_rect(c[vars], h[vars])])
    fca = FC(fc, fa)
    eca = FC(ec, ea)
    poly_list.set_facecolor(fca)
    poly_list.set_edgecolor(eca)
    poly_list.set_linewidth(ew)
    ax.add_collection(poly_list)
end

## =============================================================================
# Trajectory open loop
function trajectory_open_loop!(ax, vars, cont_sys, x0::NTuple{N, Float64}, u, nstep;
        lc = "red", lw = 1.5, mc = "black", ms = 5.0, nsub = 5) where N
    #---------------------------------------------------------------------------
    @assert length(vars) == 2 && N >= 2
    tstep = cont_sys.tstep/nsub
    Nstep = nstep*nsub + 1
    X1_vec = Vector{Float64}(undef, Nstep)
    X2_vec = Vector{Float64}(undef, Nstep)
    X1_vec[1] = x0[vars[1]]
    X2_vec[1] = x0[vars[2]]
    x = x0

    for i = 2:Nstep
        x = cont_sys.sys_map(x, u, tstep)
        X1_vec[i] = x[vars[1]]
        X2_vec[i] = x[vars[2]]
    end

    idx = 1:nsub:Nstep
    ax.plot(X1_vec, X2_vec, lw = lw, c = lc)
    ax.plot(X1_vec[idx], X2_vec[idx], ls = "None", marker = ".", ms = ms, mfc = mc, mec = mc)
end

## =============================================================================
# Images
function cell_image!(ax, vars, X_sub, U_sub, cont_sys::AB.ControlSystem{N};
        nsub = fill(5, N), fc = "blue", fa = 0.5, ec = "darkblue", ea = 1.0, ew = 1.5) where N
    #---------------------------------------------------------------------------
    @assert length(vars) == 2 && N >= 2
    fca = FC(fc, fa)
    eca = FC(ec, ea)
    verts_list = Vector{Tuple{Float64, Float64}}[]
    ns = nsub .- 1

    for x_ref in AB.enumerate_subset_ref(X_sub), u_ref in AB.enumerate_subset_ref(U_sub)
        x = AB.get_coords_by_ref(X_sub.grid_space, x_ref)
        u = AB.get_coords_by_ref(U_sub.grid_space, u_ref)
        subpos_axes = ((0:ns[i])./ns[i] .- 0.5 for i = 1:length(ns))
        subpos_iter = Iterators.product(subpos_axes...)
        x_iter = (x .+ subpos.*X_sub.grid_space.h for subpos in subpos_iter)
        Fx_iter = (cont_sys.sys_map(x, u, cont_sys.tstep) for x in x_iter)
        Fxproj_vec = [Fx[vars] for Fx in Fx_iter][:]
        hull = spatial.ConvexHull([Fx[i] for Fx in Fxproj_vec, i in 1:2])
        push!(verts_list, Fxproj_vec[hull.vertices .+ 1])
    end

    poly_list = matplotlib.collections.PolyCollection(verts_list)
    poly_list.set_facecolor(fca)
    poly_list.set_edgecolor(eca)
    poly_list.set_linewidth(ew)
    ax.add_collection(poly_list)
end

## =============================================================================
# Outer-approximation
function cell_approx!(ax, vars, X_sub, U_sub, cont_sys::AB.ControlSystem{N};
        fc = "yellow", fa = 0.5, ec = "gold", ea = 1.0, ew = 0.5) where N
    #---------------------------------------------------------------------------
    @assert length(vars) == 2 && N >= 2
    fca = FC(fc, fa)
    eca = FC(ec, ea)
    verts_list = NTuple{4, Tuple{Float64, Float64}}[]

    for x_ref in AB.enumerate_subset_ref(X_sub), u_ref in AB.enumerate_subset_ref(U_sub)
        x = AB.get_coords_by_ref(X_sub.grid_space, x_ref)
        u = AB.get_coords_by_ref(U_sub.grid_space, u_ref)
        Fx = cont_sys.sys_map(x, u, cont_sys.tstep)
        r = X_sub.grid_space.h./2 .+ cont_sys.meas_noise
        Fr = cont_sys.bound_map(r, u, cont_sys.tstep)
        Fr = Fr .+ cont_sys.meas_noise
        push!(verts_list, verts_rect(Fx[vars], Fr[vars]))
    end

    poly_list = matplotlib.collections.PolyCollection(verts_list)
    poly_list.set_facecolor(fca)
    poly_list.set_edgecolor(eca)
    poly_list.set_linewidth(ew)
    ax.add_collection(poly_list)
end

## =============================================================================
# Trajectory closed loop
function trajectory_closed_loop!(ax, vars, cont_sys, sym_model, x0, nstep;
        lc = "red", lw = 1.5, mc = "black", ms = 5.0, nsub = 5)
    #---------------------------------------------------------------------------
    X_grid = sym_model.X_grid
    U_sub = AB.NewSubSet(sym_model.U_grid)
    Y_sub = AB.NewSubSet(sym_model.Y_grid)
    for i = 1:nstep
        AB.remove_from_subset_all!(U_sub)
        AB.remove_from_subset_all!(Y_sub)
        x_ref = AB.get_ref_by_coords(X_grid, x0)
        if x_ref === X_grid.overflow_ref
            @warn("Trajectory out of domain")
            return
        end
        if ~AB.is_xref_controllable(sym_model, x_ref)
            @warn("Uncontrollable state")
            return
        end
        uref_coll = AB.get_gridspace_reftype(sym_model.U_grid)[]
        yref_coll = AB.get_gridspace_reftype(sym_model.Y_grid)[]
        AB.add_inputs_images_by_xref!(uref_coll, yref_coll, sym_model, x_ref)
        u_ref = uref_coll[1]
        u = AB.get_coords_by_ref(sym_model.U_grid, u_ref)
        Plot.trajectory_open_loop!(ax, vars, cont_sys, x0, u, 1,
            lc = lc, lw = lw, mc = mc, ms = ms, nsub = nsub)
        x0 = cont_sys.sys_map(x0, u, cont_sys.tstep)
    end
end
end  # Plot
