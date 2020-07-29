# All plotting functions

const spatial = pyimport_conda("scipy.spatial", "scipy")
FC(c, a) =  matplotlib.colors.colorConverter.to_rgba(c, alpha = a)

function verts_rect(c, h)
    return [c + [i*h[1], j*h[2]] for (i, j) in [(-1, -1), (-1, 1), (1, 1), (1, -1)]]
end

## =============================================================================
# Cells
function plot_subspace!(ax, vars, sub_space;
        fc = "red", fa = 0.5, ec = "black", ea = 1.0, ew = 1.5)
    #---------------------------------------------------------------------------
    grid_space = sub_space.grid_space
    @assert length(vars) == 2 && grid_space.dim >= 2
    fca = FC(fc, fa)
    eca = FC(ec, ea)

    pos_coll = (get_pos_by_ref(grid_space, ref) for ref in enumerate_subspace_ref(sub_space))
    verts_vec = Vector{Vector{Float64}}[]

    for pos in unique(x -> x[vars], pos_coll)
        c = get_coords_by_pos(grid_space, pos)[vars]
        push!(verts_vec, verts_rect(c, grid_space.h[vars]/2))
    end

    poly_list = matplotlib.collections.PolyCollection(verts_vec)
    poly_list.set_facecolor(fca)
    poly_list.set_edgecolor(eca)
    poly_list.set_linewidth(ew)
    ax.add_collection(poly_list)
end

## =============================================================================
# Sets
function plot_box!(ax, vars, lb, ub;
        fc = "green", fa = 0.5, ec = "black", ea = 1.0, ew = 1.5)
    #---------------------------------------------------------------------------
    @assert length(vars) == 2 && length(lb) == length(ub) >= 2
    c = (lb + ub)/2
    h = (ub - lb)/2
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
function plot_trajectory_open_loop!(ax, vars, cont_sys, x0, u, nstep;
        lc = "red", lw = 1.5, mc = "black", ms = 5.0, nsub = 5)
    #---------------------------------------------------------------------------
    @assert length(vars) == 2 && length(x0) >= 2
    tstep = cont_sys.tstep/nsub
    Nstep = nstep*nsub + 1
    X1_vec = Vector{Float64}(undef, Nstep)
    X2_vec = Vector{Float64}(undef, Nstep)
    X1_vec[1] = x0[vars[1]]
    X2_vec[1] = x0[vars[2]]
    x0 = copy(x0)

    for i = 2:Nstep
        cont_sys.sys_map!(x0, u, tstep)
        X1_vec[i] = x0[vars[1]]
        X2_vec[i] = x0[vars[2]]
    end

    idx = 1:nsub:Nstep
    ax.plot(X1_vec, X2_vec, lw = lw, c = lc)
    ax.plot(X1_vec[idx], X2_vec[idx], ls = "None", marker = ".", ms = ms, mfc = mc, mec = mc)
end

## =============================================================================
# Images
function plot_cell_image!(ax, vars, X_sub, U_sub, cont_sys;
        nsub = fill(5, X_sub.grid_space.dim), fc = "blue", fa = 0.5, ec = "darkblue", ea = 1.0, ew = 1.5)
    #---------------------------------------------------------------------------
    @assert length(vars) == 2 && X_sub.grid_space.dim >= 2
    fca = FC(fc, fa)
    eca = FC(ec, ea)
    verts_list = Vector{Vector{Float64}}[]
    ns = nsub .- 1

    for x_ref in enumerate_subspace_ref(X_sub), u_ref in enumerate_subspace_ref(U_sub)
        x = get_coords_by_ref(X_sub.grid_space, x_ref)
        u = get_coords_by_ref(U_sub.grid_space, u_ref)
        pos_ = ((0:ns[i])./ns[i] .- 0.5 for i = 1:length(ns))
        x_coll = (x + pos.*X_sub.grid_space.h for pos in Iterators.product(pos_...))
        sys_map = x -> (cont_sys.sys_map!(x, u, cont_sys.tstep); return x)
        Fx_coll = (sys_map(x) for x in x_coll)
        Fxproj_vec = [Fx[vars] for Fx in Fx_coll][:]
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
function plot_cell_approx!(ax, vars, X_sub, U_sub, cont_sys;
        fc = "yellow", fa = 0.5, ec = "gold", ea = 1.0, ew = 0.5)
    #---------------------------------------------------------------------------
    @assert length(vars) == 2 && X_sub.grid_space.dim >= 2
    fca = FC(fc, fa)
    eca = FC(ec, ea)
    verts_list = Vector{Vector{Float64}}[]

    for x_ref in enumerate_subspace_ref(X_sub), u_ref in enumerate_subspace_ref(U_sub)
        x = get_coords_by_ref(X_sub.grid_space, x_ref)
        u = get_coords_by_ref(U_sub.grid_space, u_ref)
        cont_sys.sys_map!(x, u, cont_sys.tstep)
        r = X_sub.grid_space.h/2 + cont_sys.meas_noise
        cont_sys.bound_map!(r, u, cont_sys.tstep)
        r += cont_sys.meas_noise
        push!(verts_list, verts_rect(x[vars], r[vars]))
    end

    poly_list = matplotlib.collections.PolyCollection(verts_list)
    poly_list.set_facecolor(fca)
    poly_list.set_edgecolor(eca)
    poly_list.set_linewidth(ew)
    ax.add_collection(poly_list)
end

## =============================================================================
# Trajectory closed loop
function plot_trajectory_closed_loop!(ax, vars, cont_sys, trans_map, x0, nstep;
        lc = "red", lw = 1.5, mc = "black", ms = 5.0, nsub = 5)
    #---------------------------------------------------------------------------
    X_grid = trans_map.X_grid
    x0 = copy(x0)
    for i = 1:nstep
        x_ref = get_ref_by_coords(X_grid, x0)
        if x_ref === X_grid.overflow_ref
            @warn("Trajectory out of domain")
            return
        end
        uy_ref_coll = get_transition_image(trans_map, x_ref)
        if isempty(uy_ref_coll)
            @warn("Uncontrollable state")
            return
        end
        u_ref = iterate(uy_ref_coll)[1][1]
        u = get_coords_by_ref(trans_map.U_grid, u_ref)
        plot_trajectory_open_loop!(
            ax, vars, cont_sys, x0, u, 1, lc = lc, lw = lw, mc = mc, ms = ms, nsub = nsub)
        cont_sys.sys_map!(x0, u, cont_sys.tstep)
    end
end

#=
## =============================================================================
# Sets
function plotset!(ax, vars, grid_sub, set::AbstractPolyhedron;
        rad = Inf, fc = "green", fa = 0.5, ec = "black", ea = 1.0, ew = 1.5)
    #---------------------------------------------------------------------------
    # Plots the Polyhedron set. Use rad if set is not bounded.
    @assert length(vars) == 2 && grid_space.dim >= 2
    set = Projection(set, vars)

    if rad < Inf
        H_list = constraints_list(set ∩ BallInf(zeros(2), rad))
        set = HPolyhedron(Vector{LazySets.HalfSpace{Float64,Array{Float64,1}}}(H_list))
    end

    poly_list = matplotlib.collections.PolyCollection([vertices_list(set)])
    fca = FC(fc, fa)
    eca = FC(ec, ea)
    poly_list.set_facecolor(fca)
    poly_list.set_edgecolor(eca)
    poly_list.set_linewidth(ew)
    ax.add_collection(poly_list)
end

## =============================================================================
# Images
function plotimage!(ax, vars, sym_model, iX_list, iU_list;
        nsub = fill(5, sym_model.X_space.dim), fc = "blue", fa = 0.5, ec = "darkblue", ea = 1.0, ew = 1.5)
    #---------------------------------------------------------------------------
    # Plots (convex hull) of images of cells in iX_list with inputs in iU_list
    @assert length(vars) == 2 && sym_model.X_space.dim >= 2
    if isempty(iX_list) || isempty(iU_list)
        return
    end
    iX_list_ = _preprocess_idxlist(sym_model.X_space, iX_list, soft = true)
    iU_list_ = _preprocess_idxlist(sym_model.U_space, iU_list, soft = true)
    fca = FC(fc, fa)
    eca = FC(ec, ea)

    X_space = sym_model.X_space
    U_space = sym_model.U_space
    cont_sys = sym_model.cont_sys
    idx_iter = Iterators.product(iX_list_, iU_list_)
    vertlist_list = Vector{Vector{Float64}}[]

    for (iX, iU) in idx_iter
        u = U_space.orig + U_space.cell_list[iU].*U_space.h
        nsub = nsub .- 1
        subs = [((0:nsub[i])./nsub[i] .- 0.5) for i = 1:length(nsub)]
        sub_iter = Iterators.product(subs...)
        xsub = sub -> X_space.orig + (X_space.cell_list[iX] .+ sub).*X_space.h
        sys_map = x -> cont_sys.sysbound_map(x, zeros(X_space.dim), u)[1][vars]
        F_list = [sys_map(xsub(sub)) for sub in sub_iter][:]
        convex_hull!(F_list)
        push!(vertlist_list, F_list)
    end

    poly_list = matplotlib.collections.PolyCollection(vertlist_list)
    poly_list.set_facecolor(fca)
    poly_list.set_edgecolor(eca)
    poly_list.set_linewidth(ew)
    ax.add_collection(poly_list)
end

## =============================================================================
# Overapproximation
function plotapprox!(ax, vars, sym_model, iX_list, iU_list;
        fc = "yellow", fa = 0.5, ec = "gold", ea = 1.0, ew = 0.5)
    #---------------------------------------------------------------------------
    # Plots ouverapproximation of the images of cells in iX_list with inputs in iU_list.
    # Use the "growth bound" function.
    @assert length(vars) == 2
    if isempty(iX_list) || isempty(iU_list)
        return
    end
    iX_list_ = _preprocess_idxlist(sym_model.X_space, iX_list, soft = true)
    iU_list_ = _preprocess_idxlist(sym_model.U_space, iU_list, soft = true)
    fca = FC(fc, fa)
    eca = FC(ec, ea)

    X_space = sym_model.X_space
    U_space = sym_model.U_space
    cont_sys = sym_model.cont_sys
    idx_iter = Iterators.product(iX_list_, iU_list_)
    vertlist_list = Vector{Vector{Float64}}[]

    for (iX, iU) in idx_iter
        cx = X_space.orig + X_space.cell_list[iX].*X_space.h
        u = U_space.orig + U_space.cell_list[iU].*U_space.h
        Fx, rad = cont_sys.sysbound_map(cx, X_space.h/2 + cont_sys.meas_noise, u)
        rad += cont_sys.meas_noise
        push!(vertlist_list, _rectangle_vertlist(Fx[vars], 2*rad[vars], [2, 2]))
    end

    poly_list = matplotlib.collections.PolyCollection(vertlist_list)
    poly_list.set_facecolor(fca)
    poly_list.set_edgecolor(eca)
    poly_list.set_linewidth(ew)
    ax.add_collection(poly_list)
end
=#
