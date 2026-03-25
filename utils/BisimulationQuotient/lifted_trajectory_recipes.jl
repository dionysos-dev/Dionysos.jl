using GLMakie
using GeometryBasics
using LazySets

function _vertices2d(P)
    VP = overapproximate(P, VPolygon)
    verts = vertices_list(VP)
    xs = [v[1] for v in verts]
    ys = [v[2] for v in verts]
    return xs, ys
end

_memory_to_qid(mem) = mem[2]

function make_node_zmap(bisimulation)
    nodes = sort!(collect(keys(bisimulation.part_ids)); by = x -> string(x))
    return Dict(node => Float64(i) for (i, node) in enumerate(nodes))
end

function _polygon_mesh3d(xs, ys, z)
    n = length(xs)
    pts = Point3f[(xs[i], ys[i], z) for i in 1:n]

    # triangle fan: (1,2,3), (1,3,4), ...
    faces = GLTriangleFace[]
    for i in 2:(n - 1)
        push!(faces, GLTriangleFace(1, i, i+1))
    end

    return GeometryBasics.Mesh(pts, faces)
end

function plot_lifted_bisimulation_makie!(
    ax,
    bisimulation;
    state_ids = nothing,
    node_z = make_node_zmap(bisimulation),
    color_by = :node,
    node_colors = Dict((1,) => :blue, (2,) => :orange),
    strokecolor = :black,
    strokewidth = 1.0,
    alpha = 0.6,
    show_contours = true,
)
    palette = [:red, :blue, :green, :orange, :purple, :brown, :pink, :cyan]

    ids =
        isnothing(state_ids) ? sort(collect(keys(bisimulation.states))) : collect(state_ids)

    for sid in ids
        q = bisimulation.states[sid]
        z = node_z[q.node]

        c = if color_by == :node
            get(node_colors, q.node, :gray)
        elseif color_by == :slice
            palette[mod1(q.slice, length(palette))]
        elseif color_by == :obs
            palette[mod1(q.obs + 2, length(palette))]
        elseif color_by == :state
            palette[mod1(sid, length(palette))]
        else
            color_by
        end

        for part in q.set.parts
            xs, ys = _vertices2d(part)
            msh = _polygon_mesh3d(xs, ys, z)
            Makie.mesh!(ax, msh; color = (c, alpha))
            if show_contours
                pts = Point3f[(xs[i], ys[i], z) for i in eachindex(xs)]
                Makie.lines!(
                    ax,
                    vcat(pts, pts[1:1]);
                    color = strokecolor,
                    linewidth = strokewidth,
                )
            end
        end
    end
    return ax
end

function plot_lifted_trajectory_makie!(
    ax,
    bisimulation,
    X_seq,
    M_seq;
    node_z = make_node_zmap(bisimulation),
    traj_color = :black,
    traj_linewidth = 3,
    markersize = 10,
    show_points = true,
    show_start_end = true,
)
    xs = [x[1] for x in X_seq]
    ys = [x[2] for x in X_seq]

    qids = [_memory_to_qid(m) for m in M_seq]
    zs = [node_z[bisimulation.states[qid].node] for qid in qids]

    if length(xs) == length(zs) + 1
        push!(zs, zs[end])
    elseif length(zs) == length(xs) + 1
        pop!(zs)
    elseif length(zs) != length(xs)
        error("Incompatible lengths")
    end

    pts = Point3f[(xs[i], ys[i], zs[i]) for i in eachindex(xs)]

    # trajectory line
    Makie.lines!(ax, pts; color = traj_color, linewidth = traj_linewidth)

    # intermediate points
    if show_points
        Makie.scatter!(ax, pts; color = traj_color, markersize = markersize)
    end

    # start / end markers
    if show_start_end
        Makie.scatter!(ax, [pts[1]]; color = :green, markersize = markersize * 1.5)
        Makie.scatter!(ax, [pts[end]]; color = :red, markersize = markersize * 1.5)
    end

    return ax
end
