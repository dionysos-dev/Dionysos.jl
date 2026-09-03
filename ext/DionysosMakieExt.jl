module DionysosMakieExt

# 3D "augmented" visualisation of a PCLF bisimulation quotient: quotient states are drawn as
# stacked polygons, one z-layer per automaton node, with the closed-loop trajectory augmented
# onto the layer of the node it currently sits in. Backend-agnostic — the trigger is `Makie`,
# so any backend (`GLMakie` for interaction, `CairoMakie` for figures, …) activates it.

import Dionysos
import LazySets
import Makie

const DI = Dionysos
const PCLFBQ = DI.Optim.Abstraction.PCLFBisimulationQuotient
const GB = Makie.GeometryBasics

# Assign each automaton node a distinct z-layer, ordered deterministically by name.
function _make_node_zmap(bisimulation::PCLFBQ.PCBisimulationQuotient)
    nodes = sort!(collect(keys(bisimulation.part_ids)); by = x -> string(x))
    return Dict(node => Float64(i) for (i, node) in enumerate(nodes))
end

# Planar outline of a (possibly lazy) set as x/y vertex vectors.
function _vertices2d(P)
    # A degenerate (vertexless) part has nothing to draw.
    isempty(LazySets.vertices_list(P)) && return Float64[], Float64[]
    VP = LazySets.overapproximate(P, LazySets.VPolygon)
    verts = LazySets.vertices_list(VP)
    xs = [v[1] for v in verts]
    ys = [v[2] for v in verts]
    return xs, ys
end

# A filled polygon at height `z`, triangulated as a fan (assumes a convex outline).
function _polygon_mesh3d(xs, ys, z)
    n = length(xs)
    pts = GB.Point3f[(xs[i], ys[i], z) for i in 1:n]
    faces = GB.GLTriangleFace[]
    for i in 2:(n - 1)
        push!(faces, GB.GLTriangleFace(1, i, i + 1))
    end
    return GB.Mesh(pts, faces)
end

# Append one polygon's fan triangulation to a mesh being accumulated, offsetting the face
# indices past the vertices already there.
#
# A quotient has thousands of cells, and both Makie and Plots charge per plot object rather
# than per triangle, so drawing each cell separately is what makes these figures slow. One
# mesh per colour and one NaN-separated line for every outline draw the same picture out of a
# handful of objects instead of tens of thousands.
function _append_polygon!(pts, faces, xs, ys, z)
    n = length(xs)
    n < 3 && return nothing
    off = length(pts)
    for i in 1:n
        push!(pts, GB.Point3f(xs[i], ys[i], z))
    end
    for i in 2:(n - 1)
        push!(faces, GB.GLTriangleFace(off + 1, off + i, off + i + 1))
    end
    return nothing
end

# Append a closed outline, then a NaN point so the next one starts a new stroke.
function _append_outline!(outline, xs, ys, z)
    isempty(xs) && return nothing
    for i in eachindex(xs)
        push!(outline, GB.Point3f(xs[i], ys[i], z))
    end
    push!(outline, GB.Point3f(xs[1], ys[1], z))
    push!(outline, GB.Point3f(NaN, NaN, NaN))
    return nothing
end

# The controller memory is a `(mode, quotient_state_id)` tuple.
_memory_to_qid(mem) = mem[2]

function Dionysos.plot_augmented_bisimulation!(
    ax,
    bisimulation::PCLFBQ.PCBisimulationQuotient;
    state_ids = nothing,
    node_z = _make_node_zmap(bisimulation),
    color_by = :node,
    node_colors = Dict((1,) => :blue, (2,) => :orange),
    strokecolor = :black,
    strokewidth = 1.0,
    alpha = 0.6,
    show_contours = true,
    merge_plots = true,
)
    palette = [:red, :blue, :green, :orange, :purple, :brown, :pink, :cyan]

    ids = if isnothing(state_ids)
        sort(collect(keys(bisimulation.states)))
    else
        collect(state_ids)
    end

    colour_of(sid, q) =
        if color_by == :node
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

    if !merge_plots
        for sid in ids
            q = bisimulation.states[sid]
            z = node_z[q.node]
            c = colour_of(sid, q)

            for part in q.set.array
                xs, ys = _vertices2d(part)
                length(xs) < 3 && continue
                Makie.mesh!(ax, _polygon_mesh3d(xs, ys, z); color = (c, alpha))
                if show_contours
                    pts = GB.Point3f[(xs[i], ys[i], z) for i in eachindex(xs)]
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

    # One accumulating mesh per colour, and one outline stroke for all of them: the outlines
    # share `strokecolor`, so they never need splitting.
    meshes = Dict{Any, Tuple{Vector{GB.Point3f}, Vector{GB.GLTriangleFace}}}()
    order = Any[]
    outline = GB.Point3f[]

    for sid in ids
        q = bisimulation.states[sid]
        z = node_z[q.node]
        c = colour_of(sid, q)

        pts, faces = get!(meshes, c) do
            push!(order, c)
            return (GB.Point3f[], GB.GLTriangleFace[])
        end

        for part in q.set.array
            xs, ys = _vertices2d(part)
            _append_polygon!(pts, faces, xs, ys, z)
            show_contours && _append_outline!(outline, xs, ys, z)
        end
    end

    for c in order
        pts, faces = meshes[c]
        isempty(faces) && continue
        Makie.mesh!(ax, GB.Mesh(pts, faces); color = (c, alpha))
    end

    if show_contours && !isempty(outline)
        Makie.lines!(ax, outline; color = strokecolor, linewidth = strokewidth)
    end

    return ax
end

function Dionysos.plot_augmented_trajectory!(
    ax,
    bisimulation::PCLFBQ.PCBisimulationQuotient,
    state_seq,
    memory_seq;
    node_z = _make_node_zmap(bisimulation),
    traj_color = :black,
    traj_linewidth = 3,
    markersize = 10,
    show_points = true,
    show_start_end = true,
)
    xs = [x[1] for x in state_seq]
    ys = [x[2] for x in state_seq]

    qids = [_memory_to_qid(m) for m in memory_seq]
    zs = [node_z[bisimulation.states[qid].node] for qid in qids]

    # The state and memory sequences can differ in length by one depending on whether the
    # final memory update is recorded; reconcile before zipping into 3D points.
    if length(xs) == length(zs) + 1
        push!(zs, zs[end])
    elseif length(zs) == length(xs) + 1
        pop!(zs)
    elseif length(zs) != length(xs)
        error("Incompatible state/memory sequence lengths")
    end

    pts = GB.Point3f[(xs[i], ys[i], zs[i]) for i in eachindex(xs)]

    Makie.lines!(ax, pts; color = traj_color, linewidth = traj_linewidth)

    if show_points
        Makie.scatter!(ax, pts; color = traj_color, markersize = markersize)
    end

    if show_start_end
        Makie.scatter!(ax, [pts[1]]; color = :green, markersize = markersize * 1.5)
        Makie.scatter!(ax, [pts[end]]; color = :red, markersize = markersize * 1.5)
    end

    return ax
end

end # module DionysosMakieExt
