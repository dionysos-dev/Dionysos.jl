struct Digraph{T <: Real,U}
    edges::Dict{Tuple{U,U},T}
    verts::Set{U}
end
 
function Digraph(edges::Vector{Tuple{U,U,T}}) where {T <: Real,U}
    vnames = Set{U}(v for edge in edges for v in edge[1:2])
    adjmat = Dict((edge[1], edge[2]) => edge[3] for edge in edges)
    return Digraph(adjmat, vnames)
end
 
vertices(g::Digraph) = g.verts
edges(g::Digraph)    = g.edges
 
neighbours(g::Digraph, v) = Set((b, c) for ((a, b), c) in edges(g) if a == v)
 
function dijkstrapath(g::Digraph{T,U}, source::U, dest::U) where {T, U}
    @assert source ∈ vertices(g) "$source is not a vertex in the graph"
 
    # Easy case
    if source == dest return [source], 0 end
    # Initialize variables
    inf  = typemax(T)
    dist = Dict(v => inf for v in vertices(g))
    prev = Dict(v => v   for v in vertices(g))
    dist[source] = 0
    Q = copy(vertices(g))
    neigh = Dict(v => neighbours(g, v) for v in vertices(g))
 
    # Main loop
    while !isempty(Q)
        u = reduce((x, y) -> dist[x] < dist[y] ? x : y, Q)
        pop!(Q, u)
        #if dist[u] == inf || u == dest break end
        if dist[u] == inf  break end
        for (v, cost) in neigh[u]
            alt = dist[u] + cost
            if alt < dist[v]
                dist[v] = alt
                prev[v] = u
            end
        end
    end
 
    # Return path
    rst, cost = U[], dist
    if prev[dest] == dest
        return rst, cost
    else
        while dest != source
            pushfirst!(rst, dest)
            dest = prev[dest]
        end
        pushfirst!(rst, dest)
        return rst, cost
    end
end

function reverse_graph(gc)
    return gc
end
 
# testgraph = [("a", "b", 1), ("b", "e", 2), ("a", "e", 4)]
# testgraph = [("a", "b", 7.1),  ("a", "c", 9.5),  ("a", "f", 14.5), ("b", "c", 10.5),
#              ("b", "d", 15.5), ("c", "d", 20.5), ("c", "f", 2.5),  ("d", "e", 26.5),
#              ("e", "f", 9.5),("e", "g", 9.5)]
# g = Digraph(testgraph)

# src, dst = "a", "g"
# path, cost = dijkstrapath(g, src, dst)
# println("Shortest path from $src to $dst: ", isempty(path) ? "no possible path" : join(path, " → "), " (cost $cost)")
 
# Print all possible paths
# @printf("\n%4s | %3s | %s\n", "src", "dst", "path")
# @printf("----------------\n")
# for src in vertices(g), dst in vertices(g)
#     path, cost = dijkstrapath(g, src, dst)
#     @printf("%4s | %3s | %s\n", src, dst, isempty(path) ? "no possible path" : join(path, " → ") * " ($cost)")
# end
