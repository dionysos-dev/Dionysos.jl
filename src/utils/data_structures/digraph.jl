mutable struct Digraph{T <: Real, U}
    edges::Dict{Tuple{U, U}, T}
    verts::Set{U}
end

# constructor based on transitions
function Digraph(edges::Vector{Tuple{U, U, T}}) where {T <: Real, U}
    vnames = Set{U}(v for edge in edges for v in edge[1:2])
    adjmat = Dict((edge[1], edge[2]) => edge[3] for edge in edges)
    return Digraph(adjmat, vnames)
end

function add_states!(g::Digraph{T, U}, states::Vector{U}) where {T <: Real, U}
    vnames = Set{U}(v for v in states)
    return g.verts = g.verts ∪ vnames
end

function add_transitions!(
    g::Digraph{T, U},
    edges::Vector{Tuple{U, U, T}},
) where {T <: Real, U}
    vnames = Set{U}(v for edge in edges for v in edge[1:2])
    adjmat = Dict((edge[1], edge[2]) => edge[3] for edge in edges)
    g.verts = g.verts ∪ vnames
    return g.edges = Dict(merge(g.edges, adjmat))
end

verts(g::Digraph) = g.verts
edges(g::Digraph) = g.edges

neighbours(g::Digraph, v) = Set((b, c) for ((a, b), c) in edges(g) if a == v)

function is_state(g::Digraph, s)
    return s ∈ verts(g)
end

# return 
# -the path: rst = [source p2 ... dest]
# -the cost from source to every element in the path: cost = Dict{source=>c0, p2=>c1, p3=>c2,...}
#  where c0 = 0.0, c1 is the cost from source to p2,  c2 is the cost from source to p3,...
function dijkstrapath(g::Digraph{T, U}, source::U, dest::U) where {T, U}
    @assert source ∈ verts(g) "$source is not a vertex in the graph"

    # Easy case
    if source == dest
        return [source], 0
    end
    # Initialize variables
    inf = typemax(T)
    dist = Dict(v => inf for v in verts(g))
    prev = Dict(v => v for v in verts(g))
    dist[source] = 0
    Q = copy(verts(g))
    neigh = Dict(v => neighbours(g, v) for v in verts(g))

    # Main loop
    while !isempty(Q)
        u = reduce((x, y) -> dist[x] < dist[y] ? x : y, Q)
        pop!(Q, u)
        #if dist[u] == inf || u == dest break end
        if dist[u] == inf
            break
        end
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
