using LightGraphs

mutable struct LightTopologicalGraph{GT, ET, T}
    G::GT
    Σ::Dict{ET, T}
end

function LightTopologicalGraph{T}(n::Int) where {T}
    G = LightGraphs.SimpleGraph(n)
    Σ = Dict{LightGraphs.edgetype(G), T}()
    LightTopologicalGraph(G, Σ)
end

function add_transition!(A::LightTopologicalGraph, q, r, I)
    edge = edge_object(A, q, r)
    LightGraphs.add_edge!(A.G, edge)
    A.Σ[edge] = I
end
function has_transition(A::LightTopologicalGraph, q, r)
    return LightGraphs.has_edge(A.G, edge_object(A, q, r))
end

function rem_transition!(A::LightTopologicalGraph, q, r)
    edge = edge_object(A, q, r)
    LightGraphs.rem_edge!(A.G, edge)
    delete!(A.Σ, edge)
end

function edge_object(A::LightTopologicalGraph, q, r)
    @assert 1 <= q <= LightGraphs.nv(A.G)
    @assert 1 <= r <= LightGraphs.nv(A.G)
    return LightGraphs.Edge(q, r)
end

function check_obstacles(I::AbstractPolytope,O::Vector{<:AbstractPolytope}) where T
    I_e = expansion(0.8,I)
    for obs in O
        if LazySets.issubset(I_e, obs)
            return false
        end
    end
    return true
end

function topological_graph(D::Partition{T},O) where T
    l = length(D.L)
    A = LightTopologicalGraph{T}(l)
    for i=1:l
        for j=i+1:l
            if(!is_intersection_empty(D.L[i],D.L[j]))
                I =  intersection(D.L[i],D.L[j])
                if check_obstacles(I,O)
                    add_transition!(A, i, j, I)
                end
            end
        end
    end
    return A
end

function get_transition(A::LightTopologicalGraph,i,j)
    e = edge_object(A, min(i,j), max(i,j))
    return A.Σ[e]
end

function get_path(A::LightTopologicalGraph,s::Int64,t::Int64)
    #ds = dijkstra_shortest_paths(A.G, s)
    ds = a_star(A.G, s, t)
    path = [src(e) for e in ds]
    path = [path..., t]
    return path
end
function plot_topological_graph(D::Partition,A::LightTopologicalGraph,O = [];path = [],save=nothing)
    fig = plot(aspect_ratio = 1,legend = false)
    l = length(D.L)
    plot_partition!(D)
    plot!(O,color=:black)
    for i=1:l
        plot!(Singleton(LazySets.center(D.L[i])))
        for j=i+1:l
            if(i!=j && has_transition(A, i, j))
                plot!(LineSegment(LazySets.center(D.L[i]), LazySets.center(D.L[j])),color=:blue)
            end
        end
    end
    for i=1:length(path)-1
        plot!(LineSegment(LazySets.center(D.L[path[i]]), LazySets.center(D.L[path[i+1]])),color=:red)
    end
    display(fig)
    if save != nothing
        savefig(fig, save)
    end
end

function plot_path(D::Partition,path::Vector{Int})
    fig = plot(aspect_ratio = 1,legend = false)
    for i=1:length(path)-1
        plot!(D.L[path[i]])
        plot!(LineSegment(LazySets.center(D.L[path[i]]), LazySets.center(D.L[path[i+1]])),color=:red)
    end
    plot!(D.L[path[end]])
    display(fig)
end
