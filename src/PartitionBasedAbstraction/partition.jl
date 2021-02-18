module PartitionAbstraction
using LinearAlgebra,StaticArrays,Plots,Polyhedra,LazySets
import LazySets.AbstractPolytope

include(joinpath("..", "Abstraction", "abstraction.jl"))
using .Abstraction
AB = Abstraction
include("constrainedZonotope.jl")

mutable struct Partition{T}
    X::T
    L::Vector{T}
end

function get_vertices(X::T,L::Vector{<:T}) where T
    v = []
    N = length(L)
    for i = 1:N
        for j = i+1:N
            if !is_intersection_empty(L[i],L[j])
                I = intersection(L[i],L[j])
                push!(v,vertices_list(I)...)
            end
        end
    end

    for i = 1:N
        if !is_intersection_empty(L[i],X)
            I = intersection(L[i],X)
            push!(v,vertices_list(I)...)
        end
    end
    push!(v,vertices_list(X)...)
    sort!(v)
    v2 = unique(map(x->round.(x,digits = 8),v))

    for P in L
        filter!(e->!is_interior_point(e,P, ε=1e-3),v2)
    end
    return v2
end

function expansion(α::Number,D::Partition{T}) where T
    L = [expansion(α,P) for P in D.L]
    return Partition{T}(D.X,L)
end

function plot!(L::Vector{T}) where T
    for P in L
        plot!(P)
    end
end
function plot_partition!(D::Partition)
    plot!(D.X)
    plot!(D.L)
end

function plot_partition(D::Partition;save=nothing)
    fig = plot(aspect_ratio = 1,legend = false)
    plot!(D.X)
    plot!(D.L)
    display(fig)
    if save != nothing
        savefig(fig, save)
    end
end
include("topologicalGraph.jl")
include("controlsystem.jl")
include("symbolicmodel.jl")
include("controller.jl")



end
