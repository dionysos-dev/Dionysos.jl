# ----------------------------
# ListMapping
# ----------------------------

"""
    ListMapping{N,T} <: AbstractMapping{N,T}

Finite list of concrete points (SVector{N,T}) with id=1..K.
"""
mutable struct ListMapping{N, T} <: AbstractMapping{N, T}
    id2coord::Vector{SVector{N, T}}
    coord2id::Dict{SVector{N, T}, Int}
end

get_n_state(m::ListMapping) = length(m.id2coord)
enum_states(m::ListMapping) = 1:get_n_state(m)
get_coord_by_state(m::ListMapping, q::Int) = m.id2coord[q]

function get_state_by_coord(m::ListMapping{N, T}, x) where {N, T}
    xx = x isa SVector{N, T} ? x : SVector{N, T}(x)
    return m.coord2id[xx]
end

function Base.empty!(m::ListMapping)
    empty!(m.id2coord)
    empty!(m.coord2id)
    return m
end

convert_to_list_mapping(m::ListMapping) = m

# Constructors
function ListMapping(points::AbstractVector{SVector{N, T}}) where {N, T}
    id2 = SVector{N, T}[]
    coord2 = Dict{SVector{N, T}, Int}()
    for p in points
        if !haskey(coord2, p)         # deduplicate (keep first)
            push!(id2, p)
            coord2[p] = length(id2)
        end
    end
    return ListMapping{N, T}(id2, coord2)
end

function ListMapping(points::AbstractVector{<:Tuple})
    return ListMapping([SVector(p) for p in points])
end

function ListMapping(points::AbstractVector{<:AbstractVector{T}}) where {T}
    @assert !isempty(points) "points must be non-empty"
    N = length(points[1])
    return ListMapping([SVector{N, T}(p) for p in points])
end
