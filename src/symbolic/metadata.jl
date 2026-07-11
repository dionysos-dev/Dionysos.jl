# A transition is stored and enumerated as the tuple `(target, source, symbol)`.
# Note the ordering flip vs. `add_transition!(autom, source, target, symbol)`: the
# storage order puts the target first. Use the accessors below rather than raw
# positional indexing (`t[1]`/`t[2]`/`t[3]`) at call sites.
const TransitionKey = NTuple{3, Int}

@inline transition_target(t::TransitionKey) = t[1]
@inline transition_source(t::TransitionKey) = t[2]
@inline transition_symbol(t::TransitionKey) = t[3]

abstract type AbstractTransitionMetadata end

struct NoTransitionMetadata <: AbstractTransitionMetadata end

has_metadata(::NoTransitionMetadata) = false
get_metadata(::NoTransitionMetadata, tr::TransitionKey) = nothing
add_metadata!(::NoTransitionMetadata, tr::TransitionKey, value) = nothing

abstract type AbstractTransitionInfo end

mutable struct TransitionMetadata{K, V} <: AbstractTransitionMetadata
    data::Dict{K, Vector{V}}
end

TransitionMetadata{K, V}() where {K, V} = TransitionMetadata{K, V}(Dict{K, Vector{V}}())

TransitionMetadata() = TransitionMetadata{TransitionKey, AbstractTransitionInfo}()

has_metadata(::TransitionMetadata) = true

get_metadata(md::TransitionMetadata{K, V}, tr::K) where {K, V} = get(md.data, tr, V[])

function add_metadata!(md::TransitionMetadata{K, V}, tr::K, value::V) where {K, V}
    push!(get!(md.data, tr, V[]), value)
    return nothing
end

struct ConcreteTransitionSample{X, Y} <: AbstractTransitionInfo
    x_source::X
    x_target::Y
end
