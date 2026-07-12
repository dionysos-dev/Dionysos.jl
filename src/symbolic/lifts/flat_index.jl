"""
    FlatIndex{K}

Bijection between composite keys of type `K` (e.g. an augmented-state tuple) and a
dense integer range `1:n`. Backs the flattened integer numbering used by
product/lifted symbolic models.

Both directions are O(1): `id_to_key` is a plain vector (id → key), `key_to_id` a
hash map (key → id). Absence is encoded by the id `0`, so callers can branch on
`flat_id(...) > 0` without a second lookup.
"""
struct FlatIndex{K}
    id_to_key::Vector{K}
    key_to_id::Dict{K, Int}
end

FlatIndex{K}() where {K} = FlatIndex{K}(Vector{K}(), Dict{K, Int}())

"""
    FlatIndex(keys::AbstractVector{K}) -> FlatIndex{K}

Build a `FlatIndex` from a vector of **unique** keys, numbering them `1:length(keys)`
in order.
"""
function FlatIndex(keys::AbstractVector{K}) where {K}
    n = length(keys)
    id_to_key = Vector{K}(undef, n)
    key_to_id = Dict{K, Int}()
    sizehint!(key_to_id, n)
    @inbounds for i in 1:n
        id_to_key[i] = keys[i]
        key_to_id[keys[i]] = i
    end
    return FlatIndex{K}(id_to_key, key_to_id)
end

"Number of flattened ids."
@inline n_flat(fi::FlatIndex) = length(fi.id_to_key)

"Key associated with flattened id `id` (bounds-checked by the caller)."
@inline flat_key(fi::FlatIndex, id::Int) = fi.id_to_key[id]

"Flattened id of `key`, or `0` if the key is absent."
@inline flat_id(fi::FlatIndex{K}, key::K) where {K} = get(fi.key_to_id, key, 0)
