struct IndColl{E}
    e2i::Dict{E,UInt}
    i2e::Vector{E}
end

IndColl{E}() where E = IndColl{E}(Dict{E,UInt}(), E[])

@inline function _push!(e2i::Dict{E,UInt}, i2e::Vector{E}, elem::E) where E
    if !haskey(e2i, elem)
        push!(e2i, elem => length(i2e) + 1)
        push!(i2e, elem)
    end
end

function Base.push!(coll::IndColl, elem) where E
    _push!(coll.e2i, coll.i2e, elem)
    return coll
end

function Base.append!(coll::IndColl, elems)
    for elem in elems
        _push!(coll.e2i, coll.i2e, elem)
    end
    return coll
end

get_elems(coll::IndColl) = coll.i2e
get_indexes(coll::IndColl) = eachindex(coll.i2e)
Base.length(coll::IndColl) = length(coll.i2e)
