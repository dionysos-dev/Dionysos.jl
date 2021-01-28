mutable struct SortedTupleSet{N,T<:Integer} <: AbstractSet{NTuple{N,T}}
    data::Vector{NTuple{N,T}}
    is_sorted::Bool
    function SortedTupleSet{N,T}() where {N,T}
        return new(NTuple{N,T}[], true)
    end
end

function push_new!(set::SortedTupleSet, x::Tuple)
    push!(set.data, x)
    set.is_sorted = false
    return set
end
function append_new!(set::SortedTupleSet, xs)
    append!(set.data, xs)
    set.is_sorted = false
    return set
end

function Base.empty!(set::SortedTupleSet)
    empty!(set.data)
    set.is_sorted = true
end

Base.length(set::SortedTupleSet) = length(set.data)

function ensure_sorted!(set::SortedTupleSet)
    if !set.is_sorted
        sort!(set.data)
        set.is_sorted = true
    end
end

drop_first(x::NTuple{2}) = (x[2],)
drop_first(x::NTuple{3}) = (x[2], x[3])

function fix_and_eliminate_first(set::SortedTupleSet, value)
    ensure_sorted!(set)
    idxlist = searchsorted(set.data, (value,), by = x -> x[1])
    return Base.Generator(drop_first, view(set.data, idxlist))
end

function fix_and_eliminate_tail!(output, set::SortedTupleSet, values)
    ensure_sorted!(set)
    for el in set.data
        if Base.tail(el) == values
            push!(output, first(el))
        end
    end
end
