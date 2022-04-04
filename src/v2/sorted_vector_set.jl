mutable struct SortedTupleSet{N,T} <: AbstractSet{T}
    data::Vector{T}
    is_sorted::Bool
    function SortedTupleSet{N,T}() where {N,T}
        return new(T[], true)
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

function delete!(set::SortedTupleSet, x, comparison)
    filter!(e->!comparison(e,x), set.data)
    set.is_sorted = false
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
drop_first(x::Tuple{Int,Int,Int,Float64}) = (x[2], x[3], x[4])

function fix_and_eliminate_first(set::SortedTupleSet, value)
    ensure_sorted!(set)
    idxlist = searchsorted(set.data, (value,), by = x -> x[1])
    return Base.Generator(drop_first, view(set.data, idxlist))
end

function fix_and_eliminate_tail!(output, set::SortedTupleSet, values)
    ensure_sorted!(set)
    for el in set.data
        if (el[2],el[3]) == values
            ## pourrait ajouter ici if el[4]>0 si proba non nulle, on pourrait aussi ajouter les proba, à voir
            push!(output, (first(el),el[4]))
        end
    end
end
