struct IntersectionSet
    sets::Vector{Union{Ellipsoid, HyperRectangle}}
end

function Base.in(x, iset::IntersectionSet)
    return all(x in s for s in iset.sets)
end

compare_sets(set1, set2) = get_volume(set1) < get_volume(set2)

function Base.sort!(iset::IntersectionSet, cmp::Function = compare_sets)
    return sort!(iset.sets, cmp)
end

@recipe function f(iset::IntersectionSet)
    sets = sort(iset.sets; lt = compare_sets, rev = true)
    for set in sets
        @series begin
            return set
        end
    end
end
