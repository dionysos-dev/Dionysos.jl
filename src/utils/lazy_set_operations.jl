# LazySetMinus(A, B) = A \ B
struct LazySetMinus{SA, SB}
    A::SA
    B::SB
end

function get_B(LSM::LazySetMinus)
    return LSM.B
end

get_dims(s::LazySetMinus) = get_dims(s.A)

# LazySetMinus = sets[1] ∪ sets[2] ∪ ... ∪ sets[n]
struct LazyUnionSetArray{S}
    sets::Vector{S}
end

function Base.isempty(A::LazyUnionSetArray)
    return Base.isempty(A.sets)
end

function get_sets(A::LazyUnionSetArray)
    if isempty(A)
        return []
    else
        return A.sets
    end
end

# Computes (A1 ∪ A2 ∪ ... ∪ An) ∩ B = (A1 ∩ B) ∪ (A2 ∩ B) ∪ ... ∪ (An ∩ B)
function Base.intersect(A::LazyUnionSetArray, B)
    if isempty(A)
        return LazyUnionSetArray([])
    end
    sets = typeof(A.sets[1])[]
    for set in A.sets
        push!(sets, B ∩ set)
    end
    return LazyUnionSetArray(sets)
end

# Computes A ∩ (B1 ∪ B2 ∪ ... ∪ Bn) = (A ∩ B1) ∪ (A ∩ B2) ∪ ... ∪ (A ∩ Bn)
function Base.intersect(A, B::LazyUnionSetArray)
    return B ∩ A
end

@recipe function f(set::LazyUnionSetArray; dims = [1, 2])
    dims := dims
    for elem in set.sets
        @series begin
            return elem
        end
    end
end

@recipe function f(set::LazySetMinus; dims = [1, 2])
    dims := dims
    @series begin
        return set.A
    end
    @series begin
        color := :black
        return set.B
    end
end
