# LazySetMinus(A, B) = A \ B
struct LazySetMinus{SA, SB}
    A::SA
    B::SB
end

# LazySetMinus = sets[1] ∪ sets[2] ∪ ... ∪ sets[n]
struct LazyUnionSetArray{S}
    sets::Vector{S}
end

# Computes (A1 ∪ A2 ∪ ... ∪ An) ∩ B = (A1 ∩ B) ∪ (A2 ∩ B) ∪ ... ∪ (An ∩ B)
function Base.intersect(A::LazyUnionSetArray, B)    
    sets = typeof(A.sets[1])[]
    for set in A.sets
        push!(sets, B ∩ set)
    end
    return LazyUnionSetArray(sets)
end

# Computes A ∩ (B1 ∪ B2 ∪ ... ∪ Bn) = (A ∩ B1) ∪ (A ∩ B2) ∪ ... ∪ (A ∩ Bn)
function Base.intersect(A, B::LazyUnionSetArray)    
    sets = typeof(B.sets[1])[]
    for set in B.sets
        push!(sets, A ∩ set)
    end
    return LazyUnionSetArray(sets)
end