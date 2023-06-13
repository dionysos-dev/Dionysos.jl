# LazySetMinus(A, B) = A \ B
struct LazySetMinus{SA, SB}
    A::SA
    B::SB
end

function get_B(LSM::LazySetMinus)
    return LSM.B
end

# LazySetMinus = sets[1] ∪ sets[2] ∪ ... ∪ sets[n]
struct LazyUnionSetArray{S}
    sets::Vector{S}
end

import Base.isempty

function isempty(A::LazyUnionSetArray)
    return Base.isempty(A.sets)
end

function get_sets(A::LazyUnionSetArray)
    if isempty(A)
        return []
    else
        return A.sets
    end
end

function Plots.plot!(sets::LazyUnionSetArray{S}; dims=[1,2], color=:yellow, opacity=0.2) where {S}
    for set in sets.sets
        plot!(set; dims=dims, color=color, opacity=opacity)
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

function Plots.plot!(set::LazySetMinus{SA, SB}; dims=[1,2], color=:yellow, opacity=0.2) where {SA, SB}
    plot!(set.A; dims=dims, color=color, opacity=opacity)
    plot!(set.B; dims=dims, color=:black, opacity=1.0)
end