# LazySetMinus(A, B) = A \ B
struct LazySetMinus{SA, SB}
    A::SA
    B::SB
end

function get_B(LSM::LazySetMinus)
    return LSM.B
end

get_dims(s::LazySetMinus) = get_dims(s.A)
Base.in(x, set::LazySetMinus) = x ∈ set.A && !(x ∈ set.B)

function set_in_period(
    set::LazySetMinus,
    periodic_dims::SVector{P, Int},
    periods::SVector{P, T},
    start::SVector{P, T},
) where {P, T}
    A_wrapped = set_in_period(set.A, periodic_dims, periods, start)
    B_wrapped = set_in_period(set.B, periodic_dims, periods, start)
    return LazySetMinus(A_wrapped, B_wrapped)
end

@recipe function f(
    set::LazySetMinus;
    dims = [1, 2],
    hole_color = :gray,
    hole_alpha = 1.0,
    label = "Set",
)
    dims := dims
    @series begin
        label := label
        return set.A
    end
    @series begin
        label := ""
        seriestype := :shape
        fillcolor := hole_color
        fillalpha := hole_alpha
        return set.B
    end
end

# LazySetMinus = sets[1] ∪ sets[2] ∪ ... ∪ sets[n]
struct LazyUnionSetArray{S}
    sets::Vector{S}
end

function Base.isempty(A::LazyUnionSetArray)
    return Base.isempty(A.sets)
end

Base.in(x, set::LazyUnionSetArray) = any(s -> x ∈ s, set.sets)

function get_sets(A::LazyUnionSetArray)
    if isempty(A)
        return []
    else
        return A.sets
    end
end

function set_in_period(
    A::LazyUnionSetArray,
    periodic_dims::SVector{P, Int},
    periods::SVector{P, T},
    start::SVector{P, T},
) where {T, P}
    wrapped_sets = LazyUnionSetArray([])
    for set in A.sets
        wrapped_set = set_in_period(set, periodic_dims, periods, start)
        Base.append!(wrapped_sets.sets, wrapped_set.sets)
    end
    return wrapped_sets
end

function Base.intersect(A::LazyUnionSetArray, S::LazySetMinus)
    # (⋃Ai) ∩ (B \ C)  =  ⋃(Ai ∩ (B \ C))
    sets = Any[]
    sizehint!(sets, length(A.sets))
    for ai in get_sets(A)
        push!(sets, intersect(ai, S))
    end
    return LazyUnionSetArray(sets)
end
Base.intersect(S::LazySetMinus, A::LazyUnionSetArray) = intersect(A, S)

# (A \ B) ∩ C  =  (A ∩ C) \ B
function Base.intersect(S::LazySetMinus, C)
    return LazySetMinus(intersect(S.A, C), S.B)
end

# C ∩ (A \ B)  =  (A ∩ C) \ B
Base.intersect(C, S::LazySetMinus) = intersect(S, C)

@recipe function f(set::LazyUnionSetArray; dims = [1, 2], label = "set")
    dims := dims

    first_series = true
    for elem in set.sets
        @series begin
            if first_series
                label := label
                first_series = false
            else
                label := ""
            end
            return elem
        end
    end
end
