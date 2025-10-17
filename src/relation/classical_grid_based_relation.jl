struct ClassicalGridBasedRelation <: GridBasedRelation{N, M}
    Xdom::Domain
    Udom::Domain
    grid::Grid
    input_grid::Grid
    pos2state::Dict{NTuple, Int}
    state2pos::Vector{NTuple}
end

is_deterministic(r::ClassicalGridBasedRelation)::Bool = true

get_state_domain(r::ClassicalGridBasedRelation) = r.Xdom
get_input_domain(r::ClassicalGridBasedRelation) = r.Udom




























# struct EllipsoidalGridBasedRelation <: GridBasedRelation
#     ...
# end


# struct ComposedRelation{R1<:AbstractRelation, R2<:AbstractRelation} <: AbstractRelation
#     R1::R1
#     R2::R2
# end
# get_abstract_states(R::ComposedRelation, x) = 
#     reduce(vcat, get_abstract_states(R.R2, s) for s in get_abstract_states(R.R1, x))
