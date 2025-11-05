abstract type LazySymbolicModel{N,M} <: SY.SymbolicModel{N,M} end

mutable struct LazyGridBasedSymbolicModel{N, M, S1, U} <: LazySymbolicModel{N, M}
    Xdom::S1
    state_mapping::Vector{Tuple{NTuple{N,Int}, S2}}
    input_mapping::Vector{U}
    autom::SY.AbstractAutomatonList
    pos2states::Dict{NTuple{N, Int}, Vector{Int}}
    abstract_transition_costs::Dict{Tuple{Int, Int}, Float64}
end

function add_transition_cost!(model::LazyGridBasedSymbolicModel, q::Int, s::Int, cost::Float64)
    model.abstract_transition_costs[(q, s)] = cost
end

# Constructor
function LazyGridBasedSymbolicModel(Xdom::S1, AutomatonConstructor::Function = (n, m) -> NewSortedAutomatonList(n, m)) where {S1, S2}
    N = DO.get_dim(Xdom)
    return LazyEllipsoidalModel{N, :, S1, ST.Controller}(
        Xdom,
        Tuple{NTuple{N, Int}, S2}[],
        ST.Controller[],
        AutomatonConstructor(0, 0),
        Dict{NTuple{N, Int}, Vector{Int}}(),
    )
end

# Interface methods
get_n_state(model::LazyGridBasedSymbolicModel) = length(model.state_mapping)
get_n_input(model::LazyGridBasedSymbolicModel) = length(model.input_mapping)
enum_states(model::LazyGridBasedSymbolicModel) = eachindex(model.state_mapping)
enum_inputs(model) = eachindex(model.input_mapping)
get_state_domain(model::LazyGridBasedSymbolicModel) = model.Xdom
get_input_domain(model::LazyGridBasedSymbolicModel) = nothing

function get_pos(model::LazyGridBasedSymbolicModel, q::Int)
    return model.state_mapping[q][1]
end

# function get_shape(model::LazyGridBasedSymbolicModel, q::Int)
#     return model.state_mapping[q][2]
# end

function get_center(model::LazyGridBasedSymbolicModel, q::Int)
    pos = get_pos(model, q)
    return DO.get_coord_by_pos(model.Xdom, pos)
end

function get_concrete_states(model::LazyGridBasedSymbolicModel, q::Int)
    return model.state_mapping[q][2]
end

function get_abstract_states(model::LazyGridBasedSymbolicModel, x::SVector)
    return [q for q in enum_states(model) if x ∈ get_concrete_states(model, q)]
end

function get_input_controller(model::LazyGridBasedSymbolicModel, s::Int)
    return model.input_mapping[s]
end

function get_concrete_input(model::LazyGridBasedSymbolicModel, s::Int, x::SVector)
    kappa = get_input_controller(model, s)
    return ST.get_control(kappa, x)
end

function get_states_by_pos!(model::LazyGridBasedSymbolicModel, pos::NTuple)
    return model.pos2states[pos]
end

function add_state!(model::LazyGridBasedSymbolicModel, pos, set)
    q = length(model.state_mapping) + 1
    push!(model.state_mapping, (pos, set))
    if haskey(model.pos2states, pos)
        push!(model.pos2states[pos], q)
    else
        model.pos2states[pos] = [q]
    end
    return q
end

function add_input!(model::LazyGridBasedSymbolicModel, κ::ST.Controller)
    push!(model.input_mapping, κ)
    return length(model.input_mapping)
end

# Transition interface
enum_transitions(model::LazyGridBasedSymbolicModel) = SY.enum_transitions(model.autom)
add_transition!(model::LazyGridBasedSymbolicModel, q, q′, u) = SY.add_transition!(model.autom, q, q′, u)
add_transitions!(model::LazyGridBasedSymbolicModel, translist) = SY.add_transitions!(model.autom, translist)
post(model::LazyGridBasedSymbolicModel, q, s) = SY.post(model.autom, q, s)
pre(model::LazyGridBasedSymbolicModel, q) = SY.pre(model.autom, q)
is_deterministic(model::LazyGridBasedSymbolicModel) = SY.is_deterministic(model.autom)

# Summary
function Base.show(io::IO, model::LazyGridBasedSymbolicModel)
    println(io, "LazyGridBasedSymbolicModel:")
    println(io, "  States: $(length(model.state_mapping))")
    println(io, "  Inputs: $(length(model.input_mapping))")
    println(io, "  Grid dim: $(DO.get_dim(model.Xdom))")
end

end  # module

    # x0 = get_center(model, q)
    # P = get_shape(model, q)
    # return UT.Ellipsoid(collect(P), collect(x0))