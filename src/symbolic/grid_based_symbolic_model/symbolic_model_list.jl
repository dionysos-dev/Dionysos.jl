"""
SymbolicModelList:
- XMapping / UMapping define universes and conversions
- Xset: states we build from (sources)
- Rset: allowed states as targets (superset; "relation universe allowance")
- Uset: inputs considered
"""
mutable struct SymbolicModelList{
    N,
    M,
    XM <: Mappings.AbstractMapping{N,Any},
    UM <: Mappings.AbstractMapping{M,Any},
    XS <: AbstractIdSet,
    RS <: AbstractIdSet,
    US <: AbstractIdSet,
    A,
    OS,
} <: GridBasedSymbolicModel{N,M}
    XMapping::XM
    UMapping::UM
    Xset::XS
    Rset::RS
    Uset::US
    autom::A
    original_symmodel::OS
end

# -----------------------
# State/Input enumeration
# -----------------------

get_n_state(sym::SymbolicModelList) = Mappings.get_n_state(sym.XMapping)
get_n_input(sym::SymbolicModelList) = Mappings.get_n_state(sym.UMapping)

enum_states(sym) = enum_states(sym.Xset, sym.XMapping)
enum_inputs(sym) = enum_states(sym.Uset, sym.UMapping)

is_active_state(sym,q) = contains_state(sym.Xset, sym.XMapping, q)
is_allowed_state(sym,q) = contains_state(sym.Rset, sym.XMapping, q)

enum_all_states(sym::SymbolicModelList) = Mappings.enum_states(sym.XMapping)
enum_all_inputs(sym::SymbolicModelList) = Mappings.enum_states(sym.UMapping)

enum_active_states(sym::SymbolicModelList) = enum_states(sym)
enum_allowed_states(sym::SymbolicModelList) = enum_ids(sym.Rset)

is_active_state(sym::SymbolicModelList, q::Int) = (q in sym.Xset)
is_allowed_state(sym::SymbolicModelList, q::Int) = (q in sym.Rset)

# For now, "domain" getters return sets (not grids)
get_state_domain(sym::SymbolicModelList) = sym.Xset
get_input_domain(sym::SymbolicModelList) = sym.Uset

# -----------------------
# Mapping-based conversions
# -----------------------

get_xpos_by_state(sym::SymbolicModelList, q::Int) =
    Mappings.get_pos_by_state(sym.XMapping, q)

get_state_by_xpos(sym::SymbolicModelList, pos) =
    Mappings.get_state_by_pos(sym.XMapping, pos)

is_xpos(sym::SymbolicModelList, pos) =
    Mappings.is_valid_pos(sym.XMapping, pos)

get_concrete_state(sym::SymbolicModelList, q::Int) =
    Mappings.get_coord_by_state(sym.XMapping, q)

get_abstract_state(sym::SymbolicModelList, x) =
    Mappings.get_state_by_coord(sym.XMapping, x)

get_concrete_input(sym::SymbolicModelList, u::Int) =
    Mappings.get_coord_by_state(sym.UMapping, u)

get_abstract_input(sym::SymbolicModelList, ucoord) =
    Mappings.get_state_by_coord(sym.UMapping, ucoord)

# -----------------------
# Automaton passthrough
# -----------------------

pre(sym::SymbolicModelList, target::Int) = pre(sym.autom, target)
post(sym::SymbolicModelList, source::Int, input::Int) = post(sym.autom, source, input)
enum_transitions(sym::SymbolicModelList) = enum_transitions(sym.autom)
add_transition!(sym::SymbolicModelList, q::Int, q′::Int, u::Int) = add_transition!(sym.autom, q, q′, u)
add_transitions!(sym::SymbolicModelList, translist) = add_transitions!(sym.autom, translist)
is_deterministic(sym::SymbolicModelList) = is_deterministic(sym.autom)

# -----------------------
# Transition filter logic you wanted
# -----------------------

"""
should_keep_target(sym, q′; keep_outside=false)

Keep if:
- q′ ∈ Xset
- OR (keep_outside && q′ ∈ Rset)
"""
should_keep_target(sym, q′; keep_outside=false) =
    contains_state(sym.Xset, sym.XMapping, q′) ||
    (keep_outside && contains_state(sym.Rset, sym.XMapping, q′))







# """
#     SymbolicModelList(Xdom, Udom; AutomatonConstructor)

# Constructor for a fresh (non-determinized) SymbolicModelList.
# """
# function SymbolicModelList(
#     Xdom,
#     Udom,
#     AutomatonConstructor::Function = (n, m) -> NewSortedAutomatonList(n, m),
# )
#     nx = DO.get_ncells(Xdom)
#     xint2pos = [pos for pos in DO.enum_pos(Xdom)]
#     xpos2int = Dict((pos, i) for (i, pos) in enumerate(DO.enum_pos(Xdom)))

#     customDomainList = DO.convert_to_custom_domain(Udom)
#     nu = DO.get_ncells(customDomainList)
#     uint2coord = [coord for coord in DO.enum_elems(customDomainList)]
#     ucoord2int =
#         Dict((coord, i) for (i, coord) in enumerate(DO.enum_elems(customDomainList)))
#     autom = AutomatonConstructor(nx, nu)

#     return SymbolicModelList(
#         Xdom,
#         customDomainList,
#         autom,
#         xpos2int,
#         xint2pos,
#         ucoord2int,
#         uint2coord,
#         nothing, # OS = Nothing for base models
#     )
# end

# function with_automaton(symmodel::SymbolicModelList, autom)
#     return SymbolicModelList(
#         symmodel.Xdom,
#         symmodel.Udom,
#         autom,
#         symmodel.xpos2int,
#         symmodel.xint2pos,
#         symmodel.ucoord2int,
#         symmodel.uint2coord,
#         symmodel.original_symmodel,
#     )
# end

# get_n_state(symmodel::SymbolicModelList) = Relations.get_n_state(symmodel.relation)
# enum_states(symmodel::SymbolicModelList) = Relations.enum_states(symmodel.relation)
# get_xpos_by_state(symmodel::SymbolicModelList, state) = Relations.state_to_pos(symmodel.relation, state)
# get_state_by_xpos(symmodel::SymbolicModelList, xpos) = Relations.pos_to_state(symmodel.relation, xpos)
# is_xpos(symmodel::SymbolicModelList, xpos) = Relations.is_valid_pos(symmodel.relation, xpos)
# get_state_domain(symmodel::SymbolicModelList) = symmodel.Xdom

# get_n_input(symmodel::SymbolicModelList) = length(symmodel.uint2coord)
# enum_inputs(symmodel::SymbolicModelList) = 1:get_n_input(symmodel)
# get_input_domain(symmodel::SymbolicModelList) = symmodel.Udom

# pre(symmodel::SymbolicModelList, target::Int) = pre(symmodel.autom, target)
# post(symmodel::SymbolicModelList, source::Int, input::Int) =
#     post(symmodel.autom, source, input)
# enum_transitions(symmodel::SymbolicModelList) = enum_transitions(symmodel.autom)
# add_transition!(symmodel::SymbolicModelList, q::Int, q′::Int, u::Int) =
#     add_transition!(symmodel.autom, q, q′, u)
# add_transitions!(symmodel::SymbolicModelList, translist) =
#     add_transitions!(symmodel.autom, translist)

# is_determinized(symmodel::SymbolicModelList) = !(symmodel.original_symmodel === nothing)

# """
#     is_deterministic(symmodel::SymbolicModelList) -> Bool

# Returns `true` if the symbolic model is deterministic.
# """
# is_deterministic(symmodel::SymbolicModelList) = is_deterministic(symmodel.autom)

# function get_concrete_input(
#     symmodel::SymbolicModelList{N, M, S1, S2, A, U, Nothing},
#     input::Int,
# ) where {N, M, S1, S2, A, U}
#     return symmodel.uint2coord[input]
# end

# function get_abstract_input(
#     symmodel::SymbolicModelList{N, M, S1, S2, A, U, Nothing},
#     u::U,
# ) where {N, M, S1, S2, A, U}
#     return symmodel.ucoord2int[u]
# end

# function get_concrete_input(
#     symmodel::SymbolicModelList{N, M, S1, S2, A, Tuple{Uprev, Int}, OS},
#     input::Int,
# ) where {N, M, S1, S2, A, Uprev, OS}
#     u, _ = symmodel.uint2coord[input]
#     return get_concrete_input(symmodel.original_symmodel, u)
# end

# function get_abstract_input(
#     symmodel::SymbolicModelList{N, M, S1, S2, A, Tuple{Uprev, Int}, OS},
#     u,
# ) where {N, M, S1, S2, A, Uprev, OS}
#     return get_abstract_input(symmodel.original_symmodel, u)
# end

# """
#     determinize_symbolic_model(symmodel::SymbolicModelList) -> SymbolicModelList

# Returns a determinized version of the given symbolic model by encoding each transition input as a pair `(input_symbol, target_state)`.
# """
# function determinize_symbolic_model(
#     symmodel::SymbolicModelList;
#     AutomatonConstructor::Function = (n, m) -> NewSortedAutomatonList(n, m),
# )
#     transitions = enum_transitions(symmodel)

#     U = eltype(symmodel.uint2coord)
#     new_ucoord2int = Dict{Tuple{U, Int}, Int}()
#     new_uint2coord = Tuple{U, Int}[]

#     transition_buffer = Vector{NTuple{3, Int}}()
#     for (target, source, symbol) in transitions
#         u_coord = symmodel.uint2coord[symbol]  # Get symbolic input
#         new_input = (u_coord, target)          # Determinize with symbolic input and target state

#         input_id = get!(new_ucoord2int, new_input) do
#             push!(new_uint2coord, new_input)
#             return length(new_uint2coord)
#         end

#         push!(transition_buffer, (target, source, input_id))
#     end
#     new_autom = AutomatonConstructor(get_n_state(symmodel), length(new_uint2coord))
#     add_transitions!(new_autom, transition_buffer)

#     new_symmodel = SymbolicModelList(
#         symmodel.Xdom,
#         symmodel.Udom,
#         new_autom,
#         symmodel.xpos2int,
#         symmodel.xint2pos,
#         new_ucoord2int,
#         new_uint2coord,
#         symmodel,
#     )

#     return new_symmodel
# end
