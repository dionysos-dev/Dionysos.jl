"""
    Abstract Type: SymbolicModel{N, M}

Defines a generic symbolic model interface, where:
- `N` is the state space dimension.
- `M` is the input space dimension.
"""
abstract type SymbolicModel{N, M} end

# -----------------------
# State/Input enumeration
# -----------------------

function get_state_mapping(::SymbolicModel) end
function get_input_mapping(::SymbolicModel) end

function get_state_domain(::SymbolicModel) end
function get_allowed_state_domain(::SymbolicModel) end
function get_input_domain(::SymbolicModel) end


get_n_state(sym::SymbolicModel) = MP.get_n_state(get_state_domain(sym), get_state_mapping(sym))
get_n_allowed_state(sym::SymbolicModel) = MP.get_n_state(get_allowed_state_domain(sym), get_state_mapping(sym))
get_n_input(sym::SymbolicModel) = MP.get_n_state(get_input_domain(sym), get_input_mapping(sym))

enum_states(sym::SymbolicModel) = MP.enum_states(get_state_domain(sym), get_state_mapping(sym))
enum_allowed_states(sym::SymbolicModel) = MP.enum_states(get_allowed_state_domain(sym), get_state_mapping(sym))
enum_inputs(sym::SymbolicModel) = MP.enum_states(get_input_domain(sym), get_input_mapping(sym))

is_state(sym::SymbolicModel, q::Int) = MP.contains_state(get_state_domain(sym), get_state_mapping(sym), q)
is_allowed_state(::SymbolicModel, ::Nothing) = false
is_allowed_state(sym::SymbolicModel, q::Int) = MP.contains_state(get_allowed_state_domain(sym), get_state_mapping(sym), q)
is_input(sym::SymbolicModel, q::Int) = MP.contains_state(get_input_domain(sym), get_input_mapping(sym), q)

get_state_dim(sym::SymbolicModel) = MP.get_dim(get_state_mapping(sym))
get_input_dim(sym::SymbolicModel) = MP.get_dim(get_input_mapping(sym))

get_concrete_state(sym::SymbolicModel, state) = MP.get_coord_by_state(get_state_mapping(sym), state)
get_concrete_input(sym::SymbolicModel, input) = MP.get_coord_by_state(get_input_mapping(sym), input)
get_abstract_state(sym::SymbolicModel, x) = MP.get_state_by_coord(get_state_mapping(sym), x)
get_abstract_states(sym::SymbolicModel, x) = MP.get_states_by_coord(get_state_mapping(sym), x)
is_state_cover(sym::SymbolicModel) = MP.is_state_cover(get_state_mapping(sym))

get_abstract_input(sym::SymbolicModel, u) = MP.get_state_by_coord(get_input_mapping(sym), u)

get_concrete_elem(sym::SymbolicModel, q::Int) = MP.get_elem_by_state(get_state_mapping(sym), q)

function get_automaton(::SymbolicModel) end
pre(sym::SymbolicModel, target::Int) = pre(get_automaton(sym), target)
post(sym::SymbolicModel, source::Int, input::Int) = post(get_automaton(sym), source, input)
enum_transitions(sym::SymbolicModel) = enum_transitions(get_automaton(sym))
add_transition!(sym::SymbolicModel, q::Int, q′::Int, u::Int) = add_transition!(get_automaton(sym), q, q′, u)
add_transitions!(sym::SymbolicModel, translist) = add_transitions!(get_automaton(sym), translist)
get_n_transitions(sym::SymbolicModel) = length(enum_transitions(sym))
is_deterministic(sym::SymbolicModel) = is_deterministic(get_automaton(sym))
function is_determinized(sym::SymbolicModel) end


function get_states_from_set(
    sym::SymbolicModel,
    set,
    incl_mode::MP.INCL_MODE,
)
    return MP.get_states_from_set(get_state_mapping(sym), set, incl_mode)
end

function get_states_from_set_strict(
    sym::SymbolicModel,
    set,
    incl_mode::MP.INCL_MODE,
)
    return MP.get_states_from_set_strict(get_state_mapping(sym), set, incl_mode)
end


function get_state_set_from_states(sym::SymbolicModel, states)
    return MP.stateset_from_states(get_state_mapping(sym), states)
end

"""
    determinize_symbolic_model(sym; AutomatonConstructor, convert_U_to_list)

Return a deterministic symbolic model by refining the input alphabet:
each original input `u` is replaced by a pair `(u_coord, target_state)`.

This is useful to turn nondeterministic transitions into deterministic ones
by making the target part of the symbol.
"""
function determinize_symbolic_model(
    sym::SymbolicModel{N,M};
    AutomatonConstructor::Function = (n, m) -> NewSortedAutomatonList(n, m),
) where {N,M}

    Umap  = get_input_mapping(sym)
    Xmap  = get_state_mapping(sym)
    Xset  = get_state_domain(sym)
    Rset  = get_allowed_state_domain(sym)

    trans = enum_transitions(sym)

    # infer coordinate scalar type
    u_any = first(enum_inputs(sym))
    u1 = MP.get_coord_by_state(Umap, u_any)
    T  = eltype(u1)

    TU = SVector{M+1,T} # augmented coordinate type

    new_ucoord2int = Dict{TU, Int}()
    new_uint2coord = TU[]  

    new_transitions = NTuple{3,Int}[]
    sizehint!(new_transitions, length(trans))

    for (q′, q, u) in trans
        u_coord = MP.get_coord_by_state(Umap, u)

        u_aug = TU(ntuple(i -> i <= M ? u_coord[i] : T(q′), M+1))

        u_new_id = get!(new_ucoord2int, u_aug) do
            push!(new_uint2coord, u_aug)
            length(new_uint2coord)
        end

        push!(new_transitions, (q′, q, u_new_id))
    end

    new_autom = AutomatonConstructor(get_n_state(sym), length(new_uint2coord))
    add_transitions!(new_autom, new_transitions)

    new_Umap = MP.ListMapping(new_uint2coord)

    new_sym = SymbolicModelList(
        Xmap, new_Umap;
        Xset = Xset,
        Rset = Rset,
        Uset = MP.MappingSet{M+1}(),
        automaton_constructor = (n,m)->new_autom,
        original_symmodel = sym,
        convert_U_to_list = true,
    )

    return new_sym
end
@recipe function f(
    sym::SymbolicModel;
    with_arrows = false,
    dims = [1, 2],
    value_function = nothing, 
)
    @series begin
        dims := dims
        value_function := value_function
        return ((get_state_domain(sym), get_state_mapping(sym)),)
    end
    if with_arrows
        for t in enum_transitions(sym)
            color = RGB(
                abs(0.6 * sin(t[1])),
                abs(0.6 * sin(t[1] + 2π / 3)),
                abs(0.6 * sin(t[1] - 2π / 3)),
            )
            p1 = get_concrete_state(sym, t[2])
            p2 = get_concrete_state(sym, t[1])

            @series begin
                dims := dims
                color := color
                return t[1] == t[2] ? UT.DrawPoint(p1) : UT.DrawArrow(p1, p2)
            end
        end
    end
end