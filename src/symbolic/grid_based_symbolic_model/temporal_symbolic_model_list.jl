"""
    TemporalSymbolicModelList{N, M, S1, S2, A, U, OS} <: GridBasedSymbolicModel{N, M}

A classical symbolic model where the entire domain is partitioned into grid cells. 
This version is adapted for abstractions using multiple time steps, where each input is augmented with a time step.
"""
mutable struct TemporalSymbolicModelList{
    N,
    M,
    S1 <: DO.GridDomainType{N},
    S2 <: DO.CustomList{M},
    T,
    A,
    U,
    OS,
} <: GridBasedSymbolicModel{N, M}
    Xdom::S1
    Udom::S2
    Tsteps::Vector{T}
    autom::A
    xpos2int::Dict{NTuple{N, Int}, Int}
    xint2pos::Vector{NTuple{N, Int}}

    # for easiness of implementation
    tstep2int::Dict{T, Int}                       # Maps time step to integer
    uint2coord::Vector{U}                         # Maps integer to input
    ucoord2int::Dict{U, Int}                      # Maps input to integer

    # this is the true symbolic encoding
    augmented_ucoord2int::Dict{Tuple{U, T}, Int}  # Maps (input, time) to integer
    augmented_uint2coord::Vector{Tuple{U, T}}     # Maps integer to (input, time)
    original_symmodel::OS
end

"""
    TemporalSymbolicModelList(Xdom, Udom, Tsteps; AutomatonConstructor)
Constructor for a fresh (non-determinized) TemporalSymbolicModelList.
"""
function TemporalSymbolicModelList(
    Xdom,
    Udom,
    Tsteps,
    AutomatonConstructor::Function = (n, m) -> NewSortedAutomatonList(n, m),
)
    nx = DO.get_ncells(Xdom)
    xint2pos = [pos for pos in DO.enum_pos(Xdom)]
    xpos2int = Dict((pos, i) for (i, pos) in enumerate(DO.enum_pos(Xdom)))

    customDomainList = DO.convert_to_custom_domain(Udom)

    tstep2int = Dict(t => i for (i, t) in enumerate(Tsteps))
    uint2coord = [coord for coord in DO.enum_elems(customDomainList)]
    ucoord2int = Dict((coord, i) for (i, coord) in enumerate(DO.enum_elems(customDomainList)))

    augmented_uint2coord = [(coord, t) for coord in DO.enum_elems(customDomainList) for t in Tsteps]
    augmented_ucoord2int = Dict(((coord, t), i) for (i, (coord, t)) in enumerate(augmented_uint2coord))

    nu = length(augmented_uint2coord)
    autom = AutomatonConstructor(nx, nu)  # Start with 0 inputs; will add later

    return TemporalSymbolicModelList(
        Xdom,
        customDomainList,
        Tsteps,
        autom,
        xpos2int,
        xint2pos,
        tstep2int,
        uint2coord,
        ucoord2int,
        augmented_ucoord2int,
        augmented_uint2coord,
        nothing,  # OS = Nothing for base models
    )
end

# function with_automaton(symmodel::TemporalSymbolicModelList, autom)
#     return TemporalSymbolicModelList(
#         symmodel.Xdom,
#         symmodel.Udom,
#         symmodel.Tsteps,
#         autom,
#         symmodel.xpos2int,
#         symmodel.xint2pos,
#         symmodel.augmented_ucoord2int,
#         symmodel.augmented_uint2coord,
#         symmodel.original_symmodel,
#     )
# end

get_n_state(symmodel::TemporalSymbolicModelList)      = length(symmodel.xint2pos)
get_n_input(symmodel::TemporalSymbolicModelList)      = length(symmodel.uint2coord)
get_n_aug_input(symmodel::TemporalSymbolicModelList)  = length(symmodel.augmented_uint2coord) 
get_n_time_steps(symmodel::TemporalSymbolicModelList) = length(symmodel.Tsteps)
enum_states(symmodel::TemporalSymbolicModelList)      = 1:get_n_state(symmodel)
enum_inputs(symmodel::TemporalSymbolicModelList)      = 1:get_n_input(symmodel)
enum_aug_inputs(symmodel::TemporalSymbolicModelList)  = 1:get_n_aug_input(symmodel)
enum_abstract_time_steps(symmodel::TemporalSymbolicModelList)  = 1:get_n_time_steps(symmodel)
get_state_domain(symmodel::TemporalSymbolicModelList) = symmodel.Xdom
get_input_domain(symmodel::TemporalSymbolicModelList) = symmodel.Udom

get_xpos_by_state(symmodel::TemporalSymbolicModelList, state) = symmodel.xint2pos[state]
get_state_by_xpos(symmodel::TemporalSymbolicModelList, xpos)  = symmodel.xpos2int[xpos]
is_xpos(symmodel::TemporalSymbolicModelList, xpos)            = haskey(symmodel.xpos2int, xpos)

pre(symmodel::TemporalSymbolicModelList, target::Int) = pre(symmodel.autom, target)
post(symmodel::TemporalSymbolicModelList, source::Int, input::Int) =
    post(symmodel.autom, source, input)
enum_transitions(symmodel::TemporalSymbolicModelList) = enum_transitions(symmodel.autom)
add_transition!(symmodel::TemporalSymbolicModelList, q::Int, q′::Int, u::Int) =
    add_transition!(symmodel.autom, q, q′, u)
add_transitions!(symmodel::TemporalSymbolicModelList, translist) =
    add_transitions!(symmodel.autom, translist)

is_determinized(symmodel::TemporalSymbolicModelList) = !(symmodel.original_symmodel === nothing)

"""
    is_deterministic(symmodel::SymbolicModelList) -> Bool

Returns `true` if the symbolic model is deterministic.
"""
is_deterministic(symmodel::TemporalSymbolicModelList) = is_deterministic(symmodel.autom)

function get_concrete_input(
    symmodel::TemporalSymbolicModelList{N, M, S1, S2, T, A, U, Nothing},
    input::Int,
) where {N, M, S1, S2, T, A, U}
    concrete_input = symmodel.augmented_uint2coord[input]
    return concrete_input
end

function get_abstract_input(
    symmodel::TemporalSymbolicModelList{N, M, S1, S2, T, A, U, Nothing},
    u::Tuple{U, T},
) where {N, M, S1, S2, T, A, U}
    return symmodel.augmented_ucoord2int[u]
end 

function get_concrete_true_input(
    symmodel::TemporalSymbolicModelList{N, M, S1, S2, T, A, U, Nothing},
    input::Int,
) where {N, M, S1, S2, T, A, U}
    return symmodel.uint2coord[input]
end

function get_abstract_true_input(
    symmodel::TemporalSymbolicModelList{N, M, S1, S2, T, A, U, Nothing},
    u::U,
) where {N, M, S1, S2, T, A, U}
    return symmodel.ucoord2int[u]
end

# TODO: determinized system
function get_concrete_input(
    symmodel::TemporalSymbolicModelList{N, M, S1, S2, T, A, Tuple{Uprev, Int}, OS},
    input::Int,
) where {N, M, S1, S2, T, A, Uprev, OS}
    u, _ = symmodel.uint2coord[input]
    return get_concrete_input(symmodel.original_symmodel, u)
end

function get_abstract_input(
    symmodel::TemporalSymbolicModelList{N, M, S1, S2, T, A, Tuple{Uprev, Int}, OS},
    u,
) where {N, M, S1, S2, T, A, Uprev, OS}
    return get_abstract_input(symmodel.original_symmodel, u)
end

"""
    determinize_symbolic_model(symmodel::SymbolicModelList) -> SymbolicModelList

Returns a determinized version of the given symbolic model by encoding each transition input as a pair `(input_symbol, target_state)`.
TODO
"""
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
#         new_input = (symbol, target)

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

function compute_abstract_transitions_from_rectangle!(
    symmodel::TemporalSymbolicModelList,
    reachable_set::UT.HyperRectangle,
    abstract_state,
    abstract_input,
    translist,
)
    Xdom = get_state_domain(symmodel)
    ypos_iter = DO.get_subset_pos_in_grid(Xdom, reachable_set, DO.OUTER)
    allin = true
    for ypos in ypos_iter
        if !(ypos in Xdom)
            allin = false
            break
        end
        target = get_state_by_xpos(symmodel, ypos)
        push!(translist, (target, abstract_state, abstract_input))
    end
    return allin
end

function compute_abstract_system_from_concrete_system!(
    abstract_system::TemporalSymbolicModelList,
    concrete_system_approx::Vector{ST.DiscreteTimeGrowthBound};
    verbose = false,
    update_interval = Int(1e5),
    threaded::Bool = false,
)   
    # If multithreading is not requested or only one thread is available -> sequential execution
    if !threaded || Threads.nthreads() == 1
        translist = Tuple{Int, Int, Int}[]
        growthbound_maps = [csa.growthbound_map for csa in concrete_system_approx]
        system_maps      = [ST.get_system_map(csa) for csa in concrete_system_approx]
        r = DO.get_h(DO.get_grid(get_state_domain(abstract_system))) / 2.0
        total_iterations = max(
            div(
                get_n_aug_input(abstract_system) * get_n_state(abstract_system),
                update_interval,
            ),
            1,
        )
        progress = verbose ? ProgressMeter.Progress(total_iterations) : nothing
        count = 0
        for abstract_input in enum_inputs(abstract_system)
            concrete_input = get_concrete_true_input(abstract_system, abstract_input)
            Fr = [gb(r, concrete_input) for gb in growthbound_maps]
            for abstract_state in enum_states(abstract_system)
                concrete_state = get_concrete_state(abstract_system, abstract_state)
                for abstract_timestep in enum_abstract_time_steps(abstract_system)
                    Fx = system_maps[abstract_timestep](concrete_state, concrete_input)
                    reachable_set = UT.HyperRectangle(Fx .- Fr[abstract_timestep], Fx .+ Fr[abstract_timestep])
                    Base.empty!(translist)
                    aug_abstract_input = get_abstract_input(
                        abstract_system,
                        (concrete_input, abstract_system.Tsteps[abstract_timestep]),
                    )
                    allin = compute_abstract_transitions_from_rectangle!(
                        abstract_system,
                        reachable_set,
                        abstract_state,
                        aug_abstract_input,
                        translist,
                    )
                    !allin && break # we don't consider longer time steps if not allin
                    add_transitions!(abstract_system, translist)
                    count += 1
                    verbose && count % update_interval == 0 && ProgressMeter.next!(progress)
                end
            end
        end
        verbose && ProgressMeter.finish!(progress)
        return
    end

    # ---- Multithreaded implementation ----
    # TODO
    growthbound_map = concrete_system_approx.growthbound_map
    system_map = ST.get_system_map(concrete_system_approx)
    r = DO.get_h(DO.get_grid(get_state_domain(abstract_system))) / 2.0

    inputs = collect(enum_inputs(abstract_system))
    states = collect(enum_states(abstract_system))
    Xdom = get_state_domain(abstract_system)

    input_data = Dict{Int, Tuple{Any, Any}}()
    for abstract_input in inputs
        concrete_input = get_concrete_input(abstract_system, abstract_input)
        Fr = growthbound_map(r, concrete_input)
        input_data[abstract_input] = (concrete_input, Fr)
    end

    nthreads = Threads.nthreads()

    total_work = length(inputs) * length(states)

    transitions_by_thread = [Vector{Tuple{Int, Int, Int}}() for _ in 1:nthreads]

    progress =
        verbose ? ProgressMeter.Progress(total_work ÷ max(1, update_interval ÷ 100)) :
        nothing
    progress_count = Threads.Atomic{Int}(0)

    Threads.@threads for linear_idx in 1:total_work
        tid = Threads.threadid()
        local_transitions = transitions_by_thread[tid]

        input_idx = ((linear_idx - 1) ÷ length(states)) + 1
        state_idx = ((linear_idx - 1) % length(states)) + 1

        @inbounds begin
            abstract_input = inputs[input_idx]
            abstract_state = states[state_idx]

            concrete_input, Fr = input_data[abstract_input]
            concrete_state = get_concrete_state(abstract_system, abstract_state)

            Fx = system_map(concrete_state, concrete_input)
            reachable_set = UT.HyperRectangle(Fx - Fr, Fx + Fr)

            prev_length = length(local_transitions)
            allin = compute_abstract_transitions_from_rectangle!(
                abstract_system,
                reachable_set,
                abstract_state,
                abstract_input,
                local_transitions,
            )

            if !allin
                resize!(local_transitions, prev_length)
            end
        end

        if verbose
            count_val = Threads.atomic_add!(progress_count, 1)
            if count_val % max(1, update_interval ÷ 100) == 0
                ProgressMeter.next!(progress)
            end
        end
    end

    for (_, local_transitions) in enumerate(transitions_by_thread)
        if !isempty(local_transitions)
            add_transitions!(abstract_system, local_transitions)
        end
    end

    verbose && ProgressMeter.finish!(progress)
    return
end 