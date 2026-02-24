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
is_allowed_state(sym::SymbolicModel, q::Int) = MP.contains_state(get_allowed_state_domain(sym), get_state_mapping(sym), q)
is_input(sym::SymbolicModel, q::Int) = MP.contains_state(get_input_domain(sym), get_input_mapping(sym), q)

get_state_dim(sym::SymbolicModel) = MP.get_dim(get_state_mapping(sym))
get_input_dim(sym::SymbolicModel) = MP.get_dim(get_input_mapping(sym))

get_concrete_state(sym::SymbolicModel, state) = MP.get_coord_by_state(get_state_mapping(sym), state)
get_concrete_input(sym::SymbolicModel, input) = MP.get_coord_by_state(get_input_mapping(sym), input)
get_abstract_state(sym::SymbolicModel, x) = MP.get_state_by_coord(get_state_mapping(sym), x)
get_abstract_input(sym::SymbolicModel, u) = MP.get_state_by_coord(get_input_mapping(sym), u)

get_concrete_elem(sym::SymbolicModel, q::Int) = MP.get_elem_by_state(get_state_mapping(sym), q)

function get_automaton(::SymbolicModel) end

pre(sym::SymbolicModel, target::Int) = pre(get_automaton(sym), target)
post(sym::SymbolicModel, source::Int, input::Int) = post(get_automaton(sym), source, input)
enum_transitions(sym::SymbolicModel) = enum_transitions(get_automaton(sym))
add_transition!(sym::SymbolicModel, q::Int, q′::Int, u::Int) = add_transition!(get_automaton(sym), q, q′, u)
add_transitions!(sym::SymbolicModel, translist) = add_transitions!(get_automaton(sym), translist)
is_deterministic(sym::SymbolicModel) = is_deterministic(get_automaton(sym))
get_n_transitions(sym::SymbolicModel) = length(enum_transitions(sym))

function get_states_from_set(
    sym::SymbolicModel,
    set,
    incl_mode::MP.INCL_MODE,
)
    return MP.get_states_from_set(get_state_mapping(sym), set, incl_mode)
end


"""
    GridBasedSymbolicModel{N, M} <: SymbolicModel{N, M}

An intermediate abstract type for symbolic models that rely on a grid-based discretization.
- `N`: Dimension of the state space.
- `M`: Dimension of the input space.
"""
abstract type GridBasedSymbolicModel{N, M} <: SymbolicModel{N, M} end

function compute_abstract_transitions_from_rectangle!(
    symmodel::GridBasedSymbolicModel,
    reachable_set::UT.HyperRectangle,
    abstract_state,
    abstract_input,
    translist,
)
    targets = get_states_from_set(symmodel, reachable_set, MP.OUTER)
    allin = true
    for target in targets
        if !is_allowed_state(symmodel, target)
            allin = false
            break
        end
        push!(translist, (target, abstract_state, abstract_input))
    end
    return allin
end

function compute_abstract_transitions_from_points!(
    symmodel::GridBasedSymbolicModel,
    reachable_points,
    abstract_state,
    abstract_input,
    translist,
)
    allin = true
    for y in reachable_points
        target = get_abstract_state(symmodel, y)
        if !is_allowed_state(symmodel, target)
            allin = false
            break
        end
        push!(translist, (target, abstract_state, abstract_input))
    end
    unique!(translist)
    return allin
end

function compute_abstract_system_from_concrete_system!(
    abstract_system::GridBasedSymbolicModel,
    concrete_system_approx::ST.DiscreteTimeSystemOverApproximation;
    verbose = false,
    update_interval = Int(1e5),
    progress_dt = 0.2,
    threaded::Bool = false,
)
    # If multithreading is not requested or only one thread is available -> sequential execution
    if !threaded || Threads.nthreads() == 1
        translist = Tuple{Int, Int, Int}[]
        compute_reachable_set = ST.get_over_approximation_map(concrete_system_approx)
        total_iterations = max(
            div(
                get_n_input(abstract_system) * get_n_state(abstract_system),
                update_interval,
            ),
            1,
        )
        progress = verbose ? ProgressMeter.Progress(total_iterations) : nothing
        count = 0
        for abstract_input in enum_inputs(abstract_system)
            concrete_input = get_concrete_input(abstract_system, abstract_input)
            for abstract_state in enum_states(abstract_system)
                concrete_elem = get_concrete_elem(abstract_system, abstract_state)
                reachable_set = compute_reachable_set(concrete_elem, concrete_input)
                empty!(translist)
                allin = compute_abstract_transitions_from_rectangle!(
                    abstract_system,
                    reachable_set,
                    abstract_state,
                    abstract_input,
                    translist,
                )
                allin && add_transitions!(abstract_system, translist)
                count += 1
                verbose && count % update_interval == 0 && ProgressMeter.next!(progress)
            end
        end
        verbose && ProgressMeter.finish!(progress)
        return
    end

    # ---- Multithreaded implementation ----
    compute_reachable_set = ST.get_over_approximation_map(concrete_system_approx)
    inputs = collect(enum_inputs(abstract_system))
    states = collect(enum_states(abstract_system))

    nthreads = Threads.nthreads()

    total_work = length(inputs) * length(states)

    transitions_by_thread = [Vector{Tuple{Int, Int, Int}}() for _ in 1:nthreads]

    progress =
        verbose ? ProgressMeter.Progress(total_work ÷ max(1, update_interval)) : nothing
    progress_count = Threads.Atomic{Int}(0)

    Threads.@threads for linear_idx in 1:total_work
        tid = Threads.threadid()
        local_transitions = transitions_by_thread[tid]

        input_idx = ((linear_idx - 1) ÷ length(states)) + 1
        state_idx = ((linear_idx - 1) % length(states)) + 1

        @inbounds begin
            abstract_input = inputs[input_idx]
            abstract_state = states[state_idx]

            concrete_input = get_concrete_input(abstract_system, abstract_input)
            concrete_elem = get_concrete_elem(abstract_system, abstract_state)

            reachable_set = compute_reachable_set(concrete_elem, concrete_input)

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
            if count_val % max(1, update_interval) == 0
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

function compute_abstract_system_from_concrete_system!(
    abstract_system::GridBasedSymbolicModel,
    concrete_system_approx::ST.DiscreteTimeGrowthBound;
    verbose = false,
    update_interval = Int(1e5),
    progress_dt = 0.2,
    threaded::Bool = false,
)

    # If multithreading is not requested or only one thread is available -> sequential execution
    if !threaded || Threads.nthreads() == 1
        translist = Tuple{Int, Int, Int}[]
        growthbound_map = concrete_system_approx.growthbound_map
        system_map = ST.get_system_map(concrete_system_approx)
        r = MP.get_h(MP.get_grid(get_state_mapping(abstract_system))) / 2.0
        total_iterations = max(
            div(
                get_n_input(abstract_system) * get_n_state(abstract_system),
                update_interval,
            ),
            1,
        )
        progress = verbose ? ProgressMeter.Progress(total_iterations) : nothing
        count = 0
        for abstract_input in enum_inputs(abstract_system)
            concrete_input = get_concrete_input(abstract_system, abstract_input)
            Fr = growthbound_map(r, concrete_input)
            for abstract_state in enum_states(abstract_system)
                concrete_state = get_concrete_state(abstract_system, abstract_state)
                Fx = system_map(concrete_state, concrete_input)
                reachable_set = UT.HyperRectangle(Fx - Fr, Fx + Fr)
                Base.empty!(translist)
                allin = compute_abstract_transitions_from_rectangle!(
                    abstract_system,
                    reachable_set,
                    abstract_state,
                    abstract_input,
                    translist,
                )
                allin && add_transitions!(abstract_system, translist)
                count += 1
                verbose && count % update_interval == 0 && ProgressMeter.next!(progress)
            end
        end
        verbose && ProgressMeter.finish!(progress)
        return
    end

    # ---- Multithreaded implementation ----
    growthbound_map = concrete_system_approx.growthbound_map
    system_map = ST.get_system_map(concrete_system_approx)
    r = MP.get_h(MP.get_grid(get_state_mapping(abstract_system))) / 2.0

    inputs = collect(enum_inputs(abstract_system))
    states = collect(enum_states(abstract_system))

    input_data = Dict{Int, Tuple{Any, Any}}()
    for abstract_input in inputs
        concrete_input = get_concrete_input(abstract_system, abstract_input)
        Fr = growthbound_map(r, concrete_input)
        input_data[abstract_input] = (concrete_input, Fr)
    end

    # ---- Progress bar parameters ----
    progress_dt_ns = Int(round(progress_dt * 1e9))

    total_work = length(inputs) * length(states)
    prog = verbose ? ProgressMeter.Progress(total_work) : nothing # absolute-progress bar
    global_done = Threads.Atomic{Int}(0) # global counter (updated rarely)
    last_t = time_ns() # only used by thread 1
    # ---------------------------------

    nthreads = Threads.nthreads()
    transitions_by_thread = [Vector{Tuple{Int, Int, Int}}() for _ in 1:nthreads]
    local_done = fill(0, nthreads)

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

        local_done[tid] += 1
        if local_done[tid] ≥ update_interval
            Threads.atomic_add!(global_done, local_done[tid])
            local_done[tid] = 0
        end

        # UI refresh only from thread 1 and time-limited
        if verbose && tid == 1
            t = time_ns()
            if (t - last_t) ≥ progress_dt_ns
                ProgressMeter.update!(prog, global_done[] + local_done[1])
                last_t = t
            end
        end
    end

    # flush leftovers
    for tid in 1:nthreads
        v = local_done[tid]
        if v > 0
            Threads.atomic_add!(global_done, v)
            local_done[tid] = 0
        end
    end

    if verbose
        ProgressMeter.update!(prog, global_done[])
        ProgressMeter.finish!(prog)
    end

    for (_, local_transitions) in enumerate(transitions_by_thread)
        if !isempty(local_transitions)
            add_transitions!(abstract_system, local_transitions)
        end
    end

    return
end

function compute_abstract_system_from_concrete_system!(
    abstract_system::GridBasedSymbolicModel,
    concrete_system_approx::ST.DiscreteTimeLinearized;
    verbose = false,
    update_interval = Int(1e5),
    progress_dt = 0.2,
    threaded::Bool = false,
)
    # If multithreading is not requested or only one thread is available -> sequential execution
    if !threaded || Threads.nthreads() == 1
        XMapping = get_state_mapping(abstract_system)
        N = MP.get_dim(XMapping)
        r = MP.get_h(MP.get_grid(XMapping)) / 2.0
        _H_ = SMatrix{N, N}(LA.I) .* r
        _ONE_ = ones(SVector{N})
        e = LA.norm(r, Inf)
        translist = Tuple{Int, Int, Int}[]
        error_map = concrete_system_approx.error_map
        linsys_map = concrete_system_approx.linsys_map
        total_iterations = max(
            div(
                get_n_input(abstract_system) * get_n_state(abstract_system),
                update_interval,
            ),
            1,
        )
        progress = verbose ? ProgressMeter.Progress(total_iterations) : nothing
        count = 0
        for abstract_input in enum_inputs(abstract_system)
            concrete_input = get_concrete_input(abstract_system, abstract_input)
            Fe = error_map(e, concrete_input)
            Fr = r .+ Fe
            for abstract_state in enum_states(abstract_system)
                concrete_state = get_concrete_state(abstract_system, abstract_state)
                Fx, DFx = linsys_map(concrete_state, _H_, concrete_input)
                A = inv(DFx)
                b = abs.(A) * Fr .+ 1.0
                HP = UT.CenteredPolyhedron(A, b)
                # TODO: can we improve abs.(DFx)*_ONE_?
                rad = abs.(DFx) * _ONE_ .+ Fe
                reachable_set = UT.HyperRectangle(Fx - rad, Fx + rad)

                Base.empty!(translist)

                allin = compute_abstract_transitions_from_rectangle!(
                    abstract_system,
                    reachable_set,
                    abstract_state,
                    abstract_input,
                    translist,
                )
                allin && add_transitions!(abstract_system, translist)
                count += 1
                verbose && count % update_interval == 0 && ProgressMeter.next!(progress)
            end
        end
        verbose && ProgressMeter.finish!(progress)
        return
    end

    # ---- Multithreaded implementation ----
    XMapping = get_state_mapping(abstract_system)
    N = MP.get_dim(XMapping)
    r = MP.get_h(MP.get_grid(XMapping)) / 2.0
    _H_ = SMatrix{N, N}(LA.I) .* r
    _ONE_ = ones(SVector{N})
    e = LA.norm(r, Inf)
    error_map = concrete_system_approx.error_map
    linsys_map = concrete_system_approx.linsys_map

    inputs = collect(enum_inputs(abstract_system))
    states = collect(enum_states(abstract_system))

    input_data = Dict{Int, Tuple{Any, Any, Any}}()
    for abstract_input in inputs
        concrete_input = get_concrete_input(abstract_system, abstract_input)
        Fe = error_map(e, concrete_input)
        Fr = r .+ Fe
        input_data[abstract_input] = (concrete_input, Fe, Fr)
    end

    nthreads = Threads.nthreads()

    total_work = length(inputs) * length(states)

    transitions_by_thread = [Vector{Tuple{Int, Int, Int}}() for _ in 1:nthreads]

    progress =
        verbose ? ProgressMeter.Progress(total_work ÷ max(1, update_interval)) : nothing
    progress_count = Threads.Atomic{Int}(0)

    Threads.@threads for linear_idx in 1:total_work
        tid = Threads.threadid()
        local_transitions = transitions_by_thread[tid]

        input_idx = ((linear_idx - 1) ÷ length(states)) + 1
        state_idx = ((linear_idx - 1) % length(states)) + 1

        @inbounds begin
            abstract_input = inputs[input_idx]
            abstract_state = states[state_idx]

            concrete_input, Fe, Fr = input_data[abstract_input]
            concrete_state = get_concrete_state(abstract_system, abstract_state)

            Fx, DFx = linsys_map(concrete_state, _H_, concrete_input)
            A = inv(DFx)
            b = abs.(A) * Fr .+ 1.0
            HP = UT.CenteredPolyhedron(A, b)
            # TODO: can we improve abs.(DFx)*_ONE_?
            rad = abs.(DFx) * _ONE_ .+ Fe
            reachable_set = UT.HyperRectangle(Fx - rad, Fx + rad)

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
            if count_val % max(1, update_interval) == 0
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

function compute_abstract_system_from_concrete_system!(
    abstract_system::GridBasedSymbolicModel,
    concrete_system_approx::ST.DiscreteTimeSystemUnderApproximation;
    verbose = false,
    update_interval = Int(1e5),
    progress_dt = 0.2,
    threaded::Bool = false,
)
    # If multithreading is not requested or only one thread is available -> sequential execution
    if !threaded || Threads.nthreads() == 1
        translist = Tuple{Int, Int, Int}[]
        under_approximation_map = ST.get_under_approximation_map(concrete_system_approx)
        total_iterations = max(
            div(
                get_n_input(abstract_system) * get_n_state(abstract_system),
                update_interval,
            ),
            1,
        )
        progress = verbose ? ProgressMeter.Progress(total_iterations) : nothing
        count = 0
        for abstract_input in enum_inputs(abstract_system)
            concrete_input = get_concrete_input(abstract_system, abstract_input)
            for abstract_state in enum_states(abstract_system)
                concrete_elem = get_concrete_elem(abstract_system, abstract_state)
                reachable_points = under_approximation_map(concrete_elem, concrete_input)
                empty!(translist)
                allin = compute_abstract_transitions_from_points!(
                    abstract_system,
                    reachable_points,
                    abstract_state,
                    abstract_input,
                    translist,
                )
                allin && add_transitions!(abstract_system, translist)
                count += 1
                verbose && count % update_interval == 0 && ProgressMeter.next!(progress)
            end
        end
        verbose && ProgressMeter.finish!(progress)
        return
    end

    # ---- Multithreaded implementation ----
    under_approximation_map = ST.get_under_approximation_map(concrete_system_approx)
    inputs = collect(enum_inputs(abstract_system))
    states = collect(enum_states(abstract_system))

    nthreads = Threads.nthreads()

    total_work = length(inputs) * length(states)

    transitions_by_thread = [Vector{Tuple{Int, Int, Int}}() for _ in 1:nthreads]

    progress =
        verbose ? ProgressMeter.Progress(total_work ÷ max(1, update_interval)) : nothing
    progress_count = Threads.Atomic{Int}(0)

    Threads.@threads for linear_idx in 1:total_work
        tid = Threads.threadid()
        local_transitions = transitions_by_thread[tid]

        input_idx = ((linear_idx - 1) ÷ length(states)) + 1
        state_idx = ((linear_idx - 1) % length(states)) + 1

        @inbounds begin
            abstract_input = inputs[input_idx]
            abstract_state = states[state_idx]

            concrete_input = get_concrete_input(abstract_system, abstract_input)
            concrete_elem = get_concrete_elem(abstract_system, abstract_state)

            reachable_points = under_approximation_map(concrete_elem, concrete_input)

            prev_length = length(local_transitions)
            allin = compute_abstract_transitions_from_points!(
                abstract_system,
                reachable_points,
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
            if count_val % max(1, update_interval) == 0
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

function compute_abstract_system_from_concrete_system!(
    abstract_system::GridBasedSymbolicModel,
    concrete_system_approx::ST.DiscreteTimeCenteredSimulation;
    verbose = false,
    update_interval = Int(1e5), # flush local by thread
    progress_dt = 0.2, # seconds between UI refresh
    threaded::Bool = false,
)
    # If multithreading is not requested or only one thread is available -> sequential execution
    if !threaded || Threads.nthreads() == 1
        system_map = ST.get_system_map(concrete_system_approx)
        total_iterations = max(
            div(
                get_n_input(abstract_system) * get_n_state(abstract_system),
                update_interval,
            ),
            1,
        )
        progress = verbose ? ProgressMeter.Progress(total_iterations) : nothing
        count = 0
        for abstract_input in enum_inputs(abstract_system)
            concrete_input = get_concrete_input(abstract_system, abstract_input)
            for abstract_state in enum_states(abstract_system)
                concrete_state = get_concrete_state(abstract_system, abstract_state)
                y = system_map(concrete_state, concrete_input)
                target = get_abstract_state(abstract_system, y)
                if is_allowed_state(abstract_system, target)
                    add_transitions!(
                        abstract_system,
                        ((target, abstract_state, abstract_input),),
                    )
                end
                count += 1
                verbose && count % update_interval == 0 && ProgressMeter.next!(progress)
            end
        end
        verbose && ProgressMeter.finish!(progress)
        return
    end

    # ---- Multithreaded implementation ----

    system_map = ST.get_system_map(concrete_system_approx)
    inputs = collect(enum_inputs(abstract_system))
    states = collect(enum_states(abstract_system))

    # ---- Progress bar parameters ----
    progress_dt_ns = Int(round(progress_dt * 1e9))

    total_work = length(inputs) * length(states)
    prog = verbose ? ProgressMeter.Progress(total_work) : nothing # absolute-progress bar
    global_done = Threads.Atomic{Int}(0) # global counter (updated rarely)
    last_t = time_ns() # only used by thread 1
    # ---------------------------------

    nthreads = Threads.nthreads()
    transitions_by_thread = [Vector{Tuple{Int, Int, Int}}() for _ in 1:nthreads]
    local_done = fill(0, nthreads)

    Threads.@threads for linear_idx in 1:total_work
        tid = Threads.threadid()
        local_transitions = transitions_by_thread[tid]

        input_idx = ((linear_idx - 1) ÷ length(states)) + 1
        state_idx = ((linear_idx - 1) % length(states)) + 1

        @inbounds begin
            abstract_input = inputs[input_idx]
            abstract_state = states[state_idx]

            concrete_input = get_concrete_input(abstract_system, abstract_input)
            concrete_state = get_concrete_state(abstract_system, abstract_state)

            y = system_map(concrete_state, concrete_input)
            target = get_abstract_state(abstract_system, y)
            if is_allowed_state(abstract_system, target)
                push!(local_transitions, (target, abstract_state, abstract_input))
            end
        end
        local_done[tid] += 1
        if local_done[tid] ≥ update_interval
            Threads.atomic_add!(global_done, local_done[tid])
            local_done[tid] = 0
        end

        # UI refresh only from thread 1 and time-limited
        if verbose && tid == 1
            t = time_ns()
            if (t - last_t) ≥ progress_dt_ns
                ProgressMeter.update!(prog, global_done[] + local_done[1])
                last_t = t
            end
        end
    end

    # flush leftovers
    for tid in 1:nthreads
        v = local_done[tid]
        if v > 0
            Threads.atomic_add!(global_done, v)
            local_done[tid] = 0
        end
    end

    if verbose
        ProgressMeter.update!(prog, global_done[])
        ProgressMeter.finish!(prog)
    end

    total_found_transitions = 0
    for (_, local_transitions) in enumerate(transitions_by_thread)
        if !isempty(local_transitions)
            add_transitions!(abstract_system, local_transitions)
            total_found_transitions += length(local_transitions)
        end
    end

    return
end




# get_state_grid(symmodel::GridBasedSymbolicModel) = DO.get_grid(get_state_domain(symmodel))

# function get_concrete_elem(symmodel::GridBasedSymbolicModel, state)
#     xpos = get_xpos_by_state(symmodel, state)
#     return DO.get_elem_by_pos(get_state_domain(symmodel), xpos)
# end

# @recipe function f(
#     symmodel::GridBasedSymbolicModel;
#     arrowsB = false,
#     dims = [1, 2],
#     value_function = [], # Should be a function value_function(state::Int)::Float64 or left as []
#     colormap_name = "Blues",
#     default_color = :yellow,
# )
#     # Display the cells
#     state_grid = get_state_grid(symmodel)

#     if isa(value_function, Function)
#         # 1. Extract needed parts
#         projection_map = Dict{Tuple{Int, Int}, Tuple{Float64, Any}}() # value + elem
#         for state in enum_states(symmodel)
#             v = value_function(state)
#             elem = get_concrete_elem(symmodel, state)
#             pos = get_xpos_by_state(symmodel, state)
#             x1x2 = pos[dims]
#             # Only store if it's better (or first time)
#             if haskey(projection_map, x1x2)
#                 if v < projection_map[x1x2][1]
#                     projection_map[x1x2] = (v, elem)
#                 end
#             else
#                 projection_map[x1x2] = (v, elem)
#             end
#         end
#         # 2. Determine maximum finite value (for color scaling)
#         finite_vals = filter(isfinite, getindex.(values(projection_map), 1))
#         ValueMax = isempty(finite_vals) ? 1.0 : maximum(finite_vals)
#         # 3. Setup colormap
#         cmap = Colors.colormap(colormap_name)
#         mycolorMap = UT.Colormap([0.0, ValueMax], cmap)
#         # 4. Order states by decreasing value (for proper layering in case of overlapping cells)
#         cost_ordered = sort(collect(projection_map); by = x -> -x[2][1])
#         # 5. Plot
#         @series begin
#             for (_, (value, elem)) in cost_ordered
#                 color = isfinite(value) ? UT.get_color(mycolorMap, value) : default_color

#                 @series begin
#                     color := color
#                     dims := dims
#                     label := ""
#                     return elem
#                 end
#             end
#             mycolorMap
#         end
#     else
#         @series begin
#             dims := dims
#             get_state_domain(symmodel)
#         end
#     end
#     # Display the arrows
#     if arrowsB
#         for t in enum_transitions(symmodel)
#             color = RGB(
#                 abs(0.6 * sin(t[1])),
#                 abs(0.6 * sin(t[1] + 2π / 3)),
#                 abs(0.6 * sin(t[1] - 2π / 3)),
#             )
#             p1 = DO.get_coord_by_pos(state_grid, get_xpos_by_state(symmodel, t[2]))
#             p2 = DO.get_coord_by_pos(state_grid, get_xpos_by_state(symmodel, t[1]))

#             @series begin
#                 color := color
#                 return t[1] == t[2] ? UT.DrawPoint(p1) : UT.DrawArrow(p1, p2)
#             end
#         end
#     end
# end