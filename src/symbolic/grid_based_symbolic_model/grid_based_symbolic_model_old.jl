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
    targets, allin = get_states_from_set_strict(symmodel, reachable_set, MP.OUTER)
    allin || return false
    for target in targets
        if !is_allowed_state(symmodel, target)
            return false
        end
    end
    for target in targets
        push!(translist, (target, abstract_state, abstract_input))
    end

    return true
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
        if target === nothing || !is_allowed_state(symmodel, target)
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
                if target !== nothing && is_allowed_state(abstract_system, target)
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
