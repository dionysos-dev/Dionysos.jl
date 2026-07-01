mutable struct OptimizerReachAndStayProblem{T} <: MOI.AbstractOptimizer
    # inputs
    problem::Union{Nothing, PR.ReachAndStayProblem}
    print_level::Int

    # outputs
    controller::Union{Nothing, ST.AbstractDiscreteController}
    winning_set::Any
    winning_set_complement::Any
    success::Bool
    solve_time_sec::T

    function OptimizerReachAndStayProblem{T}() where {T}
        return new{T}(
            nothing, # problem
            1,       # print_level
            nothing, # controller
            nothing, # winning_set
            nothing, # winning_set_complement
            false,   # success
            zero(T), # solve_time_sec
        )
    end
end

OptimizerReachAndStayProblem() = OptimizerReachAndStayProblem{Float64}()


MOI.is_empty(optimizer::OptimizerReachAndStayProblem) = optimizer.problem === nothing

function MOI.set(model::OptimizerReachAndStayProblem, param::MOI.RawOptimizerAttribute, value)
    return setproperty!(model, Symbol(param.name), value)
end

function MOI.get(model::OptimizerReachAndStayProblem, ::MOI.SolveTimeSec)
    return model.solve_time_sec
end

function MOI.get(model::OptimizerReachAndStayProblem, param::MOI.RawOptimizerAttribute)
    return getproperty(model, Symbol(param.name))
end


function _compute_seed(autom, Y::Set{Int}, S_cells::Set{Int})
    seed = Set{Int}()
    for source in S_cells
        for symbol in ST.enum_inputs(autom)
            destinations = ST.post(autom, source, symbol)
            if !isempty(destinations) && all(d -> d in Y, destinations)
                push!(seed, source)
                break
            end
        end
    end
    return seed
end

function _invariant_with_floor(autom, safe_cells, floor_cells)
    floor_set = Set(floor_cells)
    
    contr_tab = DiscreteControlTable(ST.get_n_state(autom))
    nstates = ST.get_n_state(autom)
    nsymbols = ST.get_n_input(autom)
    pairstable = falses(nstates, nsymbols)

    _compute_pairstable(pairstable, autom)
    nsymbolslist = vec(sum(pairstable; dims = 2))

    safeset = Set(safe_cells) ∪ floor_set

    for source in collect(safeset)
        if nsymbolslist[source] == 0 && !(source in floor_set)
            delete!(safeset, source)
        end
    end

    unsafeset = Set(1:nstates)
    setdiff!(unsafeset, safeset)

    for source in unsafeset
        for symbol in 1:nsymbols
            pairstable[source, symbol] = false
        end
    end

    nextunsafeset = Set{Int}()

    while true
        for target in unsafeset
            for (source, symbol) in ST.pre(autom, target)
                if pairstable[source, symbol]
                    pairstable[source, symbol] = false
                    nsymbolslist[source] -= 1
                    if nsymbolslist[source] == 0 && !(source in floor_set)
                        push!(nextunsafeset, source)
                    end
                end
            end
        end

        isempty(nextunsafeset) && break

        setdiff!(safeset, nextunsafeset)
        unsafeset, nextunsafeset = nextunsafeset, unsafeset
        empty!(nextunsafeset)
    end

    for source in safeset
        for symbol in 1:nsymbols
            if pairstable[source, symbol]
                add_control!(contr_tab, source, symbol)
            end
        end
    end

    controller = ST.DiscreteStaticController(safeset, contr_tab, false)
    return controller, safeset
end

# Build a reach-and-stay CONTROLLER from the winning set.  Rank-0 is the invariant
# core of T (record a STAY input there); every other winning cell records an input
# whose successors are ALL strictly closer to the core (already ranked).  Hence the
# controller actively DRIVES UP into the core and then HOLDS — not merely "stay safe".
function _build_reach_stay_controller(autom, winning, T_cells)
    nstates  = ST.get_n_state(autom)
    nsymbols = ST.get_n_input(autom)
    contr    = DiscreteControlTable(nstates)
    Wset     = winning isa Set ? winning : Set(winning)

    # rank-0: invariant core of T (the balance set), with a stay input for each cell.
    _, core = _invariant_with_floor(autom, collect(T_cells), Int[])
    coreset = core isa Set ? core : Set(core)
    ranked  = falses(nstates)
    queue   = Int[]
    for c in coreset
        ranked[c] = true
        push!(queue, c)
        for u in 1:nsymbols
            dests = ST.post(autom, c, u)
            if !isempty(dests) && all(d -> d in coreset, dests)
                set_control!(contr, c, u)
                break
            end
        end
    end

    # count not-yet-ranked successors for each (winning cell, input)
    remaining = Dict{Tuple{Int, Int}, Int}()
    for c in Wset
        ranked[c] && continue
        for u in 1:nsymbols
            dests = ST.post(autom, c, u)
            isempty(dests) || (remaining[(c, u)] = length(dests))
        end
    end

    # backward BFS: a winning cell becomes rankable once ALL successors of some input
    # are ranked (so that input drives it strictly toward the core).
    while !isempty(queue)
        t = popfirst!(queue)
        for (src, u) in ST.pre(autom, t)
            key = (src, u)
            haskey(remaining, key) || continue
            remaining[key] -= 1
            if remaining[key] == 0 && !ranked[src] && (src in Wset)
                ranked[src] = true
                set_control!(contr, src, u)
                push!(queue, src)
            end
        end
    end

    return ST.DiscreteStaticController(Set(findall(ranked)), contr, false)
end

function MOI.optimize!(optimizer::OptimizerReachAndStayProblem)
    t0 = time()

    problem = optimizer.problem
    problem === nothing && error("problem not set")

    autom = problem.system

    T_cells = Set(problem.target_set)
    S_cells = Set(problem.safe_set)
    I_cells = collect(problem.initial_set)

    optimizer.print_level >= 1 && println("compute_reach_and_stay! started")

    # Outer μY loop
    Y = Set{Int}()

    while true
        Y_old = copy(Y)

        # Compute seed = CPre(Y) ∩ S
        seed = _compute_seed(autom, Y, S_cells)

        # Inner νX: largest invariant in T ∪ seed with seed as floor
        _, X_star = _invariant_with_floor(autom, collect(T_cells), collect(seed))

        Y = X_star

        Y == Y_old && break
    end

    # Extract a controller that DRIVES toward the core then holds (not just "stay safe").
    optimizer.controller = _build_reach_stay_controller(autom, Y, T_cells)
    optimizer.winning_set = Y
    optimizer.winning_set_complement = setdiff(Set(1:ST.get_n_state(autom)), Y)
    optimizer.success = all(q -> q in Y, I_cells)

    optimizer.print_level >= 1 &&
        println("\n Reach and Stay: terminated with $(optimizer.success)")

    optimizer.solve_time_sec = time() - t0
    return
end