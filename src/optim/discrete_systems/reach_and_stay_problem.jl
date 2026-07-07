# ============================================================
# Reach-and-Stay Problem
# ============================================================

mutable struct OptimizerReachAndStayProblem{T} <: MOI.AbstractOptimizer
    problem::Union{Nothing, PR.ReachAndStayProblem}
    early_stop::Bool
    print_level::Int

    controller::Union{Nothing, ST.AbstractDiscreteController}
    winning_set::Any
    winning_set_complement::Any
    success::Bool
    solve_time_sec::T

    function OptimizerReachAndStayProblem{T}() where {T}
        return new{T}(nothing, false, 1, nothing, nothing, nothing, false, zero(T))
    end
end

OptimizerReachAndStayProblem() = OptimizerReachAndStayProblem{Float64}()

MOI.is_empty(optimizer::OptimizerReachAndStayProblem) = optimizer.problem === nothing

function MOI.set(
    model::OptimizerReachAndStayProblem,
    param::MOI.RawOptimizerAttribute,
    value,
)
    return setproperty!(model, Symbol(param.name), value)
end

function MOI.get(model::OptimizerReachAndStayProblem, ::MOI.SolveTimeSec)
    return model.solve_time_sec
end

function MOI.get(model::OptimizerReachAndStayProblem, param::MOI.RawOptimizerAttribute)
    return getproperty(model, Symbol(param.name))
end

# ------------------------------------------------------------
# Controlled predecessor:
# Pre(Y | S) = {q in S | exists u, Post(q,u) subset Y}
# ------------------------------------------------------------

function _compute_seed(autom, Y::Set{Int}, S::Set{Int})
    seed = Set{Int}()
    seen = Set{Tuple{Int, Int}}()

    for target in Y
        for (source, symbol) in ST.pre(autom, target)
            source in S || continue

            key = (source, symbol)
            key in seen && continue
            push!(seen, key)

            destinations = ST.post(autom, source, symbol)

            if !isempty(destinations) && all(d -> d in Y, destinations)
                push!(seed, source)
            end
        end
    end

    return seed
end

# ------------------------------------------------------------
# Maximal controlled invariant subset of safe_cells ∪ floor_cells
# with floor_cells protected from removal.
# ------------------------------------------------------------

function _invariant_with_floor(autom, safe_cells, floor_cells)
    nstates = ST.get_n_state(autom)
    nsymbols = ST.get_n_input(autom)

    floor_set = Set(floor_cells)
    safeset = Set(safe_cells) ∪ floor_set

    pairstable = falses(nstates, nsymbols)
    _compute_pairstable(pairstable, autom)

    outside = setdiff(Set(1:nstates), safeset)

    for source in outside
        for symbol in 1:nsymbols
            pairstable[source, symbol] = false
        end
    end

    nsymbolslist = vec(sum(pairstable; dims = 2))

    unsafeset = outside
    nextunsafeset = Set{Int}()

    while true
        for target in unsafeset
            for (source, symbol) in ST.pre(autom, target)
                if pairstable[source, symbol]
                    pairstable[source, symbol] = false
                    nsymbolslist[source] -= 1

                    if source in safeset &&
                       !(source in floor_set) &&
                       nsymbolslist[source] == 0
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

    return safeset, pairstable
end

# ------------------------------------------------------------
# Add one admissible control q -> target_set.
# ------------------------------------------------------------

function _add_one_control_to_target!(contr, autom, cells, target_set::Set{Int})
    for q in cells
        isempty(contr(q)) || continue

        for u in ST.enum_inputs(autom)
            destinations = ST.post(autom, q, u)

            if !isempty(destinations) && all(d -> d in target_set, destinations)
                add_control!(contr, q, u)
                break
            end
        end
    end

    return contr
end

# ------------------------------------------------------------
# Main solver: Algorithm (3) from the paper.
# ------------------------------------------------------------

function MOI.optimize!(optimizer::OptimizerReachAndStayProblem)
    t0 = time()

    problem = optimizer.problem
    problem === nothing && error("problem not set")

    autom = problem.system

    T = Set(problem.target_set)
    S = Set(problem.safe_set)
    I = collect(problem.initial_set)

    nstates = ST.get_n_state(autom)

    optimizer.print_level >= 1 && println("compute_reach_and_stay! started")

    contr = DiscreteControlTable(nstates)

    Y = Set{Int}()
    X_prev = Set{Int}()

    while true
        Y_old = copy(Y)

        # Inner ν fixed point:
        # X∞_{i+1} = maximal controlled invariant inside Ω ∪ Y_i
        X_star, _ = _invariant_with_floor(autom, T, Y)

        # Controls for states in Ω ∩ (X∞_{i+1} \ X∞_i)
        new_stay_cells = intersect(T, setdiff(X_star, X_prev))
        _add_one_control_to_target!(contr, autom, new_stay_cells, X_star)

        # Outer μ fixed point:
        # Y_{i+1} = Pre(X∞_{i+1})
        Y_next = _compute_seed(autom, X_star, S)

        # Controls for states in Y_{i+1} \ (Y_i ∪ Ω)
        new_reach_cells = setdiff(Y_next, union(Y, T))
        _add_one_control_to_target!(contr, autom, new_reach_cells, X_star)

        Y = Y_next
        X_prev = X_star

        if optimizer.early_stop && all(q -> q in Y, I)
            break
        end

        Y == Y_old && break
    end

    optimizer.controller = ST.DiscreteStaticController(Y, contr, false)
    optimizer.winning_set = Y
    optimizer.winning_set_complement = setdiff(Set(ST.enum_states(autom)), Y)
    optimizer.success = all(q -> q in Y, I)

    optimizer.print_level >= 1 &&
        println("\n Reach and Stay: terminated with $(optimizer.success)")

    optimizer.solve_time_sec = time() - t0
    return
end
