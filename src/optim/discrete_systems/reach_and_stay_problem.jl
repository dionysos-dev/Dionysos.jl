# ============================================================
# Reach-and-Stay Problem
# ============================================================

mutable struct OptimizerReachAndStayProblem{T} <: AbstractDionysosOptimizer
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

# RawOptimizerAttribute get/set + SolveTimeSec provided by AbstractDionysosOptimizer.

# ------------------------------------------------------------
# Controlled predecessor:
# Pre(Y | S) = {q in S | exists u, Post(q,u) subset Y}
# ------------------------------------------------------------

function _compute_seed(autom, Y::BitVector, S::BitVector, seen_pairs::BitMatrix)
    nstates = SY.get_n_state(autom)
    seed = falses(nstates)

    fill!(seen_pairs, false)

    for target in eachindex(Y)
        Y[target] || continue
        for (source, symbol) in SY.pre(autom, target)
            S[source] || continue
            seen_pairs[source, symbol] && continue

            seen_pairs[source, symbol] = true

            dests = SY.post(autom, source, symbol)

            if !isempty(dests) && all(d -> Y[d], dests)
                seed[source] = true
            end
        end
    end

    return seed
end

# ------------------------------------------------------------
# Maximal controlled invariant subset of safe_cells ∪ floor_cells.
#
# `floor_cells` are protected from removal. In Algorithm (3),
# this is used with:
#
#     safe_cells  = Ω
#     floor_cells = Y_i
#
# so the invariant is computed inside Ω ∪ Y_i.
# ------------------------------------------------------------

function _invariant_with_floor(
    autom,
    safe_cells::BitVector,
    floor_cells::BitVector,
    base_pairstable::BitMatrix,
)
    nstates = SY.get_n_state(autom)
    nsymbols = SY.get_n_input(autom)

    safeset = copy(safe_cells)
    safeset .|= floor_cells

    pairstable = copy(base_pairstable)

    unsafeset = falses(nstates)

    @inbounds for q in 1:nstates
        if !safeset[q]
            unsafeset[q] = true
            for u in 1:nsymbols
                pairstable[q, u] = false
            end
        end
    end

    nsymbolslist = _compute_nsymbolslist(pairstable)
    nextunsafeset = falses(nstates)

    while true
        fill!(nextunsafeset, false)

        @inbounds for target in 1:nstates
            unsafeset[target] || continue

            for (source, symbol) in SY.pre(autom, target)
                if pairstable[source, symbol]
                    pairstable[source, symbol] = false
                    nsymbolslist[source] -= 1

                    if safeset[source] && !floor_cells[source] && nsymbolslist[source] == 0
                        nextunsafeset[source] = true
                    end
                end
            end
        end

        has_next = false

        @inbounds for q in 1:nstates
            if nextunsafeset[q]
                safeset[q] = false
                has_next = true
            end
        end

        has_next || break

        unsafeset, nextunsafeset = nextunsafeset, unsafeset
    end

    return safeset, pairstable
end

# ------------------------------------------------------------
# Add one admissible control q -> target_set.
# ------------------------------------------------------------

function _add_one_control_to_target!(contr, autom, cells::BitVector, target_set::BitVector)
    nstates = SY.get_n_state(autom)

    @inbounds for q in 1:nstates
        cells[q] || continue
        isempty(contr(q)) || continue

        for u in SY.enum_inputs(autom)
            dests = SY.post(autom, q, u)

            if !isempty(dests) && all(d -> target_set[d], dests)
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
    nstates = SY.get_n_state(autom)
    nsymbols = SY.get_n_input(autom)

    T = _bitset_from_states(problem.target_set, nstates)
    S = _bitset_from_states(problem.safe_set, nstates)
    I = collect(problem.initial_set)

    optimizer.print_level >= 1 && println("compute_reach_and_stay! started")

    base_pairstable = _compute_base_pairstable(autom)
    seen_pairs = falses(nstates, nsymbols)

    contr = _empty_control_table(nstates)

    Y = falses(nstates)
    X_prev = falses(nstates)

    while true
        Y_old = copy(Y)

        # Inner ν fixed point:
        # X∞_{i+1} = maximal controlled invariant inside Ω ∪ Y_i
        X_star, _ = _invariant_with_floor(autom, T, Y, base_pairstable)

        # Controls for Ω ∩ (X∞_{i+1} \ X∞_i)
        new_stay_cells = X_star .& T .& .!X_prev
        _add_one_control_to_target!(contr, autom, new_stay_cells, X_star)

        # Outer μ fixed point:
        # Y_{i+1} = Pre(X∞_{i+1} | S)
        Y_next = _compute_seed(autom, X_star, S, seen_pairs)

        # Controls for Y_{i+1} \ (Y_i ∪ Ω)
        new_reach_cells = Y_next .& .!Y .& .!T
        _add_one_control_to_target!(contr, autom, new_reach_cells, X_star)

        Y = Y_next
        X_prev = X_star

        if optimizer.early_stop && all(q -> Y[q], I)
            break
        end

        Y == Y_old && break
    end

    winning_set = _set_from_bitset(Y)

    optimizer.controller = ST.DiscreteStaticController(winning_set, contr, false)
    optimizer.winning_set = winning_set
    optimizer.winning_set_complement = setdiff(Set(SY.enum_states(autom)), winning_set)
    optimizer.success = all(q -> Y[q], I)

    optimizer.print_level >= 1 &&
        println("\n Reach and Stay: terminated with $(optimizer.success)")

    optimizer.solve_time_sec = time() - t0
    return
end
