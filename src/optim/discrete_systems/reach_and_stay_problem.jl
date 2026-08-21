# ============================================================
# Reach-and-Stay Problem
# ============================================================

"""
    OptimizerReachAndStayProblem{T} <: AbstractDionysosOptimizer

Reach-and-stay synthesis on a finite automaton: given a
[`ReachAndStayProblem`](@ref Dionysos.Problem.ReachAndStayProblem) whose system is an
`AbstractAutomatonList`, compute the winning set and a **memoryless** controller.

The problem's `stay_on_first_entry` picks the algorithm. `false` is ◇□ in the literal sense and
runs the nested μ/ν fixed point of Li & Liu; `true` forbids leaving the target once entered, and
reduces to one invariance solve followed by one reachability solve.

Set `"problem"`, optionally `"early_stop"`; read back `"controller"`, `"winning_set"` and
`"winning_set_complement"`.
"""
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

# ------------------------------------------------------------
# Backward attractor: grow `won` with every state that can be forced into it in finitely many
# steps without leaving `safe`.
#
# `counter[q, u]` holds how many successors of the pair are not yet won, so a pair becomes
# usable exactly when it reaches zero and every transition is touched once over the whole
# computation: `O(E + nm)`. This replaces a per-layer `Pre(Y | S)` sweep, which re-walked all
# `E` transitions for every layer — on a domain whose reachability runs hundreds of layers
# deep, `O(depth ⋅ E)`. It is the device the reach-avoid solver already uses.
#
# A pair with no transition at all keeps its initial count of zero and is never decremented, so
# it is never enabled — the `Post(q,u) ≠ ∅` half of the controlled predecessor comes for free.
#
# Every newly won cell takes the input that won it. That is the memoryless choice Li & Liu
# require — fixed at the layer the cell first enters the winning set, never revised — and it is
# free, the enabling pair being known already. `no_control` names cells that may be won but
# must not be given one here, because they are the caller's to assign.
#
# `seen` is what guards the counters, and it is deliberately not `won`: a counter may only be
# decremented once per transition, ever. A cell can be won here in one round and then join the
# invariant core in the next, and relaxing it a second time would drop some predecessor's count
# to zero while a successor of that pair was still unwon — a strictly larger, unsound winning
# set. Callers keep `seen` across rounds and seed only cells it does not already hold.
#
# `frontier` is consumed. `stop` is polled once per layer, for `early_stop`.
# ------------------------------------------------------------

function _attract!(
    contr,
    autom,
    won::BitVector,
    seen::BitVector,
    frontier::Vector{Int},
    safe::BitVector,
    counter,
    no_control::BitVector,
    stop,
)
    next = Int[]

    while !isempty(frontier)
        empty!(next)

        for target in frontier
            for (source, symbol) in SY.pre(autom, target)
                seen[source] && continue
                safe[source] || continue

                if decrease_counter!(counter, source, symbol) == 0
                    seen[source] = true
                    won[source] = true
                    push!(next, source)
                    no_control[source] || add_control!(contr, source, symbol)
                end
            end
        end

        frontier, next = next, frontier
        stop !== nothing && stop() && break
    end

    return won
end

# ------------------------------------------------------------
# Maximal controlled invariant subset of safe_cells ∪ floor_cells.
#
# `floor_cells` are protected from removal. In Algorithm 3 of
# Li & Liu, "Robustly Complete Synthesis of Memoryless Controllers for
# Nonlinear Systems with Reach-and-Stay Specifications" (IEEE TAC, 2020),
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

    # Same dead-state seeding as `compute_largest_invariant_set`. Floor cells stay protected:
    # a pair with no transition is never enabled by `_attract!` either, so a floor cell — which
    # had to be won through one — always has a pair.
    @inbounds for q in 1:nstates
        if safeset[q] && !floor_cells[q] && nsymbolslist[q] == 0
            safeset[q] = false
            unsafeset[q] = true
        end
    end

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
# ◇□ (the default): Algorithm 3 of Li & Liu, "Robustly Complete Synthesis of Memoryless
# Controllers for Nonlinear Systems with Reach-and-Stay Specifications" (IEEE TAC, 2020).
#
# The nested μ/ν fixed point. The inner ν recomputes the invariant inside `Ω ∪ Y_i` — the
# target *together with* what is already winning — so a target cell may be certified by
# leaving the target into won territory and coming back. Finitely many such departures are
# allowed, which is exactly what ◇□ means and why this winning set is the larger one.
#
# The controller comes out *memoryless* — the paper's title — because of when inputs are
# chosen, not because of the fixed point itself. Each cell is given one input at the moment it
# is first won and never revised: `_add_one_control_to_target!` for a target cell joining the
# invariant core, `_attract!` for a cell reached on the way in.
# ------------------------------------------------------------

function _solve_eventually_always(autom, T, S, I, early_stop)
    nstates = SY.get_n_state(autom)

    base_pairstable = _compute_base_pairstable(autom)
    counter = _counter_dense(autom)

    contr = _empty_control_table(nstates)

    # `Y` accumulates rather than being rebuilt per round: `X∞` grows monotonically, so `Pre`
    # of it does too, and the attractor below only ever adds. Keeping one `counter` across
    # rounds is what makes the whole outer loop cost `O(E + nm)` in total instead of per round.
    Y = falses(nstates)
    X_prev = falses(nstates)
    seen = falses(nstates)

    while true
        # Inner ν fixed point:
        # X∞_{i+1} = maximal controlled invariant inside Ω ∪ Y_i
        X_star, _ = _invariant_with_floor(autom, T, Y, base_pairstable)

        # Controls for Ω ∩ (X∞_{i+1} \ X∞_i)
        new_stay_cells = X_star .& T .& .!X_prev
        _add_one_control_to_target!(contr, autom, new_stay_cells, X_star)

        # Outer μ fixed point. Taking the whole attractor of `X∞_{i+1}` rather than a single
        # `Pre` layer reaches the same limit — at the fixed point `Y ⊆ X∞`, so one `Pre` layer
        # of `X∞` already contains the attractor — but in a handful of rounds instead of one
        # per reachability layer.
        # Seed with whatever `X∞` has gained, and only that — see `_attract!` on why a cell may
        # be relaxed at most once. A core cell inside the safe set is winning by definition; it
        # already carries an input, from `new_stay_cells` above if it is a target, from the
        # round that first reached it otherwise.
        won_before = count(Y)

        frontier = Int[]
        for q in 1:nstates
            (X_star[q] && !seen[q]) || continue
            seen[q] = true
            push!(frontier, q)
            S[q] && (Y[q] = true)
        end
        isempty(frontier) && break

        # Ω is the caller's to control: a target cell earns its input from `new_stay_cells`
        # above once it joins the invariant core, never a reach input from the attractor.
        _attract!(contr, autom, Y, seen, frontier, S, counter, T, nothing)

        X_prev = X_star

        if early_stop && all(q -> Y[q], I)
            break
        end

        count(Y) == won_before && break
    end

    return Y, contr
end

# ------------------------------------------------------------
# "Stay from the first entry": safety, then reachability.
#
# The invariant core `X∞` is the maximal controlled invariant subset of the target **alone**
# — no floor, so it is computed once and never widened by what becomes winning later. The run
# is then steered to `X∞` and stays there from the moment it arrives.
#
# Target cells outside `X∞` are removed from the region the reachability may use. Such a cell
# is *in* the target but can only continue by leaving it, so routing through it would break
# the very property this variant exists to enforce — and counting it as winning would promise
# something the controller cannot deliver.
# ------------------------------------------------------------

function _solve_stay_on_first_entry(autom, T, S, I, early_stop)
    nstates = SY.get_n_state(autom)

    base_pairstable = _compute_base_pairstable(autom)
    contr = _empty_control_table(nstates)

    X_inf, _ = _invariant_with_floor(autom, T, falses(nstates), base_pairstable)
    _add_one_control_to_target!(contr, autom, X_inf, X_inf)

    # Reachability may not pass through the target outside its invariant core.
    S_reach = S .& .!(T .& .!X_inf)

    # The core is fixed here — computed once, never widened — so the whole reachability is one
    # attractor rather than a fixed point that re-sweeps the automaton per layer. The core is
    # won from the start, which is also what keeps it from being handed a reach input over the
    # stay input it already has.
    won = copy(X_inf)
    stop = early_stop ? () -> all(q -> won[q], I) : nothing

    # One pass, so `seen` and `won` start out the same; the core is both already winning and
    # already relaxed.
    _attract!(
        contr,
        autom,
        won,
        copy(X_inf),
        findall(X_inf),
        S_reach,
        _counter_dense(autom),
        falses(nstates),
        stop,
    )

    # The core is winning too: it is where the run ends up and stays.
    return won, contr
end

function MOI.optimize!(optimizer::OptimizerReachAndStayProblem)
    t0 = time()

    problem = optimizer.problem
    problem === nothing && error("problem not set")

    autom = problem.system
    nstates = SY.get_n_state(autom)

    T = _bitset_from_states(problem.target_set, nstates)
    S = _bitset_from_states(problem.safe_set, nstates)
    I = collect(problem.initial_set)

    optimizer.print_level >= 1 && println("compute_reach_and_stay! started")

    solve =
        problem.stay_on_first_entry ? _solve_stay_on_first_entry : _solve_eventually_always
    Y, contr = solve(autom, T, S, I, optimizer.early_stop)

    winning_set = _set_from_bitset(Y)

    optimizer.controller = ST.DiscreteStaticController(winning_set, contr, false)
    optimizer.winning_set = winning_set
    optimizer.winning_set_complement = setdiff(Set(SY.enum_states(autom)), winning_set)
    optimizer.success = covers_initial_set(q -> Y[q], I)

    optimizer.print_level >= 1 &&
        println("\n Reach and Stay: terminated with $(optimizer.success)")

    optimizer.solve_time_sec = time() - t0
    return
end
