# ============================================================
# Safety Control
# ============================================================

"""
    OptimizerSafetyProblem{T} <: AbstractDionysosOptimizer

Safety synthesis on a finite automaton: given a [`SafetyProblem`](@ref Dionysos.Problem.SafetyProblem)
whose system is an `AbstractAutomatonList`, compute the **maximal controlled-invariant subset** of
the safe set.

Set `"problem"`; read back `"controller"`, `"invariant_set"` and `"invariant_set_complement"`.
The controller keeps *every* surviving input at each state, not one, so downstream code is free
to choose among them.
"""
mutable struct OptimizerSafetyProblem{T} <: AbstractDionysosOptimizer
    # inputs
    problem::Union{Nothing, PR.SafetyProblem}
    print_level::Int

    # outputs
    controller::Union{Nothing, ST.AbstractDiscreteController}
    invariant_set::Any
    invariant_set_complement::Any
    success::Bool
    solve_time_sec::T

    function OptimizerSafetyProblem{T}() where {T}
        return new{T}(
            nothing, # problem
            1,       # print_level
            nothing, # controller
            nothing, # invariant_set
            nothing, # invariant_set_complement
            false,   # success
            zero(T), # solve_time_sec
        )
    end
end

OptimizerSafetyProblem() = OptimizerSafetyProblem{Float64}()

MOI.is_empty(optimizer::OptimizerSafetyProblem) = optimizer.problem === nothing

function compute_largest_invariant_set(autom::SY.AbstractAutomatonList, safelist)
    nstates = SY.get_n_state(autom)
    nsymbols = SY.get_n_input(autom)

    safeset = _bitset_from_states(safelist, nstates)

    pairstable = _compute_base_pairstable(autom)

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

    # The loop below removes a state only when a pair being *disabled* drops its count to
    # zero, and a count that starts at zero never drops. A state every input drives out of the
    # domain would therefore stay invariant forever, predecessors included. Seed it here.
    @inbounds for q in 1:nstates
        if safeset[q] && nsymbolslist[q] == 0
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

                    if safeset[source] && nsymbolslist[source] == 0
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

    contr_tab = DiscreteControlTable(nstates)

    @inbounds for q in 1:nstates
        safeset[q] || continue

        for u in 1:nsymbols
            if pairstable[q, u]
                add_control!(contr_tab, q, u)
            end
        end
    end

    inv_set = _set_from_bitset(safeset)

    safe_bits_original = _bitset_from_states(safelist, nstates)
    invc_bits = safe_bits_original .& .!safeset
    invc_set = _set_from_bitset(invc_bits)

    controller = ST.DiscreteStaticController(inv_set, contr_tab, false)

    return controller, inv_set, invc_set
end

function MOI.optimize!(optimizer::OptimizerSafetyProblem)
    t0 = time()

    problem = optimizer.problem
    problem === nothing && error("problem not set")

    autom = problem.system

    optimizer.print_level >= 1 && println("compute_controller_safe! started")

    controller, inv_set, invc_set = compute_largest_invariant_set(autom, problem.safe_set)

    optimizer.controller = controller
    optimizer.invariant_set = inv_set
    optimizer.invariant_set_complement = invc_set

    optimizer.success = covers_initial_set(q -> q in inv_set, problem.initial_set)

    optimizer.print_level >= 1 && println("\n Safety: terminated with $(optimizer.success)")

    optimizer.solve_time_sec = time() - t0
    return
end
