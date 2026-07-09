# ============================================================
# Safety Control
# ============================================================

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

# RawOptimizerAttribute get/set + SolveTimeSec provided by AbstractDionysosOptimizer.

function compute_largest_invariant_set(autom::ST.AbstractAutomatonList, safelist)
    nstates = ST.get_n_state(autom)
    nsymbols = ST.get_n_input(autom)

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
    nextunsafeset = falses(nstates)

    while true
        fill!(nextunsafeset, false)

        @inbounds for target in 1:nstates
            unsafeset[target] || continue

            for (source, symbol) in ST.pre(autom, target)
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

    optimizer.success = all(q -> q in inv_set, problem.initial_set)

    optimizer.print_level >= 1 && println("\n Safety: terminated with $(optimizer.success)")

    optimizer.solve_time_sec = time() - t0
    return
end
