# ============================================================
# Safety Control
# ============================================================

mutable struct OptimizerSafetyProblem{T} <: MOI.AbstractOptimizer
    # inputs
    problem::Union{Nothing, PR.SafetyProblem}
    print_level::Int

    # outputs
    controller::Union{Nothing, MS.ConstrainedBlackBoxMap}
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

function MOI.set(model::OptimizerSafetyProblem, param::MOI.RawOptimizerAttribute, value)
    return setproperty!(model, Symbol(param.name), value)
end

function MOI.get(model::OptimizerSafetyProblem, ::MOI.SolveTimeSec)
    return model.solve_time_sec
end

function MOI.get(model::OptimizerSafetyProblem, param::MOI.RawOptimizerAttribute)
    return getproperty(model, Symbol(param.name))
end

function MOI.optimize!(optimizer::OptimizerSafetyProblem)
    t0 = time()

    problem = optimizer.problem
    problem === nothing && error("problem not set")

    autom = problem.system

    optimizer.print_level >= 1 && println("compute_controller_safe! started")

    controller, inv_set, invc_set =
        compute_largest_invariant_set(autom, optimizer.problem.safe_set)

    optimizer.controller = controller
    optimizer.invariant_set = inv_set
    optimizer.invariant_set_complement = invc_set

    # success check: initial_set ⊆ invariant_set
    optimizer.success = all(q -> q in inv_set, problem.initial_set)

    optimizer.print_level >= 1 && println("\n Safety: terminated with $(optimizer.success)")

    optimizer.solve_time_sec = time() - t0
    return
end

function compute_largest_invariant_set(autom::ST.AbstractAutomatonList, safelist;)
    contr_tab = SymbolicControlTable(ST.get_n_state(autom))
    nstates = ST.get_n_state(autom)
    nsymbols = ST.get_n_input(autom)
    pairstable = [false for i in 1:nstates, j in 1:nsymbols]

    _compute_pairstable(pairstable, autom)
    nsymbolslist = sum(pairstable; dims = 2)

    # Remove unsafe states
    safeset = Set(safelist)
    for source in safeset
        if nsymbolslist[source] == 0
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

    # Iterate until convergence
    while true
        for target in unsafeset
            for soursymb in ST.pre(autom, target)
                if pairstable[soursymb[1], soursymb[2]]
                    pairstable[soursymb[1], soursymb[2]] = false
                    nsymbolslist[soursymb[1]] -= 1
                    if nsymbolslist[soursymb[1]] == 0
                        push!(nextunsafeset, soursymb[1])
                    end
                end
            end
        end

        if isempty(nextunsafeset)
            break
        end

        setdiff!(safeset, nextunsafeset)
        unsafeset, nextunsafeset = nextunsafeset, unsafeset
        empty!(nextunsafeset)
    end

    # Populate controller
    for source in safeset
        for symbol in 1:nsymbols
            if pairstable[source, symbol]
                set_control!(contr_tab, source, symbol)
            end
        end
    end
    unsafeset = setdiff(Set(safelist), safeset)
    controller = to_ms_controller(contr_tab)
    return controller, safeset, unsafeset
end

function _compute_pairstable(pairstable, autom)
    for target in ST.enum_states(autom)
        for soursymb in ST.pre(autom, target)
            pairstable[soursymb[1], soursymb[2]] = true
        end
    end
end
