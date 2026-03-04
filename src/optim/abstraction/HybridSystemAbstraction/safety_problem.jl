mutable struct OptimizerSafetyProblem{T} <: MOI.AbstractOptimizer
    # Inputs
    concrete_problem::Union{Nothing, PR.SafetyProblem}
    abstract_system::Union{Nothing, SY.TimedHybridSymbolicModel}

    # Specific parameters
    optimizer_kwargs_dict::Union{Nothing, Any}

    print_level::Int

    # Outputs
    abstract_problem::Union{Nothing, PR.SafetyProblem}
    abstract_controller::Union{Nothing, MS.AbstractMap}
    abstract_problem_time_sec::T
    success::Bool


    function OptimizerSafetyProblem{T}() where {T}
        return new{T}(
            nothing,
            nothing,
            nothing,
            1,
            nothing,
            nothing,
            0.0,
            false,
        )
    end
end

OptimizerSafetyProblem() = OptimizerSafetyProblem{Float64}()

MOI.is_empty(optimizer::OptimizerSafetyProblem) =
    optimizer.concrete_problem === nothing

function MOI.set(
    model::OptimizerSafetyProblem,
    param::MOI.RawOptimizerAttribute,
    value,
)
    return setproperty!(model, Symbol(param.name), value)
end

function MOI.get(model::OptimizerSafetyProblem, ::MOI.SolveTimeSec)
    return model.abstract_problem_time_sec
end

function MOI.get(model::OptimizerSafetyProblem, param::MOI.RawOptimizerAttribute)
    return getproperty(model, Symbol(param.name))
end

function reset!(model::OptimizerSafetyProblem)
    model.abstract_problem = nothing
    model.abstract_system = nothing
    model.abstract_controller = nothing
    model.abstract_problem_time_sec = 0.0
    model.success = false
    return model
end

function MOI.optimize!(optimizer::OptimizerSafetyProblem)
    t_ref = time()

    optimizer.abstract_system === nothing &&
        error("Abstract system is not defined. Ensure abstraction is computed first.")
    optimizer.concrete_problem === nothing &&
        error("Concrete problem is not defined.")

    abstract_system = optimizer.abstract_system

    # Build abstract problem
    abstract_problem = build_abstract_problem(optimizer.concrete_problem, abstract_system)
    optimizer.abstract_problem = abstract_problem

    optimizer.print_level >= 1 && println("compute_controller_reachability! started")

    abstract_controller, invariant_set_symbols, _ = compute_largest_invariant_set_timed_hybrid(
        abstract_problem.system.symbolic_automaton,
        collect(abstract_problem.safe_set),
    )
    optimizer.abstract_controller = abstract_controller
    optimizer.success = ⊆(abstract_problem.initial_set, invariant_set_symbols)

    if optimizer.success
        optimizer.print_level >= 1 && println("✅ Safety problem is solvable: initial set is safe-controllable")
    else
        optimizer.print_level >= 1 && println("⚠️ Warning: initial set is only partially safe-controllable")
    end
    optimizer.print_level >= 1 && println("\n Reachability: terminated with $(optimizer.success)")
    optimizer.print_level >= 1 && println("Controllable set size: $(length(invariant_set_symbols))")

    optimizer.abstract_problem_time_sec = time() - t_ref
    return
end

function build_abstract_problem(
    concrete_problem::PR.SafetyProblem,
    abstract_system::SY.TimedHybridSymbolicModel,
)
    concrete_initial_state = concrete_problem.initial_set
    abstract_initial_set = [
        SY.get_abstract_state(
            abstract_system,
            concrete_initial_state,
        ),
    ]

    concrete_safe_set = concrete_problem.safe_set
    abstract_safe_set = SY.get_states_from_set(
        abstract_system,
        concrete_safe_set[1], # state
        concrete_safe_set[2], # time 
        concrete_safe_set[3]; # mode 
        domain = MP.OUTER,
    )

    return PR.SafetyProblem(
        abstract_system,
        abstract_initial_set,
        abstract_safe_set,
        concrete_problem.time,
    )
end

function safe(concrete_problem::PR.SafetyProblem, aug_state)
    (x, t, k) = aug_state
    (Xs_safe, Ts_safe, Ns_safe) = concrete_problem.safe_set
    idx = findfirst(==(k), Ns_safe)
    if isnothing(idx)
        return false
    end
    X_set = Xs_safe[idx]
    T_set = Ts_safe[idx]
    in_X = x ∈ X_set
    in_T = t ≥ T_set.lb[1] && t ≤ T_set.ub[1]

    return in_X && in_T
end

# ================================================================
# Specialized invariant set computation for timed hybrid systems
# ================================================================
mutable struct SymbolicControlTable
    U::Vector{Vector{Int}}
end
SymbolicControlTable(nstates::Int) = SymbolicControlTable([Int[] for _ in 1:nstates])

function ensure_states!(C::SymbolicControlTable, nstates::Int)
    n = length(C.U)
    nstates <= n && return
    resize!(C.U, nstates)
    for i in (n + 1):nstates
        C.U[i] = Int[]
    end
end

add_control!(C::SymbolicControlTable, qs::Int, u::Int) =
    (ensure_states!(C, qs); push!(C.U[qs], u))
is_defined(C::SymbolicControlTable, qs::Int) = !isempty(C.U[qs])

struct PredicateDomain{F}
    ;
    pred::F;
end
Base.in(x, X::PredicateDomain) = X.pred(x)

function to_ms_abstract_controller(C::SymbolicControlTable)
    qfun = (qs::Int) -> C.U[qs]
    X = PredicateDomain((qs::Int) -> is_defined(C, qs))
    return MS.ConstrainedBlackBoxMap(1, 1, qfun, X)  # (qs)->Vector{Int}
end

function compute_largest_invariant_set_timed_hybrid(autom, safelist)
    nstates = SY.get_n_state(autom)
    nsymbols = SY.get_n_input(autom)

    pairstable = [false for i in 1:nstates, j in 1:nsymbols]

    for target in SY.enum_states(autom)
        for (source, symbol) in SY.pre(autom, target)
            pairstable[source, symbol] = true
        end
    end

    nsymbolslist = sum(pairstable; dims = 2)

    safeset = Set(safelist)
    terminal_states = Set{Int}()
    for source in safeset
        if nsymbolslist[source] == 0
            push!(terminal_states, source)
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
    iteration = 0
    controller = SymbolicControlTable(nstates)

    while true
        iteration += 1
        truly_unsafe = setdiff(unsafeset, terminal_states)

        for target in truly_unsafe
            for (source, symbol) in SY.pre(autom, target)
                if pairstable[source, symbol]
                    pairstable[source, symbol] = false
                    nsymbolslist[source] -= 1

                    if nsymbolslist[source] == 0 && !(source in terminal_states)
                        push!(nextunsafeset, source)
                    end
                end
            end
        end

        if isempty(nextunsafeset)
            break
        end

        setdiff!(safeset, nextunsafeset)
        union!(unsafeset, nextunsafeset)
        nextunsafeset = Set{Int}()
    end

    for source in safeset
        for symbol in 1:nsymbols
            if pairstable[source, symbol]
                add_control!(controller, source, symbol)
            end
        end
    end

    invariant_set_complement = setdiff(Set(safelist), safeset)
    return to_ms_abstract_controller(controller), safeset, invariant_set_complement
end