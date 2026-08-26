const MP = DI.Mapping
const SY = DI.Symbolic

# ============================================================
# Quotient automaton wrapper
# ============================================================

struct QuotientAutomaton{QT} <: SY.AbstractAutomatonList
    quotient::QT
    qids::Vector{Int}
    id2idx::Dict{Int, Int}
    post_tab::Vector{Vector{Vector{Int}}}
    pre_tab::Vector{Vector{Tuple{Int, Int}}}
    ninput::Int
end

SY.get_n_state(Q::QuotientAutomaton) = length(Q.qids)
SY.get_n_input(Q::QuotientAutomaton) = Q.ninput
SY.post(Q::QuotientAutomaton, s::Int, u::Int) = Q.post_tab[s][u]
SY.pre(Q::QuotientAutomaton, t::Int) = Q.pre_tab[t]

function SY.enum_transitions(Q::QuotientAutomaton)
    trans = Tuple{Int, Int, Int}[]
    for s in 1:SY.get_n_state(Q), u in 1:SY.get_n_input(Q)
        for s2 in SY.post(Q, s, u)
            push!(trans, (s2, s, u))
        end
    end
    return trans
end

function QuotientAutomaton(T::PCBisimulationQuotient)
    qids = sort(collect(keys(T.states)))
    id2idx = Dict(qid => i for (i, qid) in enumerate(qids))

    modes = Int[]
    for q in values(T.states), (m, _) in q.next
        push!(modes, m)
    end
    ninput = isempty(modes) ? 0 : maximum(modes)

    n = length(qids)
    post_tab = [[Int[] for _ in 1:ninput] for _ in 1:n]
    pre_tab = [Tuple{Int, Int}[] for _ in 1:n]

    for (i, qid) in enumerate(qids)
        q = T.states[qid]
        for (m, dst_qid) in q.next
            haskey(id2idx, dst_qid) || continue
            j = id2idx[dst_qid]
            push!(post_tab[i][m], j)
            push!(pre_tab[j], (i, m))
        end
    end

    for i in 1:n, u in 1:ninput
        post_tab[i][u] = unique(post_tab[i][u])
    end

    return QuotientAutomaton(T, qids, id2idx, post_tab, pre_tab, ninput)
end

from_autom_to_bis_state(Q::QuotientAutomaton, qs::Int) = Q.qids[qs]

function from_autom_to_bis_states(Q::QuotientAutomaton, qs_list)
    return [from_autom_to_bis_state(Q, qs) for qs in qs_list]
end

function from_autom_to_bis_value_function(Q::QuotientAutomaton, V::AbstractDict)
    out = Dict{Int, Float64}()
    for (qs, v) in V
        out[from_autom_to_bis_state(Q, qs)] = Float64(v)
    end
    return out
end

function from_autom_to_bis_value_function(Q::QuotientAutomaton, V::AbstractVector)
    out = Dict{Int, Float64}()
    for qs in eachindex(V)
        out[from_autom_to_bis_state(Q, qs)] = Float64(V[qs])
    end
    return out
end

function get_states_from_set(Q::QuotientAutomaton, X0)
    X0h = UT._as_hpolytope(X0)
    out = Int[]
    for (i, qid) in enumerate(Q.qids)
        q = Q.quotient.states[qid]
        UT.is_disjoint(q.set, X0h) || push!(out, i)
    end
    return out
end

# ============================================================
# Optimizer
# ============================================================

mutable struct OptimizerCoSafeLTLOnQuotient{T} <: OP.AbstractDionysosOptimizer
    # inputs
    concrete_problem::Union{Nothing, PR.CoSafeLTLProblem}
    bisimulation_quotient::Any
    ap_to_obs::Dict{Symbol, Int}
    early_stop::Bool
    sparse_input::Bool
    print_level::Int

    # outputs / internals
    quotient_automaton::Any
    abstract_optimizer::Union{Nothing, OPDS.OptimizerCoSafeLTLProblem}
    abstract_problem::Union{Nothing, PR.CoSafeLTLProblem}
    abstract_controller::Union{Nothing, ST.AbstractDiscreteController}
    controllable_set::Union{Nothing, Vector{Int}}
    uncontrollable_set::Union{Nothing, Vector{Int}}
    value_fun_tab::Any
    success::Bool
    solve_time_sec::T

    function OptimizerCoSafeLTLOnQuotient{T}() where {T}
        return new{T}(
            nothing,
            nothing,
            Dict{Symbol, Int}(),
            true,
            false,
            1,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            false,
            zero(T),
        )
    end
end

OptimizerCoSafeLTLOnQuotient() = OptimizerCoSafeLTLOnQuotient{Float64}()

MOI.is_empty(opt::OptimizerCoSafeLTLOnQuotient) = opt.concrete_problem === nothing

function MOI.optimize!(optimizer::OptimizerCoSafeLTLOnQuotient)
    t0 = time()

    optimizer.concrete_problem === nothing && error("concrete_problem not set")
    optimizer.bisimulation_quotient === nothing && error("bisimulation_quotient not set")
    isempty(optimizer.ap_to_obs) && error("ap_to_obs not set")

    concrete_problem = optimizer.concrete_problem
    T = optimizer.bisimulation_quotient

    Q = QuotientAutomaton(T)
    optimizer.quotient_automaton = Q

    abstract_problem = build_abstract_problem(concrete_problem, Q, optimizer.ap_to_obs)
    optimizer.abstract_problem = abstract_problem

    abstract_optimizer = MOI.instantiate(OPDS.OptimizerCoSafeLTLProblem)
    MOI.set(abstract_optimizer, MOI.RawOptimizerAttribute("problem"), abstract_problem)
    MOI.set(
        abstract_optimizer,
        MOI.RawOptimizerAttribute("early_stop"),
        optimizer.early_stop,
    )
    MOI.set(
        abstract_optimizer,
        MOI.RawOptimizerAttribute("sparse_input"),
        optimizer.sparse_input,
    )
    MOI.set(
        abstract_optimizer,
        MOI.RawOptimizerAttribute("print_level"),
        optimizer.print_level,
    )

    MOI.optimize!(abstract_optimizer)

    optimizer.abstract_optimizer = abstract_optimizer
    optimizer.abstract_controller = abstract_optimizer.controller

    optimizer.controllable_set =
        sort(from_autom_to_bis_states(Q, abstract_optimizer.controllable_set))
    optimizer.uncontrollable_set =
        sort(from_autom_to_bis_states(Q, abstract_optimizer.uncontrollable_set))
    optimizer.value_fun_tab =
        from_autom_to_bis_value_function(Q, abstract_optimizer.value_fun_tab)

    optimizer.success = success(abstract_optimizer, abstract_problem.initial_set)
    optimizer.print_level >= 1 && println("Initial set controllable: ", optimizer.success)
    optimizer.solve_time_sec = time() - t0

    return
end

function build_abstract_problem(
    concrete_problem::PR.CoSafeLTLProblem,
    Q::QuotientAutomaton,
    ap_to_obs::Dict{Symbol, Int};
    atol::Float64 = 0.0,
)
    init_states = get_states_from_set(Q, concrete_problem.initial_set)
    isempty(init_states) && error("Initial set does not intersect any quotient state.")

    lab_abs = quotient_labeling_from_obs(Q, ap_to_obs)

    return PR.CoSafeLTLProblem(
        Q,
        init_states,
        concrete_problem.spec,
        lab_abs,
        concrete_problem.ap_semantics,
    )
end

# ============================================================
# Success condition
# ============================================================

"""
    success(abstract_optimizer, abstract_initial_set) -> Bool

Whether every abstract state the problem starts from ended up controllable.

The sub-solver reports success against the states it was *constructed* from, and under
`early_stop = false` those are every state of the quotient: the controller is deliberately
built over the whole domain so that the controllable set can be inspected. Its flag
therefore answers whether the entire domain is controllable, which a single state violating
the specification outright makes false, and says nothing about the problem posed. This asks
the question the caller means instead — is the initial set covered?
"""
function success(abstract_optimizer, abstract_initial_set)
    controllable = Set(abstract_optimizer.controllable_set)
    return OPDS.covers_initial_set(q -> q in controllable, abstract_initial_set)
end

# ============================================================
# Build quotient AP labeling from obs
# ============================================================

function quotient_labeling_from_obs(Q::QuotientAutomaton, ap_to_obs::Dict{Symbol, Int})
    lab = Dict{Symbol, Vector{Int}}()
    for (ap, obsval) in ap_to_obs
        lab[ap] =
            Int[i for i in 1:length(Q.qids) if Q.quotient.states[Q.qids[i]].obs == obsval]
    end
    return lab
end

# ============================================================
# Concretize the controller
# ============================================================

# Locating the quotient state a concrete point sits in, cheapest first.
#
# The closed loop asks this at every step, and the three are tried in order: the successors of
# the current state, then the rest of its node, then the whole quotient. The first tier answers
# almost always -- it is the transition structure doing its job -- and each fallback is both
# wider and more expensive, the last being linear in the number of states.
#
# The fallbacks are not defensive padding. The partition does not quite cover the domain: every
# cut insets by `atol`, so thin gaps run along the shell boundaries (see
# `OptimizerBisimulationQuotient`). A trajectory that lands in one belongs to no state at all,
# and a trajectory that drifts across a boundary lands in a state that is not among the
# predicted successors. The wider searches are what keep simulation going in both cases.
function _find_qid_in_node(T::PCBisimulationQuotient, node, x)
    for qid in get(T.part_ids, node, Int[])
        q = T.states[qid]
        if x ∈ q.set
            return qid
        end
    end
    return nothing
end

function _find_successor_qid(T::PCBisimulationQuotient, qid::Int, x)
    q = T.states[qid]
    for (_, dst_qid) in q.next
        dst = T.states[dst_qid]
        if x ∈ dst.set
            return dst_qid
        end
    end
    return nothing
end

function _find_qid_same_node(T::PCBisimulationQuotient, qid::Int, x)
    node = T.states[qid].node
    return _find_qid_in_node(T, node, x)
end

function _find_qid_global(T::PCBisimulationQuotient, x)
    for (qid, q) in T.states
        if x ∈ q.set
            return qid
        end
    end
    return nothing
end

function solve_concrete_problem_lifted(
    Q::QuotientAutomaton,
    abstract_controller::ST.AbstractDiscreteController,
)
    T = Q.quotient

    qid_to_qs(qid::Int) = Q.id2idx[qid]
    qs_to_qid(qs::Int) = Q.qids[qs]

    function pick_symbol(us)
        us === nothing && return nothing
        if us isa AbstractVector
            isempty(us) && return nothing
            return first(us)
        end
        return us
    end

    h_conc = function (mem, x)
        qa, qid = mem

        haskey(T.states, qid) || return nothing
        q = T.states[qid]

        qid_use = if x ∈ q.set
            qid
        else
            qid_same = _find_qid_same_node(T, qid, x)
            isnothing(qid_same) ? nothing : qid_same
        end

        qid_use === nothing && return nothing

        qs_dense = qid_to_qs(qid_use)
        us = ST.output_control(abstract_controller, qa, qs_dense)
        return us
    end

    g_conc = function (mem, x_for_update)
        qa, qid = mem

        haskey(T.states, qid) || return mem

        qid_next = _find_successor_qid(T, qid, x_for_update)

        if isnothing(qid_next)
            qid_next = _find_qid_same_node(T, qid, x_for_update)
        end

        if isnothing(qid_next)
            qid_next = _find_qid_global(T, x_for_update)
        end

        isnothing(qid_next) && return mem

        qs_dense_next = qid_to_qs(qid_next)
        qa_next = ST.update_state(abstract_controller, qa, qs_dense_next)

        return (qa_next, qid_next)
    end

    is_defined_memx = function (memx)
        mem, x = memx
        qa, qid = mem

        haskey(T.states, qid) || return false

        qid_use = if x ∈ T.states[qid].set
            qid
        else
            _find_qid_same_node(T, qid, x)
        end

        isnothing(qid_use) && return false

        qs_dense = qid_to_qs(qid_use)
        u = ST.output_control(abstract_controller, qa, qs_dense)

        return u !== nothing
    end

    X_memx = ST.PredicateDomain(is_defined_memx)

    x0_abs = ST.initial_state(abstract_controller)

    return ST.DiscreteDynamicController(x0_abs, X_memx, g_conc, h_conc, false)
end

function solve_concrete_problem(opt::OptimizerCoSafeLTLOnQuotient)
    Q = opt.quotient_automaton
    Cabs = opt.abstract_controller

    Q === nothing && error("No quotient_automaton available.")
    Cabs === nothing && error("No abstract_controller available.")

    return solve_concrete_problem_lifted(Q, Cabs)
end

function initial_lifted_controller_memory(opt::OptimizerCoSafeLTLOnQuotient, x0)
    Q = opt.quotient_automaton
    absopt = opt.abstract_optimizer
    absprob = opt.abstract_problem

    Q === nothing && error("No quotient_automaton available.")
    absopt === nothing && error("No abstract_optimizer available.")
    absprob === nothing && error("No abstract_problem available.")

    T = Q.quotient
    product_automaton_opt = absopt.product_automaton_optimizer
    P = product_automaton_opt.problem.system
    W = product_automaton_opt.controllable_set

    labeling =
        absprob.labeling isa Function ? absprob.labeling :
        OPDS.labeling_function_from_state_sets(absprob.labeling)

    spec = absprob.spec
    spec isa OPDS.AbstractSpecStepper ||
        error("Unsupported spec type $(typeof(spec)); expected AbstractSpecStepper")

    candidates = Tuple{Int, Int}[]
    for (qs, qid) in enumerate(Q.qids)
        if x0 ∈ T.states[qid].set
            push!(candidates, (qs, qid))
        end
    end

    isempty(candidates) &&
        error("Initial concrete state is not contained in any quotient cell.")

    for (qs, qid) in candidates
        qa_init = OPDS.step(spec, OPDS.init_state(spec), labeling(qs))
        p0 = get(P.pid, (qs, qa_init), nothing)
        if p0 !== nothing && p0 in W
            return (qa_init, qid)
        end
    end

    return error("Initial concrete state belongs to quotient cells, but none is winning.")
end

function initial_controller_memory(opt::OptimizerCoSafeLTLOnQuotient, x0)
    return initial_lifted_controller_memory(opt, x0)
end

function simulate_closed_loop(
    f::HybridSystems.HybridSystem,
    controller::ST.AbstractController,
    x0,
    mem0;
    N::Int = 20,
    update_on_next::Bool = true,
)
    A = f.resetmaps

    xs = Vector{typeof(x0)}()
    us = Int[]
    mems = Any[]

    x = x0
    mem = mem0

    push!(xs, x)
    push!(mems, mem)

    for _ in 1:N
        u = ST.output_control(controller, mem, x)
        u === nothing && break

        if u isa AbstractVector
            isempty(u) && break
            u = first(u)
        end

        push!(us, u)

        Ak = A[Int(u)]
        xnext = if Ak isa AbstractMatrix
            Ak * x
        elseif hasproperty(Ak, :A)
            getproperty(Ak, :A) * x
        else
            error("Cannot simulate reset map of type $(typeof(Ak))")
        end

        x_for_update = update_on_next ? xnext : x
        mem = ST.update_state(controller, mem, x_for_update)

        x = xnext
        push!(xs, x)
        push!(mems, mem)
    end

    return (X = xs, U = us, M = mems)
end
