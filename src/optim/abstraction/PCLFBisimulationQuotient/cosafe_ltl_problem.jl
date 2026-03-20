const MP = DI.Mapping
const SY = DI.Symbolic

import MathematicalSystems as MS

# ============================================================
# Quotient automaton wrapper
# ============================================================

struct QuotientAutomaton{QT} <: SY.AbstractAutomatonList{0, 0}
    quotient::QT
    qids::Vector{Int}                      # dense index -> quotient state id
    id2idx::Dict{Int, Int}                 # quotient state id -> dense index
    post_tab::Vector{Vector{Vector{Int}}}  # post_tab[s][u]
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

mutable struct OptimizerCoSafeLTLOnQuotient{T} <: MOI.AbstractOptimizer
    # inputs
    concrete_problem::Union{Nothing, PR.CoSafeLTLProblem}
    bisimulation_quotient::Any
    ap_to_obs::Dict{Symbol, Int}
    sparse_input::Bool
    print_level::Int

    # outputs / internals
    quotient_automaton::Any
    abstract_optimizer::Union{Nothing, SY.OptimizerCoSafeLTLProblem}
    abstract_problem::Union{Nothing, PR.CoSafeLTLProblem}
    abstract_controller::Union{Nothing, MS.SystemWithOutput}
    qa0::Union{Nothing, Int}
    success::Bool
    solve_time_sec::T

    function OptimizerCoSafeLTLOnQuotient{T}() where {T}
        return new{T}(
            nothing,                    # concrete_problem
            nothing,                    # bisimulation_quotient
            Dict{Symbol, Int}(),        # ap_to_obs
            false,                      # sparse_input
            1,                          # print_level
            nothing,                    # quotient_automaton
            nothing,                    # abstract_optimizer
            nothing,                    # abstract_problem
            nothing,                    # abstract_controller
            nothing,                    # qa0
            false,                      # success
            zero(T),                    # solve_time_sec
        )
    end
end

OptimizerCoSafeLTLOnQuotient() = OptimizerCoSafeLTLOnQuotient{Float64}()

MOI.is_empty(opt::OptimizerCoSafeLTLOnQuotient) = opt.concrete_problem === nothing

function MOI.set(opt::OptimizerCoSafeLTLOnQuotient, p::MOI.RawOptimizerAttribute, v)
    return setproperty!(opt, Symbol(p.name), v)
end

function MOI.get(opt::OptimizerCoSafeLTLOnQuotient, p::MOI.RawOptimizerAttribute)
    return getproperty(opt, Symbol(p.name))
end

MOI.get(opt::OptimizerCoSafeLTLOnQuotient, ::MOI.SolveTimeSec) = opt.solve_time_sec

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

    abstract_optimizer = MOI.instantiate(SY.OptimizerCoSafeLTLProblem)
    MOI.set(abstract_optimizer, MOI.RawOptimizerAttribute("problem"), abstract_problem)
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
    optimizer.qa0 = abstract_optimizer.qa0

    # For each concrete initial state, there exists at least one winning lifted representative
    optimizer.success = success(abstract_optimizer, concrete_problem.initial_set)
    optimizer.print_level >= 1 &&
        println("Success of concrete problem: ", optimizer.success)
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
        concrete_problem.strict_spot,
    )
end

# ============================================================
# Success condition
# ============================================================

# assuming concrete_initial_set is a singleton, success if:
# there exists a winning lifted representative for the concrete initial point,
function success(abstract_optimizer, concrete_initial_set)
    return any(p -> p in abstract_optimizer.controllable_set, abstract_optimizer.init_set)
end

# ============================================================
# Lift initial set onto quotient states
# ============================================================

function get_states_from_set(Q::QuotientAutomaton, X0)
    X0h = _as_hpolytope(X0)
    out = Int[]
    for (i, qid) in enumerate(Q.qids)
        q = Q.quotient.states[qid]
        I = set_intersection(q.set, X0h)
        is_nonempty_set(I) && push!(out, i)
    end
    return out
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

struct PredicateDomain{F}
    pred::F
end
Base.in(x, X::PredicateDomain) = X.pred(x)

function _find_qid_in_node(T::PCBisimulationQuotient, node, x)
    for qid in get(T.part_ids, node, Int[])
        q = T.states[qid]
        if x ∈ q.set
            return qid
        end
    end
    return nothing
end

# Efficiency: find among successors first
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

# Fallback: global search in same node
function _find_qid_same_node(T::PCBisimulationQuotient, qid::Int, x)
    node = T.states[qid].node
    return _find_qid_in_node(T, node, x)
end

# Fallback: global search (slow)
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
    abstract_controller::MS.SystemWithOutput,
)
    T = Q.quotient

    # abstract controller callables
    h_abs = abstract_controller.outputmap.h      # expects ((qa, qs_dense))
    g_abs = MS.mapping(abstract_controller.s)    # expects (qa, qs_dense)

    # ---------------------------------------------------------
    # Conversion helpers
    # ---------------------------------------------------------
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

    # ---------------------------------------------------------
    # Concrete output map: ((qa,qid), x) -> u
    # ---------------------------------------------------------
    h_conc = function (mem, x)
        qa, qid = mem

        haskey(T.states, qid) || return nothing
        q = T.states[qid]
        # Prefer the current cell if x still belongs to it
        qid_use = if x ∈ q.set
            qid
        else
            # fallback: search in same node partition
            qid_same = _find_qid_same_node(T, qid, x)
            isnothing(qid_same) ? nothing : qid_same
        end

        qid_use === nothing && return nothing

        qs_dense = qid_to_qs(qid_use)
        us = h_abs((qa, qs_dense))
        u_sym = pick_symbol(us)
        u_sym === nothing && return nothing

        # for the quotient, mode index is already the concrete switched input
        return u_sym
    end

    # ---------------------------------------------------------
    # Memory update: ((qa,qid), x_for_update) -> (qa_next, qid_next)
    # ---------------------------------------------------------
    g_conc = function (mem, x_for_update)
        qa, qid = mem

        haskey(T.states, qid) || return mem

        # First try successor cells
        qid_next = _find_successor_qid(T, qid, x_for_update)

        # Then fallback to same node
        if isnothing(qid_next)
            qid_next = _find_qid_same_node(T, qid, x_for_update)
        end

        # Then fallback globally
        if isnothing(qid_next)
            qid_next = _find_qid_global(T, x_for_update)
        end

        # If still nothing, keep memory unchanged
        isnothing(qid_next) && return mem

        qs_dense_next = qid_to_qs(qid_next)
        qa_next = g_abs(qa, qs_dense_next)

        return (qa_next, qid_next)
    end

    # ---------------------------------------------------------
    # Domain predicate
    # ---------------------------------------------------------
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
        us = h_abs((qa, qs_dense))
        u_sym = pick_symbol(us)

        return u_sym !== nothing
    end
    X_memx = PredicateDomain(is_defined_memx)

    # ---------------------------------------------------------
    # Dimensions
    # ---------------------------------------------------------
    first_q = first(values(T.states))
    nx = dim(first_q.set)
    nu = 1

    # memory is a 2-tuple (qa,qid)
    outmap = MS.ConstrainedBlackBoxMap(2, nu, memx -> begin
        mem, x = memx
        h_conc(mem, x)
    end, X_memx)

    memsys = MS.BlackBoxControlDiscreteSystem(
        (mem, x_for_update) -> g_conc(mem, x_for_update),
        2,
        nx,
    )

    return MS.SystemWithOutput(memsys, outmap)
end

function solve_concrete_problem(opt::OptimizerCoSafeLTLOnQuotient)
    Q = opt.quotient_automaton
    Cabs = opt.abstract_controller

    Q === nothing && error("No quotient_automaton available.")
    Cabs === nothing && error("No abstract_controller available.")

    return solve_concrete_problem_lifted(Q, Cabs)
end

import Spot
function initial_concrete_controller_memory(opt::OptimizerCoSafeLTLOnQuotient, x0)
    Q = opt.quotient_automaton
    absopt = opt.abstract_optimizer
    absprob = opt.abstract_problem

    Q === nothing && error("No quotient_automaton available.")
    absopt === nothing && error("No abstract_optimizer available.")
    absprob === nothing && error("No abstract_problem available.")

    T = Q.quotient
    P = absopt.product_autom
    W = absopt.controllable_set

    labeling =
        absprob.labeling isa Function ? absprob.labeling :
        SY.labeling_function_from_state_sets(absprob.labeling)

    spec0 = absprob.spec
    spec = if spec0 isa Spot.SpotFormula
        SY.spot_stepper(spec0)
    elseif spec0 isa SY.AbstractSpecStepper
        spec0
    else
        error("Unsupported spec type $(typeof(spec0))")
    end

    candidates = Tuple{Int, Int}[]  # (qs_dense, qid)
    for (qs, qid) in enumerate(Q.qids)
        if x0 ∈ T.states[qid].set
            push!(candidates, (qs, qid))
        end
    end

    isempty(candidates) &&
        error("Initial concrete state is not contained in any quotient cell.")

    for (qs, qid) in candidates
        qa_init = SY.step(spec, SY.init_state(spec), labeling(qs))
        p0 = get(P.pid, (qs, qa_init), nothing)
        if p0 !== nothing && p0 in W
            return (qa_init, qid)
        end
    end

    return error("Initial concrete state belongs to quotient cells, but none is winning.")
end

function initial_controller_memory(opt::OptimizerCoSafeLTLOnQuotient, x0)
    return initial_concrete_controller_memory(opt, x0)
end

function simulate_closed_loop(
    f::HybridSystems.HybridSystem,
    controller::MS.SystemWithOutput,
    x0,
    mem0;
    N::Int = 20,
    update_on_next::Bool = true,
)
    # plant matrices / reset maps
    A = f.resetmaps

    # controller callables
    h = controller.outputmap.h
    g = MS.mapping(controller.s)

    xs = Vector{typeof(x0)}()
    us = Int[]
    mems = Any[]

    x = x0
    mem = mem0

    push!(xs, x)
    push!(mems, mem)

    for k in 1:N
        # controller output
        u = h((mem, x))
        u === nothing && break

        # if controller returns a vector, pick first
        if u isa AbstractVector
            isempty(u) && break
            u = first(u)
        end

        push!(us, u)

        # plant update
        Ak = A[Int(u)]
        xnext = if Ak isa AbstractMatrix
            Ak * x
        elseif hasproperty(Ak, :A)
            getproperty(Ak, :A) * x
        else
            error("Cannot simulate reset map of type $(typeof(Ak))")
        end

        # controller memory update
        x_for_update = update_on_next ? xnext : x
        mem = g(mem, x_for_update)

        x = xnext
        push!(xs, x)
        push!(mems, mem)
    end

    return (X = xs, U = us, M = mems)
end
