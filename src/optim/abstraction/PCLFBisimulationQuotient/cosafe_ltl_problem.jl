mutable struct OptimizerCoSafeLTLOnQuotient{T} <: MOI.AbstractOptimizer
    # -----------------------------
    # inputs
    # -----------------------------
    concrete_problem::Union{Nothing, PR.CoSafeLTLProblem}
    bisimulation_quotient::Any

    # map AP => observation label in quotient
    # example: Dict(:p => 1, :q => 2, :terminal => -1)
    ap_to_obs::Dict{Symbol, Int}

    # semantics for initial-set lifting onto quotient states
    initial_set_semantics::Any

    # -----------------------------
    # outputs / internals
    # -----------------------------
    quotient_automaton::Any
    abstract_problem::Union{Nothing, PR.CoSafeLTLProblem}
    product_system::Any
    abstract_controller_product::Union{Nothing, MS.ConstrainedBlackBoxMap}
    abstract_controller::Union{Nothing, MS.SystemWithOutput}
    qa0::Any
    controllable_set::Any
    uncontrollable_set::Any
    value_fun_tab::Any
    abstract_problem_time_sec::T
    success::Bool
    print_level::Int

    function OptimizerCoSafeLTLOnQuotient{T}() where {T}
        return new{T}(
            nothing,                    # concrete_problem
            nothing,                    # bisimulation_quotient
            Dict{Symbol, Int}(),        # ap_to_obs
            MP.OUTER,                   # initial_set_semantics
            nothing,                    # quotient_automaton
            nothing,                    # abstract_problem
            nothing,                    # product_system
            nothing,                    # abstract_controller_product
            nothing,                    # abstract_controller
            nothing,                    # qa0
            nothing,                    # controllable_set
            nothing,                    # uncontrollable_set
            nothing,                    # value_fun_tab
            zero(T),                    # solve time
            false,                      # success
            1,                          # print_level
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

MOI.get(opt::OptimizerCoSafeLTLOnQuotient, ::MOI.SolveTimeSec) = opt.abstract_problem_time_sec





# ============================================================
# Quotient automaton wrapper
# ============================================================

struct QuotientAutomaton{QT} <: SY.AbstractAutomatonList{0, 0}
    quotient::QT
    qids::Vector{Int}                      # dense index i -> original quotient state id
    id2idx::Dict{Int, Int}                 # original quotient state id -> dense index
    post_tab::Vector{Vector{Vector{Int}}}  # post_tab[s][u] = successors
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
    for q in values(T.states)
        for (m, _) in q.next
            push!(modes, m)
        end
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

# ============================================================
# AP coverage check for quotient mode
# ============================================================

collect_aps(φ::Spot.SpotFormula) = [Symbol(ap) for ap in Spot.atomic_prop_collect(φ)]

function _check_ap_coverage_on_quotient(
    concrete_problem::PR.CoSafeLTLProblem,
    ap_to_obs::Dict{Symbol, Int},
)
    if concrete_problem.spec isa Spot.SpotFormula
        aps = Set(collect_aps(concrete_problem.spec))
        provided = Set(keys(ap_to_obs))
        missing = setdiff(aps, provided)
        if !isempty(missing)
            error("Missing AP -> obs mapping in `ap_to_obs` for: $(collect(missing)).")
        end
    end
    return nothing
end

# ============================================================
# Lift initial set onto quotient states
# ============================================================

function quotient_states_from_set(
    Q::QuotientAutomaton,
    X0
)
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

function quotient_lab_abs_dense(
    Q::QuotientAutomaton,
    ap_to_obs::Dict{Symbol, Int},
)
    lab = Dict{Symbol, Vector{Int}}()

    for (ap, obsval) in ap_to_obs
        lab[ap] = Int[
            i for i in 1:length(Q.qids)
            if Q.quotient.states[Q.qids[i]].obs == obsval
        ]
    end

    return lab
end

# ============================================================
# Labeling function for product construction
# ============================================================

function labeling_function_from_quotient_obs(
    Q::QuotientAutomaton,
    ap_to_obs::Dict{Symbol, Int},
)
    obs_to_aps = Dict{Int, Vector{Symbol}}()
    for (ap, ob) in ap_to_obs
        get!(obs_to_aps, ob, Symbol[])
        push!(obs_to_aps[ob], ap)
    end

    return function (qs::Int)
        qid = Q.qids[qs]
        q = Q.quotient.states[qid]
        return Tuple(get(obs_to_aps, q.obs, Symbol[]))
    end
end

# ============================================================
# Build abstract problem on quotient
# ============================================================

function build_abstract_problem_on_quotient(
    concrete_problem::PR.CoSafeLTLProblem,
    Q::QuotientAutomaton,
    ap_to_obs::Dict{Symbol, Int};
    initial_set_semantics = MP.OUTER,
    atol::Float64 = 0.0,
)
    _check_ap_coverage_on_quotient(concrete_problem, ap_to_obs)

    init_states = quotient_states_from_set(
        Q,
        concrete_problem.initial_set,
        initial_set_semantics;
        atol = atol,
    )

    isempty(init_states) && error("Initial set does not intersect any quotient state.")

    lab_abs = quotient_lab_abs_dense(Q, ap_to_obs)

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
# Main optimize!
# ============================================================

function MOI.optimize!(optimizer::OptimizerCoSafeLTLOnQuotient)
    t0 = time()

    optimizer.bisimulation_quotient === nothing &&
        error("bisimulation_quotient not set")
    optimizer.concrete_problem === nothing &&
        error("concrete_problem not set")
    isempty(optimizer.ap_to_obs) &&
        error("ap_to_obs not set")

    concrete_problem = optimizer.concrete_problem
    T = optimizer.bisimulation_quotient

    # (1) wrap quotient as automaton
    Q = QuotientAutomaton(T)
    optimizer.quotient_automaton = Q

    # (2) build abstract problem on quotient
    abstract_problem = build_abstract_problem_on_quotient(
        concrete_problem,
        Q,
        optimizer.ap_to_obs;
        initial_set_semantics = optimizer.initial_set_semantics,
    )
    optimizer.abstract_problem = abstract_problem

    if optimizer.print_level > 0
        println("Quotient states           : ", SY.get_n_state(Q))
        println("Quotient inputs           : ", SY.get_n_input(Q))
        println("Initial abstract states   : ", length(abstract_problem.initial_set))
        for (ap, states) in abstract_problem.labeling
            println("AP = $ap, #states = $(length(states))")
        end
    end

    # (3) build spec stepper
    spec = if abstract_problem.spec isa Spot.SpotFormula
        spot_stepper(abstract_problem.spec)
    elseif abstract_problem.spec isa AbstractSpecStepper
        abstract_problem.spec
    else
        error("spec must be SpotFormula or AbstractSpecStepper")
    end

    # (4) build quotient labeling function
    abstract_labeling = labeling_function_from_quotient_obs(Q, optimizer.ap_to_obs)

    # (5) build product automaton
    P = build_product_automaton(
        Q,
        spec,
        abstract_labeling;
        initial_set = abstract_problem.initial_set,
    )
    optimizer.product_system = P

    # (6) accepting product states
    accQ = accepting_states(spec)
    target_set = [p for p in 1:SY.get_n_state(P) if P.rev[p][2] in accQ]
    isempty(target_set) && error("Empty target_set in product automaton.")

    # (7) initial product states
    init_abs = abstract_problem.initial_set
    init_prod = Int[]
    for qs in init_abs
        if haskey(P.pid, (qs, init_state(spec)))
            push!(init_prod, P.pid[(qs, init_state(spec))])
        end
    end
    isempty(init_prod) && error("No reachable initial product states.")

    # (8) solve reachability / worst-case cost
    controller_product_automaton, controllableP, uncontrollableP, V =
        SY.compute_worst_case_uniform_cost_controller(
            P,
            target_set;
            initial_set = init_prod,
        )

    # (9) build finite-memory controller on quotient states
    abstract_controller, qa0 = build_abstract_fm_controller_ms(
        abstract_problem.labeling,
        spec,
        controller_product_automaton,
        P.pid,
    )

    optimizer.abstract_controller_product = controller_product_automaton
    optimizer.abstract_controller = abstract_controller
    optimizer.qa0 = qa0
    optimizer.controllable_set = controllableP
    optimizer.uncontrollable_set = uncontrollableP
    optimizer.value_fun_tab = V
    optimizer.success = all(p -> p in controllableP, init_prod)
    optimizer.abstract_problem_time_sec = time() - t0

    if optimizer.print_level > 0
        println("Product states            : ", SY.get_n_state(P))
        println("Initial product states    : ", length(init_prod))
        println("Target product states     : ", length(target_set))
        println("Success                   : ", optimizer.success)
    end

    return
end