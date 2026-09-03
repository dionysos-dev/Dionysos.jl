const MP = DI.Mapping
const SY = DI.Symbolic

# ============================================================
# Quotient automaton wrapper
# ============================================================

# The abstract input alphabet is the *edge choice* `(mode, destination node)`, not the bare
# mode. On a deterministic graph the two coincide up to renaming and nothing changes. On a
# nondeterministic path-complete graph — an incomplete one necessarily is — a mode has several
# destination nodes, and that choice is the controller's own (it is the controller announcing
# its future), not the environment's: an alphabet of bare modes would fold it into the
# ∀-successor position of the fixed points and synthesis would be needlessly conservative. The
# concrete boundary translates back: actuation uses `mode_of_input`, the announcement steers the
# controller's memory.
struct QuotientAutomaton{QT} <: SY.AbstractAutomatonList
    quotient::QT
    qids::Vector{Int}
    id2idx::Dict{Int, Int}
    inputs::Vector{Tuple{Int, Any}}          # symbol -> (mode, destination node)
    input_index::Dict{Tuple{Int, Any}, Int}
    post_tab::Vector{Vector{Vector{Int}}}
    pre_tab::Vector{Vector{Tuple{Int, Int}}}
    ninput::Int
end

SY.get_n_state(Q::QuotientAutomaton) = length(Q.qids)
SY.get_n_input(Q::QuotientAutomaton) = Q.ninput
SY.post(Q::QuotientAutomaton, s::Int, u::Int) = Q.post_tab[s][u]
SY.pre(Q::QuotientAutomaton, t::Int) = Q.pre_tab[t]

"The concrete mode actuated by abstract input symbol `u`."
mode_of_input(Q::QuotientAutomaton, u::Int) = Q.inputs[u][1]

"The distinct concrete modes appearing in the alphabet, sorted."
modes_of(Q::QuotientAutomaton) = sort(unique(first.(Q.inputs)))

function SY.enum_transitions(Q::QuotientAutomaton)
    trans = Tuple{Int, Int, Int}[]
    for s in 1:SY.get_n_state(Q), u in 1:SY.get_n_input(Q)
        for s2 in SY.post(Q, s, u)
            push!(trans, (s2, s, u))
        end
    end
    return trans
end

function QuotientAutomaton(quotient::PCBisimulationQuotient)
    qids = sort(collect(keys(quotient.states)))
    id2idx = Dict(qid => i for (i, qid) in enumerate(qids))

    pairs = Set{Tuple{Int, Any}}()
    for q in values(quotient.states), (m, dst_qid) in q.next
        haskey(quotient.states, dst_qid) || continue
        push!(pairs, (m, quotient.states[dst_qid].node))
    end
    inputs = Vector{Tuple{Int, Any}}(sort!(collect(pairs); by = string))
    input_index = Dict{Tuple{Int, Any}, Int}(p => k for (k, p) in enumerate(inputs))
    ninput = length(inputs)

    n = length(qids)
    post_tab = [[Int[] for _ in 1:ninput] for _ in 1:n]
    pre_tab = [Tuple{Int, Int}[] for _ in 1:n]

    for (i, qid) in enumerate(qids)
        q = quotient.states[qid]
        for (m, dst_qid) in q.next
            haskey(id2idx, dst_qid) || continue
            j = id2idx[dst_qid]
            k = input_index[(m, quotient.states[dst_qid].node)]
            push!(post_tab[i][k], j)
            push!(pre_tab[j], (i, k))
        end
    end

    for i in 1:n, u in 1:ninput
        post_tab[i][u] = unique(post_tab[i][u])
    end

    return QuotientAutomaton(
        quotient,
        qids,
        id2idx,
        inputs,
        input_index,
        post_tab,
        pre_tab,
        ninput,
    )
end

from_autom_to_bis_state(Q::QuotientAutomaton, qs::Int) = Q.qids[qs]

# Indices past `length(Q.qids)` are the pessimistic completion's sink (`SY.complete_with_sink`),
# which corresponds to no quotient state; every translation back skips it.
_is_bis_state(Q::QuotientAutomaton, qs::Int) = qs <= length(Q.qids)

function from_autom_to_bis_states(Q::QuotientAutomaton, qs_list)
    return [from_autom_to_bis_state(Q, qs) for qs in qs_list if _is_bis_state(Q, qs)]
end

function from_autom_to_bis_value_function(Q::QuotientAutomaton, V::AbstractDict)
    out = Dict{Int, Float64}()
    for (qs, v) in V
        _is_bis_state(Q, qs) || continue
        out[from_autom_to_bis_state(Q, qs)] = Float64(v)
    end
    return out
end

function from_autom_to_bis_value_function(Q::QuotientAutomaton, V::AbstractVector)
    out = Dict{Int, Float64}()
    for qs in eachindex(V)
        _is_bis_state(Q, qs) || continue
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

"""
    OptimizerCoSafeLTLOnQuotient{T} <: Dionysos.Optim.AbstractDionysosOptimizer

Co-safe LTL synthesis on a bisimulation quotient built by [`OptimizerBisimulationQuotient`](@ref).

Set `concrete_problem` (a `PR.CoSafeLTLProblem`), `bisimulation_quotient`, and `ap_to_obs`
mapping each atomic proposition to the observation label the quotient carries for it. The
quotient's sparse state ids are packed into a dense automaton, the product with the
specification is solved by `OPDS.OptimizerCoSafeLTLProblem`, and the winning states come back
as `controllable_set`.

`early_stop` decides what synthesis is asked for, and so what `success` means: `true` restricts
the product's initial states to those meeting the problem's initial set, so `success` reports
whether *that* set is controllable; `false` seeds the product with the whole domain, so
`success` reports whether *every* state is. A run that reaches its own initial set but not the
whole domain therefore ends `success = false` with a perfectly usable controller.
"""
mutable struct OptimizerCoSafeLTLOnQuotient{T} <: OP.AbstractDionysosOptimizer
    # inputs
    concrete_problem::Union{Nothing, PR.CoSafeLTLProblem}
    bisimulation_quotient::Any
    ap_to_obs::Dict{Symbol, Int}
    early_stop::Bool
    sparse_input::Bool
    print_level::Int
    # Optional polyhedral backend (e.g. CDDLib.Library()). When set, a folded run measures how
    # much of the slice family the quotient's cells cover and stores `uncovered_fraction` — the
    # atol-erosion caveat every verification result owes; see `covered_fraction`.
    coverage_backend::Any

    # outputs / internals
    quotient_automaton::Any
    abstract_optimizer::Union{Nothing, OPDS.OptimizerCoSafeLTLProblem}
    abstract_problem::Union{Nothing, PR.CoSafeLTLProblem}
    abstract_controller::Union{Nothing, ST.AbstractDiscreteController}
    controllable_set::Union{Nothing, Vector{Int}}
    uncontrollable_set::Union{Nothing, Vector{Int}}
    value_fun_tab::Any
    # Whether the modes were the environment's (`ST.environment_input`), i.e. the quotient was
    # folded and the run answered the universal question. A folded run has no controller.
    environment_folded::Bool
    # 1 − covered_fraction of the quotient, measured when `coverage_backend` is set on a folded
    # run; `nothing` when not measured.
    uncovered_fraction::Union{Nothing, Float64}
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
            nothing, # coverage_backend
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            false,
            nothing, # uncovered_fraction
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
    quotient = optimizer.bisimulation_quotient

    Q = QuotientAutomaton(quotient)
    optimizer.quotient_automaton = Q

    abstract_problem = build_abstract_problem(concrete_problem, Q, optimizer.ap_to_obs)

    # Who owns the switching signal decides the question. `ControlledSwitching` (or any system
    # declaring no environment) leaves the modes as the controller's alphabet — synthesis, as
    # always. `AutonomousSwitching` hands them to the environment: the quotient is folded so the
    # unchanged solver computes `∀` over the modes, after the pessimistic completion that keeps a
    # missing environment move from being a vacuous win.
    environment = ST.environment_input(concrete_problem.system)
    optimizer.environment_folded = environment !== nothing
    if optimizer.environment_folded
        length(environment) == length(modes_of(Q)) || error(
            "The environment owns $(length(environment)) modes but the quotient carries " *
            "$(length(modes_of(Q))); the quotient was built from a different system.",
        )
        # Missing behaviour is judged per MODE: the alphabet's edge symbols are finer than the
        # environment's choices, and a layer that never carried some announcement is a
        # structural non-edge, not a dropped behaviour.
        modes = modes_of(Q)
        mode_group =
            [findfirst(==(mode_of_input(Q, k)), modes) for k in 1:SY.get_n_input(Q)]
        completed, _, ncompleted = SY.complete_with_sink(Q; groups = mode_group)
        if ncompleted > 0 && optimizer.print_level >= 1
            println(
                "Pessimistic completion: $ncompleted (state, mode) pairs had no successor " *
                "and were routed to a losing sink.",
            )
        end
        # The atol-erosion caveat: the cells cover slightly less than the slices they were
        # carved from, and a point the quotient does not cover is a point this universal answer
        # says nothing about. Measured only on request — it costs a volume computation.
        if optimizer.coverage_backend !== nothing
            optimizer.uncovered_fraction =
                1.0 - covered_fraction(quotient; backend = optimizer.coverage_backend)
            optimizer.print_level >= 1 && println(
                "Coverage: the quotient leaves ",
                round(100 * optimizer.uncovered_fraction; digits = 3),
                "% of the slice family uncovered; the verified set says nothing there.",
            )
        end
        abstract_problem =
            PR.remake(abstract_problem; system = SY.FoldedAutomaton(completed))
    end
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
#
# How far the search is allowed to go depends on what is being asked, and the asymmetry is
# deliberate. Tracking where the state went may fall all the way back to a global search, since a
# wrong guess only costs accuracy. Emitting a control, and deciding whether one exists, stops at
# the current node: the control was synthesised for a particular state, and handing back the one
# belonging to some unrelated cell that happens to contain the point would not be sound. A point
# reachable only by the global search is therefore trackable but not controllable, which is the
# honest answer rather than a gap.
function _find_qid_in_node(quotient::PCBisimulationQuotient, node, x)
    for qid in get(quotient.part_ids, node, Int[])
        q = quotient.states[qid]
        if x ∈ q.set
            return qid
        end
    end
    return nothing
end

function _find_successor_qid(quotient::PCBisimulationQuotient, qid::Int, x)
    q = quotient.states[qid]
    for (_, dst_qid) in q.next
        dst = quotient.states[dst_qid]
        if x ∈ dst.set
            return dst_qid
        end
    end
    return nothing
end

function _find_qid_same_node(quotient::PCBisimulationQuotient, qid::Int, x)
    node = quotient.states[qid].node
    return _find_qid_in_node(quotient, node, x)
end

function _find_qid_global(quotient::PCBisimulationQuotient, x)
    for (qid, q) in quotient.states
        if x ∈ q.set
            return qid
        end
    end
    return nothing
end

"""
    build_concrete_controller(automaton, abstract_controller)

Lower a controller synthesised on the quotient automaton onto concrete states.

Takes the automaton and controller directly; [`solve_concrete_problem`](@ref) is the same
thing reached through the optimizer.
"""
function build_concrete_controller(
    Q::QuotientAutomaton,
    abstract_controller::ST.AbstractDiscreteController,
)
    quotient = Q.quotient

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

        haskey(quotient.states, qid) || return nothing
        q = quotient.states[qid]

        qid_use = if x ∈ q.set
            qid
        else
            qid_same = _find_qid_same_node(quotient, qid, x)
            isnothing(qid_same) ? nothing : qid_same
        end

        qid_use === nothing && return nothing

        qs_dense = qid_to_qs(qid_use)
        us = ST.output_control(abstract_controller, qa, qs_dense)
        # The abstract symbols are edge choices (mode, announced node); the plant is actuated
        # with the mode alone — the announcement lives on in the controller's memory.
        us === nothing && return nothing
        syms = us isa AbstractVector ? us : [us]
        isempty(syms) && return nothing
        return [mode_of_input(Q, s) for s in syms]
    end

    g_conc = function (mem, x_for_update)
        qa, qid = mem

        haskey(quotient.states, qid) || return mem

        # Under a nondeterministic announcement the next point lies in one cell of *several*
        # layers, and only the announced one continues the winning strategy — so among the
        # successors that contain the point, prefer one where the abstract controller stays
        # live, and only then fall back to the first geometric match.
        candidates = Int[]
        for (_, dst) in quotient.states[qid].next
            haskey(quotient.states, dst) || continue
            if x_for_update ∈ quotient.states[dst].set && !(dst in candidates)
                push!(candidates, dst)
            end
        end
        qid_next = nothing
        for cand in candidates
            qa_cand = ST.update_state(abstract_controller, qa, qid_to_qs(cand))
            out = ST.output_control(abstract_controller, qa_cand, qid_to_qs(cand))
            if out !== nothing && !(out isa AbstractVector && isempty(out))
                qid_next = cand
                break
            end
        end
        if qid_next === nothing && !isempty(candidates)
            qid_next = first(candidates)
        end

        if isnothing(qid_next)
            qid_next = _find_qid_same_node(quotient, qid, x_for_update)
        end

        if isnothing(qid_next)
            qid_next = _find_qid_global(quotient, x_for_update)
        end

        isnothing(qid_next) && return mem

        qs_dense_next = qid_to_qs(qid_next)
        qa_next = ST.update_state(abstract_controller, qa, qs_dense_next)

        return (qa_next, qid_next)
    end

    is_defined_memx = function (memx)
        mem, x = memx
        qa, qid = mem

        haskey(quotient.states, qid) || return false

        qid_use = if x ∈ quotient.states[qid].set
            qid
        else
            _find_qid_same_node(quotient, qid, x)
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

"""
    solve_concrete_problem(optimizer) -> ST.DiscreteDynamicController

The concrete controller for a solved [`OptimizerCoSafeLTLOnQuotient`](@ref).

This is the entry point callers want after `MOI.optimize!`: it unpacks the quotient automaton
and the controller synthesised on it and hands both to [`build_concrete_controller`](@ref).
Errors if the optimizer has not been run.
"""
function solve_concrete_problem(opt::OptimizerCoSafeLTLOnQuotient)
    opt.environment_folded && error(
        "The modes were the environment's (the system declares autonomous switching), so this " *
        "run answered the universal question: every switching sequence from the controllable " *
        "set satisfies the specification. A verification returns that set, not a controller — " *
        "read the \"controllable_set\" attribute. For synthesis, declare the switching " *
        "controlled (`ST.with_switching(system, HybridSystems.ControlledSwitching())`).",
    )

    Q = opt.quotient_automaton
    Cabs = opt.abstract_controller

    Q === nothing && error("No quotient_automaton available.")
    Cabs === nothing && error("No abstract_controller available.")

    return build_concrete_controller(Q, Cabs)
end

"""
    initial_controller_memory(optimizer, x0) -> (spec_state, quotient_state)

The memory a freshly started controller carries at `x0`.

`x0` may sit in more than one quotient state -- the shells share their boundaries -- so every
candidate is tried and the first whose product state is winning is taken. Failing to be in any
cell and being in cells that are all losing are different failures and are reported as such.
"""
function initial_controller_memory(opt::OptimizerCoSafeLTLOnQuotient, x0)
    Q = opt.quotient_automaton
    absopt = opt.abstract_optimizer
    absprob = opt.abstract_problem

    Q === nothing && error("No quotient_automaton available.")
    absopt === nothing && error("No abstract_optimizer available.")
    absprob === nothing && error("No abstract_problem available.")

    quotient = Q.quotient
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
        if x0 ∈ quotient.states[qid].set
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

"""
    verification_counterexample(opt, x0; nmax = 200)

The environment's winning play from `x0`, for a folded (verification) run in which `x0` was not
verified: a switching word no controller survives, as evidence rather than a bare "no".

The extraction walks the product of the folded quotient with the specification monitor. A state
outside the verified set always has, by the fixed point's complement, a successor outside it, and
a walk that follows such successors never reaches an accepting product state — so the emitted
word never completes a good prefix, and the run it induces violates the specification. The walk
ends when a product state repeats (a lasso: prefix + cycle, the environment repeating the cycle
forever) or when it enters the pessimistic completion's sink, where the abstraction has no
successor for some mode and the word beyond it is unmodelled.

Returns `(; modes, qids, lasso_start, entered_sink, X)`:

- `modes` — the switching word, one mode per step;
- `qids` — the quotient state ids visited (the sink, having no quotient state, ends the list);
- `lasso_start` — index into `modes` where the repeating cycle begins, or `0` if the walk ended
  in the sink;
- `entered_sink` — whether it did; such a counterexample shows where the abstraction's coverage
  ends rather than a concrete violation, and `X` stops at that point;
- `X` — the concrete replay of the word from `x0`, exact for a switched linear system since the
  quotient is a bisimulation.
"""
function verification_counterexample(opt::OptimizerCoSafeLTLOnQuotient, x0; nmax::Int = 200)
    opt.environment_folded || error(
        "Counterexamples are extracted from a verification run; this run synthesised — its " *
        "evidence is the controller, not an environment play.",
    )

    Q = opt.quotient_automaton
    absopt = opt.abstract_optimizer
    absprob = opt.abstract_problem
    quotient = Q.quotient

    product_automaton_opt = absopt.product_automaton_optimizer
    P = product_automaton_opt.problem.system
    W = product_automaton_opt.controllable_set

    labeling =
        absprob.labeling isa Function ? absprob.labeling :
        OPDS.labeling_function_from_state_sets(absprob.labeling)
    spec = absprob.spec

    # The completed quotient automaton underneath the fold: its per-symbol successors name the
    # environment's move behind every folded edge, and the symbol's mode is the word letter.
    completed = SY.inner(absprob.system)
    nsymbols = SY.get_n_input(completed)
    nquotient = length(Q.qids)

    # The product state x0 starts from. Boundary points sit in several cells; any unverified one
    # carries a counterexample, so the first is taken.
    p0 = nothing
    for (qs, qid) in enumerate(Q.qids)
        x0 ∈ quotient.states[qid].set || continue
        qa = OPDS.step(spec, OPDS.init_state(spec), labeling(qs))
        p = get(P.pid, (qs, qa), nothing)
        if p !== nothing && !(p in W)
            p0 = p
            break
        end
    end
    p0 === nothing && error(
        "Every quotient cell containing x0 is verified; there is no counterexample to extract.",
    )

    modes = Int[]
    qids = Int[quotient.states[Q.qids[P.rev[p0][1]]].id]
    seen = Dict{Int, Int}(p0 => 1)
    lasso_start = 0
    entered_sink = false

    p = p0
    for _ in 1:nmax
        qs = P.rev[p][1]
        # The complement of the attractor: some successor also avoids the verified set.
        p2 = nothing
        for cand in SY.post(P, p, 1)
            if !(cand in W)
                p2 = cand
                break
            end
        end
        p2 === nothing && error(
            "No unverified successor from an unverified product state; the fixed point and " *
            "the product disagree, which indicates a bug rather than a property of the system.",
        )

        qs2 = P.rev[p2][1]
        k = findfirst(k -> qs2 in SY.post(completed, qs, k), 1:nsymbols)
        push!(modes, mode_of_input(Q, k))

        if qs2 > nquotient
            entered_sink = true
            break
        end
        push!(qids, Q.qids[qs2])

        if haskey(seen, p2)
            lasso_start = seen[p2]
            break
        end
        seen[p2] = length(modes) + 1
        p = p2
    end

    # Concrete replay: the word applied to x0. Exactness is the bisimulation's promise; the
    # trajectory is returned so the violation can be inspected on the real system.
    A = UT.mode_matrices(opt.concrete_problem.system)
    x = collect(Float64, x0)
    X = [copy(x)]
    for m in modes
        x = A[m] * x
        push!(X, copy(x))
    end

    return (; modes, qids, lasso_start, entered_sink, X)
end

function simulate_closed_loop(
    f::HybridSystems.HybridSystem,
    controller::ST.AbstractController,
    x0,
    mem0;
    N::Int = 20,
    update_on_next::Bool = true,
)
    # one extraction for the whole run, and the same one the quotient was built from -- the
    # mode matrices were being unpacked again on every step, by a second copy of the rule
    A = UT.mode_matrices(f)

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

        xnext = A[Int(u)] * x

        x_for_update = update_on_next ? xnext : x
        mem = ST.update_state(controller, mem, x_for_update)

        x = xnext
        push!(xs, x)
        push!(mems, mem)
    end

    return (X = xs, U = us, M = mems)
end
