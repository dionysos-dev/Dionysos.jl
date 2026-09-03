# ============================================================
# CoSafe LTL Control
# ============================================================

"""
    OptimizerCoSafeLTLProblem{T} <: AbstractDionysosOptimizer

Co-safe LTL synthesis on a finite automaton: given a
[`CoSafeLTLProblem`](@ref Dionysos.Problem.CoSafeLTLProblem) whose system is an
`AbstractAutomatonList`, build the synchronous product of the abstraction with a deterministic
monitor over the atomic propositions, then run the reachability solver on that product with the
monitor's accepting states as target.

The monitor is an interface, not a fixed type: any object answering `step`, an initial state and
a set of accepting states will do. Loading [Spot](https://spot.lre.epita.fr/) supplies one by
translating an LTL formula.

Because the strategy depends on the monitor state and not on the abstract state alone, the
returned controller is **dynamic** — tabulated at construction time, so it stays plain
serializable data even when the monitor is a closure.

Set `"problem"`, optionally `"early_stop"`, `"sparse_input"`; read back `"controller"`,
`"controllable_set"` and `"uncontrollable_set"`.
"""
mutable struct OptimizerCoSafeLTLProblem{T} <: AbstractDionysosOptimizer
    # inputs
    problem::Union{Nothing, PR.CoSafeLTLProblem}
    early_stop::Bool
    sparse_input::Bool
    print_level::Int

    # outputs / internals
    controller::Union{Nothing, ST.AbstractDiscreteController}
    update_on_next::Bool
    controllable_set::Any
    uncontrollable_set::Any
    value_fun_tab::Any
    product_automaton_optimizer::Any
    success::Bool
    solve_time_sec::T

    function OptimizerCoSafeLTLProblem{T}() where {T}
        return new{T}(
            nothing,   # problem
            true,      # early_stop
            false,     # sparse_input
            1,         # print_level
            nothing,   # controller
            true,      # update_on_next
            nothing,   # controllable_set
            nothing,   # uncontrollable_set
            nothing,   # value_fun_tab
            nothing,   # product_automaton_optimizer
            false,     # success
            zero(T),   # solve_time_sec
        )
    end
end

OptimizerCoSafeLTLProblem() = OptimizerCoSafeLTLProblem{Float64}()

MOI.is_empty(opt::OptimizerCoSafeLTLProblem) = opt.problem === nothing

function MOI.optimize!(optimizer::OptimizerCoSafeLTLProblem)
    t0 = time()

    problem = optimizer.problem

    problem === nothing && error("problem not set")
    _check_ap_coverage(problem)

    if optimizer.print_level > 0 && problem.labeling isa AbstractDict
        for (ap, states) in problem.labeling
            println("AP=$ap  #states=$(length(states))")
        end
    end

    # (1) Finite system
    autom = problem.system

    # (2) Specification automaton / monitor
    spec = if problem.spec isa AbstractSpecStepper
        problem.spec
    else
        error("problem.spec must be an AbstractSpecStepper")
    end

    labeling = if problem.labeling isa Function
        problem.labeling
    elseif problem.labeling isa AbstractDict
        labeling_function_from_state_sets(problem.labeling)
    else
        error(
            "problem.labeling must be either a labeling function or Dict{Symbol,<:AbstractVector{Int}}",
        )
    end

    # (3) Product automaton
    construction_states =
        optimizer.early_stop ? problem.initial_set : collect(1:SY.get_n_state(autom))
    product_autom =
        build_product_automaton(autom, spec, labeling; initial_set = construction_states)

    # (4) Solve reachability problem
    product_initial_set = Int[]
    for qs in construction_states
        ap0 = labeling(qs)
        qa_init = step(spec, init_state(spec), ap0)
        p0 = get(product_autom.pid, (qs, qa_init), nothing)
        p0 === nothing && continue
        push!(product_initial_set, p0)
    end
    isempty(product_initial_set) && error("Empty product initial set.")

    # The alphabet the acceptance is judged against is what the system can emit, not the free
    # powerset of atomic propositions — the labels of every abstract state, whether or not the
    # product reached it, since an extension can.
    used_labels = unique(labeling(qs) for qs in 1:SY.get_n_state(autom))
    accQ = accepting_states(spec, used_labels)
    product_target_set =
        [p for p in 1:SY.get_n_state(product_autom) if product_autom.rev[p][2] in accQ]
    isempty(product_target_set) &&
        error("Empty product target_set (AP mismatch or acceptance not found).")

    product_automaton_problem = PR.OptimalControlProblem(
        product_autom,
        product_initial_set,
        product_target_set,
        nothing,  # state_cost
        nothing,  # transition_cost
        0,
    )

    product_automaton_optimizer = MOI.instantiate(OptimizerOptimalControlProblem)
    MOI.set(
        product_automaton_optimizer,
        MOI.RawOptimizerAttribute("problem"),
        product_automaton_problem,
    )
    MOI.set(
        product_automaton_optimizer,
        MOI.RawOptimizerAttribute("early_stop"),
        optimizer.early_stop,
    )
    MOI.set(
        product_automaton_optimizer,
        MOI.RawOptimizerAttribute("sparse_input"),
        optimizer.sparse_input,
    )
    MOI.set(
        product_automaton_optimizer,
        MOI.RawOptimizerAttribute("print_level"),
        optimizer.print_level,
    )

    MOI.optimize!(product_automaton_optimizer)

    optimizer.product_automaton_optimizer = product_automaton_optimizer
    product_controller = product_automaton_optimizer.controller
    product_controllable_set = product_automaton_optimizer.controllable_set
    product_value_fun_tab = product_automaton_optimizer.value_fun_tab
    success = product_automaton_optimizer.success

    problem.labeling isa AbstractDict ||
        error("build_dynamic_controller currently requires dictionary labeling.")

    # (5) Wrap product controller into finite-memory controller on autom
    controller = build_dynamic_controller(
        problem.labeling,
        spec,
        product_controller,
        product_autom.pid,
    )

    optimizer.controller = controller
    optimizer.update_on_next = true
    optimizer.controllable_set = project_initial_memory_controllable_set(
        product_autom,
        product_controllable_set,
        spec,
        labeling,
        SY.get_n_state(autom),
    )
    optimizer.uncontrollable_set = project_initial_memory_uncontrollable_set(
        product_autom,
        product_controllable_set,
        spec,
        labeling,
        SY.get_n_state(autom),
    )
    optimizer.value_fun_tab = project_initial_memory_value_function(
        product_autom,
        product_value_fun_tab,
        spec,
        labeling,
        SY.get_n_state(autom),
    )
    optimizer.success = success
    optimizer.print_level >= 1 && println("Success: ", success)
    optimizer.solve_time_sec = time() - t0

    return
end

# ============================================================
# Helpers
# ============================================================

function _check_ap_coverage(problem)
    if !(problem.labeling isa AbstractDict)
        return nothing
    end

    if problem.spec isa FunctionMonitor && !isempty(problem.spec.aps)
        provided = Set(keys(problem.labeling))
        missing = setdiff(Set(problem.spec.aps), provided)
        if !isempty(missing)
            error("Missing AP sets in problem.labeling for: $(collect(missing)).")
        end
    end

    return nothing
end

# ============================================================
# Spec interface (user-defined monitor)
# ============================================================

abstract type AbstractSpecStepper end

# Step a spec automaton/monitor: qa_next = step(spec, qa, ap_tuple)
# Must return an Int state id.
step(::AbstractSpecStepper, ::Int, ::Tuple{Vararg{Symbol}}) = error("step not implemented")

# Optional: return the initial spec state qa0.
init_state(::AbstractSpecStepper) = 1

# Optional: return the set of "done/accepting" spec states for co-safe reduction.
accepting_states(::AbstractSpecStepper) = error("accepting_states not implemented")

# A state is "done" only relative to what the system can say: a good prefix must survive every
# extension over the labels the product actually emits, not over the free powerset of atomic
# propositions — mutually exclusive observations, for instance, never emit conjunctions. Steppers
# that can exploit the restriction (the spot-backed one does) override this two-argument form; the
# fallback ignores the alphabet.
accepting_states(spec::AbstractSpecStepper, used_labels) = accepting_states(spec)

# --------------------------
# User-defined deterministic monitor
# --------------------------
struct FunctionMonitor <: AbstractSpecStepper
    qa0::Int
    acc::Set{Int}
    stepfun::Function              # (qa::Int, ap::Tuple) -> qa2::Int
    aps::Vector{Symbol}
end

FunctionMonitor(qa0::Int, acc::Set{Int}, stepfun::Function) =
    FunctionMonitor(qa0, acc, stepfun, Symbol[])

step(M::FunctionMonitor, qa::Int, ap::Tuple{Vararg{Symbol}}) = M.stepfun(qa, ap)
init_state(M::FunctionMonitor) = M.qa0
accepting_states(M::FunctionMonitor) = M.acc

function labeling_function_from_state_sets(lab_abs::Dict{Symbol, <:AbstractVector{Int}})
    nmax = 0
    for states in values(lab_abs)
        isempty(states) || (nmax = max(nmax, maximum(states)))
    end

    lab_bits = Dict{Symbol, BitVector}()
    for (ap, states) in lab_abs
        bits = falses(nmax)
        for q in states
            bits[q] = true
        end
        lab_bits[ap] = bits
    end

    aps = collect(keys(lab_bits))

    return function (qs::Int)
        true_aps = Symbol[]
        for ap in aps
            bits = lab_bits[ap]
            if qs <= length(bits) && bits[qs]
                push!(true_aps, ap)
            end
        end
        return Tuple(true_aps)
    end
end

# ============================================================
# Product automaton (System × SpecStepper)
# ============================================================

struct ProductAutomaton{SYS, LAB, STEP} <: SY.AbstractAutomatonList
    sys::SYS
    labeling::LAB
    spec::STEP

    pid::Dict{Tuple{Int, Int}, Int}
    rev::Vector{Tuple{Int, Int}}

    postmap::Vector{Vector{Int}}              # pair_id(p,u) -> successors
    premap::Vector{Vector{Tuple{Int, Int}}}   # target -> (source, input)

    ninput::Int
end

_product_pair_id(P::ProductAutomaton, p::Int, u::Int) = (p - 1) * P.ninput + u

SY.get_n_state(P::ProductAutomaton) = length(P.rev)
SY.get_n_input(P::ProductAutomaton) = P.ninput
SY.pre(P::ProductAutomaton, t::Int) = P.premap[t]
SY.post(P::ProductAutomaton, p::Int, u::Int) = P.postmap[_product_pair_id(P, p, u)]

SY.enum_transitions(P::ProductAutomaton) = begin
    trans = Tuple{Int, Int, Int}[]
    for p in 1:SY.get_n_state(P)
        for u in 1:SY.get_n_input(P)
            for p2 in SY.post(P, p, u)
                push!(trans, (p2, p, u))
            end
        end
    end
    return trans
end

# build the product automaton only from the reachable states starting from the initial set:
function build_product_automaton(
    sys::SY.AbstractAutomatonList,
    spec::AbstractSpecStepper,
    labeling;
    initial_set = 1:SY.get_n_state(sys),
)
    nS = SY.get_n_state(sys)
    nU = SY.get_n_input(sys)
    qa0 = init_state(spec)

    labels = Vector{Any}(undef, nS)
    for qs in 1:nS
        labels[qs] = labeling(qs)
    end

    pid = Dict{Tuple{Int, Int}, Int}()
    rev = Tuple{Int, Int}[]

    postmap = Vector{Vector{Int}}()
    premap = Vector{Vector{Tuple{Int, Int}}}()

    seen_succ = Int[]
    discovered = Bool[]

    function add_product_state!(qs::Int, qa::Int)
        push!(rev, (qs, qa))

        for _ in 1:nU
            push!(postmap, Int[])
        end

        push!(premap, Tuple{Int, Int}[])
        push!(seen_succ, 0)
        push!(discovered, false)

        return length(rev)
    end

    function getpid!(qs::Int, qa::Int)
        return get!(pid, (qs, qa)) do
            return add_product_state!(qs, qa)
        end
    end

    work = Int[]

    for qs in initial_set
        qa_init = step(spec, qa0, labels[qs])
        p0 = getpid!(qs, qa_init)

        if !discovered[p0]
            discovered[p0] = true
            push!(work, p0)
        end
    end

    stamp = 0

    i = 1
    while i <= length(work)
        p = work[i]
        i += 1

        qs, qa = rev[p]

        for u in 1:nU
            pair_id = (p - 1) * nU + u
            succs = postmap[pair_id]
            empty!(succs)

            stamp += 1

            for qs2 in SY.post(sys, qs, u)
                qa2 = step(spec, qa, labels[qs2])
                p2 = getpid!(qs2, qa2)

                if seen_succ[p2] != stamp
                    seen_succ[p2] = stamp
                    push!(succs, p2)
                    push!(premap[p2], (p, u))
                end

                if !discovered[p2]
                    discovered[p2] = true
                    push!(work, p2)
                end
            end
        end
    end

    return ProductAutomaton(sys, labeling, spec, pid, rev, postmap, premap, nU)
end

# Tabulate the closed-loop memory controller: the spec stepper (which may wrap a
# non-serializable monitor) is evaluated here, once per (memory state, label),
# so the returned controller is plain data (JLD2-serializable).
function build_dynamic_controller(
    lab_abs::Dict{Symbol, Vector{Int}},
    spec::AbstractSpecStepper,
    contrP::ST.AbstractDiscreteController,
    pid::Dict{Tuple{Int, Int}, Int},
)
    nmax = 0
    for states in values(lab_abs)
        isempty(states) || (nmax = max(nmax, maximum(states)))
    end

    lab_bits = Dict{Symbol, BitVector}()

    for (ap, states) in lab_abs
        bits = falses(nmax)
        for q in states
            bits[q] = true
        end
        lab_bits[ap] = bits
    end

    aps = collect(keys(lab_bits))
    qa0 = init_state(spec)

    labeling_qs = function (qs::Int)
        trues = Symbol[]

        for ap in aps
            bits = lab_bits[ap]
            if qs <= length(bits) && bits[qs]
                push!(trues, ap)
            end
        end

        return Tuple(trues)
    end

    # Intern the distinct label tuples of 1:nmax (plus the empty label used for
    # any state outside the labeled range).
    label_ids = Dict{Tuple{Vararg{Symbol}}, Int}(() => 1)
    label_tuples = Tuple{Vararg{Symbol}}[()]
    label_of_state = Vector{Int}(undef, nmax)
    for qs in 1:nmax
        lab = labeling_qs(qs)
        label_of_state[qs] = get!(label_ids, lab) do
            push!(label_tuples, lab)
            return length(label_tuples)
        end
    end
    default_label = label_ids[()]

    # Tabulate the spec automaton over the memory states reachable in the
    # product plus every state they can step to.
    qas = Set{Int}([qa0])
    for (_, qa) in keys(pid)
        push!(qas, qa)
    end
    step_map = Dict{Tuple{Int, Int}, Int}()
    work = collect(qas)
    while !isempty(work)
        qa = pop!(work)
        for (lid, lab) in enumerate(label_tuples)
            qa2 = step(spec, qa, lab)
            step_map[(qa, lid)] = qa2
            if !(qa2 in qas)
                push!(qas, qa2)
                push!(work, qa2)
            end
        end
    end

    return ST.AutomatonMemoryController(
        qa0,
        label_of_state,
        default_label,
        step_map,
        pid,
        contrP,
    )
end

function project_initial_memory_controllable_set(
    product_autom,
    product_controllable_set,
    spec,
    labeling,
    nsys,
)
    Wsys = Set{Int}()
    Wprod = Set(product_controllable_set)
    qa0 = init_state(spec)

    for qs in 1:nsys
        qa_init = step(spec, qa0, labeling(qs))
        p = get(product_autom.pid, (qs, qa_init), nothing)
        if p !== nothing && p in Wprod
            push!(Wsys, qs)
        end
    end
    return Wsys
end

function project_initial_memory_uncontrollable_set(
    product_autom,
    product_controllable_set,
    spec,
    labeling,
    nsys,
)
    Wsys = project_initial_memory_controllable_set(
        product_autom,
        product_controllable_set,
        spec,
        labeling,
        nsys,
    )
    Usys = Set(1:nsys)
    for qs in Wsys
        delete!(Usys, qs)
    end
    return Usys
end

function project_initial_memory_value_function(
    product_autom,
    product_value_fun_tab,
    spec,
    labeling,
    nsys::Int,
)
    Vsys = Dict{Int, Float64}()
    qa0 = init_state(spec)

    for qs in 1:nsys
        qa_init = step(spec, qa0, labeling(qs))
        p = get(product_autom.pid, (qs, qa_init), nothing)
        p === nothing && continue

        vp = product_value_fun_tab[p]
        vp isa Real || continue

        Vsys[qs] = Float64(vp)
    end

    return Vsys
end
