# ============================================================
# CoSafe LTL Control
# ============================================================

import Spot

mutable struct OptimizerCoSafeLTLProblem{T} <: MOI.AbstractOptimizer
    # inputs
    problem::Union{Nothing, PR.CoSafeLTLProblem}
    early_stop::Bool
    sparse_input::Bool
    print_level::Int

    # outputs / internals
    controller::Union{Nothing, MS.SystemWithOutput}
    qa0::Any
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
            nothing,   # qa0
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

function MOI.set(opt::OptimizerCoSafeLTLProblem, p::MOI.RawOptimizerAttribute, v)
    return setproperty!(opt, Symbol(p.name), v)
end

function MOI.get(opt::OptimizerCoSafeLTLProblem, p::MOI.RawOptimizerAttribute)
    return getproperty(opt, Symbol(p.name))
end

MOI.get(opt::OptimizerCoSafeLTLProblem, ::MOI.SolveTimeSec) = opt.solve_time_sec

function MOI.optimize!(optimizer::OptimizerCoSafeLTLProblem)
    t0 = time()

    problem = optimizer.problem

    problem === nothing && error("problem not set")
    _check_ap_coverage(problem)

    if optimizer.print_level > 0
        for (ap, states) in problem.labeling
            println("AP=$ap  #states=$(length(states))")
        end
    end

    # (1) Finite system
    autom = problem.system

    # (2) Specification automaton
    spec = if problem.spec isa Spot.SpotFormula
        spot_stepper(problem.spec)
    elseif problem.spec isa AbstractSpecStepper
        problem.spec
    else
        error("spec must be SpotFormula or AbstractSpecStepper")
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
        optimizer.early_stop ? problem.initial_set : collect(1:get_n_state(autom))
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

    accQ = accepting_states(spec)
    labels_seen = Set{Any}()
    for qs in 1:get_n_state(autom)
        push!(labels_seen, labeling(qs))
    end
    product_target_set =
        [p for p in 1:get_n_state(product_autom) if product_autom.rev[p][2] in accQ]
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
    product_uncontrollable_set = product_automaton_optimizer.uncontrollable_set
    product_value_fun_tab = product_automaton_optimizer.value_fun_tab
    success = product_automaton_optimizer.success

    problem.labeling isa AbstractDict ||
        error("build_fm_controller_ms currently requires dictionary labeling.")

    # (5) Wrap product controller into finite-memory controller on autom
    controller, qa0 = build_fm_controller_ms(
        problem.labeling,
        spec,
        product_controller,
        product_autom.pid,
    )

    optimizer.controller = controller
    optimizer.qa0 = qa0
    optimizer.update_on_next = true
    optimizer.controllable_set = project_initial_memory_controllable_set(
        product_autom,
        product_controllable_set,
        spec,
        labeling,
        get_n_state(autom),
    )
    optimizer.uncontrollable_set = project_initial_memory_uncontrollable_set(
        product_autom,
        product_controllable_set,
        spec,
        labeling,
        get_n_state(autom),
    )
    optimizer.value_fun_tab = project_initial_memory_value_function(
        product_autom,
        product_value_fun_tab,
        spec,
        labeling,
        get_n_state(autom),
    )
    optimizer.success = success
    optimizer.print_level >= 1 && println("Success: ", success)
    optimizer.solve_time_sec = time() - t0

    return
end

# ============================================================
# Helpers
# ============================================================

# Small check: if the spec is a Spot formula, ensure all APs mentioned in the formula
# have a set provided in `concrete_problem.labeling`.
# Extra keys in the dict are fine.)
collect_aps(φ::Spot.SpotFormula) = [Symbol(ap) for ap in Spot.atomic_prop_collect(φ)]

function _check_ap_coverage(problem)
    if !(problem.spec isa Spot.SpotFormula)
        return nothing
    end

    problem.labeling isa AbstractDict || return nothing

    aps = Set(collect_aps(problem.spec))
    provided = Set(keys(problem.labeling))
    missing = setdiff(aps, provided)
    if !isempty(missing)
        error("Missing AP sets in problem.labeling for: $(collect(missing)).")
    end
    return nothing
end

# ============================================================
# Spec interface (either Spot DRA or user-defined monitor)
# ============================================================

abstract type AbstractSpecStepper end

# Step a spec automaton/monitor: qa_next = step(spec, qa, ap_tuple)
# Must return an Int state id.
step(::AbstractSpecStepper, ::Int, ::Tuple{Vararg{Symbol}}) = error("step not implemented")

# Optional: return the initial spec state qa0.
init_state(::AbstractSpecStepper) = 1

# Optional: return the set of "done/accepting" spec states for co-safe reduction.
# If not provided, we will attempt to compute it (Spot) or require the user to provide.
accepting_states(::AbstractSpecStepper) = error("accepting_states not implemented")

# --------------------------
# User-defined deterministic monitor
# --------------------------
struct FunctionMonitor <: AbstractSpecStepper
    qa0::Int
    acc::Set{Int}
    stepfun::Function  # (qa::Int, ap::Tuple) -> qa2::Int
end

step(M::FunctionMonitor, qa::Int, ap::Tuple{Vararg{Symbol}}) = M.stepfun(qa, ap)
init_state(M::FunctionMonitor) = M.qa0
accepting_states(M::FunctionMonitor) = M.acc

# --------------------------
# Spot-based stepper (DRA used as a monitor)
# --------------------------
struct SpotDRAstepper <: AbstractSpecStepper
    φ::Spot.SpotFormula
    dra::Any
    qa0::Int
    qa_dead::Int
    aps::Vector{Symbol}  # AP universe used for done-state detection
end

init_state(S::SpotDRAstepper) = S.qa0

@inline function _nextstate_int(dra, qa::Int, ap::Tuple, qa_dead::Int)
    qa2 = Spot.nextstate(dra, qa, ap)
    return qa2 === nothing ? qa_dead : qa2
end

step(S::SpotDRAstepper, qa::Int, ap::Tuple{Vararg{Symbol}}) =
    (qa == S.qa_dead) ? S.qa_dead : _nextstate_int(S.dra, qa, ap, S.qa_dead)

# Heuristic done-state detector for Spot DRA used as a co-safe monitor.
# A state q is considered "done" if for all valuations v:
#   nextstate(dra, q, v) is either q OR nothing.
# This is a practical criterion that works well for many co-safe task specs.
function _all_valuations(aps::Vector{Symbol})
    n = length(aps)
    vals = Vector{Tuple{Vararg{Symbol}}}(undef, 1 << n)
    for mask in 0:((1 << n) - 1)
        trues = Symbol[]
        for (i, a) in enumerate(aps)
            if (mask >> (i-1)) & 1 == 1
                push!(trues, a)
            end
        end
        vals[mask + 1] = Tuple(trues)
    end
    return vals
end

function _cosafe_done_states_dra(dra; aps::Vector{Symbol}, q0::Int = 1)
    vals = _all_valuations(aps)

    # reachable spec states (ignore `nothing`)
    reachable = Set{Int}([q0])
    queue = [q0]
    while !isempty(queue)
        q = popfirst!(queue)
        for v in vals
            q2 = Spot.nextstate(dra, q, v)
            q2 === nothing && continue
            if !(q2 in reachable)
                push!(reachable, q2)
                push!(queue, q2)
            end
        end
    end

    done = Set{Int}()
    for q in reachable
        ok = true
        for v in vals
            q2 = Spot.nextstate(dra, q, v)
            if !(q2 === nothing || q2 == q)
                ok = false
                break
            end
        end
        ok && push!(done, q)
    end
    return done, reachable
end

accepting_states(S::SpotDRAstepper) = begin
    doneQ, _ = _cosafe_done_states_dra(S.dra; aps = S.aps, q0 = S.qa0)
    isempty(doneQ) && error(
        "Could not find any 'done' spec states with the co-safe heuristic. " *
        "The formula may not be co-safe for this pipeline, or the AP set is mismatched. " *
        "Try providing a FunctionMonitor with explicit accepting states.",
    )
    doneQ
end

function spot_stepper(
    φ::Spot.SpotFormula;
    qa0::Int = 1,
    qa_dead::Int = 0,
    aps::Union{Nothing, Vector{Symbol}} = nothing,
)
    aps_use = aps === nothing ? collect_aps(φ) : aps
    try
        dra = Spot.DeterministicRabinAutomata(φ)
        return SpotDRAstepper(φ, dra, qa0, qa_dead, aps_use)
    catch
        error(
            "Spot could not provide a DeterministicRabinAutomata with nextstate(). " *
            "Provide a FunctionMonitor wrapper.",
        )
    end
end

function labeling_function_from_state_sets(lab_abs::Dict{Symbol, <:AbstractVector{Int}})
    # speed: BitSet membership is cheap
    lab_bits = Dict{Symbol, BitSet}()
    for (ap, states) in lab_abs
        lab_bits[ap] = BitSet(states)
    end
    aps = collect(keys(lab_bits))

    return function (qs::Int)
        true_aps = Symbol[]
        for ap in aps
            (qs in lab_bits[ap]) && push!(true_aps, ap)
        end
        return Tuple(true_aps)
    end
end

# ============================================================
# Product automaton (System × SpecStepper)
# ============================================================

struct ProductAutomaton{SYS, LAB, STEP} <: AbstractAutomatonList{0, 0}
    sys::SYS
    labeling::LAB
    spec::STEP

    pid::Dict{Tuple{Int, Int}, Int}
    rev::Vector{Tuple{Int, Int}}
    post_tab::Vector{Vector{Vector{Int}}}
    pre_tab::Vector{Vector{Tuple{Int, Int}}}
    ninput::Int
end

get_n_state(P::ProductAutomaton) = length(P.rev)
get_n_input(P::ProductAutomaton) = P.ninput
pre(P::ProductAutomaton, t::Int) = P.pre_tab[t]
post(P::ProductAutomaton, s::Int, u::Int) = P.post_tab[s][u]

function enum_transitions(P::ProductAutomaton)
    trans = Tuple{Int, Int, Int}[]
    for q in 1:get_n_state(P), u in 1:get_n_input(P)
        for q2 in post(P, q, u)
            push!(trans, (q2, q, u))
        end
    end
    return trans
end

# build the product automaton only from the reachable states starting from the initial set:
function build_product_automaton(
    sys::AbstractAutomatonList,
    spec::AbstractSpecStepper,
    labeling;
    initial_set = 1:get_n_state(sys),
)
    nU = get_n_input(sys)
    qa0 = init_state(spec)

    pid = Dict{Tuple{Int, Int}, Int}()
    rev = Tuple{Int, Int}[]

    getpid(qs::Int, qa::Int) =
        get!(pid, (qs, qa)) do
            push!(rev, (qs, qa))
            return length(rev)
        end

    # BFS
    work = Int[]
    inqueue = BitSet()
    for qs in initial_set
        ap0 = labeling(qs)
        qa_init = step(spec, qa0, ap0)
        p0 = getpid(qs, qa_init)
        push!(work, p0)
        push!(inqueue, p0)
    end

    i = 1
    while i <= length(work)
        p = work[i];
        i += 1
        (qs, qa) = rev[p]
        for u in 1:nU
            for qs2 in post(sys, qs, u)
                ap = labeling(qs2)
                qa2 = step(spec, qa, ap)
                p2 = getpid(qs2, qa2)
                if !(p2 in inqueue)
                    push!(work, p2);
                    push!(inqueue, p2)
                end
            end
        end
    end

    nP = length(rev)
    post_tab = [[Int[] for _ in 1:nU] for _ in 1:nP]
    pre_tab = [Tuple{Int, Int}[] for _ in 1:nP]

    for p in 1:nP
        (qs, qa) = rev[p]
        for u in 1:nU
            succs = Int[]
            for qs2 in post(sys, qs, u)
                ap = labeling(qs2)
                qa2 = step(spec, qa, ap)
                p2 = pid[(qs2, qa2)]
                push!(succs, p2)
                push!(pre_tab[p2], (p, u))
            end
            post_tab[p][u] = unique(succs)
        end
    end
    return ProductAutomaton(sys, labeling, spec, pid, rev, post_tab, pre_tab, nU)
end

function build_fm_controller_ms(
    lab_abs::Dict{Symbol, Vector{Int}},
    spec::AbstractSpecStepper,
    contrP::MS.AbstractMap,
    pid::Dict{Tuple{Int, Int}, Int},
)
    # fast labeling_qs(qs) -> Tuple{Symbol,...}
    lab_bits = Dict{Symbol, BitSet}((ap => BitSet(states)) for (ap, states) in lab_abs)
    aps = collect(keys(lab_bits))

    labeling_qs = function (qs::Int)
        trues = Symbol[]
        for ap in aps
            (qs in lab_bits[ap]) && push!(trues, ap)
        end
        return Tuple(trues)
    end

    qa0 = init_state(spec)

    kP = contrP.h

    # output map: (qa, qs) -> u_sym (or Vector) or nothing
    h = function (qa::Int, qs::Int)
        p = get(pid, (qs, qa), nothing)
        p === nothing && return nothing
        return kP(p)
    end

    # memory update: (qa, qs_for_update) -> qa_next
    g = function (qa::Int, qs_for_update::Int)
        return step(spec, qa, labeling_qs(qs_for_update))
    end

    # Domain for output map: defined iff pid exists and product-controller defined
    isdef_qaqs = function (qaqs)
        qa, qs = qaqs
        p = get(pid, (qs, qa), nothing)
        p === nothing && return false
        out = kP(p)
        out === nothing && return false
        (out isa AbstractVector) && return !isempty(out)
        return true
    end
    X_qaqs = PredicateDomain(isdef_qaqs)

    # Build MS objects
    # outputmap expects input (qa,qs); we pack into a 2-tuple input for ConstrainedBlackBoxMap
    outmap = MS.ConstrainedBlackBoxMap(2, 1, qaqs -> begin
        qa, qs = qaqs
        h(qa, qs)
    end, X_qaqs)

    memsys = MS.BlackBoxDiscreteSystem((qa, qs_for_update) -> g(qa, qs_for_update), 1)

    return MS.SystemWithOutput(memsys, outmap), qa0
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
