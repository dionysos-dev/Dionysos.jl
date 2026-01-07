# ============================================================
# Optimizer
# ============================================================

mutable struct OptimizerCoSafeLTLProblem{T} <: MOI.AbstractOptimizer
    # inputs
    concrete_problem::Union{Nothing, PR.CoSafeLTLProblem}
    abstract_system::Union{Nothing, SY.SymbolicModelList}

    # outputs / internals
    abstract_problem::Union{Nothing, PR.CoSafeLTLProblem}
    product_system::Any
    abstract_controller_product::Union{Nothing, MS.ConstrainedBlackBoxMap}
    abstract_controller::Union{Nothing, MS.SystemWithOutput}
    qa0::Any # abstract/concrete controller initial state
    update_on_next::Bool
    controllable_set::Any
    uncontrollable_set::Any
    value_fun_tab::Any
    abstract_problem_time_sec::T
    success::Bool
    print_level::Int
    sparse_input::Bool

    function OptimizerCoSafeLTLProblem{T}() where {T}
        return new{T}(
            nothing, nothing,
            nothing,
            nothing,
            nothing, nothing, nothing,
            true,
            nothing, nothing,
            nothing,
            0.0,
            false,
            1,
            false,
        )
    end
end

OptimizerCoSafeLTLProblem() = OptimizerCoSafeLTLProblem{Float64}()

MOI.is_empty(opt::OptimizerCoSafeLTLProblem) = opt.concrete_problem === nothing

function MOI.set(opt::OptimizerCoSafeLTLProblem, p::MOI.RawOptimizerAttribute, v)
    return setproperty!(opt, Symbol(p.name), v)
end
function MOI.get(opt::OptimizerCoSafeLTLProblem, p::MOI.RawOptimizerAttribute)
    return getproperty(opt, Symbol(p.name))
end
MOI.get(opt::OptimizerCoSafeLTLProblem, ::MOI.SolveTimeSec) = opt.abstract_problem_time_sec

# ============================================================
# Lifting: concrete CoSafeLTLProblem -> abstract CoSafeLTLProblem
# concrete_problem.labeling :: Dict{Symbol, LazySet}
# abstract_problem.labeling  :: Dict{Symbol, Vector{Int}}
# ============================================================

"""
Small check: if the spec is a Spot formula, ensure all APs mentioned in the formula
have a set provided in `concrete_problem.labeling`.
(Extra keys in the dict are fine.)
"""
collect_aps(φ::Spot.SpotFormula) = [Symbol(ap) for ap in atomic_prop_collect(φ)]

function _check_ap_coverage(concrete_problem)
    if concrete_problem.spec isa Spot.SpotFormula
        aps = Set(collect_aps(concrete_problem.spec))
        provided = Set(keys(concrete_problem.labeling))
        missing = setdiff(aps, provided)
        if !isempty(missing)
            error("Missing AP sets in concrete_problem.labeling for: $(collect(missing)).")
        end
    end
    return nothing
end

function build_abstract_problem(
    concrete_problem::PR.CoSafeLTLProblem,
    abstract_system::SY.SymbolicModelList,
)
    _check_ap_coverage(concrete_problem)

    # lift initial set (conservative for initial condition is usually OUTER)
    init_states = SY.get_states_from_set(
        abstract_system,
        concrete_problem.initial_set,
        DO.OUTER,
    )

    # lift each AP set to a set of symbolic states
    lab_abs = Dict{Symbol, Vector{Int}}()
    for (ap, setX) in concrete_problem.labeling
        sem = get(concrete_problem.ap_semantics, ap, DO.INNER)
        lab_abs[ap] = SY.get_states_from_set(abstract_system, setX, sem)
    end

    # NOTE: adapt this constructor call to match your actual PR.CoSafeLTLProblem definition
    return PR.CoSafeLTLProblem(
        abstract_system,
        init_states,
        concrete_problem.spec,
        lab_abs,
        concrete_problem.ap_semantics,
        concrete_problem.strict_spot,
    )
end

function labeling_function_from_state_sets(lab_abs::Dict{Symbol,<:AbstractVector{Int}})
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
# optimize!
# ============================================================

function MOI.optimize!(optimizer::OptimizerCoSafeLTLProblem)
    t0 = time()

    optimizer.abstract_system === nothing && error("abstract_system not set")
    optimizer.concrete_problem === nothing && error("concrete_problem not set")

    concrete_problem = optimizer.concrete_problem
    abs_sys = optimizer.abstract_system

    # (1) Lift to abstract problem (same struct, labeling becomes Dict{Symbol,Vector{Int}})
    abstract_problem = build_abstract_problem(concrete_problem, abs_sys)
    for (ap, states) in abstract_problem.labeling
        println("AP=$ap  #states=$(length(states))")
    end
    optimizer.abstract_problem = abstract_problem

    # (2) Build spec stepper
    spec = if abstract_problem.spec isa Spot.SpotFormula
        spot_stepper(abstract_problem.spec)
    elseif abstract_problem.spec isa AbstractSpecStepper
        abstract_problem.spec
    else
        error("spec must be SpotFormula or AbstractSpecStepper")
    end

    # (3) Build abstract labeling function (qs -> Tuple{Symbol,...})
    abstract_labeling = labeling_function_from_state_sets(abstract_problem.labeling)

    # (4) Product automaton on the *abstract* automaton
    P = build_product_automaton(
        abs_sys.autom,
        spec,
        abstract_labeling;
        initial_set = abstract_problem.initial_set,
    )
    optimizer.product_system = P

    # (5) Targets from accepting memory states
    accQ = accepting_states(spec)
    target_set = [p for p in 1:SY.get_n_state(P) if P.rev[p][2] in accQ]
    isempty(target_set) && error("Empty target_set (AP mismatch or acceptance not found).")

    # (6) Initial product states
    init_states = abstract_problem.initial_set
    init_prod = [P.pid[(qs, init_state(spec))] for qs in init_states]

    controller_product_automaton, controllableP, uncontrollableP, V =
        SY.compute_worst_case_uniform_cost_controller(P, target_set; initial_set = init_prod)

    # (7) Wrap product controller into finite-memory controller on abstract system
    abstract_controller, qa0 = build_abstract_fm_controller_ms(abstract_problem.labeling, spec, controller_product_automaton, P.pid)

    optimizer.abstract_controller_product = controller_product_automaton
    optimizer.abstract_controller = abstract_controller
    optimizer.qa0 = qa0
    optimizer.controllable_set = controllableP
    optimizer.uncontrollable_set = uncontrollableP
    optimizer.value_fun_tab = V
    optimizer.success = all(p -> p in controllableP, init_prod)

    optimizer.abstract_problem_time_sec = time() - t0
    return
end


# ============================================================
# Everything below is your SpecStepper + ProductAutomaton code
# (kept as-is, only tiny fixes if needed)
# ============================================================

# ============================================================
# Spec interface (either Spot DRA or user-defined monitor)
# ============================================================

abstract type AbstractSpecStepper end

"""
Step a spec automaton/monitor:
  qa_next = step(spec, qa, ap_tuple)

Must return an Int state id.
"""
step(::AbstractSpecStepper, ::Int, ::Tuple{Vararg{Symbol}}) = error("step not implemented")

"""
Optional: return the initial spec state qa0.
"""
init_state(::AbstractSpecStepper) = 1

"""
Optional: return the set of "done/accepting" spec states for co-safe reduction.
If not provided, we will attempt to compute it (Spot) or require the user to provide.
"""
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
    dra
    qa0::Int
    qa_dead::Int
    aps::Vector{Symbol}  # AP universe used for done-state detection
end

init_state(S::SpotDRAstepper) = S.qa0

@inline function _nextstate_int(dra, qa::Int, ap::Tuple, qa_dead::Int)
    qa2 = nextstate(dra, qa, ap)
    return qa2 === nothing ? qa_dead : qa2
end

step(S::SpotDRAstepper, qa::Int, ap::Tuple{Vararg{Symbol}}) =
    (qa == S.qa_dead) ? S.qa_dead : _nextstate_int(S.dra, qa, ap, S.qa_dead)

"""
Heuristic done-state detector for Spot DRA used as a co-safe monitor.

A state q is considered "done" if for all valuations v:
  nextstate(dra, q, v) is either q OR nothing.

This is a practical criterion that works well for many co-safe task specs.
"""
function _all_valuations(aps::Vector{Symbol})
    n = length(aps)
    vals = Vector{Tuple{Vararg{Symbol}}}(undef, 1 << n)
    for mask in 0:(1<<n)-1
        trues = Symbol[]
        for (i, a) in enumerate(aps)
            if (mask >> (i-1)) & 1 == 1
                push!(trues, a)
            end
        end
        vals[mask+1] = Tuple(trues)
    end
    return vals
end

function _cosafe_done_states_dra(dra; aps::Vector{Symbol}, q0::Int=1)
    vals = _all_valuations(aps)

    # reachable spec states (ignore `nothing`)
    reachable = Set{Int}([q0])
    queue = [q0]
    while !isempty(queue)
        q = popfirst!(queue)
        for v in vals
            q2 = nextstate(dra, q, v)
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
            q2 = nextstate(dra, q, v)
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
    doneQ, _ = _cosafe_done_states_dra(S.dra; aps=S.aps, q0=S.qa0)
    isempty(doneQ) && error(
        "Could not find any 'done' spec states with the co-safe heuristic. " *
        "The formula may not be co-safe for this pipeline, or the AP set is mismatched. " *
        "Try providing a FunctionMonitor with explicit accepting states."
    )
    doneQ
end

collect_aps(φ::Spot.SpotFormula) = [Symbol(ap) for ap in atomic_prop_collect(φ)]

function spot_stepper(φ::Spot.SpotFormula; qa0::Int=1, qa_dead::Int=0, aps::Union{Nothing,Vector{Symbol}}=nothing)
    aps_use = aps === nothing ? collect_aps(φ) : aps
    try
        dra = DeterministicRabinAutomata(φ)
        return SpotDRAstepper(φ, dra, qa0, qa_dead, aps_use)
    catch
        error("Spot could not provide a DeterministicRabinAutomata with nextstate(). " *
              "Provide a FunctionMonitor wrapper.")
    end
end


# ============================================================
# Product automaton (System × SpecStepper)
# ============================================================

struct ProductAutomaton{SYS,LAB,STEP} <: SY.AbstractAutomatonList{0,0}
    sys::SYS
    labeling::LAB
    spec::STEP

    pid::Dict{Tuple{Int,Int},Int}
    rev::Vector{Tuple{Int,Int}}
    post_tab::Vector{Vector{Vector{Int}}}
    pre_tab::Vector{Vector{Tuple{Int,Int}}}
    ninput::Int
end

SY.get_n_state(P::ProductAutomaton) = length(P.rev)
SY.get_n_input(P::ProductAutomaton) = P.ninput
SY.pre(P::ProductAutomaton, t::Int) = P.pre_tab[t]
SY.post(P::ProductAutomaton, s::Int, u::Int) = P.post_tab[s][u]

function SY.enum_transitions(P::ProductAutomaton)
    trans = Tuple{Int,Int,Int}[]
    for q in 1:SY.get_n_state(P), u in 1:SY.get_n_input(P)
        for q2 in SY.post(P, q, u)
            push!(trans, (q2, q, u))
        end
    end
    return trans
end

function build_product_automaton(sys::SY.AbstractAutomatonList, spec::AbstractSpecStepper, labeling;
    initial_set = 1:SY.get_n_state(sys),
)
    nU = SY.get_n_input(sys)
    qa0 = init_state(spec)

    pid = Dict{Tuple{Int,Int},Int}()
    rev = Tuple{Int,Int}[]

    getpid(qs::Int, qa::Int) = get!(pid, (qs, qa)) do
        push!(rev, (qs, qa))
        length(rev)
    end

    # BFS
    work = Int[]
    inqueue = BitSet()
    for qs in initial_set
        p0 = getpid(qs, qa0)
        push!(work, p0); push!(inqueue, p0)
    end

    i = 1
    while i <= length(work)
        p = work[i]; i += 1
        (qs, qa) = rev[p]
        for u in 1:nU
            for qs2 in SY.post(sys, qs, u)
                ap  = labeling(qs2)
                qa2 = step(spec, qa, ap)
                p2  = getpid(qs2, qa2)
                if !(p2 in inqueue)
                    push!(work, p2); push!(inqueue, p2)
                end
            end
        end
    end

    nP = length(rev)
    post_tab = [[Int[] for _ in 1:nU] for _ in 1:nP]
    pre_tab  = [Tuple{Int,Int}[] for _ in 1:nP]

    for p in 1:nP
        (qs, qa) = rev[p]
        for u in 1:nU
            succs = Int[]
            for qs2 in SY.post(sys, qs, u)
                ap  = labeling(qs2)
                qa2 = step(spec, qa, ap)
                p2  = pid[(qs2, qa2)]
                push!(succs, p2)
                push!(pre_tab[p2], (p, u))
            end
            post_tab[p][u] = unique(succs)
        end
    end
    return ProductAutomaton(sys, labeling, spec, pid, rev, post_tab, pre_tab, nU)
end

"""
Build abstract finite-memory controller as MS.SystemWithOutput.

Inputs:
- lab_abs :: Dict{Symbol, Vector{Int}}   (abstract labeling as state sets)
- spec    :: AbstractSpecStepper         (your step(spec, qa, ap_tuple))
- contrP  :: MS.AbstractMap              (controller on product states p -> u_sym or Vector{Int})
- pid     :: Dict{Tuple{Int,Int},Int}    ((qs,qa) -> p)

Output:
- Cabs :: MS.SystemWithOutput where:
    outputmap: (qa, qs) -> u_sym (or Vector{Int}) or nothing
    memsys:    (qa, qs_for_update) -> qa_next
- qa0  :: Int (initial abstract control state)
"""
function build_abstract_fm_controller_ms(
    lab_abs::Dict{Symbol, Vector{Int}},
    spec::AbstractSpecStepper,
    contrP::MS.AbstractMap,
    pid::Dict{Tuple{Int,Int},Int},
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