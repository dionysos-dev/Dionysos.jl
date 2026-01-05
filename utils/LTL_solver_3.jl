using Plots
using Spot
using Dionysos
const DI = Dionysos
const SY = DI.Symbolic
const ST = DI.System

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


# ============================================================
# Finite-memory controller on the ORIGINAL system
# (wrap product controller + spec memory)
# ============================================================

struct FiniteMemoryLTLController{SYS,LAB,STEP,CP}
    sys::SYS
    labeling::LAB
    spec::STEP
    contrP::CP
    pid::Dict{Tuple{Int,Int},Int}
end

"""
Get one allowed control from a SymbolicControllerList for a given product state p.
This is robust to different controller APIs: it tries ST.get_control, then falls back.
"""
function _get_one_control(contrP, p::Int)
    # Try common getter
    if hasmethod(ST.get_control, Tuple{typeof(contrP), Int})
        return ST.get_control(contrP, p)
    end
    # Fallback: if controls are stored as (state,input) pairs in `contrP.data`
    if hasproperty(contrP, :data)
        data = getproperty(contrP, :data)
        # data likely contains tuples (p,u)
        for (pp, u) in data
            if pp == p
                return u
            end
        end
    end
    error("Don't know how to extract a control from $(typeof(contrP)). Please adapt _get_one_control.")
end

function control(ctrl::FiniteMemoryLTLController, qs::Int, qa::Int)
    p = get(ctrl.pid, (qs, qa), nothing)
    p === nothing && return nothing
    return _get_one_control(ctrl.contrP, p)
end

@inline function update_memory(ctrl::FiniteMemoryLTLController, qa::Int, qs_next::Int)
    ap = ctrl.labeling(qs_next)
    return step(ctrl.spec, qa, ap)
end

function controller_for_original_system(P::ProductAutomaton, contrP)
    return FiniteMemoryLTLController(P.sys, P.labeling, P.spec, contrP, P.pid)
end

function rollout_original(sys, ctrl::FiniteMemoryLTLController; qs0::Int=1, horizon::Int=10, succ_rule::Symbol=:first)
    qs = qs0
    qa = init_state(ctrl.spec)

    println("Rollout:")
    for k in 1:horizon
        u = control(ctrl, qs, qa)
        if u === nothing
            println("  step=$k  qs=$qs  qa=$qa  -> no control (losing)")
            return
        end
        println("  step=$k  qs=$qs  qa=$qa  -> u=$u")

        succ = SY.post(sys, qs, u)
        isempty(succ) && (println("  dead-end at qs=$qs u=$u"); return)
        qs_next = (succ_rule == :first) ? first(succ) : succ[rand(1:length(succ))]

        qa = update_memory(ctrl, qa, qs_next)
        qs = qs_next
    end
end


# ============================================================
# One high-level function: solve co-safe LTL on an abstract TS
# ============================================================

"""
Solve a co-safe LTL problem on a finite nondeterministic TS using a product + reachability solver.

You can provide either:
- spec = Spot formula (Spot.SpotFormula)  -> Spot-based stepper (may fail for some formulas)
- spec = FunctionMonitor(...)            -> user-defined deterministic monitor (always works)

Returns:
- P: product automaton
- contrP: product controller
- ctrl: finite-memory controller on original system
- V: value function on product
"""
function solve_cosafe_ltl(sys::SY.AbstractAutomatonList, spec_in, labeling;
    initial_set = 1:SY.get_n_state(sys),
    strict_spot::Bool=false,
)
    # Build spec stepper
    spec = if spec_in isa Spot.SpotFormula
        φ = spec_in
        # Optional conservative check
        if strict_spot && !(is_reachability(φ) || is_constrained_reachability(φ))
            error("Spot formula not in recognized reachability/constrained-reachability fragment (strict_spot=true).")
        end
        spot_stepper(φ)  # may throw if Spot can't determinize
    elseif spec_in isa AbstractSpecStepper
        spec_in
    else
        error("spec_in must be a SpotFormula or an AbstractSpecStepper (e.g., FunctionMonitor).")
    end

    # Product
    P = build_product_automaton(sys, spec, labeling; initial_set=initial_set)

    # Targets from accepting spec states
    accQ = accepting_states(spec)
    target_set = [p for p in 1:SY.get_n_state(P) if P.rev[p][2] in accQ]
    isempty(target_set) && error("Empty target_set. Check labeling/AP consistency or provide explicit monitor acceptance.")

    # Initial product ids
    init_prod = [P.pid[(qs, init_state(spec))] for qs in initial_set]

    contrP, controllableP, uncontrollableP, V =
        SY.compute_worst_case_uniform_cost_controller(P, target_set; initial_set=init_prod)

    ctrl = controller_for_original_system(P, contrP)
    return P, contrP, ctrl, V
end






struct ToyAutomatonDanger <: SY.AbstractAutomatonList{5,2}
    post_tab::Vector{Vector{Vector{Int}}}
    pre_tab::Vector{Vector{Tuple{Int,Int}}}
end

function ToyAutomatonDanger()
    nS, nU = 5, 2
    post_tab = [[Int[] for _ in 1:nU] for _ in 1:nS]
    pre_tab  = [Tuple{Int,Int}[] for _ in 1:nS]

    add!(q, u, q2) = (push!(post_tab[q][u], q2); push!(pre_tab[q2], (q, u)))

    # Same base structure:
    add!(1, 1, 2); add!(1, 2, 4);
    add!(2, 1, 3); add!(2, 2, 1)

    # Key change: from g1 (state 3), input 2 can go to danger (2) or obstacle (5)
    add!(3, 1, 5)
    add!(3, 2, 2); 

    add!(4, 1, 4); add!(4, 2, 4)
    add!(5, 1, 1); # add!(5, 1, 2)

    return ToyAutomatonDanger(post_tab, pre_tab)
end



SY.get_n_state(::ToyAutomatonDanger) = 5
SY.get_n_input(::ToyAutomatonDanger) = 2
SY.post(A::ToyAutomatonDanger, q::Int, u::Int) = A.post_tab[q][u]
SY.pre(A::ToyAutomatonDanger, q::Int) = A.pre_tab[q]
SY.enum_transitions(A::ToyAutomatonDanger) =
    [(q2, q, u) for q in 1:5, u in 1:2 for q2 in SY.post(A, q, u)]


struct ProductAutomatonGeneric{SYS,LAB} <: SY.AbstractAutomatonList{0,0}
    sys::SYS
    labeling::LAB
    pid::Dict{Tuple{Int,Int},Int}
    rev::Vector{Tuple{Int,Int}}
    post_tab::Vector{Vector{Vector{Int}}}
    pre_tab::Vector{Vector{Tuple{Int,Int}}}
    ninput::Int
end

SY.get_n_state(P::ProductAutomatonGeneric) = length(P.rev)
SY.get_n_input(P::ProductAutomatonGeneric) = P.ninput
SY.pre(P::ProductAutomatonGeneric, t::Int) = P.pre_tab[t]
SY.post(P::ProductAutomatonGeneric, s::Int, u::Int) = P.post_tab[s][u]


"""
Deterministic monitor for a specific co-safe pattern:
  G(!obs) ∧ F( g1 ∧ G(!danger) ∧ F g2 )

States:
  0 = dead
  1 = before g1
  2 = after g1, waiting for g2 (danger forbidden)
  3 = done
"""
struct MonitorG1NoDangerUntilG2 end

@inline function mon_next(::MonitorG1NoDangerUntilG2, q::Int, ap::Tuple{Vararg{Symbol}})
    obs    = (:obs in ap)
    g1     = (:g1 in ap)
    g2     = (:g2 in ap)
    danger = (:danger in ap)

    obs && return 0                 # always forbidden

    q == 0 && return 0              # dead absorbs
    q == 3 && return 3              # done absorbs

    if q == 1
        # before g1
        if g1
            # reaching g1 starts phase2; if g2 is already true, we can finish immediately
            return g2 ? 3 : 2
        else
            return 1
        end
    end

    @assert q == 2
    # after g1: danger forbidden until g2
    if g2
        return 3
    elseif danger
        return 0
    else
        return 2
    end
end

sys = ToyAutomatonDanger()

##  Example monitor: accepting state is 3, initial is 1
# function labeling_toy_danger(qs::Int)
#     ap = Symbol[]
#     qs == 3 && push!(ap, :g1)
#     qs == 1 && push!(ap, :g2)
#     qs == 4 && push!(ap, :obs)
#     qs == 2 && push!(ap, :danger)   # <-- danger region
#     return Tuple(ap)
# end
# mon = FunctionMonitor(
#     1,                      # qa0
#     Set([3]),               # accepting spec states
#     (qa, ap) -> mon_next(MonitorG1NoDangerUntilG2(), qa, ap),
# )

# P, contrP, ctrl, V = solve_cosafe_ltl(sys, mon, labeling_toy_danger; initial_set=[1])
# rollout_original(sys, ctrl; qs0=1, horizon=10)




## Example SPOT:
φ = ltl"!obs U g2"
φ = ltl"G(!obs) & F(g1 & F(g2))"
# φ = ltl"G(!obs) & F( g1 & ( (!danger) U g2 ) )"

function labeling_toy(qs::Int)
    ap = Symbol[]
    qs == 3 && push!(ap, :g1)
    qs == 1 && push!(ap, :g2)
    qs == 4 && push!(ap, :obs)
    # qs == 2 && push!(ap, :danger)   # <-- danger region
    return Tuple(ap)
end
P, contrP, ctrl, V = solve_cosafe_ltl(sys, φ, labeling_toy; initial_set=[1], strict_spot=false)
rollout_original(sys, ctrl; qs0=1, horizon=10)

# using TikzPictures
# translator = LTLTranslator(buchi=true, deterministic=true, state_based_acceptance=true)
# surveillance_aut = Spot.translate(translator, φ)
# pic = plot_automata(surveillance_aut)
# save(PDF("test"), pic)

