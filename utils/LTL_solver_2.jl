using Plots
using Spot
using Dionysos
const DI = Dionysos
const SY = DI.Symbolic
const ST = DI.System

# ============================================================
# Utilities
# ============================================================

"""
Safe wrapper around Spot.nextstate:
- returns an Int state id
- maps `nothing` to `qa_dead`
"""
@inline function nextstate_int(dra, q::Int, v::Tuple, q_dead::Int)
    q2 = nextstate(dra, q, v)
    return q2 === nothing ? q_dead : q2
end

"""
All valuations over a set of atomic propositions `aps`.
Each valuation is a Tuple{Vararg{Symbol}} of props that are true.
"""
function all_valuations(aps::Vector{Symbol})
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

"""
Heuristic "accepting sink" detector for Spot.jl DRA used as a co-safe monitor.

A state q is considered "done" if for all valuations v:
  nextstate(dra, q, v) is either q OR nothing.

This avoids the too-strict requirement δ(q,v)=q for all v and accommodates Spot.jl's `nothing`.
"""
function cosafe_done_states_dra(dra; aps::Vector{Symbol}, q0::Int=1)
    vals = all_valuations(aps)

    # BFS over reachable spec states (including `nothing` transitions by ignoring them)
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


# ============================================================
# Product automaton TS × Spec-automaton (co-safe monitor)
# ============================================================

"""
Product automaton of:
- Dionysos nondet controlled TS (AbstractAutomatonList)
- Spot deterministic automaton (we only need nextstate)

States are (qs, qa) pairs:
- qs: system abstract state (Int)
- qa: spec automaton state (Int, with qa_dead=0)
"""
struct ProductAutomaton{SYS,DRA,LAB} <: SY.AbstractAutomatonList{0,0}
    sys::SYS
    dra::DRA
    labeling::LAB

    qa0::Int
    qa_dead::Int

    # bijection (qs,qa) ↔ pid
    pid::Dict{Tuple{Int,Int},Int}
    rev::Vector{Tuple{Int,Int}}

    # transitions
    post_tab::Vector{Vector{Vector{Int}}}          # post_tab[p][u] -> p′ list
    pre_tab::Vector{Vector{Tuple{Int,Int}}}        # pre_tab[p′] -> (p,u) list

    ninput::Int
end

# Required interface for AbstractAutomatonList (IMPORTANT: define in Dionysos.Symbolic namespace)
SY.get_n_state(P::ProductAutomaton) = length(P.rev)
SY.get_n_input(P::ProductAutomaton) = P.ninput
SY.pre(P::ProductAutomaton, target::Int) = P.pre_tab[target]
SY.post(P::ProductAutomaton, source::Int, symbol::Int) = P.post_tab[source][symbol]

function SY.enum_transitions(P::ProductAutomaton)
    trans = Tuple{Int,Int,Int}[]
    for q in 1:SY.get_n_state(P), u in 1:SY.get_n_input(P)
        for q2 in SY.post(P, q, u)
            push!(trans, (q2, q, u))
        end
    end
    return trans
end

"""
Build product automaton.

By default, we explore only from `initial_set` of system states.
- initial_set should be an iterable of system states (Int).
"""
function build_product_automaton(
    sys::SY.AbstractAutomatonList,
    dra,
    labeling;
    initial_set = 1:SY.get_n_state(sys),
    qa0::Int = 1,
    qa_dead::Int = 0,
)
    nS = SY.get_n_state(sys)
    nU = SY.get_n_input(sys)

    pid = Dict{Tuple{Int,Int},Int}()
    rev = Tuple{Int,Int}[]

    getpid(qs::Int, qa::Int) = get!(pid, (qs, qa)) do
        push!(rev, (qs, qa))
        length(rev)
    end

    # --- reachability exploration (BFS over product states) ---
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
                ap = labeling(qs2)
                qa2 = (qa == qa_dead) ? qa_dead : nextstate_int(dra, qa, ap, qa_dead)
                p2 = getpid(qs2, qa2)

                if !(p2 in inqueue)
                    push!(work, p2); push!(inqueue, p2)
                end
            end
        end
    end

    nP = length(rev)
    post_tab = [[Int[] for _ in 1:nU] for _ in 1:nP]
    pre_tab  = [Tuple{Int,Int}[] for _ in 1:nP]

    # --- fill adjacency lists ---
    for p in 1:nP
        (qs, qa) = rev[p]
        for u in 1:nU
            succs = Int[]
            for qs2 in SY.post(sys, qs, u)
                ap = labeling(qs2)
                qa2 = (qa == qa_dead) ? qa_dead : nextstate_int(dra, qa, ap, qa_dead)
                p2 = pid[(qs2, qa2)]
                push!(succs, p2)
                push!(pre_tab[p2], (p, u))
            end
            post_tab[p][u] = unique(succs)
        end
    end

    return ProductAutomaton(sys, dra, labeling, qa0, qa_dead, pid, rev, post_tab, pre_tab, nU)
end

# ============================================================
# Controller
# ============================================================

struct FiniteMemoryLTLController{SYS,DRA,LAB,CP}
    sys::SYS
    dra::DRA
    labeling::LAB
    contrP::CP                 # controller on product states (p -> u)
    pid::Dict{Tuple{Int,Int},Int}  # (qs,qa) -> p

    qa0::Int
    qa_dead::Int
end

"""
Return control input u for the ORIGINAL system state qs, given current memory qa.
Returns `nothing` if no control is defined (losing state).
"""
function control(ctrl::FiniteMemoryLTLController, qs::Int, qa::Int)
    p = get(ctrl.pid, (qs, qa), nothing)
    p === nothing && return nothing
    u = ST.get_control(ctrl.contrP, p)
    return u
end

@inline function update_memory(ctrl::FiniteMemoryLTLController, qa::Int, qs_next::Int, step_spec)
    ap = ctrl.labeling(qs_next)
    return (qa == ctrl.qa_dead) ? ctrl.qa_dead : (step_spec(ctrl.dra, qa, ap) === nothing ? ctrl.qa_dead : step_spec(ctrl.dra, qa, ap))
end

function controller_for_original_system(sys, dra, labeling, P, contrP; qa0=1, qa_dead=0)
    return FiniteMemoryLTLController(sys, dra, labeling, contrP, P.pid, qa0, qa_dead)
end

"""
Simulate a finite-memory LTL controller on the ORIGINAL system.

Arguments:
- sys: your original AbstractAutomatonList
- ctrl: FiniteMemoryLTLController
- qs0: initial system state
- horizon: number of steps
- succ_rule: how the environment picks a successor among post(sys,qs,u)
    - default: :first (pick first successor)
"""
function rollout_original(sys, ctrl; qs0::Int=1, horizon::Int=10, succ_rule::Symbol=:first,step_spec)
    qs = qs0
    qa = ctrl.qa0

    println("Rollout:")
    for k in 1:horizon
        u = control(ctrl, qs, qa)
        if u === nothing
            println("  step=$k  qs=$qs  qa=$qa  -> no control (losing)")
            return
        end
        println("  step=$k  qs=$qs  qa=$qa  -> u=$u")
        succ = SY.post(sys, qs, u)
        if isempty(succ)
            println("  dead-end transition at qs=$qs with u=$u")
            return
        end
        qs_next = (succ_rule == :first) ? first(succ) : succ[rand(1:length(succ))]
        qa = update_memory(ctrl, qa, qs_next, step_spec)
        qs = qs_next
    end
end


# ============================================================
# Toy nondeterministic controlled automaton (Dionysos interface)
# ============================================================

struct ToyAutomaton <: SY.AbstractAutomatonList{5,2}
    post_tab::Vector{Vector{Vector{Int}}}          # post_tab[q][u] = [q′...]
    pre_tab::Vector{Vector{Tuple{Int,Int}}}        # pre_tab[q′] = [(q,u)...]
end

function ToyAutomaton()
    nS, nU = 5, 2
    post_tab = [[Int[] for _ in 1:nU] for _ in 1:nS]
    pre_tab  = [Tuple{Int,Int}[] for _ in 1:nS]

    add!(q, u, q2) = (push!(post_tab[q][u], q2); push!(pre_tab[q2], (q, u)))

    # transitions
    add!(1, 1, 2); add!(1, 2, 5)
    add!(2, 1, 3); add!(2, 2, 4); add!(2, 2, 5)
    add!(3, 1, 4); add!(3, 2, 5)
    add!(4, 1, 4); add!(4, 2, 3)  # goal2 absorbing
    add!(5, 1, 2); add!(5, 2, 3)  # obstacle absorbing

    return ToyAutomaton(post_tab, pre_tab)
end

# Required interface for AbstractAutomatonList
SY.get_n_state(::ToyAutomaton) = 5
SY.get_n_input(::ToyAutomaton) = 2
SY.post(A::ToyAutomaton, q::Int, u::Int) = A.post_tab[q][u]
SY.pre(A::ToyAutomaton, q::Int) = A.pre_tab[q]
SY.enum_transitions(A::ToyAutomaton) =
    [(q2, q, u) for q in 1:5, u in 1:2 for q2 in SY.post(A, q, u)]


# ============================================================
# Demo: co-safe LTL "avoid obs, reach g1 then g2"
# ============================================================

function is_supported_cosafe_fragment(φ::Spot.SpotFormula)
    return is_reachability(φ) || is_constrained_reachability(φ)
end

"""
Check whether φ is acceptable for the co-safe LTL pipeline.

Modes:
- strict=true: require Spot's syntactic classification (reachability / constrained reachability)
- strict=false: allow other formulas but require that our monitor-based acceptance
  heuristic finds some target states (behavioral check).
"""
function require_cosafe_for_dionysos!(φ::Spot.SpotFormula; strict::Bool=true)
    if is_supported_cosafe_fragment(φ)
        return true
    end
    if strict
        error("Formula is not in the supported co-safe fragment (reachability / constrained reachability) per Spot.jl. " *
              "Either simplify the spec, or run with strict=false to attempt the heuristic co-safe pipeline.")
    end
    return false
end

function collect_aps(φ::Spot.SpotFormula)
    return [Symbol(ap) for ap in atomic_prop_collect(φ)]
end
function ensure_nonempty_target_set!(target_set::Vector{Int})
    isempty(target_set) && error("Empty target_set on the product. Check labeling/AP consistency or formula validity.")
end


# # function labeling_toy(q::Int)
# #     q == 3 && return (:g1,)
# #     q == 4 && return (:g2,)
# #     q == 2 && return (:obs,)
# #     return ()
# # end

# goal1set = [3]
# goal2set = [4]
# obstset = [2]
# function labeling(qs::Int)
#     ap = Symbol[]
#     qs in goal1set && push!(ap, :g1)
#     qs in goal2set && push!(ap, :g2)
#     qs in obstset  && push!(ap, :obs)
#     return Tuple(ap)
# end

# φ = ltl"G(!obs) & F(g1 & F(g2))"

# aps = collect_aps(φ)
# println(require_cosafe_for_dionysos!(φ; strict=false))   # allow beyond reach-avoid if needed

# dra = DeterministicRabinAutomata(φ)

# init_sys = [1]
# P = build_product_automaton(sys, dra, labeling; initial_set=init_sys)

# doneQ, reachableQ = cosafe_done_states_dra(dra; aps=aps, q0=1)
# target_set = [p for p in 1:SY.get_n_state(P) if P.rev[p][2] in doneQ]
# ensure_nonempty_target_set!(target_set)

# # solve reachability on product
# init_prod = [P.pid[(qs, P.qa0)] for qs in init_sys]
# contrP, controllableP, uncontrollableP, V = SY.compute_worst_case_uniform_cost_controller(P, target_set; initial_set=init_prod)


# # 3) System
# sys = ToyAutomaton()

# # 4) Product (start from initial system state 1 for a clean test)
# initial_set=[1]
# P = build_product_automaton(sys, dra, labeling_toy; initial_set=initial_set)

# # 5) Target set in the product = all states with qa == q_accept
# # target_set = [p for p in 1:SY.get_n_state(P) if P.rev[p][2] == q_accept]

# # 6) Solve (your existing worst-case reachability solver)
# init_prod = [P.pid[(qs, P.qa0)] for qs in initial_set]
# contrP, controllableP, uncontrollableP, V = SY.compute_worst_case_uniform_cost_controller(P, target_set; initial_set=init_prod)

# @info "Product size" nP=SY.get_n_state(P) nU=SY.get_n_input(P)
# @info "Winning states count" length(controllableP)

# # ============================================================
# # (Optional) Quick visualization of value function on product
# # ============================================================

# # Plot values for product states grouped by system state
# function plot_product_values(P::ProductAutomaton, V::Vector{Float64})
#     qs_list = [P.rev[p][1] for p in 1:length(P.rev)]
#     scatter(qs_list, V;
#         xlabel="System state qs",
#         ylabel="Value (steps-to-target; Inf=losing)",
#         title="Co-safe LTL on product: value vs system state",
#         legend=false,
#     )
# end

# display(plot_product_values(P, V))
# println(contrP)
# println(controllableP)
# println(uncontrollableP)
# println(V)


# println("target_set = ", target_set)
# println("target product states: ")
# for p in target_set
#     println("p=$p  (qs, qa)=", P.rev[p])
# end


# ctrl = controller_for_original_system(sys, dra, labeling_toy, P, contrP; qa0=1, qa_dead=0)
# rollout_original(sys, ctrl; qs0=1, horizon=10, succ_rule=:first)



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
    add!(1, 1, 2);
    add!(2, 1, 3); add!(2, 2, 5)

    # Key change: from g1 (state 3), input 2 can go to danger (2) or obstacle (5)
    add!(3, 1, 1)
    add!(3, 2, 2); 

    add!(4, 1, 4); add!(4, 2, 4)
    add!(5, 1, 1); add!(5, 2, 2)

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

function build_product_automaton_generic(
    sys::SY.AbstractAutomatonList,
    step_spec::Function;                      # (qa, ap_tuple) -> qa2
    labeling,
    initial_set = 1:SY.get_n_state(sys),
    qa0::Int = 1,
)
    nU = SY.get_n_input(sys)

    pid = Dict{Tuple{Int,Int},Int}()
    rev = Tuple{Int,Int}[]

    getpid(qs::Int, qa::Int) = get!(pid, (qs, qa)) do
        push!(rev, (qs, qa))
        length(rev)
    end

    # BFS on product
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
                ap = labeling(qs2)
                qa2 = step_spec(qa, ap)
                p2 = getpid(qs2, qa2)
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
                ap = labeling(qs2)
                qa2 = step_spec(qa, ap)
                p2 = pid[(qs2, qa2)]
                push!(succs, p2)
                push!(pre_tab[p2], (p, u))
            end
            post_tab[p][u] = unique(succs)
        end
    end

    # You can reuse your existing ProductAutomaton type; just store rev/pid/post/pre/ninput.
    return ProductAutomatonGeneric(sys, labeling, pid, rev, post_tab, pre_tab, nU)
end


function labeling_toy_danger(qs::Int)
    ap = Symbol[]
    qs == 3 && push!(ap, :g1)
    qs == 1 && push!(ap, :g2)
    qs == 4 && push!(ap, :obs)
    qs == 2 && push!(ap, :danger)   # <-- danger region
    return Tuple(ap)
end


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

    println(q)
    println(ap)
    println(g1)
    println(g2)
    println(danger)
    println(obs)
    println()
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
mon = MonitorG1NoDangerUntilG2()
step_spec = (qa, ap) -> mon_next(mon, qa, ap)

P = build_product_automaton_generic(sys, step_spec; labeling=labeling_toy_danger, initial_set=[1], qa0=1)

target_set = [p for p in 1:SY.get_n_state(P) if P.rev[p][2] == 3]  # qa==3 is done

init_prod = [P.pid[(1, 1)]]   # (qs=1, qa0=1)
contrP, controllableP, uncontrollableP, V =
    SY.compute_worst_case_uniform_cost_controller(P, target_set; initial_set=init_prod)

println(contrP)
println(controllableP)
println(uncontrollableP)
println(V)

# ctrl = controller_for_original_system(sys, dra, labeling_toy, P, contrP; qa0=1, qa_dead=0)
rollout_original(sys, ctrl; qs0=1, horizon=10, succ_rule=:first,step_spec)




