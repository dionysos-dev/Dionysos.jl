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
@inline function nextstate_int(dra, qa::Int, ap::Tuple, qa_dead::Int)
    qa2 = nextstate(dra, qa, ap)
    return qa2 === nothing ? qa_dead : qa2
end

"""
Run automaton on a sequence of valuations (each valuation is a tuple of Symbols).
Returns the visited automaton states.
"""
function run_spec(dra, q0::Int, vals::Vector{<:Tuple})
    qs = Int[q0]
    q = q0
    for v in vals
        q = nextstate(dra, q, v)
        push!(qs, q)
    end
    return qs
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

@inline function update_memory(ctrl::FiniteMemoryLTLController, qa::Int, qs_next::Int)
    ap = ctrl.labeling(qs_next)
    return (qa == ctrl.qa_dead) ? ctrl.qa_dead : (nextstate(ctrl.dra, qa, ap) === nothing ? ctrl.qa_dead : nextstate(ctrl.dra, qa, ap))
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
function rollout_original(sys, ctrl; qs0::Int=1, horizon::Int=10, succ_rule::Symbol=:first)
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
        qa = update_memory(ctrl, qa, qs_next)
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

"""
Atomic proposition labeling for toy system states.
Returns a tuple of Symbols that are true.
"""
function labeling_toy(q::Int)
    q == 3 && return (:g1,)
    q == 4 && return (:g2,)
    q == 2 && return (:obs,)
    return ()
end


# ============================================================
# Demo: co-safe LTL "avoid obs, reach g1 then g2"
# ============================================================

# 1) Spec
φ = ltl"G(!obs) & F(g1 & F(g2))"
dra = DeterministicRabinAutomata(φ)

# 2) Discover an "accepting" automaton state for this co-safe task
#    (for co-safe tasks, the accepting state is typically a sink after success)
vals_good = [(), (:g1,), (), (:g2,)]
q_trace = run_spec(dra, 1, vals_good)
q_accept = last(q_trace)
@info "Automaton trace" q_trace q_accept

# 3) System
sys = ToyAutomaton()

# 4) Product (start from initial system state 1 for a clean test)
initial_set=[1]
P = build_product_automaton(sys, dra, labeling_toy; initial_set=initial_set)

# 5) Target set in the product = all states with qa == q_accept
target_set = [p for p in 1:SY.get_n_state(P) if P.rev[p][2] == q_accept]

# 6) Solve (your existing worst-case reachability solver)
init_prod = [P.pid[(qs, P.qa0)] for qs in initial_set]
contrP, controllableP, uncontrollableP, V = SY.compute_worst_case_uniform_cost_controller(P, target_set; initial_set=init_prod)

@info "Product size" nP=SY.get_n_state(P) nU=SY.get_n_input(P)
@info "Winning states count" length(controllableP)

# ============================================================
# (Optional) Quick visualization of value function on product
# ============================================================

# Plot values for product states grouped by system state
function plot_product_values(P::ProductAutomaton, V::Vector{Float64})
    qs_list = [P.rev[p][1] for p in 1:length(P.rev)]
    scatter(qs_list, V;
        xlabel="System state qs",
        ylabel="Value (steps-to-target; Inf=losing)",
        title="Co-safe LTL on product: value vs system state",
        legend=false,
    )
end

display(plot_product_values(P, V))
println(contrP)
println(controllableP)
println(uncontrollableP)
println(V)


println("target_set = ", target_set)
println("target product states: ")
for p in target_set
    println("p=$p  (qs, qa)=", P.rev[p])
end


ctrl = controller_for_original_system(sys, dra, labeling_toy, P, contrP; qa0=1, qa_dead=0)
rollout_original(sys, ctrl; qs0=1, horizon=10, succ_rule=:first)







