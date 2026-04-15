# ============================================================
# CoSafe LTL Control
# ============================================================

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
        optimizer.early_stop ? problem.initial_set : collect(1:ST.get_n_state(autom))
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
    product_target_set =
        [p for p in 1:ST.get_n_state(product_autom) if product_autom.rev[p][2] in accQ]
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
        ST.get_n_state(autom),
    )
    optimizer.uncontrollable_set = project_initial_memory_uncontrollable_set(
        product_autom,
        product_controllable_set,
        spec,
        labeling,
        ST.get_n_state(autom),
    )
    optimizer.value_fun_tab = project_initial_memory_value_function(
        product_autom,
        product_value_fun_tab,
        spec,
        labeling,
        ST.get_n_state(autom),
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

struct ProductAutomaton{SYS, LAB, STEP} <: ST.AbstractAutomatonList{0, 0}
    sys::SYS
    labeling::LAB
    spec::STEP

    pid::Dict{Tuple{Int, Int}, Int}
    rev::Vector{Tuple{Int, Int}}
    post_tab::Vector{Vector{Vector{Int}}}
    pre_tab::Vector{Vector{Tuple{Int, Int}}}
    ninput::Int
end

ST.get_n_state(P::ProductAutomaton) = length(P.rev)
ST.get_n_input(P::ProductAutomaton) = P.ninput
ST.pre(P::ProductAutomaton, t::Int) = P.pre_tab[t]
ST.post(P::ProductAutomaton, s::Int, u::Int) = P.post_tab[s][u]

function enum_transitions(P::ProductAutomaton)
    trans = Tuple{Int, Int, Int}[]
    for q in 1:ST.get_n_state(P), u in 1:ST.get_n_input(P)
        for q2 in ST.post(P, q, u)
            push!(trans, (q2, q, u))
        end
    end
    return trans
end

# build the product automaton only from the reachable states starting from the initial set:
function build_product_automaton(
    sys::ST.AbstractAutomatonList,
    spec::AbstractSpecStepper,
    labeling;
    initial_set = 1:ST.get_n_state(sys),
)
    nU = ST.get_n_input(sys)
    qa0 = init_state(spec)

    pid = Dict{Tuple{Int, Int}, Int}()
    rev = Tuple{Int, Int}[]

    getpid(qs::Int, qa::Int) =
        get!(pid, (qs, qa)) do
            push!(rev, (qs, qa))
            return length(rev)
        end

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
        p = work[i]
        i += 1
        (qs, qa) = rev[p]
        for u in 1:nU
            for qs2 in ST.post(sys, qs, u)
                ap = labeling(qs2)
                qa2 = step(spec, qa, ap)
                p2 = getpid(qs2, qa2)
                if !(p2 in inqueue)
                    push!(work, p2)
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
            for qs2 in ST.post(sys, qs, u)
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

    h = function (qa::Int, qs::Int)
        p = get(pid, (qs, qa), nothing)
        p === nothing && return nothing
        return kP(p)
    end

    g = function (qa::Int, qs_for_update::Int)
        return step(spec, qa, labeling_qs(qs_for_update))
    end

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
