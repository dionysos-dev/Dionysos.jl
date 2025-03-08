export SampleBasedAbstraction

module SampleBasedAbstraction

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const PR = DI.Problem

using JuMP
using DataStructures
using Profile

"""
    Optimizer{T} <: MOI.AbstractOptimizer

TODO
"""
mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    concrete_problem::Union{Nothing, PR.ProblemType}
    abstract_problem::Union{Nothing, PR.OptimalControlProblem, PR.SafetyProblem}
    abstract_system::Union{Nothing, SY.SymbolicModelList}
    abstract_controller::Union{Nothing, UT.SortedTupleSet{3, NTuple{3, Int}}}
    value_function::Union{Nothing, DefaultDict{Int, Float64}}
    concrete_controller::Any
    state_grid::Union{Nothing, DO.Grid}
    input_grid::Union{Nothing, DO.Grid}
    solve_time_sec::T
    function Optimizer{T}() where {T}
        return new{T}(nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, 0.0)
    end
end
Optimizer() = Optimizer{Float64}()

MOI.is_empty(optimizer::Optimizer) = optimizer.concrete_problem === nothing

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    return setproperty!(model, Symbol(param.name), value)
end
function MOI.get(model::Optimizer, ::MOI.SolveTimeSec)
    return model.solve_time_sec
end
function MOI.get(model::Optimizer, param::MOI.RawOptimizerAttribute)
    return getproperty(model, Symbol(param.name))
end

function build_abstraction(concrete_system, state_grid::DO.Grid, input_grid::DO.Grid)
    Xfull = DO.DomainList(state_grid)
    DO.add_set!(Xfull, concrete_system.X, DO.INNER)
    Ufull = DO.DomainList(input_grid)
    DO.add_set!(Ufull, concrete_system.U, DO.CENTER)
    abstract_system = SY.NewSymbolicModelListList(Xfull, Ufull)
    @time SY.compute_symmodel_from_data2!(abstract_system, concrete_system.f)
    return abstract_system
end


function build_abstract_problem(
    concrete_problem::PR.OptimalControlProblem,
    abstract_system::SY.SymbolicModelList,
)
    state_grid = abstract_system.Xdom.grid
    Xinit = DO.DomainList(state_grid)
    DO.add_subset!(Xinit, abstract_system.Xdom, concrete_problem.initial_set, DO.OUTER)
    Xtarget = DO.DomainList(state_grid)
    DO.add_subset!(Xtarget, abstract_system.Xdom, concrete_problem.target_set, DO.INNER)
    init_list = [SY.get_state_by_xpos(abstract_system, pos) for pos in DO.enum_pos(Xinit)]
    target_list =
        [SY.get_state_by_xpos(abstract_system, pos) for pos in DO.enum_pos(Xtarget)]
    return PR.OptimalControlProblem(
        abstract_system,
        init_list,
        target_list,
        concrete_problem.state_cost, # TODO this is the continuous cost, not the abstraction
        concrete_problem.transition_cost, # TODO this is the continuous cost, not the abstraction
        concrete_problem.time, # TODO this is the continuous time, not the number of transition
    )
end

function build_abstract_problem(
    concrete_problem::PR.SafetyProblem,
    abstract_system::SY.SymbolicModelList,
)
    state_grid = abstract_system.Xdom.grid
    Xinit = DO.DomainList(state_grid)
    DO.add_subset!(Xinit, abstract_system.Xdom, concrete_problem.initial_set, DO.OUTER)
    Xsafe = DO.DomainList(state_grid)
    DO.add_subset!(Xsafe, abstract_system.Xdom, concrete_problem.safe_set, DO.INNER)
    init_list = [SY.get_state_by_xpos(abstract_system, pos) for pos in DO.enum_pos(Xinit)]
    safe_list = [SY.get_state_by_xpos(abstract_system, pos) for pos in DO.enum_pos(Xsafe)]
    return PR.SafetyProblem(
        abstract_system,
        init_list,
        safe_list,
        concrete_problem.time, # TODO this is the continuous time, not the number of transition
    )
end

function solve_abstract_problem(abstract_system, abstract_problem::PR.OptimalControlProblem)
    abstract_controller = NewControllerList()
    @time compute_controller_reach!(
        abstract_controller,
        abstract_system,
        abstract_problem.system.autom,
        abstract_problem.initial_set,
        abstract_problem.target_set,
    )
    return abstract_controller
end

function solve_abstract_problem(abstract_system, abstract_problem::PR.SafetyProblem)
    abstract_controller = NewControllerList()
    compute_controller_safe!(
        abstract_controller,
        abstract_problem.system.autom,
        abstract_problem.initial_set,
        abstract_problem.safe_set,
    )
    return abstract_controller
end

function solve_value_function(abstract_controller, abstract_system, abstract_problem::PR.OptimalControlProblem)
    value_function = DefaultDict{Int, Float64}(typemax(Int))
    println("compute_value_function! started")
    @time compute_value_function!(
        value_function,
        abstract_controller,
        abstract_system,
        abstract_problem.system.autom,
        abstract_problem.initial_set,
        abstract_problem.target_set,
    )
    println("compute_value_function! terminated")
    return value_function
end

function solve_concrete_problem(abstract_system, abstract_controller)
    function concrete_controller(x; param = false)
        xpos = DO.get_pos_by_coord(abstract_system.Xdom.grid, x)
        if !(xpos ∈ abstract_system.Xdom)
            @warn("State out of domain")
            return nothing
        end
        source = SY.get_state_by_xpos(abstract_system, xpos)
        symbollist = UT.fix_and_eliminate_first(abstract_controller, source)
        if isempty(symbollist)
            @warn("Uncontrollable state")
            return nothing
        end
        if param
            temp = rand(collect(symbollist))
            symbol, tstep = temp[1], temp[2]
        else
            temp = first(symbollist)
            symbol, tstep = temp[1], temp[2] # symbol, timestep
        end
        upos = SY.get_upos_by_symbol(abstract_system, symbol) 
        u = DO.get_coord_by_pos(abstract_system.Udom.grid, upos)
        return u, tstep
    end
end

function MOI.optimize!(optimizer::Optimizer)
    t_ref = time()

    # Build the abstraction
    abstract_system = build_abstraction(
        optimizer.concrete_problem.system,
        optimizer.state_grid,
        optimizer.input_grid,
    )
    optimizer.abstract_system = abstract_system

    # Build the abstract problem
    abstract_problem = build_abstract_problem(optimizer.concrete_problem, abstract_system)
    optimizer.abstract_problem = abstract_problem

    # Solve the abstract problem
    abstract_controller = solve_abstract_problem(abstract_system, abstract_problem)
    optimizer.abstract_controller = abstract_controller
    
    # Solve the concrete problem
    optimizer.concrete_controller =
        solve_concrete_problem(abstract_system, abstract_controller)

    value_function = solve_value_function(abstract_controller, abstract_system, abstract_problem)
    optimizer.value_function = value_function

    optimizer.solve_time_sec = time() - t_ref
    return
end

NewControllerList() = UT.SortedTupleSet{3, NTuple{3, Int}}() # pour chaque source un input : rajouter l'info du timestep

function compute_value_function!(g, contr, abstract_system, autom, initlist, targetlist::Vector{Int})
    # data
    num_targets_unreachable = zeros(Int, autom.nstates, autom.nsymbols,100)
    _compute_num_targets_unreachable(num_targets_unreachable, autom)
    init_set = BitSet(initlist)
    target_set = BitSet(targetlist)
    current_targets = copy(targetlist)
    next_targets = Int[]
    # compute the value function
    for t in target_set
        g[t] = 0.0
    end
    num_init_unreachable = length(init_set)
    post = Int[]
    while !isempty(current_targets) && !iszero(num_init_unreachable)
        empty!(next_targets)
        for target in current_targets
            for (source, symbol, p) in SY.pre(autom, target)   # pre = predecesseur
                if !(source in target_set) &&                           
                iszero(num_targets_unreachable[source, symbol, p+1] -= 1)
                    #println(collect(UT.fix_and_eliminate_first(contr, source)))
                    symbollist = UT.fix_and_eliminate_first(contr, source)
                    if isempty(symbollist)
                        continue
                    end
                    u_policy, t_step = collect(symbollist)[1]
                    if u_policy == symbol && t_step == p
                        push!(target_set, source)
                        push!(next_targets, source)
                        if source in init_set
                            num_init_unreachable -= 1
                        end
                        empty!(post)
                        SY.compute_post!(post, autom, source, symbol, p)
                        gmax = maximum(g[q] for q in post)
                        time = 0.2 + p * 0.1 
                        # u_symb = SY.get_upos_by_symbol(abstract_system, symbol)
                        # true_u = DO.get_coord_by_pos(abstract_system.Udom.grid, u_symb)
                        cost = time
                        g[source] = cost + gmax
                    end
                end
            end
        end
        current_targets, next_targets = next_targets, current_targets
    end
end 

function _compute_num_targets_unreachable(num_targets_unreachable, autom)
    # for each target, we count how many times each source/symbol pair can reach it in 1 step
    for target in 1:(autom.nstates)
        for soursymb in SY.pre(autom, target)
            num_targets_unreachable[soursymb[1], soursymb[2], soursymb[3]+1] += 1
        end
    end
end

function _compute_controller_reach!(
    contr,
    autom,
    init_set,
    target_set,
    num_targets_unreachable,
    current_targets,
    next_targets,
)
    num_init_unreachable = length(init_set)
    while !isempty(current_targets) && !iszero(num_init_unreachable)
        empty!(next_targets)
        for target in current_targets
            for (source, symbol, p) in SY.pre(autom, target)   # pre = predecesseur
                if !(source in target_set) &&                           
                   iszero(num_targets_unreachable[source, symbol, p+1] -= 1) # si c'est 0, c'est déterministe
                    push!(target_set, source)
                    push!(next_targets, source)
                    # controleur mis à jour 
                    UT.push_new!(contr, (source, symbol, p)) 
                    if source in init_set
                        num_init_unreachable -= 1
                    end
                end
            end
        end
        current_targets, next_targets = next_targets, current_targets
    end
    return iszero(num_init_unreachable)
end

function _data(contr, autom, initlist, targetlist)
    num_targets_unreachable = zeros(Int, autom.nstates, autom.nsymbols, 100)
    _compute_num_targets_unreachable(num_targets_unreachable, autom)
    initset = BitSet(initlist)
    targetset = BitSet(targetlist)
    current_targets = copy(targetlist)
    next_targets = Int[]
    return initset, targetset, num_targets_unreachable, current_targets, next_targets
end

function _compute_controller_reach_with_cost!(
    contr,
    abstract_system,
    autom,
    init_set,                # set of initial states
    target_set,              # set of target states
    num_targets_unreachable, # see _compute_num_targets_unreachable
    N, 
    g, 
    f,
)     
    println("cost considered")
    num_init_unreachable = length(init_set)
    M = Int[]
    controlDict = Dict{Int, Tuple{Int, Int}}()   # controlDict = {source : (symbol, timestep)} (the policy)
    for t in target_set
        g[t] = 0
        enqueue!(N, (t, (-1,-1)), 0)
    end
    post = Int[]
    iter = 0
    while !isempty(N) && !iszero(num_init_unreachable)
        (q, input_opt) = dequeue!(N)
        u_opt, t_opt = input_opt
        if !(q in M) 
            push!(M, q)
            if q in init_set
                num_init_unreachable -= 1
            end
            iter += 1
            if iter % 500 == 0 
                println("states in M $iter out of $(autom.nstates)")
                println("   initial set reachability : $(length(init_set) - num_init_unreachable)")
            end
            if u_opt != -1
                g[q] = f[q,input_opt]
                controlDict[q] = (u_opt, t_opt)
            end
            for (qminus, symbol, p) in SY.pre(autom, q) 
                if !(qminus in M) && iszero(num_targets_unreachable[qminus, symbol, p+1] -= 1)
                    empty!(post)
                    SY.compute_post!(post, autom, qminus, symbol, p)
                    gmax = maximum(g[qprime] for qprime in post)
                    time = 0.2 + p * 0.1 #0.1 * 1.5 ^ p
                    # u_symb = SY.get_upos_by_symbol(abstract_system, symbol)
                    # true_u = DO.get_coord_by_pos(abstract_system.Udom.grid, u_symb)
                    cost = time #* true_u[1]^2
                    f[qminus, (symbol,p)] = cost + gmax
                    enqueue!(N, (qminus, (symbol,p)), f[qminus, (symbol,p)])
                end
            end
        end
    end
    println("number of init unreachable $num_init_unreachable out of $(length(init_set))") 
    for q in keys(controlDict)
        UT.push_new!(contr, (q, controlDict[q][1], controlDict[q][2]))
    end
    return iszero(num_init_unreachable)
end

function _data_cost(contr, autom, initlist, targetlist)
    num_targets_unreachable = zeros(Int, autom.nstates, autom.nsymbols, 100)
    _compute_num_targets_unreachable(num_targets_unreachable, autom)
    initset = BitSet(initlist)
    targetset = BitSet(targetlist)
    N = PriorityQueue{Tuple{Int,Tuple{Int,Int}}, Float64}()            
    g = DefaultDict{Int, Float64}(typemax(Int))             # g = {source : cost} (belmann value function)
    f = DefaultDict{Tuple{Int,Tuple{Int,Int}}, Float64}(typemax(Int))  # f = {(source, input) : cost} (belmann value function)
    return initset, targetset, num_targets_unreachable, N, g, f
end

function compute_controller_reach!(contr, abstract_system, autom, initlist, targetlist::Vector{Int})
    println("compute_controller_reach! started")
    # TODO: try to infer whether num_targets_unreachable is sparse or not,
    # and if sparse, use a dictionary instead
    cost = true
    if cost 
        if !_compute_controller_reach_with_cost!(
            contr,
            abstract_system,
            autom,
            _data_cost(contr, autom, initlist, targetlist)...,
        )
            println("\ncompute_controller_reach_with_cost! terminated without covering init set")
            # ProgressMeter.finish!(prog)
            return
        end
    else 
        if !_compute_controller_reach!(
            contr,
            autom,
            _data(contr, autom, initlist, targetlist)...,
        )
            println("\ncompute_controller_reach! terminated without covering init set")
            # ProgressMeter.finish!(prog)
            return
        end
    end
    # ProgressMeter.finish!(prog)
    return println("\ncompute_controller_reach! terminated with success")
end

function _compute_pairstable(pairstable, timesteps, autom)
    for target in 1:(autom.nstates)
        for soursymb in SY.pre(autom, target)
            pairstable[soursymb[1], soursymb[2]] = true
            timesteps[soursymb[1], soursymb[2]] = soursymb[3]
        end
    end
end

function compute_controller_safe!(contr, autom, initlist, safelist)
    println("compute_controller_safe! started")
    nstates = autom.nstates
    nsymbols = autom.nsymbols
    pairstable = [false for i in 1:nstates, j in 1:nsymbols]
    timesteps = Dict{Tuple{Int, Int}, Int}()
    _compute_pairstable(pairstable, timesteps, autom)
    nsymbolslist = sum(pairstable; dims = 2)
    safeset = Set(safelist)
    for source in safeset
        if nsymbolslist[source] == 0
            delete!(safeset, source)
        end
    end
    unsafeset = Set(1:nstates)
    setdiff!(unsafeset, safeset)
    for source in unsafeset
        for symbol in 1:nsymbols
            pairstable[source, symbol] = false
        end
    end 
    nextunsafeset = Set{Int}()

    # prog = ProgressUnknown("# iterations computing controller:")
    while true
        # ProgressMeter.next!(prog)
        for target in unsafeset
            for soursymb in SY.pre(autom, target)
                if pairstable[soursymb[1], soursymb[2]]
                    pairstable[soursymb[1], soursymb[2]] = false
                    nsymbolslist[soursymb[1]] -= 1
                    if nsymbolslist[soursymb[1]] == 0
                        push!(nextunsafeset, soursymb[1])
                    end
                end
            end
        end
        if isempty(nextunsafeset)
            break
        end
        setdiff!(safeset, nextunsafeset)
        unsafeset, nextunsafeset = nextunsafeset, unsafeset
        empty!(nextunsafeset)
    end
    # ProgressMeter.finish!(prog)
    for source in safeset
        for symbol in 1:nsymbols
            if pairstable[source, symbol] # source symb timestep
                UT.push_new!(contr, (source, symbol, timesteps[source, symbol]))   # push (source, symbol, timestep)
            end
        end
    end
    l = length(safeset)
    println("number of safe cells $l")
    if ⊆(initlist, safeset)
        println("\ncompute_controller_safe! terminated with success")
    else
        println("\ncompute_controller_safe! terminated without covering init set")
    end
end

end
