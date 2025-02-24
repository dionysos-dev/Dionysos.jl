export EllipsoidsAbstraction

module EllipsoidsAbstraction
using LinearAlgebra, JuMP

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const PR = DI.Problem

"""
    Optimizer{T} <: MOI.AbstractOptimizer

Abstraction-based solver for which the domain is covered with ellipsoidal cells, independently of the control task.
"""
mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    concrete_problem::Union{Nothing, PR.OptimalControlProblem}
    abstract_problem::Union{Nothing, PR.OptimalControlProblem}
    abstract_system::Union{Nothing, SY.SymbolicModelList}
    abstract_controller::Union{Nothing, UT.SortedTupleSet{2, NTuple{2, Int}}}
    concrete_controller::Any
    abstract_lyap_fun::Union{Nothing, Any}
    concrete_lyap_fun::Union{Nothing, Any}
    state_grid::Union{Nothing, DO.GridEllipsoidalRectangular}
    lyap::Union{Nothing, Any}
    transitionCost::Union{Nothing, Dict}
    transitionCont::Union{Nothing, Dict}
    sdp_solver::Union{Nothing, MOI.OptimizerWithAttributes}
    solve_time_sec::T
    function Optimizer{T}() where {T}
        return new{T}(
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            0.0,
        )
    end
end

Optimizer() = Optimizer{Float64}()

MOI.is_empty(optimizer::Optimizer) = optimizer.concrete_problem === nothing

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    if param.name == "concrete_problem"
        if !(value isa PR.OptimalControlProblem)
            throw(MOI.UnsupportedAttribute(param, "$(typeof(value)) not supported"))
        end
    end
    return setproperty!(model, Symbol(param.name), value)
end

function MOI.get(model::Optimizer, ::MOI.SolveTimeSec)
    return model.solve_time_sec
end
function MOI.get(model::Optimizer, param::MOI.RawOptimizerAttribute)
    return getproperty(model, Symbol(param.name))
end

function build_abstraction(
    concrete_problem,
    state_grid::DO.GridEllipsoidalRectangular,
    opt_sdp::MOI.OptimizerWithAttributes,
)
    concrete_system = concrete_problem.system
    # The state space
    domainX = DO.DomainList(state_grid) # state space
    X = concrete_system.ext[:X]

    DO.add_set!(domainX, X, DO.INNER)

    # The input space 
    domainU = domainX
    DO.add_set!(domainU, X, DO.INNER)

    # The symbolic model for the state-feedback abstraction
    abstract_system = SY.NewSymbolicModelListList(domainX, domainU)
    empty!(abstract_system.autom)

    # Now let us define the L matrix defining the stage cost $\mathcal{J}(x,u) = ||L \cdot [x; u ; 1]||^2_2$
    Q_aug = UT.get_full_psd_matrix(concrete_problem.transition_cost[1][1])
    eigen_Q = eigen(Q_aug)
    L = (sqrt.(eigen_Q.values) .* (eigen_Q.vectors'))'

    # ## Building the abstraction
    # We then initialize the dictionaries for saving the cost and the controller associated with each transition in the abstraction
    transitionCont = Dict() # dictionary with controller associated each transition
    transitionCost = Dict() # dictionary with cost of each transition
    # And finally build the state-feedback abstraction 
    U = concrete_system.ext[:U]
    W = concrete_system.ext[:W]
    t = @elapsed SY.compute_symmodel_from_hybridcontrolsystem!(
        abstract_system,
        transitionCont,
        transitionCost,
        concrete_system,
        W,
        L,
        U,
        opt_sdp,
    )
    println("Abstraction created in $t seconds with $(length(transitionCost)) transitions")
    return abstract_system, transitionCont, transitionCost
end

function build_abstract_problem(
    concrete_problem::PR.OptimalControlProblem,
    abstract_system::SY.SymbolicModelList,
)
    state_grid = abstract_system.Xdom.grid
    Xinit = DO.DomainList(state_grid) # set of initial cells
    DO.add_coord!(Xinit, concrete_problem.initial_set)

    Xfinal = DO.DomainList(state_grid) # set of target cells
    DO.add_coord!(Xfinal, Vector(concrete_problem.target_set))

    initlist = [SY.get_state_by_xpos(abstract_system, pos) for pos in DO.enum_pos(Xinit)]
    finallist = [SY.get_state_by_xpos(abstract_system, pos) for pos in DO.enum_pos(Xfinal)]

    return PR.OptimalControlProblem(
        abstract_system,
        initlist,
        finallist,
        concrete_problem.state_cost,
        concrete_problem.transition_cost,
        concrete_problem.time,
    )
end

# Before synthesizing our controller let us define this auxiliary function to transform the transition cost from Dict to Vector{Tuple{String, String, Int64}}
function transFdict(transitionCost)
    testgraph = Vector{Tuple{Int64, Int64, Float64}}()
    for (key, value) in transitionCost
        push!(testgraph, ((key[1]), (key[2]), value))
    end
    return testgraph
end

# We perform the synthesis of the discrete controller through Dijkstra's algorithm to the reversed graph associated to the abstraction
function solve_abstract_problem(abstract_problem, transitionCost)
    testgraph = transFdict(transitionCost) # uses said function
    src, dst = abstract_problem.initial_set[1], abstract_problem.target_set[1]

    rev_graph = [(t[2], t[1], t[3]) for t in testgraph] # we applied dijkstra to the reversed graph

    gc = UT.Digraph(rev_graph)
    t = @elapsed rev_path, lyap_fun = UT.dijkstrapath(gc, dst, src) # gets optimal path
    path = reverse(rev_path)
    # println("Shortest path from $src to $dst found in $t seconds:\n ", isempty(path) ? "no possible path" : join(path, " → "), " (cost $(cost[dst]))")
    abstract_controller = UT.SortedTupleSet{2, NTuple{2, Int}}()
    for l in 1:(length(path) - 1)
        new_action = (path[l], path[l + 1])
        UT.push_new!(abstract_controller, new_action)
    end
    println("Abstract controller computed with: ", path[end] == dst ? "succes" : "failure")
    return abstract_controller, lyap_fun
end

function solve_concrete_problem(abstract_system, abstract_controller, transitionCont)
    state_grid = abstract_system.Xdom.grid
    function concrete_controller(x)
        currState = SY.get_states_by_xpos(
            abstract_system,
            DO.crop_to_domain(abstract_system.Xdom, DO.get_all_pos_by_coord(state_grid, x)),
        )
        next_action = nothing
        for action in abstract_controller.data
            if (action[1] ∩ currState) ≠ []
                next_action = action
            end
        end
        c_eval = ST.get_c_eval(transitionCont[next_action])
        return c_eval(x)
    end
    return concrete_controller
end

function build_abstract_lyap_fun(lyap)
    return abstract_lyap_fun(state) = lyap[state]
end

function build_concrete_lyap_fun(abstract_system, abstract_lyap_fun)
    state_grid = abstract_system.Xdom.grid
    function concrete_lyap_fun(x)
        l_state = SY.get_states_by_xpos(
            abstract_system,
            DO.crop_to_domain(abstract_system.Xdom, DO.get_all_pos_by_coord(state_grid, x)),
        )
        lyap_min = Inf
        for state in l_state
            lyap = abstract_lyap_fun(state)
            lyap < lyap_min ? lyap_min = lyap : lyap_min = lyap_min
        end
        return lyap_min
    end
    return concrete_lyap_fun
end

function MOI.optimize!(optimizer::Optimizer)
    t_ref = time()

    concrete_problem = optimizer.concrete_problem
    state_grid = optimizer.state_grid
    # Build the abstraction
    abstract_system, transitionCont, transitionCost =
        build_abstraction(concrete_problem, state_grid, optimizer.sdp_solver)
    optimizer.abstract_system = abstract_system
    optimizer.transitionCont = transitionCont
    optimizer.transitionCost = transitionCost
    # Build the abstract problem
    abstract_problem = build_abstract_problem(concrete_problem, abstract_system)
    optimizer.abstract_problem = abstract_problem
    # Solve the abstract problem
    abstract_controller, lyap = solve_abstract_problem(abstract_problem, transitionCost)
    optimizer.abstract_controller = abstract_controller
    optimizer.lyap = lyap
    optimizer.abstract_lyap_fun = build_abstract_lyap_fun(lyap)
    # Solve the concrete problem
    optimizer.concrete_controller =
        solve_concrete_problem(abstract_system, abstract_controller, transitionCont)
    optimizer.concrete_lyap_fun =
        build_concrete_lyap_fun(abstract_system, optimizer.abstract_lyap_fun)

    optimizer.solve_time_sec = time() - t_ref
    return
end

end
