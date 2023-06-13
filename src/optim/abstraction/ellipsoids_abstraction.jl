export EllipsoidsAbstraction

module EllipsoidsAbstraction

using LinearAlgebra
using JuMP

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const CO = DI.Control
const SY = DI.Symbolic
const PR = DI.Problem

mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    state_grid::Union{Nothing, DO.GridEllipsoidalRectangular}
    problem::Union{Nothing, PR.OptimalControlProblem}
    symmodel::Union{Nothing, SY.SymbolicModelList}
    transitionCost::Union{Nothing, Dict}
    transitionCont::Union{Nothing, Dict}
    abstract_problem::Union{Nothing, PR.OptimalControlProblem}
    abstract_controller::Union{Nothing, UT.SortedTupleSet{2,NTuple{2,Int}}}
    concrete_controller
    lyap_fun::Union{Nothing, Any}
    ip_solver::Union{Nothing, MOI.OptimizerWithAttributes}
    sdp_solver::Union{Nothing, MOI.OptimizerWithAttributes}
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
        )
    end
end

Optimizer() = Optimizer{Float64}()

MOI.is_empty(optimizer::Optimizer) = optimizer.problem === nothing

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    if param.name == "problem"
        if !(value isa PR.OptimalControlProblem)
            throw(MOI.UnsupportedAttribute(param, "$(typeof(value)) not supported"))
        end
    end
    setproperty!(model, Symbol(param.name), value)
end
function MOI.get(model::Optimizer, param::MOI.RawOptimizerAttribute)
    getproperty(model, Symbol(param.name))
end

function get_guaranteed_cost(optimizer::Optimizer, x)
    symmodel = optimizer.symmodel
    state_grid = symmodel.Xdom.grid
    l_state = SY.get_all_states_by_xpos(symmodel, DO.crop_to_domain(symmodel.Xdom, DO.get_all_pos_by_coord(state_grid , x)))
    lyap_fun = optimizer.lyap_fun
    lyap_min = Inf
    for state in l_state
        lyap = lyap_fun[state]
        lyap<lyap_min ? lyap_min = lyap : lyap_min=lyap_min
    end
    return lyap_min
end

function build_abstraction(
    problem,
    state_grid::DO.GridEllipsoidalRectangular,
    opt_sdp::MOI.OptimizerWithAttributes,
    opt_ip::MOI.OptimizerWithAttributes
)
    system = problem.system    
    # The state space
    domainX = DO.DomainList(state_grid); # state space
    X = system.ext[:X]
    DO.add_set!(domainX, X, DO.INNER) 
    
    # The input space 
    domainU = domainX; 
    DO.add_set!(domainU, state_grid.rect, DO.INNER) 

    # The symbolic model for the state-feedback abstraction
    symmodel = SY.NewSymbolicModelListList(domainX, domainU);
    empty!(symmodel.autom)

    # Now let us define the L matrix defining the stage cost $\mathcal{J}(x,u) = ||L \cdot [x; u ; 1]||^2_2$
    Q_aug = CO.get_full_psd_matrix(problem.transition_cost[1][1])
    eigen_Q = eigen(Q_aug);
    L = (sqrt.(eigen_Q.values).*(eigen_Q.vectors'))'

    # ## Building the abstraction
    # We then initialize the dictionaries for saving the cost and the controller associated with each transition in the abstraction
    transitionCont = Dict() # dictionary with controller associated each transition
    transitionCost = Dict() # dictionary with cost of each transition
    # And finally build the state-feedback abstraction 
    U = system.ext[:U]
    W = system.ext[:W]
    t = @elapsed SY.compute_symmodel_from_hybridcontrolsystem!(symmodel, transitionCont, transitionCost, system, W, L, U, opt_sdp, opt_ip);
    println("Abstraction created in $t seconds with $(length(transitionCost)) transitions")
    symmodel, transitionCont, transitionCost
end

function build_abstract_problem(
    problem::PR.OptimalControlProblem,
    symmodel::SY.SymbolicModelList
)
    state_grid = symmodel.Xdom.grid
    Xinit = DO.DomainList(state_grid) # set of initial cells
    DO.add_coord!(Xinit, problem.initial_set)

    Xfinal = DO.DomainList(state_grid) # set of target cells
    DO.add_coord!(Xfinal, Vector(problem.target_set))

    initlist = [SY.get_state_by_xpos(symmodel, pos) for pos in DO.enum_pos(Xinit)]; 
    finallist = [SY.get_state_by_xpos(symmodel, pos) for pos in DO.enum_pos(Xfinal)];

    return PR.OptimalControlProblem(
        symmodel,
        initlist,
        finallist,
        problem.state_cost, 
        problem.transition_cost, 
        problem.time,
    )
end


# Before synthesizing our controller let us define this auxiliary function to transform the transition cost from Dict to Vector{Tuple{String, String, Int64}}
function transFdict(transitionCost)
    testgraph = Vector{Tuple{Int64, Int64, Float64}}()
    for (key, value) in transitionCost
        push!(testgraph,((key[1]),(key[2]), value))
    end
    return testgraph
end

# We perform the synthesis of the discrete controller through Dijkstra's algorithm to the reversed graph associated to the abstraction
function solve_abstract_problem(abstract_problem, transitionCost)
    testgraph = transFdict(transitionCost); # uses said function
    src, dst = abstract_problem.initial_set[1], abstract_problem.target_set[1]

    rev_graph = [(t[2],t[1],t[3]) for t in testgraph] # we applied dijkstra to the reversed graph

    gc = UT.Digraph(rev_graph) 
    t = @elapsed rev_path, lyap_fun = UT.dijkstrapath(gc, dst, src) # gets optimal path
    path = reverse(rev_path)
    # println("Shortest path from $src to $dst found in $t seconds:\n ", isempty(path) ? "no possible path" : join(path, " → "), " (cost $(cost[dst]))")
    abstract_controller = CO.NewControllerList();
    for l = 1:length(path)-1 
        new_action = (path[l], path[l+1])
        UT.push_new!(abstract_controller, new_action)
    end
    println("Abstract controller computed with: ", path[end]==dst ? "succes" : "failure")
    return abstract_controller, lyap_fun
end

function solve_concrete_problem(symmodel, abstract_controller, transitionCont)
    state_grid = symmodel.Xdom.grid
    function concrete_controller(x)
        currState = SY.get_all_states_by_xpos(symmodel, DO.crop_to_domain(symmodel.Xdom, DO.get_all_pos_by_coord(state_grid, x)))
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

function MOI.optimize!(optimizer::Optimizer)
    concrete_problem = optimizer.problem
    state_grid = optimizer.state_grid
    # design the abstraction
    symmodel, transitionCont, transitionCost = build_abstraction(concrete_problem, state_grid, optimizer.sdp_solver, optimizer.ip_solver)
    optimizer.symmodel = symmodel
    optimizer.transitionCont = transitionCont
    optimizer.transitionCost = transitionCost
    # design the abstract problem
    abstract_problem = build_abstract_problem(concrete_problem, symmodel)
    optimizer.abstract_problem = abstract_problem
    # solve the abstract problem
    abstract_controller, lyap_fun = solve_abstract_problem(abstract_problem, transitionCost)
    optimizer.abstract_controller = abstract_controller
    optimizer.lyap_fun = lyap_fun
    # solve the concrete problem
    optimizer.concrete_controller = solve_concrete_problem(symmodel, abstract_controller, transitionCont)
    return 
end












end