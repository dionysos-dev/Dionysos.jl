export EllipsoidsAbstraction

module EllipsoidsAbstraction

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const CO = DI.Control
const SY = DI.Symbolic
const PR = DI.Problem

using JuMP
using LinearAlgebra

mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    state_grid::Union{Nothing, DO.GridEllipsoidalRectangular}
    problem::Union{Nothing, PR.OptimalControlProblem}
    symmodel::Union{Nothing, SY.SymbolicModelList}
    transitionCost::Union{Nothing, Dict}
    transitionCont::Union{Nothing, Dict}
    controller::Union{Nothing,UT.SortedTupleSet{2,NTuple{2,Int}}}
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

function build_abstraction(
    problem,
    state_grid::DO.GridEllipsoidalRectangular,
    opt_sdp::MOI.OptimizerWithAttributes,
    opt_ip::MOI.OptimizerWithAttributes
)
    system = problem.system
    # Instantiate $\mathcal{X}_d$ as `domainX`
    
    domainX = DO.DomainList(state_grid); # state space
    DO.add_set!(domainX, state_grid.rect, DO.OUTER)

    # the input space 

    domainU = domainX; # input space
    DO.add_set!(domainU, state_grid.rect, DO.OUTER)

    #and the symbolic model for the state-feedback abstraction $\tilde{\mathcal{S}}$
    symmodel = SY.NewSymbolicModelListList(domainX, domainU);
    empty!(symmodel.autom)

    
    # Now let us define the L matrix defining the stage cost $\mathcal{J}(x,u) = ||L \cdot [x; u ; 1]||^2_2$
    Q_aug = CO.get_full_psd_matrix(problem.transition_cost[1][1])
    eigen_Q = eigen(Q_aug);
    L = (sqrt.(eigen_Q.values).*(eigen_Q.vectors'))'

    # ## Building the abstraction

    # We then initialize the dictionaries for saving the cost and the controller associated with each transition in the abstraction

    transitionCost = Dict()  #dictionary with cost of each transition
    transitionCont = Dict() #dictionary with controller associated each transition


    # and finally build the state-feedback abstraction 
 
    U = system.ext[:U]
    W = system.ext[:W]
    t = @elapsed SY.compute_symmodel_from_hybridcontrolsystem!(symmodel,transitionCost, transitionCont, system, W, L, U, opt_sdp, opt_ip);
 
    # println("Abstraction created in $t seconds with $(length(transitionCost)) transitions")
    symmodel, transitionCost, transitionCont
end
 
function MOI.optimize!(optimizer::Optimizer)
    problem = optimizer.problem
    system = problem.system
    state_grid = optimizer.state_grid
    symmodel, transitionCost, transitionCont = build_abstraction(problem, state_grid, optimizer.sdp_solver, optimizer.ip_solver)
    optimizer.symmodel = symmodel
    optimizer.transitionCost = transitionCost
    optimizer.transitionCont = transitionCont
    
    # Now let us prepare to synthesize our controller. Define the specifications
    

    Xinit = DO.DomainList(state_grid) # set of initial cells
    DO.add_coord!(Xinit, problem.initial_set)

    Xfinal = DO.DomainList(state_grid) # set of target cells
    DO.add_coord!(Xfinal, Vector(problem.target_set))

    Xobstacles = DO.DomainList(state_grid) # set of obstacle cells
    for o in system.ext[:obstacles]
        DO.add_set!(Xobstacles, o, DO.OUTER) 
    end

    initlist = [SY.get_state_by_xpos(symmodel, pos) for pos in DO.enum_pos(Xinit)]; 
    finallist = [SY.get_state_by_xpos(symmodel, pos) for pos in DO.enum_pos(Xfinal)];
    obstaclelist = [SY.get_state_by_xpos(symmodel, pos) for pos in DO.enum_pos(Xobstacles)];


    # ## Control synthesis

    # Before synthesizing our controller let us define this auxiliary function to transform the transition cost from Dict to Vector{Tuple{String, String, Int64}} and use it. It also removes states associated with obstacles in $\mathcal{O}$

    function transFdict(transitionCost)
        testgraph = Vector{Tuple{Int64, Int64, Float64}}()
        for (key, value) in transitionCost
            
            if key[2] ∉ obstaclelist
                push!(testgraph,((key[1]),(key[2]),value))
            end
        end
        return testgraph
    end

    testgraph = transFdict(transitionCost); # uses said function


    # We perform the synthesis of the discrete controller through Dijkstra's algorithm to the reversed graph associated to the abstraction $\tilde{\mathcal{S}}$

    controller = CO.NewControllerList();

    src, dst = initlist[1], finallist[1] # initial and goal sets

    rev_graph = [(t[2],t[1],t[3]) for t in testgraph] # we applied dijkstra to the reversed graph

    gc = UT.Digraph(rev_graph) 
    t = @elapsed rev_path, lyap_fun = UT.dijkstrapath(gc, dst, src) # gets optimal path
    path = reverse(rev_path)


    # println("Shortest path from $src to $dst found in $t seconds:\n ", isempty(path) ? "no possible path" : join(path, " → "), " (cost $(cost[dst]))")

    # We create the list of control actions along the optimal path
    for l = 1:length(path)-1 
        new_action = (path[l], path[l+1])
        UT.push_new!(controller, new_action)
    end
    optimizer.controller = controller
    optimizer.lyap_fun = lyap_fun
    return 
end

end