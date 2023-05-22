export Abstraction

module Abstraction

using JuMP
using LinearAlgebra

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const CO = DI.Control
const SY = DI.Symbolic
const PR = DI.Problem

mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    state_grid::Union{Nothing, DO.Grid}
    input_grid::Union{Nothing, DO.Grid}
    problem::Union{Nothing, PR.OptimalControlProblem, PR.SafetyProblem}
    symmodel::Union{Nothing, SY.SymbolicModelList}
    abstract_problem::Union{Nothing, PR.OptimalControlProblem, PR.SafetyProblem}
    abstract_controller::Union{Nothing, UT.SortedTupleSet{2,NTuple{2,Int}}}
    controller
    function Optimizer{T}() where {T}
        return new{T}(
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing
        )
    end
end
Optimizer() = Optimizer{Float64}()

MOI.is_empty(optimizer::Optimizer) = optimizer.problem === nothing

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    setproperty!(model, Symbol(param.name), value)
end
function MOI.get(model::Optimizer, param::MOI.RawOptimizerAttribute)
    getproperty(model, Symbol(param.name))
end


function build_abstraction(
    system,
    state_grid::DO.Grid,
    input_grid::DO.Grid,
)
    Xfull = DO.DomainList(state_grid)
    DO.add_set!(Xfull, system.X, DO.OUTER)
    Ufull = DO.DomainList(input_grid)
    DO.add_set!(Ufull, system.U, DO.CENTER)
    symmodel = SY.NewSymbolicModelListList(Xfull, Ufull)
    @time SY.compute_symmodel_from_controlsystem!(
        symmodel,
        system.f,
    )
    return symmodel
end


function build_abstract_problem(
    problem::PR.OptimalControlProblem,
    symmodel::SY.SymbolicModelList,
)
    state_grid = symmodel.Xdom.grid
    Xinit = DO.DomainList(state_grid)
    DO.add_subset!(Xinit, symmodel.Xdom, problem.initial_set, DO.OUTER)
    Xtarget = DO.DomainList(state_grid)
    DO.add_subset!(Xtarget, symmodel.Xdom, problem.target_set, DO.OUTER)
    return PR.OptimalControlProblem(
        symmodel,
        Xinit,
        Xtarget,
        problem.state_cost, # TODO this is the continuous cost, not the abstraction
        problem.transition_cost, # TODO this is the continuous cost, not the abstraction
        problem.time, # TODO this is the continuous time, not the number of transition
    )
end

function build_abstract_problem(
    problem::PR.SafetyProblem,
    symmodel::SY.SymbolicModelList,
)
    state_grid = symmodel.Xdom.grid
    Xinit = DO.DomainList(state_grid)
    DO.add_subset!(Xinit, symmodel.Xdom, problem.initial_set, DO.INNER)
    Xsafe = DO.DomainList(state_grid)
    DO.add_subset!(Xsafe, symmodel.Xdom, problem.safe_set, DO.INNER)
    return PR.SafetyProblem(
        symmodel,
        Xinit,
        Xsafe,
        problem.time, # TODO this is the continuous time, not the number of transition
    )
end

function _compute_controller!(controller, discrete_problem::PR.OptimalControlProblem)
    init_list = [SY.get_state_by_xpos(discrete_problem.system, pos) for pos in DO.enum_pos(discrete_problem.initial_set)]
    target_list = [SY.get_state_by_xpos(discrete_problem.system, pos) for pos in DO.enum_pos(discrete_problem.target_set)]
    CO.compute_controller_reach!(
        controller,
        discrete_problem.system.autom,
        init_list,
        target_list,
    )
    return
end

function _compute_controller!(controller, discrete_problem::PR.SafetyProblem)
    init_list = [SY.get_state_by_xpos(discrete_problem.system, pos) for pos in DO.enum_pos(discrete_problem.initial_set)]
    safe_list = [SY.get_state_by_xpos(discrete_problem.system, pos) for pos in DO.enum_pos(discrete_problem.safe_set)]
    CO.compute_controller_safe!(
        controller,
        discrete_problem.system.autom,
        init_list,
        safe_list,
    )
    return
end

function MOI.optimize!(optimizer::Optimizer)
    # design the abstraction
    symmodel = build_abstraction(optimizer.problem.system, optimizer.state_grid, optimizer.input_grid)
    optimizer.symmodel = symmodel
    # design the abstract problem
    abstract_problem = build_abstract_problem(optimizer.problem, symmodel)
    optimizer.abstract_problem = abstract_problem
    # solve the abstract problem
    optimizer.abstract_controller = CO.NewControllerList()
    @time _compute_controller!(
        optimizer.abstract_controller,
        abstract_problem,
    )
    # refine the abstract controllerto the concrete controller
    function controller(x;param=false)
        xpos = DO.get_pos_by_coord(symmodel.Xdom.grid, x)
        if !(xpos ∈ symmodel.Xdom)
            @warn("State out of domain")
            return nothing
        end
        source = SY.get_state_by_xpos(symmodel, xpos)
        symbollist = UT.fix_and_eliminate_first(optimizer.abstract_controller, source)
        if isempty(symbollist)
            @warn("Uncontrollable state")
            return nothing
        end
        if param
            symbol = rand(collect(symbollist))[1]
        else
            symbol = first(symbollist)[1]
        end
        upos = SY.get_upos_by_symbol(symmodel, symbol)
        u = DO.get_coord_by_pos(symmodel.Udom.grid, upos)
        return u
    end
    optimizer.controller = controller
    return 
end


mutable struct OptimizerEllipsoids{T} <: MOI.AbstractOptimizer
    state_grid::Union{Nothing, DO.GridEllipsoidalRectangular}
    problem::Union{Nothing, PR.OptimalControlProblem}
    symmodel::Union{Nothing, SY.SymbolicModelList}
    transitionCost::Union{Nothing, Dict}
    transitionKappa::Union{Nothing, Dict}
    controller::Union{Nothing,UT.SortedTupleSet{2,NTuple{2,Int}}}
    lyap_fun::Union{Nothing, Any}
    ip_solver::Union{Nothing, MOI.OptimizerWithAttributes}
    sdp_solver::Union{Nothing, MOI.OptimizerWithAttributes}
    function OptimizerEllipsoids{T}() where {T}
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
OptimizerEllipsoids() = OptimizerEllipsoids{Float64}()

MOI.is_empty(optimizer::OptimizerEllipsoids) = optimizer.problem === nothing

function MOI.set(model::OptimizerEllipsoids, param::MOI.RawOptimizerAttribute, value)
    setproperty!(model, Symbol(param.name), value)
end
function MOI.get(model::OptimizerEllipsoids, param::MOI.RawOptimizerAttribute)
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
    transitionKappa = Dict() #dictionary with controller associated each transition


    # and finally build the state-feedback abstraction 
 
    U = system.ext[:U]
    W = system.ext[:W]
    t = @elapsed SY.compute_symmodel_from_hybridcontrolsystem!(symmodel,transitionCost, transitionKappa, system, W, L, U, opt_sdp, opt_ip);
 
    # println("Abstraction created in $t seconds with $(length(transitionCost)) transitions")
    symmodel, transitionCost, transitionKappa
end
 
function MOI.optimize!(optimizer::OptimizerEllipsoids)
    problem = optimizer.problem
    system = problem.system
    state_grid = optimizer.state_grid
    symmodel, transitionCost, transitionKappa = build_abstraction(problem, state_grid, optimizer.sdp_solver, optimizer.ip_solver)
    optimizer.symmodel = symmodel
    optimizer.transitionCost = transitionCost
    optimizer.transitionKappa = transitionKappa
    
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












mutable struct OptimizerLazyEllipsoids <: MOI.AbstractOptimizer
    problem #::Union{Nothing, PR.OptimalControlProblem}
    distance
    rand_state
    new_conf
    keep
    stop_crit
    RRTstar::Bool
    compute_transition
    maxIter::Int    
    maxδx
    maxδu
    λ
    sdp_opt
    k1
    k2
    tree::Union{Nothing, UT.Tree} #later we could create a symmodel with overelapping cells

    # sdp_solver::Union{Nothing, MOI.OptimizerWithAttributes}
end

function build_OptimizerLazyEllipsoids(problem, distance, rand_state, new_conf, keep, stop_crit, RRTstar, compute_transition, maxIter, maxδx, maxδu, λ, sdp_opt, k1, k2)
    return OptimizerLazyEllipsoids(problem, distance, rand_state, new_conf, keep, stop_crit, RRTstar, compute_transition, maxIter, maxδx, maxδu, λ, sdp_opt, k1, k2, nothing)
end

OptimizerLazyEllipsoids() = OptimizerLazyEllipsoids{Float64}()

MOI.is_empty(optimizer::OptimizerLazyEllipsoids) = optimizer.problem === nothing

function MOI.set(model::OptimizerLazyEllipsoids, param::MOI.RawOptimizerAttribute, value)
    setproperty!(model, Symbol(param.name), value)
end
function MOI.get(model::OptimizerLazyEllipsoids, param::MOI.RawOptimizerAttribute)
    getproperty(model, Symbol(param.name))
end
 
function MOI.optimize!(optimizer::OptimizerLazyEllipsoids)
    problem = optimizer.problem
    Einit = problem.initial_set
    Etarget = problem.target_set

    distance = optimizer.distance
    rand_state = optimizer.rand_state
    new_conf = optimizer.new_conf
    keep = optimizer.keep
    stop_crit = optimizer.stop_crit
    maxIter = optimizer.maxIter
    RRTstar = optimizer.RRTstar
    compute_transition = optimizer.compute_transition
    k1 = optimizer.k1
    k2 = optimizer.k2
    tree = UT.RRT(Etarget, Einit, distance, rand_state, new_conf, keep, stop_crit, optimizer; maxIter=maxIter, RRTstar=RRTstar, compute_transition, k1=k1, k2=k2)
    
    optimizer.tree = tree 
    return 
end

end
