using ..AlternatingSimulation
const AS = AlternatingSimulation

using ..Lazy_abstraction
const LA = Lazy_abstraction

struct RectanglularObstaclesDomain{VT} <: AbstractSet{VT}
    X::AB.HyperRectangle{VT}
    O::Vector{AB.HyperRectangle{VT}}
end
function Base.in(pos, dom::RectanglularObstaclesDomain)
    return !mapreduce(Base.Fix1(in, pos), |, dom.O, init=!in(pos, dom.X))
end
function build_dom(hx)
    X, O = build_domain()
    grid = D.build_grid_in_rec(X, hx)
    periodic = [3,4]
    periods = [2*π,2*π]
    T0 = [-π,-π]
    d = D.RectanglularObstacles(X, O)
    return X, D.GeneralDomainList(grid, d; periodic, periods, T0)
end

function h(node::LA.S.Node,problem::LA.LazyAbstraction)
    source = node.state.source
    symmodel = problem.symmodel
    xpos = AB.get_xpos_by_state(symmodel, source)
    x = AB.get_coord_by_pos(symmodel.Xdom.grid, xpos)

    heuristic = problem.heuristic_data
    symmodel2 = heuristic.symmodel
    xpos2 = AB.get_pos_by_coord(symmodel2.Xdom.grid, x)
    source2 = AB.get_state_by_xpos(symmodel2, xpos2)
    return heuristic.dists[source2]
end

function build_heuristic_data(X,contsys,Udom,_I_)
    # build the alternating simulation
    hx_coarse = [20.0, 20.0, 1.2, 2*π]
    hx_medium = [5.0, 5.0, 0.4, 0.4]
    hx_fine = [1.0, 1.0, 0.05, 0.05]
    periodic = [3,4]
    periods = [2*π,2*π]
    T0 = [-π,-π]
    Xdom = D.GeneralDomainList(hx_coarse;periodic=periodic,periods=periods,T0=T0)
    AB.add_set!(Xdom, X, AB.OUTER)
    symmodel = AB.NewSymbolicModelListList(Xdom, Udom)
    problem = AS.symmodelProblem(symmodel,contsys,compute_reachable_set,minimum_transition_cost,AS.get_possible_transitions_2)
    autom = AS.build_alternating_simulation(problem)
    symmodel = AB.with_automaton(symmodel, autom)
    # build the heurisitic
    initlist = U.get_symbol(symmodel,_I_,AB.OUTER)
    heuristic_data = AS.build_heuristic(symmodel,initlist)
    println("Heuristic ended")
    #fig = plot(aspect_ratio = 1,legend = false)
    #UT.plot_domain!(heuristic_data.symmodel.Xdom)
    #AS.plot_heuristic!(heuristic_data)
    #display(fig)
    return heuristic_data
end


function test_astart()
    hx_coarse = [20.0, 20.0, 1.2, 2*π]
    hx_medium = [5.0, 5.0, 0.4, 0.4]
    hx_fine = [1.0, 1.0, 0.05, 0.05]

    X, Xdom = build_dom(hx_coarse)
    Udom = build_input()
    # build system
    contsys = build_system(X)
    symmodel = D._SymbolicModel(Xdom, Udom)
    # control problem
    x0 = SVector(30.0, 5.0,π/2.0,0.0)
    RI = SVector(1.0, 1.0, 0.05, 0.05)
    _I_ = AB.HyperRectangle(SVector(30.0, 5.0,π/2.0,0.0)-RI, SVector(30.0, 5.0,π/2.0,0.0)+RI)
    println("INIT------------------")
    initlist = U.get_symbol(symmodel,_I_,AB.OUTER)
    @show initlist
#    println("TARGET----------------")
#    _T_ = AB.HyperRectangle(SVector(25.0, 12.0,π/2.0-0.1,-0.1), SVector(35.0, 19.0,π/2.0+0.1,0.1))
#    targetlist = U.get_symbol(symmodel,_T_,AB.INNER)
#    # Heuristic data
#    fig = plot(aspect_ratio = 1,legend = false)
#    heuristic_data = build_heuristic_data(X,contsys,Udom,_I_)
#
#    # Lazy Abstraction implementation
#    time = @elapsed begin
#        problem,success = LA.compute_controller(symmodel, contsys, initlist, targetlist, transition_cost, pre_image, post_image, h, heuristic_data=heuristic_data)
#        contr = problem.contr
#        @show success
#    end
#    println("total time: lazy abstraction + controller: ", time)
#    LA.plot_result!(problem,x0=x0)
#    display(fig)
end
