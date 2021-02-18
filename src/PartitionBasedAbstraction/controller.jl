
function Base.isempty(set::AB.SortedTupleSet)
    return isempty(set.data)
end

function localise(partition::Partition,Q::AbstractPolytope)
    for (i,P) in enumerate(partition.L)
        if issubset(Q,P)
            return i
        end
    end
    @warn("Out of domain")
end

function control(X0::AbstractPolytope,Xf::AbstractPolytope,partition::Partition,topological_graph::LightTopologicalGraph,contsys::ControlSystem{N},obstacles,_U_,hu,hx_l) where N
    # We need the H-representation for `is_intersection_empty` so we compute it here
    # to avoid having to compute it many times
    h_obstacles = LazySets.tohrep.(obstacles)
    p0 = localise(partition,X0)
    pf = localise(partition,Xf)
    path = get_path(topological_graph,p0,pf)
    control_list = []
    _I_ = X0
    time_abstraction = []
    time_controller = []
    for i=1:length(path)
        #temporary
        #parameters of the abstraction, could decide here a different state-input space discretization
        hx = SVector{N}(hx_l[i])
        push!(time_abstraction,@elapsed symmodel = build_abstraction(partition.L[path[i]],hx,_U_,hu,contsys,h_obstacles))
        #local inputs
        initlist = get_symbols(_I_,symmodel,AB.OUTER)
        #local targets
        if i==length(path)
            _T_ = Xf
        else
            # donc là pour le moement, je change de controller des que j'atteint une contraction de l'intersection
            _T_ = get_transition(topological_graph,path[i],path[i+1])#expansion(0.7,get_transition(topological_graph,path[i],path[i+1]))
        end
        targetlist = get_symbols(_T_,symmodel,AB.INNER)
        #print_symbolset(symmodel,initlist,targetlist)

        contr = AB.NewControllerList()
        push!(time_controller,@elapsed AB.compute_controller_reach!(contr, symmodel.autom, initlist, targetlist))
        if isempty(contr)
            @warn("Controller empty")
            #do something: refinement of abstraction and/or change
            #the path in topological graph
        else
            push!(control_list,(symmodel,contr,initlist,targetlist)) #or I could send _I_ and _T_
            _I_ = _T_
        end
    end
    println()
    println("Elapsed time")
    println("Abstractions : ",time_abstraction)
    println("Controllers  : ",time_controller)
    println("Total Abstractions : ",sum(time_abstraction))
    println("Total Controllers  : ",sum(time_controller))
    return control_list
end

#execute a given number of control input
function trajectory(contsys::ControlSystem{N,T}, symmodel::AB.SymbolicModelList{N}, contr::AB.SortedTupleSet{2,Int}, x0, nstep; randchoose = false, display = false) where{N,T}
    traj = []
    for i in 1:nstep
        push!(traj,x0)
        xpos = AB.get_pos_by_coord(symmodel.Xdom.grid, x0)
        if !(xpos ∈ symmodel.Xdom)
            @warn("Trajectory out of domain")
            return
        end
        source = AB.get_state_by_xpos(symmodel, xpos)
        println(source)
        symbollist = AB.fix_and_eliminate_first(contr, source)
        if isempty(symbollist)
            @warn("Uncontrollable state")
            return
        end
        if randchoose
            symbol = rand(collect(symbollist))[1]
        else
            symbol = first(symbollist)[1]
        end
        if(display)
            plot!(Singleton(x0),color = :red)
        end
        upos = AB.get_upos_by_symbol(symmodel, symbol)
        u = AB.get_coord_by_pos(symmodel.Udom.grid, upos)
        x0 = contsys.sys_map(x0, u, contsys.tstep)
    end
    return traj
end

#execute the control input until reach a target point
function trajectory_reach(contsys::ControlSystem{N,T}, symmodel::AB.SymbolicModelList{N}, contr::AB.SortedTupleSet{2,Int}, x0, targetlist; randchoose = false) where {N,T}
    traj = []
    while true
        push!(traj,x0)
        xpos = AB.get_pos_by_coord(symmodel.Xdom.grid, x0)
        if !(xpos ∈ symmodel.Xdom)
            @warn("Trajectory out of domain")
            return (traj,false)
        end
        source = AB.get_state_by_xpos(symmodel, xpos)
        if source ∈ targetlist
            break
        end
        symbollist = AB.fix_and_eliminate_first(contr, source)
        if isempty(symbollist)
            @warn("Uncontrollable state")
            return (traj,false)
        end
        if randchoose
            symbol = rand(collect(symbollist))[1]
        else
            symbol = first(symbollist)[1]
        end

        upos = AB.get_upos_by_symbol(symmodel, symbol)
        u = AB.get_coord_by_pos(symmodel.Udom.grid, upos)
        x0 = contsys.sys_map(x0, u, contsys.tstep)
    end
    return (traj,true)
end

#execute the control input until reach a target point
function trajectory_reach_full(contsys::ControlSystem{N,T},control_list,x0; randchoose = false) where {N,T}
    traj = []
    for (symmodel, contr, initlist, targetlist) in control_list
        (trajectory,success) = trajectory_reach(contsys, symmodel, contr, x0, targetlist; randchoose)
        push!(traj,trajectory)
        x0 = trajectory[end]
        if !success
            break
        end
    end
    return traj
end

function rectangle(c,r)
    Shape(c[1].-r[1] .+ [0,2*r[1],2*r[1],0], c[2].-r[2] .+ [0,0,2*r[2],2*r[2]])
end

function print_trajectory!(control_data,traj;full=false) where {N,T}
    (symmodel,contr,initlist,targetlist) = control_data
    domain = symmodel.Xdom
    grid = domain.grid
    h = grid.h
    if full == true
        for pos in domain.elems
            center = AB.get_coord_by_pos(grid, pos)
            plot!(rectangle(center,h./2), opacity=.2,color=:blue)
        end
        print_list_state(symmodel,targetlist,:red)
    end
    for i=1:length(traj)-1
        plot!([traj[i][1],traj[i+1][1]], [traj[i][2],traj[i+1][2]],color =:red,linewidth = 2)
        if i>1
            scatter!([traj[i][1]],[traj[i][2]],color =:red,markersize=2)
        end
    end
    scatter!([traj[1][1]],[traj[1][2]],color =:green,markersize=3)
    scatter!([traj[end][1]],[traj[end][2]],color =:yellow,markersize=3)
end

function print_trajectory_full(control_list,traj;full=false,partition=nothing,obstacles=[],X0=nothing,Xf=nothing,save = nothing)
    fig = plot(aspect_ratio = 1,legend = false)
    if partition != nothing
        plot!(partition.L)
    end
    if X0 != nothing
        plot!(X0)
    end
    if Xf != nothing
        plot!(Xf)
    end
    for obs in obstacles
        plot!(obs,color=:black)
    end
    if full
        print_list_state(control_list[1][1],control_list[1][3],:green)
    end
    for i=1:length(traj)
        print_trajectory!(control_list[i],traj[i];full=full)
    end
    display(fig)
    if save != nothing
        savefig(fig, save)
    end
end

function print_symbolset(symmodel,contr,source)
    symbollist = AB.fix_and_eliminate_first(contr, source)

    upos = AB.get_upos_by_symbol(symmodel, symbol)
    u = AB.get_coord_by_pos(symmodel.Udom.grid, upos)
    x0 = contsys.sys_map(x0, u, contsys.tstep)
    plot!(Singleton(x0),color=:yellow)
end


function print_list_state(symmodel,list,col)
    for state in list
        pos = AB.get_xpos_by_state(symmodel, state)
        grid = symmodel.Xdom.grid
        h = grid.h
        center = AB.get_coord_by_pos(grid, pos)
        plot!(rectangle(center,h./2), opacity=.2,color=col)
    end
end
