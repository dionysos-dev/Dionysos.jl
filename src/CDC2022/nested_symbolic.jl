
struct Param
    # MC::Bool
    N::Int
    # post_image::Function
    # pre_image::Function
end

mutable struct NestedSymbolicModel{B,A,N,M} #<: DO.SymbolicModel{N,M}
    Xdom::NestedDomain
    Udom::B
    autom::A
    xpos2int::Dict{Tuple{Int,NTuple{N,Int}},Int}
    xint2pos::Vector{Tuple{Int,NTuple{N,Int}}}
    uint2input::Vector{M}
    active::Vector{Bool}
    mc#::MarkovChain
    param::Param
end


function NewNestedSymbolicModel(dom::DO.GeneralDomainList, Udom, param)
    N = DO.get_dim(dom)
    Xdom = NestedDomain(dom)
    nu = DO.get_ncells(Udom)
    uint2pos = [input for input in DO.enum_pos(Udom)]
    #upos2int = Dict((pos, i) for (i, pos) in enumerate(DO.enum_pos(Udom)))
    symmodel = NestedSymbolicModel(
        Xdom,
        Udom,
        NewProbaAutomaton(0, nu),
        Dict{Tuple{Int,NTuple{N,Int}},Int}(),
        Tuple{Int,NTuple{N,Int}}[],
        uint2pos,
        Bool[],
        nothing,
        param
    )
end


function with_automaton(symmodel::NestedSymbolicModel, autom)
    return NestedSymbolicModel(
        symmodel.Xdom,
        symmodel.Udom,
        autom,
        symmodel.xpos2int,
        symmodel.xint2pos,
        symmodel.upos2int,
        symmodel.uint2pos,
    )
end


function get_ncells(symmodel::NestedSymbolicModel)
    return length(findall(symmodel.active))
end

function get_all_ncells(symmodel::NestedSymbolicModel)
    return length(symmodel.xint2pos)
end

function get_cells(symmodel::NestedSymbolicModel)
    return findall(symmodel.active)
end

function get_state_by_xpos(
    symmodel::NestedSymbolicModel,pos,l
) where {N,M,T}
    dom = symmodel.Xdom.domains[l]
    pos = DO.set_in_period_pos(dom,pos)
    id = get(symmodel.xpos2int, (l,pos), nothing)
    created = false
    if id === nothing
        if pos in dom
            push!(symmodel.active,true)
            created = true
            push!(symmodel.xint2pos, (l,pos))
            id = length(symmodel.xint2pos)
            symmodel.xpos2int[(l,pos)] = id
            i = HybridSystems.add_state!(symmodel.autom)#AB.HybridSystems.add_state!(symmodel.autom)
            @assert i == id
        else
            error("$pos is not in state domain $(symmodel.Xdom)")
        end
    end
    return id::Int
end


function get_xpos_by_state(symmodel::NestedSymbolicModel, state::Int)
    return symmodel.xint2pos[state]
end

function get_state_by_coord(symmodel, coord)
    Ndomain = symmodel.Xdom
    pos, l =  DO.get_pos_by_coord(Ndomain, coord)
    if is_pos(Ndomain, pos, l)
        s = get_state_by_xpos(symmodel, pos, l)
        if in(symmodel, s)
            return s
        end
    else
        error("outside domain !")
        return 0
    end
    return 0
end

function get_coord_by_state(symmodel::NestedSymbolicModel, state)
    (l,pos) = get_xpos_by_state(symmodel, state)
    return DO.get_coord_by_pos(symmodel.Xdom,l,pos)
end

function get_upos_by_symbol(symmodel::NestedSymbolicModel, symbol::Int)
    return symmodel.uint2pos[symbol]
end

function get_symbol_by_upos(symmodel::NestedSymbolicModel, upos)
    return symmodel.upos2int[upos]
end

function get_volume(symmodel,state)
    (l,pos) = get_xpos_by_state(symmodel, state)
    return DO.get_volume(symmodel.Xdom.domains[l].grid)
end

function delete_state!(symmodel::NestedSymbolicModel,s::Int)
    #should we delete in probaautomaton the transition related to node s, no should be done before deleting the node.
    symmodel.active[s] = false
end

function Base.in(symmodel::NestedSymbolicModel, s::Int)
    return symmodel.active[s]
end

function delete_transitions_post!(symmodel, source, symbol)
    delete_transition_post!(symmodel.autom,source,symbol)
end

function cut_cell!(symmodel, s)
    l,pos = get_xpos_by_state(symmodel, s)
    subpos = cut_pos!(symmodel.Xdom, pos, l)
    delete_state!(symmodel, s)
    subcells = Int[]
    for spos in subpos
        subcell = get_state_by_xpos(symmodel, spos, l+1)
        push!(subcells, subcell)
    end
    return subcells
end

function get_transitions_pre(symmodel::NestedSymbolicModel, s::Int)
    active = symmodel.active
    translist = pre(symmodel.autom, s)
    list = []
    for e in translist
        if active[e[1]]
            push!(list,e)
        end
    end
    return list
end

function get_transitions_post(symmodel::NestedSymbolicModel, s::Int)
    active = symmodel.active
    translist = post(symmodel.autom, s)
    list = []
    for e in translist
        if e[1] == 0 || active[e[1]]
            push!(list,e)
        end
    end
    return list
end

function update_ingoing_transitions_MC!(symmodel,sys,s,subcells)
    autom = symmodel.autom
    for (source,symbol,proba) in get_transitions_pre(symmodel, s)
        if source != s
            delete_transitions_post!(symmodel, source, symbol)
            u = DO.enum_pos(symmodel.Udom)[symbol]
            transitions = get_transitions_MC(symmodel,sys,source,symbol,u)
            add_transitions!(symmodel.autom, transitions)
        end
    end
end

function update_outgoing_transitions_MC!(symmodel,sys,s,subcells)
    autom = symmodel.autom
    for source in subcells
        for (symbol,u) in enumerate(DO.enum_pos(symmodel.Udom))
            transitions = get_transitions_MC(symmodel,sys,source,symbol,u)
            add_transitions!(autom, transitions)
        end
    end
end

function split_node!(symmodel,sys,s)
    subcells = cut_cell!(symmodel, s)
    update_ingoing_transitions_MC!(symmodel,sys,s,subcells)
    update_outgoing_transitions_MC!(symmodel,sys,s,subcells)
end

function Plots.plot(symmodel::NestedSymbolicModel;dims=[1,2],annotate=false)
    fig = plot(aspect_ratio = 1,legend = false)
    plot!(symmodel;dims=dims,annotate=annotate)
    display(fig)
end

function Plots.plot!(symmodel::NestedSymbolicModel;dims=[1,2],annotate=false)
    colors = :yellow
    for value in symmodel.xint2pos
        (l,pos) = value
        dom = symmodel.Xdom.domains[l]
        grid = dom.grid
        if annotate
           s = get_state_by_xpos(symmodel, pos, l)
           center = DO.get_coord_by_pos(dom.grid, pos)
           annotate!([(center[dims[1]], center[dims[2]], text(s, :red))])
        end
        DO.plot_elem!(grid,pos, opacity=.9,color=colors)
    end
end

function circleShape(h,k,r)
    angle = LinRange(0, 2*pi, 500)
    h .+ r*sin.(angle),k.+r*cos.(angle)
end

function plot_automaton(symmodel,bool=false,src=0)
    fig = plot(aspect_ratio = 1,legend = false)
    plot!(symmodel)
    autom = symmodel.autom
    for s in get_cells(symmodel)
        (l,xpos) = get_xpos_by_state(symmodel, s)
        grid = symmodel.Xdom.domains[l].grid
        x = DO.get_coord_by_pos(grid, xpos)
        r = DO.get_h(grid)/2.0
        plot!(circleShape(x[1],x[2],r[1]/4.0), seriestype = [:shape], lw=1.0,c=:blue,fillalpha=1.0,linecolor=:black)
    end
    for target in get_cells(symmodel)
        (ly,ypos) = get_xpos_by_state(symmodel, target)
        gridy = symmodel.Xdom.domains[ly].grid
        y = DO.get_coord_by_pos(gridy, ypos)
        for (source,symbol,proba) in get_transitions_pre(symmodel, target)
            if !bool || source == src
                (lx,xpos) = get_xpos_by_state(symmodel, source)
                gridx = symmodel.Xdom.domains[lx].grid
                x = DO.get_coord_by_pos(gridx, xpos)
                if source != target
                    d = y-x
                    quiver!([x[1]],[x[2]],quiver=([d[1]],[d[2]]),color=:black,linewidth=2,title=src)
                else
                    h = DO.get_h(gridx)
                    plot!(circleShape(x[1]+h[1]/12.0,x[2]+h[1]/12.0,h[1]/12.0), seriestype = [:shape], lw=2.0,c=:blue,fillalpha=0.0,linecolor=:black)
                end
            end
        end
    end
    display(fig)
end

function plot_trajectory!(sys,x0,l;color=:red,dims=[1,2])
    x = copy(x0)
    u = SVector(0.0,0.0)
    for i=1:l
        Fx = sys.sys_map(x, u, sys.tstep)
        plot!([x[dims[1]],Fx[dims[1]]], [x[dims[2]],Fx[dims[2]]],color = color,linewidth = 2)
        if i>1
            scatter!([x[dims[1]]],[x[dims[2]]],color =:grey,markersize=1) #:yellow
        end
        x = Fx
    end
    scatter!([x0[dims[1]]],[x0[dims[2]]],color =:green,markersize=2)
    scatter!([x[dims[1]]],[x[dims[2]]],color =:red,markersize=2)
    scatter!([0.0],[0.0],color =:red,markersize=2)
end

function get_transitions(symmodel,sys,source,symbol,u)
    return get_transitions_MC(symmodel,sys,source,symbol,u)
end

function compute_transition!(symmodel,sys,source)
    for (symbol,u) in enumerate(DO.enum_pos(symmodel.Udom))
        transitions = get_transitions(symmodel,sys,source,symbol,u)
        add_transitions!(symmodel.autom, transitions)
    end
end

function compute_symbolic_full_domain!(symmodel, sys)
    for pos in DO.enum_pos(symmodel.Xdom.domains[1])
        source = get_state_by_xpos(symmodel,pos,1)
        compute_transition!(symmodel,sys,source)
    end
end









function update_MC!(symmodel::NestedSymbolicModel)
    #symmodel.mc = build_Markov_Chain(get_cells(symmodel),symmodel.autom.transitions)
    symmodel.mc = build_Markov_Chain(symmodel)#build_Markov_Chain(symmodel)
end

function plot_shannon_entropy(sys,symmodel::NestedSymbolicModel;dims=[1,2],xlims=[],ylims=[])
    fig = plot(aspect_ratio = 1,legend = false,title="outgoing shannon-entropy",xlims=xlims,ylims=ylims)
    plot_shannon_entropy!(symmodel;dims=dims)
    x0 = SVector(2.0,2.0)
    plot_trajectory!(sys,x0,200)
    l = 150
    plot_trajectory!(sys,SVector(4.0,0.8),l;color=:green)
    plot_trajectory!(sys,SVector(-4.0,-0.8),l;color=:green)
    display(fig)
end

function plot_shannon_entropy!(symmodel::NestedSymbolicModel;dims=[1,2],xlims=[],ylims=[])
    mc = symmodel.mc
    max_val = max(mc.shannon_entropy...)
    for s in get_cells(symmodel)
        (l,pos) = get_xpos_by_state(symmodel, s)
        x = get_coord_by_state(symmodel,s)
        if abs(x[1])<10 && x[2]<20 && x[2]>-8.0
            dom = symmodel.Xdom.domains[l]
            val = get_shannon_entropy(mc, s)
            opacity = val/max_val
            if opacity > 0.01
                DO.plot_elem!(dom.grid, pos, dims=dims,opacity=opacity*2.0,color=:yellow)
            end
        end
    end
end

function plot_steady_state(sys,symmodel::NestedSymbolicModel;dims=[1,2],fact=1.0,tol=0.0,xlims=[],ylims=[],x0=nothing)
    fig = plot(aspect_ratio = 1,legend = false,title="steady-state probability",xlims=xlims,ylims=ylims,gridlinewidth=0)# xlims=[-5.0,5.0]
    plot_steady_state!(symmodel;dims=dims,fact=fact,tol=tol)
    if x0 != nothing
        plot_trajectory!(sys,x0,300,dims=dims) #150
    end
    display(fig)
end

function plot_steady_state!(symmodel::NestedSymbolicModel;dims=[1,2],fact=1.0,tol=0.0,color=:yellow)
    mc = symmodel.mc
    max_val = max(mc.steady_state...)
    for s in get_cells(symmodel)
        (l,pos) = get_xpos_by_state(symmodel, s)
        dom = symmodel.Xdom.domains[l]
        val = get_steady_state(mc, s)
        opa = val/max_val
        if val>fact
            opa = 1.0
        end
        if opa > tol
            DO.plot_elem!(dom.grid, pos, dims=dims,opacity=opa*12.0,color=color) #*12.0
        end
    end
end

function get_steady_state_region_threshold(symmodel, th::Float64)
    mc = symmodel.mc
    max_val = max(mc.steady_state...)
    region = Dict{Int, Any}()
    for s in get_cells(symmodel)
        (l,pos) = get_xpos_by_state(symmodel, s)
        dom = symmodel.Xdom.domains[l]
        val = get_steady_state(mc, s)
        opacity = val/max_val
        if val>th
            region[s] = true
        end
    end
    return region
end
function fixed_point!(sys,symmodel,region)
    u = SVector(0.0,0.0)
    previous = 0
    next = 1
    while previous != next
        previous = length(region)
        for source in keys(region)
            #transitions = get_transitions_MC(symmodel,sys,source,1,u;tstep=sys.tstep*40.0)
            N = 500
            (l,xpos) = get_xpos_by_state(symmodel, source)
            points = sample_elem(symmodel.Xdom.domains[l].grid, xpos, N)
            targetlist = Int[]
            for x in points
                l = 25
                for i=1:l
                    x = sys.sys_map(x, u, sys.tstep)
                end
                target = get_state_by_coord(symmodel,x)
                push!(targetlist,target)
            end
            for target in targetlist
                if !haskey(region,target)
                    pop!(region,source)
                    println("LKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK")
                    println(get_coord_by_state(symmodel,source))
                    println(source)
                    #region[source] = false
                    break
                end
            end
        end
        next = length(region)
    end
end

function plot_region(symmodel,sys,region)
    fig = plot(aspect_ratio = 1,legend = false,title="region")
    plot_region!(symmodel,region)
    x0 = SVector(2.0,2.0)
    plot_trajectory!(sys,x0,200)
    display(fig)
end

function plot_region!(symmodel,region)
    for s in keys(region)
        println(s)
        (l,pos) = get_xpos_by_state(symmodel, s)
        dom = symmodel.Xdom.domains[l]
        DO.plot_elem!(dom.grid, pos,opacity=1,color=:yellow)
    end
end

function plot_entropy(sys,symmodel::NestedSymbolicModel;dims=[1,2])
    fig = plot(aspect_ratio = 1,legend = false,title="entropy")
    plot_entropy!(symmodel;dims=dims)
    x0 = SVector(2.0,2.0)
    plot_trajectory!(sys,x0,200)
    display(fig)
end

function plot_entropy!(symmodel::NestedSymbolicModel;dims=[1,2])
    mc = symmodel.mc
    max_val = max(mc.entropy...)
    for s in  get_cells(symmodel)
        (l,pos) = get_xpos_by_state(symmodel, s)
        dom = symmodel.Xdom.domains[l]
        val = get_entropy(mc, s)
        opacity = val/max_val
        if opacity > 0.01
            DO.plot_elem!(dom.grid, pos, dims=dims,opacity=opacity,color=:yellow)
        end
    end
end

function plot_SCC(sys,symmodel::NestedSymbolicModel;dims=[1,2])
    mc = symmodel.mc
    fig = plot(aspect_ratio = 1,legend = false,title="SCC recurrent")
    for (i,class) in enumerate(mc.SCC[1])
        if mc.SCC[2][i]
            for (j,idx) in enumerate(class)
                s = mc.symmodel_from_mc[idx]
                (l,pos) = get_xpos_by_state(symmodel, s)
                dom = symmodel.Xdom.domains[l]
                DO.plot_elem!(dom.grid, pos, dims=dims,color=i)
            end
        end
    end
    x0 = SVector(2.0,2.0)
    plot_trajectory!(sys,x0,200)
    display(fig)
    fig = plot(aspect_ratio = 1,legend = false,title="SCC non recurrent ")
    for (i,class) in enumerate(mc.SCC[1])
        if !mc.SCC[2][i]
            for (j,idx) in enumerate(class)
                s = mc.symmodel_from_mc[idx]
                (l,pos) = get_xpos_by_state(symmodel, s)
                dom = symmodel.Xdom.domains[l]
                DO.plot_elem!(dom.grid, pos, dims=dims,color=i)
            end
        end
    end
    x0 = SVector(2.0,2.0)
    plot_trajectory!(sys,x0,200)
    display(fig)
end


function plot_vector_field(symmodel::NestedSymbolicModel,f;dims=[1,2])
    fig = plot(aspect_ratio = 1,legend = false)
    plot!(symmodel,dims=dims)
    for (i,value) in enumerate(symmodel.xint2pos)
        if i%1 == 0
            (l,pos) = value
            dom = symmodel.Xdom.domains[l]
            grid = dom.grid
            center = DO.get_coord_by_pos(grid, pos)
            d = f(center,0.0)
            d =d/norm(d,2)
            quiver!([center[1]],[center[2]],quiver=([d[1]],[d[2]]),color=:black,linewidth=1)
        end
    end
    display(fig)
end

function plot_map(symmodel::NestedSymbolicModel,f;dims=[1,2])
    fig = plot(aspect_ratio = 1,legend = false)
    tab = []
    val_max = -Inf
    for (i,value) in enumerate(symmodel.xint2pos)
        (l,pos) = value
        dom = symmodel.Xdom.domains[l]
        grid = dom.grid
        center = DO.get_coord_by_pos(grid, pos)
        val = f(center,0.0)
        val_max = max(val_max,val)
        push!(tab,(pos,val))
    end
    grid = symmodel.Xdom.domains[1].grid
    for (pos,val) in tab
        opacity = val/val_max
        # println()
        # println(val)
        # println(opacity)
        DO.plot_elem!(grid, pos, opacity=opacity,color=:yellow)
    end
    # x0 = SVector(2.0,2.0)
    # plot_trajectory!(sys,x0,200)
    display(fig)
end


function plot_Jacobian(sys,symmodel::NestedSymbolicModel,Jacobian;dims=[1,2])
    fig = plot(aspect_ratio = 1,legend = false,title="Stability")
    for (i,value) in enumerate(symmodel.xint2pos)
        (l,pos) = value
        dom = symmodel.Xdom.domains[l]
        grid = dom.grid
        center = DO.get_coord_by_pos(grid, pos)
        J = Jacobian(center,0.0)
        eigenvalues = eigvals(Array(J))
        r1 = real(eigenvalues[1])
        r2 = real(eigenvalues[2])
        if r1<0 && r2<0
            DO.plot_elem!(grid, pos, opacity=1.0,color=:green)
        elseif r1>0 && r2>0
            DO.plot_elem!(grid, pos, opacity=1.0,color=:red)
        else
            DO.plot_elem!(grid, pos, opacity=1.0,color=:yellow)
        end
    end
    x0 = SVector(2.0,2.0)
    plot_trajectory!(sys,x0,200)
    display(fig)
end
