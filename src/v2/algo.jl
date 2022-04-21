

mutable struct AnalysisProblem{SM,C,PoI,PoII}#PoI<:Function
    sys::C
    symmodel::SM
    post_image::PoI # function to compute the list of potential post-image of a cell for a given input
    pre_image::PoII # function to compute the list of potential post-image of a cell for a given input
    GB::Bool
end

function NewAnalysisProblem(symmodel::MultiSymbolicModel,sys::AB.ControlSystem{N,T},post_image,pre_image,GB) where {N,T}
    return AnalysisProblem(sys,symmodel,post_image,pre_image,GB)
end


#########################################################################

function sample_from_rec(rec,N)
    L = []
    n = length(rec.lb)
    for i in 1:n
        unifdist = Uniform(rec.lb[i], rec.ub[i])
        push!(L, Distributions.rand(unifdist, N))
    end
    points = []
    for i in 1:N
        vec = Float64[L[j][i] for j in 1:n]
        push!(points,SVector{n,Float64}(vec))
    end
    return points
end

function get_cell(symmodel,x)
    Xdom = symmodel.Xdom
    for l in 1:Xdom.levels
        xpos = AB.get_pos_by_coord(Xdom, l, x)
        if is_pos(Xdom,xpos,l)
            s = get_state_by_xpos(symmodel, xpos, l)
            if in(symmodel,s)
                return s
            end
        else
            #println("outside domain !")
            return 0
        end
    end
    return 0
end

function count_occurences(tab)
    sort!(tab)
    n = length(tab)
    occurences = Tuple{Int,Int}[]
    val = tab[1]
    count = 1
    for i in 2:n
        if tab[i] == val
            count+=1
        else
            push!(occurences,(val,count))
            val = tab[i]
            count = 1
        end
    end
    push!(occurences,(val,count))
    return occurences
end

function sample_elem(grid::AB.GridFree,xpos,N::Int)
    x = AB.get_coord_by_pos(grid, xpos)
    r = grid.h/2
    rec = AB.HyperRectangle(x .- r, x .+ r)
    return sample_from_rec(rec,N)
end

function sample_elem(Dgrid::AB.DeformedGrid,xpos,N::Int)
    points = sample_elem(Dgrid.grid, xpos, N)
    return [Dgrid.f(x) for x in points]
end

# Monte Carlo sampling
function get_transitions_MC(symmodel,sys,source,symbol,u)
    N = 500
    (l,xpos) = get_xpos_by_state(symmodel, source)
    points = sample_elem(symmodel.Xdom.domains[l].grid, xpos, N)
    targetlist = Int[]
    for x in points
        Fx = sys.sys_map(x, u, sys.tstep)
        target = get_cell(symmodel,Fx)
        push!(targetlist,target)
    end
    occurences = count_occurences(targetlist)
    probaList = [(target,occ/N) for (target,occ) in occurences]
    transitions = []
    for (target,proba) in probaList
        push!(transitions, (target, source, symbol, proba))
    end
    return transitions
end

function update_ingoing_transitions_MC!(prob,s,subcells)
    symmodel = prob.symmodel
    autom = symmodel.autom
    sys = prob.sys
    for (source,symbol,proba) in get_transitions_pre(symmodel, s)
        if source != s
            delete_transitions_post!(symmodel, source, symbol)
            u = AB.enum_pos(symmodel.Udom)[symbol]
            transitions = get_transitions_MC(symmodel,sys,source,symbol,u)
            add_transitions!(symmodel.autom, transitions)
        end
    end
end

function update_outgoing_transitions_MC!(prob,s,subcells)
    symmodel = prob.symmodel
    autom = symmodel.autom
    sys = prob.sys
    for source in subcells
        for (symbol,u) in enumerate(AB.enum_pos(symmodel.Udom))
            transitions = get_transitions_MC(symmodel,sys,source,symbol,u)
            add_transitions!(autom, transitions)
        end
    end
end
#########################################################################
# Growth-bound function

function fraction_volume(rec1::AB.HyperRectangle,rec2::AB.HyperRectangle)
    recI = AB.intersect(rec1,rec2)
    vI = AB.volume(recI)
    vrec = AB.volume(rec1)
    return vI/vrec
end

function compute_proba(symmodel,sys,source,target,symbol)
    dom = symmodel.Xdom.domains[1]
    rec1 = get_post_rec(symmodel, sys, source, symbol)
    rec1list = D.set_rec_in_period(dom.periodic,dom.periods,dom.T0,rec1)
    rec2 = get_rec(symmodel, target)
    v_I = 0.0
    v_post = 0.0
    for rec in rec1list
        recI = AB.intersect(rec,rec2)
        v_I += AB.volume(recI)
        v_post += AB.volume(rec)
    end
    return v_I/v_post
end

function get_post_rec(symmodel, sys, source, symbol)
    u = AB.enum_pos(symmodel.Udom)[symbol]
    (l,xpos) = get_xpos_by_state(symmodel, source)
    dom = symmodel.Xdom.domains[l]
    grid = dom.grid
    x = AB.get_coord_by_pos(grid, xpos)
    tstep = sys.tstep
    Fx = sys.sys_map(x, u, tstep)
    r = dom.grid.h/2.0 + sys.measnoise
    Fr = r
    return AB.HyperRectangle(Fx .- Fr, Fx .+ Fr)
end

function get_transitions_GB(prob,symmodel,sys,source,symbol,u)
    (l,xpos) = get_xpos_by_state(symmodel, source)
    over_approx = prob.post_image(symmodel,sys,xpos,l,u)
    transitions = Tuple{Int,Int,Int,Float64}[]
    for target in over_approx
        proba = compute_proba(symmodel,sys, source, target, symbol)
        push!(transitions, (target, source, symbol, proba))
    end
    return transitions
end


function update_ingoing_transitions_GB!(prob,s,subcells)
    symmodel = prob.symmodel
    autom = symmodel.autom
    for (source,symbol,proba) in get_transitions_pre(symmodel,s)
        if source != s
            # get post rec of source but does'nt take into account the periods
            rec = get_post_rec(symmodel, prob.sys, source, symbol)
            # rec in the levels of subcells
            (l,xpos) = get_xpos_by_state(symmodel, s)
            dom = symmodel.Xdom.domains[l+1]

            list_rect_post = D.set_rec_in_period(dom.periodic,dom.periods,dom.T0,rec)

            # rec of subcells
            lbI = Tuple([p*2 for p in xpos])
            ubI = Tuple([p*2+1 for p in xpos])
            rect_subcells = AB.HyperRectangle(lbI,ubI) #this rec is in the domain even if periodic
            # intersection of pos in subcells level
            over_approx = Int[]
            for rec in list_rect_post
                rect_post = AB.get_pos_lims_outer(dom.grid, rec)
                rectI = intersect(rect_post, rect_subcells)
                ypos_iter = Iterators.product(AB._ranges(rectI)...)
                for ypos in ypos_iter
                    ypos = D.set_in_period_pos(dom,ypos)
                    target = get_state_by_xpos(symmodel, ypos, l+1)
                    push!(over_approx, target)
                end
            end
            sort!(over_approx)
            unique!(over_approx)
            targetlist = Tuple{Int,Int,Int,Float64}[]
            for target in over_approx
                proba = compute_proba(symmodel, prob.sys, source, target, symbol)
                push!(targetlist, (target, source, symbol, proba))
            end
            add_transitions!(symmodel.autom, targetlist)
        end
    end
end

function update_outgoing_transitions_GB!(prob,s,subcells)
    symmodel = prob.symmodel
    autom = symmodel.autom
    targetlist = get_transitions_post(symmodel, s)
    (l,pos) = get_xpos_by_state(symmodel, s)
    dom = symmodel.Xdom.domains[l+1]
    rect_subcells = get_subpos(pos)
    # compute transitions from subcells to previous outgoing cells
    for xpos in Iterators.product(AB._ranges(rect_subcells)...)
        source = get_state_by_xpos(symmodel, xpos, l+1)
        rec_post = get_post_rec(symmodel, prob.sys, source, 1)
        list_rect_post = D.set_rec_in_period(dom.periodic,dom.periods,dom.T0,rec_post)
        over_approx = Int[]
        for (target,symbol,proba) in targetlist
            (ly,ypos) = get_xpos_by_state(symmodel, target)
            ydom = symmodel.Xdom.domains[ly]
            yrec = AB.HyperRectangle(ypos,ypos)
            for rec in list_rect_post
                rect_post = AB.get_pos_lims_outer(ydom.grid, rec)
                if AB.isintersect(rect_post,yrec)
                    push!(over_approx, target)
                    break
                end
            end
        end
        # compute the transitions from subcells to subcells
        rectI = AB.get_pos_lims_outer(dom.grid,rec_post)
        ypos_iter = Iterators.product(AB._ranges(rectI)...)
        for rec in list_rect_post
            rect_post = AB.get_pos_lims_outer(dom.grid, rec)
            rectI = intersect(rect_post, rect_subcells)
            ypos_iter = Iterators.product(AB._ranges(rectI)...)
            for ypos in ypos_iter
                ypos = D.set_in_period_pos(dom,ypos)
                target = get_state_by_xpos(symmodel, ypos, l+1)
                push!(over_approx, target)
            end
        end
        # add transitions
        sort!(over_approx)
        unique!(over_approx)
        targetlist_source = Tuple{Int,Int,Int,Float64}[]
        for target in over_approx
            proba = compute_proba(symmodel, prob.sys, source, target, 1)
            push!(targetlist_source, (target, source, 1, proba))
        end
        add_transitions!(symmodel.autom, targetlist_source)
    end
end


#########################################################################
function get_rec(symmodel, source)
    (l,xpos) = get_xpos_by_state(symmodel, source)
    grid = symmodel.Xdom.domains[l].grid
    x = AB.get_coord_by_pos(grid, xpos)
    r = grid.h/2
    return AB.HyperRectangle(x .- r, x .+ r)
end

function get_subpos(pos)
    lbI = Tuple([p*2 for p in pos])
    ubI = Tuple([p*2+1 for p in pos])
    rectI = AB.HyperRectangle(lbI,ubI)
    return rectI
end

function in_obstacle(prob,s)
    l,xpos = get_xpos_by_state(prob.symmodel, s)
    grid = prob.symmodel.Xdom.domains[l].grid
    x = AB.get_coord_by_pos(grid, xpos)
    r = grid.h/2.0
    rec = AB.HyperRectangle(x .- r, x .+ r)
    return AB.isintersect(prob.obstacle,rec)
end

function cut_cell!(prob::AnalysisProblem, s)
    l,pos = get_xpos_by_state(prob.symmodel, s)
    symmodel = prob.symmodel
    if prob.symmodel.Xdom.levels < l+1
        add_sub_dom!(symmodel.Xdom)
    end
    dom = symmodel.Xdom.domains[l+1]
    delete_state!(symmodel,s)
    subpos = Iterators.product(AB._ranges(get_subpos(pos))...)
    subcells = Int[]
    for pos in subpos
        subcell = get_state_by_xpos(symmodel,pos,l+1)
        push!(subcells,subcell)
        #push!(prob.isobstacle, in_obstacle(prob,subcell))
    end
    return subcells
end

function get_transitions(prob,symmodel,sys,source,symbol,u)
    if prob.GB
        return get_transitions_GB(prob,symmodel,sys,source,symbol,u)
    else
        return get_transitions_MC(symmodel,sys,source,symbol,u)
    end
end

function compute_transition!(prob,symmodel,sys,source)
    for (symbol,u) in enumerate(AB.enum_pos(symmodel.Udom))
        transitions = get_transitions(prob,symmodel,sys,source,symbol,u)
        add_transitions!(symmodel.autom, transitions)
    end
end

function compute_symbolic_full_domain!(prob,symmodel,sys)
    for pos in AB.enum_pos(symmodel.Xdom.domains[1])
        source = get_state_by_xpos(symmodel,pos,1)
        #push!(prob.isobstacle, in_obstacle(prob,source))
        compute_transition!(prob,symmodel,sys,source)
    end
end

function update_ingoing_transitions!(prob,s,subcells)
    if prob.GB
        update_ingoing_transitions_GB!(prob,s,subcells)
    else
        update_ingoing_transitions_MC!(prob,s,subcells)
    end
end

function update_outgoing_transitions!(prob,s,subcells)
    if prob.GB
        update_outgoing_transitions_GB!(prob,s,subcells)
    else
        update_outgoing_transitions_MC!(prob,s,subcells)
    end
end


function get_cell_to_cut(prob)
    return get_highest_entropies(prob.symmodel)
end

function split_node(prob,s)
    subcells = cut_cell!(prob, s)
    update_ingoing_transitions!(prob,s,subcells)
    update_outgoing_transitions!(prob,s,subcells)
    #print_automaton(prob.symmodel.autom)
end


function get_volume_safe_set(prob)
    v = 0.0
    symmodel = prob.symmodel
    for (source, symbol) in prob.safe
        v += AB.volume(get_rec(symmodel, source))
    end
    return v
end


function plot_some(symmodel,source,cells)
    fig = plot(aspect_ratio = 1,legend = false,title="steady-state probability")
    plot_steady_state!(symmodel,fact=0.0)
    Xdom = symmodel.Xdom
    for cell in cells
        if cell != 0
            (l,xpos) = get_xpos_by_state(symmodel, cell)
            dom = Xdom.domains[l]
            AB.plot_elem!(dom.grid, xpos,opacity=0.5,color=:green)
        end
    end
    (l,xpos) = get_xpos_by_state(symmodel, source)
    dom = Xdom.domains[l]
    AB.plot_elem!(dom.grid, xpos,opacity=1.0,color=:red)
    display(fig)
end

function solve!(prob::AnalysisProblem)
    symmodel = prob.symmodel
    sys = prob.sys
    println("Symbolic model")
    compute_symbolic_full_domain!(prob,symmodel,sys)
    println("Markov chain")
    update_MC!(prob.symmodel)
    #print_automaton(prob.symmodel.autom)

    #plot_automaton(prob)
    #plot(prob.symmodel)
    for k=1:0

        # s = get_cell_to_cut(prob)
        # println(s)
        # #println(prob.symmodel.mc)
        # split_node(prob,s)
        #plot(prob.symmodel)
        plot_steady_state(sys,prob.symmodel)
        plot_steady_state(sys,prob.symmodel,bool=true)
        plot_shannon_entropy(sys,prob.symmodel)
        #plot_entropy(sys,prob.symmodel)
        #plot_SCC(sys,prob.symmodel)
        #plot_shannon_entropy(prob.symmodel,300)
        #plot_automaton(prob)
    end
end

# pairstable_ij est true si etant en etat i et appliquant input j, je n'arrive pas dans le unsafe set
function _compute_pairstable(pairstable, autom, safelist)
    for target in safelist
        for soursymb in pre(autom, target)
            if soursymb[1] in safelist
                pairstable[soursymb[1], soursymb[2]] = true
            end
        end
    end
end

#safelist: list des cells blue
function compute_controller_safe!(prob, autom, initlist, safelist)
    activelist = findall(prob.active)

    isobstacle = prob.isobstacle
    safelist = filter(e->!isobstacle[e], activelist)
    println("compute_controller_safe! started")
    nstates = autom.nstates
    nsymbols = autom.nsymbols
    pairstable = [false for i in 1:nstates, j in 1:nsymbols]
    _compute_pairstable(pairstable, autom, safelist)
    nsymbolslist = sum(pairstable, dims = 2)
    safeset = Set(safelist)
    for source in safeset
        if nsymbolslist[source] == 0
            delete!(safeset, source)
        end
    end
    unsafeset = Set(activelist)
    setdiff!(unsafeset, safeset)
    for source in unsafeset
        for symbol in 1:nsymbols
            pairstable[source, symbol] = false
        end
    end
    nextunsafeset = Set{Int}()
    while true
        for target in unsafeset
            for soursymb in pre(autom, target)
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
    for source in safeset
        for symbol in 1:nsymbols
            if pairstable[source, symbol]
                prob.safe[source] = symbol
            end
        end
    end
end

function circleShape(h,k,r)
    angle = LinRange(0, 2*pi, 500)
    h .+ r*sin.(angle),k.+r*cos.(angle)
end

function plot_rotatedRectangle(rec,point,θ;dims=[1,2])
    d1 = dims[1]; d2 = dims[2]
    vertices = [SVector(rec.lb[d1],rec.lb[d2]), SVector(rec.lb[d1],rec.ub[d2]), SVector(rec.ub[d1],rec.lb[d2]), SVector(rec.ub[d1],rec.ub[d2])]
    v_r = [AB.rotate((x-point),θ)+point for x in vertices]
    plot!(VPolytope(v_r))
end

function plot_automaton(prob::AnalysisProblem;bool=false,src=0)
    fig = plot(aspect_ratio = 1,legend = false)
    symmodel = prob.symmodel
    plot!(symmodel)
    autom = symmodel.autom
    for s in get_cells(symmodel)
        (l,xpos) = get_xpos_by_state(symmodel, s)
        grid = symmodel.Xdom.domains[l].grid
        x = AB.get_coord_by_pos(grid, xpos)
        r = AB.get_h(grid)/2.0
        plot!(circleShape(x[1],x[2],r[1]/4.0), seriestype = [:shape], lw=1.0,c=:blue,fillalpha=1.0,linecolor=:black)
        # if prob.isobstacle[s]
        #     plot!(circleShape(x[1],x[2],r[1]/4.0), seriestype = [:shape], lw=1.0,c=:red,fillalpha=1.0,linecolor=:black)
        # if haskey(prob.safe,s)
        #     plot!(circleShape(x[1],x[2],r[1]/4.0), seriestype = [:shape], lw=1.0,c=:green,fillalpha=1.0,linecolor=:black)
        # else
        #     plot!(circleShape(x[1],x[2],r[1]/4.0), seriestype = [:shape], lw=1.0,c=:blue,fillalpha=1.0,linecolor=:black)
        # end
    end
    for target in get_cells(symmodel)
        (ly,ypos) = get_xpos_by_state(symmodel, target)
        gridy = symmodel.Xdom.domains[ly].grid
        y = AB.get_coord_by_pos(gridy, ypos)
        for (source,symbol,proba) in get_transitions_pre(symmodel, target)
            if !bool || source == src
                (lx,xpos) = get_xpos_by_state(symmodel, source)
                gridx = symmodel.Xdom.domains[lx].grid
                x = AB.get_coord_by_pos(gridx, xpos)
                if source != target
                    d = y-x
                    quiver!([x[1]],[x[2]],quiver=([d[1]],[d[2]]),color=:black,linewidth=2,title=src)
                else
                    h = AB.get_h(gridx)
                    plot!(circleShape(x[1]+h[1]/12.0,x[2]+h[1]/12.0,h[1]/12.0), seriestype = [:shape], lw=2.0,c=:blue,fillalpha=0.0,linecolor=:black)
                end
            end
        end
    end
    display(fig)
end
