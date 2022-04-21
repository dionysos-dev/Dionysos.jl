



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

# function get_cell(symmodel,x)
#     Xdom = symmodel.Xdom
#     for l in 1:Xdom.levels
#         xpos = AB.get_pos_by_coord(Xdom, l, x)
#         if is_pos(Xdom,xpos,l)
#             s = get_state_by_xpos(symmodel, xpos, l)
#             if in(symmodel,s)
#                 return s
#             end
#         else
#             println("outside domain !")
#             return 0
#         end
#     end
#     return 0
# end

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
function get_transitions_MC(symmodel,sys,source,symbol,u;tstep=sys.tstep)
    N = 5000
    (l,xpos) = get_xpos_by_state(symmodel, source)
    points = sample_elem(symmodel.Xdom.domains[l].grid, xpos, N)
    targetlist = Int[]
    for x in points
        Fx = sys.sys_map(x, u, tstep)
        target = get_state_by_coord(symmodel,Fx)
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

function update_ingoing_transitions_MC!(symmodel,sys,s,subcells)
    autom = symmodel.autom
    for (source,symbol,proba) in get_transitions_pre(symmodel, s)
        if source != s
            delete_transitions_post!(symmodel, source, symbol)
            u = AB.enum_pos(symmodel.Udom)[symbol]
            transitions = get_transitions_MC(symmodel,sys,source,symbol,u)
            add_transitions!(symmodel.autom, transitions)
        end
    end
end

function update_outgoing_transitions_MC!(symmodel,sys,s,subcells)
    autom = symmodel.autom
    for source in subcells
        for (symbol,u) in enumerate(AB.enum_pos(symmodel.Udom))
            transitions = get_transitions_MC(symmodel,sys,source,symbol,u)
            add_transitions!(autom, transitions)
        end
    end
end



##########################################################################

function get_subpos(pos)
    lbI = Tuple([p*2 for p in pos])
    ubI = Tuple([p*2+1 for p in pos])
    rectI = AB.HyperRectangle(lbI,ubI)
    return rectI
end


function cut_cell!(symmodel, s)
    l,pos = get_xpos_by_state(symmodel, s)
    if symmodel.Xdom.levels < l+1
        add_sub_dom!(symmodel.Xdom)
    end
    dom = symmodel.Xdom.domains[l+1]
    delete_state!(symmodel,s)
    subpos = Iterators.product(AB._ranges(get_subpos(pos))...)
    subcells = Int[]
    for pos in subpos
        subcell = get_state_by_xpos(symmodel,pos,l+1)
        push!(subcells,subcell)
    end
    return subcells
end

function get_transitions(symmodel,sys,source,symbol,u)
    return get_transitions_MC(symmodel,sys,source,symbol,u)

end

function compute_transition!(symmodel,sys,source)
    for (symbol,u) in enumerate(AB.enum_pos(symmodel.Udom))
        transitions = get_transitions(symmodel,sys,source,symbol,u)
        add_transitions!(symmodel.autom, transitions)
    end
end

function compute_symbolic_full_domain!(symmodel,sys)
    for pos in AB.enum_pos(symmodel.Xdom.domains[1])
        source = get_state_by_xpos(symmodel,pos,1)
        compute_transition!(symmodel,sys,source)
    end
end


function split_node!(symmodel,sys,s)
    subcells = cut_cell!(symmodel, s)
    update_ingoing_transitions_MC!(symmodel,sys,s,subcells)
    update_outgoing_transitions_MC!(symmodel,sys,s,subcells)
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
        x = AB.get_coord_by_pos(grid, xpos)
        r = AB.get_h(grid)/2.0
        plot!(circleShape(x[1],x[2],r[1]/4.0), seriestype = [:shape], lw=1.0,c=:blue,fillalpha=1.0,linecolor=:black)
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
