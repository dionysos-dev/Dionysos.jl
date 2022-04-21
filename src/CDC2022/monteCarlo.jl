
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

function sample_elem(grid::DO.GridFree, xpos, N::Int)
    x = DO.get_coord_by_pos(grid, xpos)
    r = grid.h/2
    rec = UT.HyperRectangle(x .- r, x .+ r)
    return sample_from_rec(rec,N)
end

function sample_elem(Dgrid::DO.DeformedGrid, xpos, N::Int)
    points = sample_elem(Dgrid.grid, xpos, N)
    return [Dgrid.f(x) for x in points]
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


# Monte Carlo sampling     N = 5000
function get_transitions_MC(symmodel, sys, source, symbol, u; N=symmodel.param.N, tstep=sys.tstep)
    (l,xpos) = get_xpos_by_state(symmodel, source)
    grid = get_grid(symmodel.Xdom, l)
    points = sample_elem(grid, xpos, N)
    targetlist = Int[]
    for x in points
        Fx = sys.sys_map(x, u, tstep)
        target = get_state_by_coord(symmodel,Fx)
        push!(targetlist,target)
    end
    occurences = count_occurences(targetlist)
    probaList = [(target,occ/N) for (target,occ) in occurences]
    transitions = Tuple{Int64,Int64,Int64,Float64}[]
    for (target,proba) in probaList
        push!(transitions, (target, source, symbol, proba))
    end
    return transitions
end
