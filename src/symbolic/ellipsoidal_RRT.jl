using LinearAlgebra, IntervalArithmetic, Plots
using ..Utils


function rectangle(c,r)
    Shape(c[1].-r[1] .+ [0,2*r[1],2*r[1],0], c[2].-r[2] .+ [0,0,2*r[2],2*r[2]])
end

function plot_box!(X::IntervalBox;dims=[1,2], color=:red, opacity=1.0)
    xmin = Vector(map(x-> x.lo, X.v))
    xmax = Vector(map(x-> x.hi, X.v))
    c = (xmin + xmax)./2
    h = xmax-xmin
    plot!(rectangle(c[dims],h./2), opacity=opacity, color=color, legend=false)
end

function sample_box(X::IntervalBox)
    return Vector(map(x-> x.lo + (x.hi-x.lo)*rand(), X.v))
end

function sample_ellipsoid(E;N=500)
    box = UT.get_min_bounding_box(E)
    points = [sample_box(box) for i in 1:N]
    filter!(x->x∈E, points)
    return points
end

# data-driven check
function check_controller(E1,kappa,E2,f_eval,Ts;wnew = zeros(n_w))
    samples = sample_ellipsoid(E1;N=500)
    for x in samples
        unew = kappa*[x-E1.c;1]
        xnew = f_eval(x, unew, wnew, Ts)
        if !(xnew ∈ E2)
            return false
        end
    end
    return true
end

# return the nodes containing x in incresing order of costs
function get_nodes_from_x(rrt, x; earlyStop=false)
    stage = rrt.treeLeaves
    nodes = []
    while !isempty(stage)
        for node in stage
            if x ∈ node.state
                push!(nodes,node)
                if earlyStop
                    return nodes
                end
            end
        end
        stage = filter(x -> x!==nothing, unique(map(x-> x.parent, stage)))
    end
    sort!(nodes, by=UT.compare, rev=false)
    return nodes
end

function check_covered(rrt, x)
    if isempty(get_nodes_from_x(rrt, x))
        return false
    else 
        return true
    end
end

function simulate(rrt, f_eval, Ts, x)
    EF = rrt.treeRoot.state
    nodes = get_nodes_from_x(rrt, x)
    currNode = nodes[1]
    if currNode != nothing
        trajx = [x]
        trajE = [currNode.state]
        while !(x ∈ EF)
            kappa = currNode.action
            unew = kappa*[x-currNode.state.c;1]
            wnew = zeros(2)
            x = f_eval(x, unew, wnew, Ts)
            currNode = currNode.parent
            if !(x ∈ currNode.state)
                println("ERROR")
                break
            end
            push!(trajx, x)
            push!(trajE, currNode.state)
        end
        return trajx, trajE
    else
        println("point not covered by the abstraction")
    end
end

function plot_traj!(trajx, trajE)
    for E in trajE
        UT.plotE!(E, color=:blue)
    end
    for i in 1:length(trajx)-1
        UT.plot_arrow!(trajx[i], trajx[i+1])
    end
end


# PSEUDO_CODE
# function buildRRT(problem, qinit, qtarget)
#     rrt = UT.RRT(qinit)
#     while !stop(rrt, problem, qtarget)
#         qrand = get_random_state(rrt, problem)
#         closestNodes, dists = findNClosestNode(rrt, qrand, distance, N=1)  
#         for closeNode in closestNodes
#             qnew, action, cost = get_new_state(rrt, closeNode, qrand, problem)
#             if add_state(qnew, problem)
#                 add_node!(rrt, qnew, closeNode, action, closeNode.path_cost+cost)
#             end
#         end
#     end
#     return rrt
# end


# problem : E0, EF, obstacles, S
# system : f_eval, X, U, Ub, Ts, fT, x, u, w
# affine approximation param : maxRadius, maxΔu, ΔX, ΔU, ΔW
# algorithm param : sdp_opt, distance, get_random_state, get_new_state, keep
function build_lazy_ellipsoidal_abstraction(E0, EF, obstacles, S, f_eval, X, U, Ub, Ts, fT, x, u, w, maxRadius, maxΔu, ΔX, ΔU, ΔW, sdp_opt, distance, get_random_state, get_new_state, keep; maxIter=100)
    println("START")
    rrt = UT.RRT(EF)
    newEllipsoids = [EF]
    bestDist = UT.centerDistance(E0,EF)
    while !any(map(E->(E0 ∈ E),newEllipsoids)) && maxIter>0
        print("Iterations2Go:\t")
        println(maxIter)
        Erand = get_random_state(rrt,X,E0)
        closestNodes, dists = UT.findNClosestNode(rrt, Erand, distance, N=1)  
        newStates = []
        for Eclose in closestNodes
            Enew, kappa, cost = get_new_state(f_eval, Ts, Eclose, Erand, U,S, Ub, maxRadius, maxΔu, sdp_opt, fT, x, u, w)
            if Enew != nothing
                push!(newStates, (Enew, kappa, cost, Eclose))
            end
        end
        newStates = keep(rrt, newStates, E0, obstacles) 
        for data in newStates
            Enew, kappa, cost, Eclose = data
            UT.add_node!(rrt, Enew, Eclose, kappa, Eclose.path_cost + cost)
            bestDist = min(bestDist, UT.centerDistance(E0, Enew))
        end
        newEllipsoids = [s[1] for s in newStates]
        print("\tClosest Dist: ")
        println(bestDist)
        maxIter-=1
    end
    return rrt
end