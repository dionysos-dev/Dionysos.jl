using LinearAlgebra, IntervalArithmetic, Plots
using ..Utils

UT = Utils

#later we will create a structure for the lazy ellispoidal based abstraction containing 
# -Tree 
# -Graph keeping all transitions computed (even tose delete from rewiring).

" return the nodes containing x in increasing order of costs "
function get_nodes_from_x(tree::UT.Tree, x; earlyStop=false)
    stage = tree.leaves
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

function check_covered(tree::UT.Tree, x)
    if isempty(get_nodes_from_x(tree, x, earlyStop=true))
        return false
    else 
        return true
    end
end

function simulate(tree::UT.Tree, f_eval, Ts, x)
    EF = tree.root.state
    nodes = get_nodes_from_x(tree, x)
    currNode = nodes[1]
    if currNode !== nothing
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

function plot_traj!(trajx, trajE; color=:black)
    for E in trajE
        UT.plotE!(E, color=:blue)
    end
    for i in 1:length(trajx)-1
        UT.plot_arrow!(trajx[i], trajx[i+1], color=color)
    end
end



