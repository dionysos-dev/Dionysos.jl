using LinearAlgebra
using ..Utils

"
    Basic generic RRT algorihtm
    SI : initial state (this will be the root of the tree) ;
    SF : target state that we try to reach ;
    distance : function that defines a metric between the states ;
    rand_state : function that returns a random candidate state ;
    new_conf : function that returns a reachable state, with the action and the cost ;
    keep : function to filter the node that we want to add during an iteration;
    stop_crit : stop criteria
"
function RRT(SI, SF, distance, rand_state, new_conf, keep, stop_crit, data; maxIter=100, RRTstar=false, compute_transition)
    tree = UT.Tree(SI)
    bestDist = distance(SI, SF)
    LNnew = [tree.root]
    while !stop_crit(tree, LNnew, SI, SF, distance, data) && maxIter>0
        Srand = rand_state(tree, SI, SF, distance)
        LNnear = UT.kNearestNeighbors(tree, Srand, distance)
        # generate state we want to add in one iteration
        LSACnew = [] 
        for Nnear in LNnear
            Snew, action, cost = new_conf(tree, Nnear, Srand, data)
            Snew===nothing ?  nothing : push!(LSACnew, (Snew, action, cost, Nnear))
        end
        LSACnew = keep(LSACnew, tree, SI, SF, distance, data)
        # add the new nodes
        LNnew = []
        for (Snew, action, cost, Nnear) in LSACnew
            push!(LNnew, UT.add_node!(tree, Snew, Nnear, action, cost))
            bestDist = min(bestDist, distance(SF, Snew))
        end
        # RRT*
        if RRTstar
            for Nnew in LNnew
                LNclose, ~ = UT.kNearestNeighbors(tree, Nnew, distance, k=1)
                for Nclose in LNclose
                    ans, action, cost = compute_transition(Nclose.state, Nnew.state, data)
                    if ans
                        if cost + Nnew.path_cost  < Nclose.path_cost
                            UT.rewire(tree, Nclose, Nnew, action, cost)
                        end
                    end
                end
            end
        end
        print("\tClosest Dist: ")
        println(bestDist)
        maxIter-=1
    end
    return tree
end
