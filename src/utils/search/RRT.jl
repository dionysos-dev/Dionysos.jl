"
    Basic generic RRT algorihtm
    SI : initial state (this will be the root of the tree) ;
    SF : target state that we try to reach ;
    distance : function that defines a metric between the states ;
    rand_state : function that returns a random candidate state ;
    new_conf : function that returns a reachable state, with the action and the cost ;
    keep : function to filter the node that we want to add during an iteration;
    stop_crit : stop criteria ;
    RRTstar : boolean to use RRT* ;
    compute_transition : to compute transition between to given state (if RRTstar is true).
"
function RRT(
    SI,
    SF,
    distance,
    rand_state,
    new_conf,
    keep,
    stop_crit,
    data;
    maxIter = 100,
    RRTstar = false,
    compute_transition,
    k1 = 1,
    k2 = 1,
)
    tree = Tree(SI)
    bestDist = distance(SI, SF)
    LNnew = [tree.root]
    while !stop_crit(tree, LNnew, SI, SF, distance, data) && maxIter > 0
        print("Iterations2Go:\t")
        println(maxIter)
        Srand = rand_state(tree, SI, SF, distance, data)
        LNnear, ~ = kNearestNeighbors(tree, Srand, distance; k = k1)
        # generate state we want to add in one iteration
        LSACnew = []
        for Nnear in LNnear
            Snew, action, cost = new_conf(tree, Nnear, Srand, data)
            Snew === nothing ? nothing : push!(LSACnew, (Snew, action, cost, Nnear))
        end
        LSACnew = keep(tree, LSACnew, SI, SF, distance, data)
        # add the new nodes
        LNnew = []
        for (Snew, action, cost, Nnear) in LSACnew
            push!(LNnew, add_node!(tree, Snew, Nnear, action, cost))
            bestDist = min(bestDist, distance(SF, Snew))
        end
        # RRT*
        if RRTstar
            for Nnew in LNnew
                LNclose, ~ = kNearestNeighbors(tree, Nnew, distance; k = k2)
                for Nclose in LNclose
                    ans, action, cost = compute_transition(Nclose.state, Nnew.state, data)
                    if ans
                        if cost + Nnew.path_cost < Nclose.path_cost
                            println("REWIRE")
                            rewire(tree, Nclose, Nnew, action, cost)
                        end
                    end
                end
            end
        end
        print("\tClosest Dist: ")
        println(bestDist)
        maxIter -= 1
    end
    return tree
end
