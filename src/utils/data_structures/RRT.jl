

mutable struct RRT
    treeRoot::Node
    treeLeaves::Vector{Node}
    nNodes::Int
end

function RRT(state)
    treeRoot = Node(state)
    treeLeaves = [treeRoot] 
    return RRT(treeRoot,treeLeaves,1)
end

function collectNodes(rrt::RRT)
    leaves = copy(rrt.treeLeaves)
    allNodes = copy(leaves)
    while !isempty(leaves)
        parents = filter(x -> x!==nothing, unique(map(x-> x.parent, leaves)))
        append!(allNodes, parents)
        leaves = parents
    end
    return allNodes 
end 

function collectStates(rrt::RRT)
    allNodes = collectNodes(rrt)
    return [node.state for node in allNodes]
end 

function get_nLeaves(rrt::RRT)
    return length(rrt.treeLeaves)
end

function get_nNodes(rrt::RRT)
    return rrt.nNodes
end

function findNmin(tab,N)
    idx = sortperm(tab)
    Nidx = idx[1:min(N,length(tab))]
    return tab[Nidx], Nidx
end

function findNClosestNode(rrt::RRT, state, distance; N=1) 
    allNodes = collectNodes(rrt)
    dists = map(e-> e===nothing ? Inf : distance(e.state, state), allNodes)
    d, idx = findNmin(dists,N)
    return allNodes[idx], d
end


function add_node!(rrt::RRT, state, parent, action, path_cost)
    push!(rrt.treeLeaves, Node(state; parent=parent, action=action, path_cost=path_cost))
    setdiff!(rrt.treeLeaves, [parent])
    rrt.nNodes = rrt.nNodes + 1
end

function add_closest_node!(rrt::RRT, state, distance, get_action)
    closestNode, dists = findNClosestNode(rrt,state, distance) 
    parent = closestNode[1]
    action = get_action(state, parent.state)
    add_node!(rrt, state, parent, action, dists[1] + parent.path_cost)
end

function get_path(rrt::RRT, node::Node)
    return reverse!(path(node))
end

function compare(node)
    return node.path_cost
end

# assuming a positive cost function
function get_max_Node(rrt::RRT)
    sortedLeaves = sort(rrt.treeLeaves, by=compare, rev=true)
    return sortedLeaves[1]
end

function get_max_path_cost(rrt::RRT)
    return get_max_Node(rrt).path_cost
end

function get_min_Node(rrt::RRT)
    return rrt.treeRoot
end

function get_min_path_cost(rrt::RRT)
    return get_min_Node(rrt).path_cost
end

function print_data(rrt::RRT)
    println()
    println("Number of nodes  : ", get_nNodes(rrt))
    println("Number of leaves : ", get_nLeaves(rrt))
    println("Maximum value    : ", get_max_Node(rrt).path_cost)
    println("Minimal value    : ", get_min_Node(rrt).path_cost)
    println()
end


# should just replace plotE! by generic plot! of the state, to get fully generic plot RRT function
function plot_RRT!(rrt::RRT)
    # create a Colormap
    vmin = get_min_path_cost(rrt)
    vmax = get_max_path_cost(rrt)
    if vmin==vmax
        vmax += 1.0e-6
    end
    colorMap = Colormap([vmin,vmax], Colors.colormap("Blues"))
    # plot the nodes of the tree
    allNodes = collectNodes(rrt)
    sort!(allNodes,  by=compare, rev=true)
    for node in allNodes
        plotE!(node.state;color=get_color(colorMap,node.path_cost))
        plot_colorBar!(colorMap)
    end
    # plot edges of the tree
    leaves = copy(rrt.treeLeaves)
    while !isempty(leaves)
        for leave in leaves
            if leave.parent!=nothing
                plot_arrow!(leave.state.c,leave.parent.state.c)
            end
        end
        parents = filter(x -> x!==nothing, unique(map(x-> x.parent, leaves)))
        leaves = parents
    end
end

function plot_path!(rrt, node)
    path = get_path(rrt, node)
    # create a Colormap
    vmin = path[end].path_cost
    vmax = path[1].path_cost
    colorMap = Colormap([vmin,vmax], Colors.colormap("Blues"))
    # plot the nodes of the tree
    sortedPath = sort(path,  by=compare, rev=true)
    for node in sortedPath
        plotE!(node.state;color=get_color(colorMap,node.path_cost))
    end
    plot_colorBar!(colorMap)
    # plot edges of the tree
    for i in 1:length(path)-1
        plot_arrow!(path[i].state.c,path[i+1].state.c)
    end
end