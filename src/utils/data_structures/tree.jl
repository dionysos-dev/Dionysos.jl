"
cost: cost to reach its parent
"
mutable struct NodeT{S}
    state::S
    parent::Union{Nothing, NodeT{S}}
    action::Any
    cost::Float64
    path_cost::Float64
    depth::Int
    children::Vector{NodeT{S}}
end

function NodeT(
    state::S;
    parent = nothing,
    action = nothing,
    cost = 0.0,
    path_cost = 0.0,
    children = NodeT{S}[],
) where {S}
    depth = parent !== nothing ? parent.depth + 1 : 0
    return NodeT(state, parent, action, cost, path_cost, depth, children)
end

get_state(node::NodeT) = node.state
get_parent(node::NodeT) = node.parent
get_action(node::NodeT) = node.action
get_cost(node::NodeT) = node.cost
get_path_cost(node::NodeT) = node.path_cost

"
Tree structure with
- cost for transitions (the cost function a non-negative function);
- an underlying metric between the states that are encapsulated in the nodes of the tree
"
mutable struct Tree
    root::NodeT
    leaves::Vector{NodeT}
    nNodes::Int
end

function Tree(state)
    root = NodeT(state)
    leaves = [root]
    return Tree(root, leaves, 1)
end

function get_nLeaves(tree::Tree)
    return length(tree.leaves)
end

function get_nNodes(tree::Tree)
    return tree.nNodes
end

function is_leave(node::NodeT)
    return length(node.children) == 0
end

function add_child!(tree::Tree, parent::NodeT, child::NodeT)
    if is_leave(parent)
        setdiff!(tree.leaves, [parent])
    end
    return push!(parent.children, child)
end

function delete_child!(tree::Tree, parent::NodeT, child::NodeT)
    setdiff!(parent.children, [child])
    if is_leave(parent)
        push!(tree.leaves, parent)
    end
end

"add a node as a leave"
function add_node!(
    tree::Tree,
    state,
    parent,
    action,
    cost;
    path_cost = parent.path_cost + cost,
)
    newNode =
        NodeT(state; parent = parent, action = action, cost = cost, path_cost = path_cost)
    add_child!(tree, parent, newNode)
    push!(tree.leaves, newNode)
    tree.nNodes = tree.nNodes + 1
    return newNode
end

"assuming that the path_cost of a node has changed (and its depth), 
we should propagate the new cost to its children"
function propagate_cost_to_leaves(node::NodeT)
    for child in node.children
        child.path_cost = node.path_cost + child.cost
        child.depth = node.depth + 1
        propagate_cost_to_leaves(child)
    end
end

# change the parent of the node to newParent
function rewire(tree::Tree, node::NodeT, newParent::NodeT, action, cost::Float64)
    delete_child!(tree, node.parent, node)
    add_child!(tree, newParent, node)
    node.parent = newParent
    node.action = action
    node.cost = cost
    node.path_cost = cost + newParent.path_cost
    node.depth = newParent.depth + 1
    return propagate_cost_to_leaves(node)
end

function collect_children!(node::NodeT, nodeAccumulator)
    for child in node.children
        push!(nodeAccumulator, child)
        collect_children!(child, nodeAccumulator)
    end
end

"Return a list with all the children of node"
function collect_children(node::NodeT)
    allNodes = []
    collect_children!(node, allNodes)
    return allNodes
end

# Return a list with node and all its children
function collect_nodes(node::NodeT)
    allNodes = collect_children(node)
    push!(allNodes, node)
    return allNodes
end

# Return a list with all the children of node
function collect_nodes(tree::Tree)
    return collect_nodes(tree.root)
end

function get_nodes(tree::Tree, state, compare)
    function explore!(node, nodeAccumulator)
        if node !== nothing && compare(node.state, state)
            push!(nodeAccumulator, node)
        end
        for child in node.children
            explore!(child, nodeAccumulator)
        end
    end
    nodes = []
    explore!(tree.root, nodes)
    return nodes
end

function collect_states(tree::Tree)
    allNodes = collect_nodes(tree)
    return [node.state for node in allNodes]
end

"Create a list of nodes from the root to this node."
function path(node::NodeT)
    x, result = node, [node]
    while x.parent !== nothing
        push!(result, x.parent)
        x = x.parent
    end
    return reverse!(result)
end

"return the path from node to the root of the tree"
function get_path(node::NodeT)
    return reverse!(path(node))
end

function compare(node::NodeT)
    return node.path_cost
end

# assuming a positive cost function
function get_max_Node(tree::Tree)
    sortedLeaves = sort(tree.leaves; by = compare, rev = true)
    return sortedLeaves[1]
end

function get_max_path_cost(tree::Tree)
    return get_max_Node(tree).path_cost
end

function get_min_Node(tree::Tree)
    return tree.root
end

function get_min_path_cost(tree::Tree)
    return get_min_Node(tree).path_cost
end

####### tree with underlying distance between states #######

function findkmin(tab, N)
    idx = sortperm(tab)
    Nidx = idx[1:min(N, length(tab))]
    return tab[Nidx], Nidx
end

# when you give a node, you return the k nearest neighbors, except those on the path from the node to the root
function kNearestNeighbors(tree::Tree, node::NodeT, distance; k = 1)
    allNodes = collect_nodes(tree)
    path = get_path(node)

    pertinentNodes = filter(e -> !(e ∈ path), allNodes)
    dists = map(e -> e === nothing ? Inf : distance(e.state, node.state), pertinentNodes)

    d, idx = findkmin(dists, k)
    return pertinentNodes[idx], d
end

function kNearestNeighbors(tree::Tree, state, distance; k = 1)
    allNodes = collect_nodes(tree)
    dists = map(e -> e === nothing ? Inf : distance(e.state, state), allNodes)
    d, idx = findkmin(dists, k)
    return allNodes[idx], d
end

# add a node in tree whose the parent'state is the closest to state
function add_closest_node!(tree::Tree, state, distance, get_action)
    closestNode, dists = kNearestNeighbors(tree, state, distance)
    parent = closestNode[1]
    action, cost = get_action(state, parent.state)
    newNode = add_node!(tree, state, parent, action, cost)
    return newNode
end

function Base.show(io::IO, tree::Tree)
    println(io, "Number of nodes  : ", get_nNodes(tree))
    println(io, "Number of leaves : ", get_nLeaves(tree))
    println(io, "Minimal value    : ", get_min_Node(tree).path_cost)
    return println(io, "Maximum value    : ", get_max_Node(tree).path_cost)
end

@recipe function f(node::NodeT; pathB = false, cost = true)
    path = get_path(node)
    # create a Colormap
    vmin = path[end].path_cost
    vmax = path[1].path_cost
    colorMap = Colormap([vmin, vmax], Colors.colormap("Blues"))
    pathB ? path = sort(path; by = compare, rev = true) : path = [node]
    sortedPath = sort(path; by = compare, rev = true)
    for node in sortedPath
        @series begin
            cost ? color := get_color(colorMap, node.path_cost) : color := :yellow
            return node.state
        end
    end
    if cost
        @series begin
            colorMap
        end
    end
    for i in 1:(length(path) - 1)
        @series begin
            DrawArrow(path[i].state.c, path[i + 1].state.c)
        end
    end
end

@recipe function f(tree::Tree; arrowsB = true, cost = true)
    # create a Colormap
    vmin = get_min_path_cost(tree)
    vmax = get_max_path_cost(tree)
    colorMap = Colormap([vmin, vmax], Colors.colormap("Blues"))
    # plot the nodes of the tree
    allNodes = collect_nodes(tree)
    sort!(allNodes; by = compare, rev = true)
    for node in allNodes
        @series begin
            cost ? color := get_color(colorMap, node.path_cost) : color := :yellow
            return node.state
        end
    end
    if cost
        @series begin
            colorMap
        end
    end
    if arrowsB
        leaves = copy(tree.leaves)
        while !isempty(leaves)
            for leave in leaves
                if leave.parent !== nothing
                    @series begin
                        return DrawArrow(leave.state.c, leave.parent.state.c)
                    end
                end
            end
            parents = filter(x -> x !== nothing, unique(map(x -> x.parent, leaves)))
            leaves = parents
        end
    end
end
