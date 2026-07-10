"""
    NodeT{S, A}

A node of a [`Tree`](@ref): a `state` of type `S`, the `action` of type `A` that reached it
from its parent (`nothing` for the root), and the cost bookkeeping used by RRT/RRT\\*.

`A` defaults to `Any` (the action type is generally only known lazily while the tree grows);
instantiate `Tree{S, A}` with a concrete `A` when the action type is known up front.
"""
mutable struct NodeT{S, A}
    state::S
    parent::Union{Nothing, NodeT{S, A}}
    action::A
    cost::Float64
    path_cost::Float64
    depth::Int
    children::Vector{NodeT{S, A}}
end

function NodeT{S, A}(
    state;
    parent = nothing,
    action = nothing,
    cost = 0.0,
    path_cost = 0.0,
    children = NodeT{S, A}[],
) where {S, A}
    depth = parent !== nothing ? parent.depth + 1 : 0
    return NodeT{S, A}(state, parent, action, cost, path_cost, depth, children)
end

# Convenience constructor: infer the state type, default the action type to `Any`.
NodeT(state::S; kwargs...) where {S} = NodeT{S, Any}(state; kwargs...)

get_state(node::NodeT) = node.state
get_parent(node::NodeT) = node.parent
get_action(node::NodeT) = node.action
get_cost(node::NodeT) = node.cost
get_path_cost(node::NodeT) = node.path_cost

"""
    Tree{S, A}

A rooted tree of [`NodeT`](@ref)`{S, A}` nodes, tracking its leaves and node count. Built
incrementally by [`add_node!`](@ref); the search backbone of the RRT-based abstractions.
"""
mutable struct Tree{S, A}
    root::NodeT{S, A}
    leaves::Vector{NodeT{S, A}}
    n_nodes::Int
end

function Tree(state::S) where {S}
    root = NodeT{S, Any}(state)
    return Tree{S, Any}(root, NodeT{S, Any}[root], 1)
end

function get_n_leaves(tree::Tree)
    return length(tree.leaves)
end

function get_n_nodes(tree::Tree)
    return tree.n_nodes
end

function is_leaf(node::NodeT)
    return length(node.children) == 0
end

function add_child!(tree::Tree, parent::NodeT, child::NodeT)
    if is_leaf(parent)
        setdiff!(tree.leaves, [parent])
    end
    return push!(parent.children, child)
end

function delete_child!(tree::Tree, parent::NodeT, child::NodeT)
    setdiff!(parent.children, [child])
    if is_leaf(parent)
        push!(tree.leaves, parent)
    end
end

"""
    add_node!(tree, state, parent, action, cost; path_cost = parent.path_cost + cost)

Create a node for `state`, attach it to `parent` with the `action` and edge `cost` that
reached it, and register it as a leaf of `tree`. Returns the new node.
"""
function add_node!(
    tree::Tree{S, A},
    state,
    parent,
    action,
    cost;
    path_cost = parent.path_cost + cost,
) where {S, A}
    new_node = NodeT{S, A}(
        state;
        parent = parent,
        action = action,
        cost = cost,
        path_cost = path_cost,
    )
    add_child!(tree, parent, new_node)
    push!(tree.leaves, new_node)
    tree.n_nodes = tree.n_nodes + 1
    return new_node
end

function propagate_cost_to_leaves(node::NodeT)
    for child in node.children
        child.path_cost = node.path_cost + child.cost
        child.depth = node.depth + 1
        propagate_cost_to_leaves(child)
    end
end

# change the parent of the node to new_parent
function rewire(tree::Tree, node::NodeT, new_parent::NodeT, action, cost::Float64)
    delete_child!(tree, node.parent, node)
    add_child!(tree, new_parent, node)
    node.parent = new_parent
    node.action = action
    node.cost = cost
    node.path_cost = cost + new_parent.path_cost
    node.depth = new_parent.depth + 1
    return propagate_cost_to_leaves(node)
end

function collect_children!(node::NodeT, node_accumulator)
    for child in node.children
        push!(node_accumulator, child)
        collect_children!(child, node_accumulator)
    end
end

function collect_children(node::NodeT{S, A}) where {S, A}
    all_nodes = NodeT{S, A}[]
    collect_children!(node, all_nodes)
    return all_nodes
end

# Return a list with node and all its children
function collect_nodes(node::NodeT{S, A}) where {S, A}
    all_nodes = collect_children(node)
    push!(all_nodes, node)
    return all_nodes
end

# Return a list with all the children of node
function collect_nodes(tree::Tree)
    return collect_nodes(tree.root)
end

function get_nodes(tree::Tree{S, A}, state, compare) where {S, A}
    function explore!(node, node_accumulator)
        if node !== nothing && compare(node.state, state)
            push!(node_accumulator, node)
        end
        for child in node.children
            explore!(child, node_accumulator)
        end
    end
    nodes = NodeT{S, A}[]
    explore!(tree.root, nodes)
    return nodes
end

function collect_states(tree::Tree)
    all_nodes = collect_nodes(tree)
    return [node.state for node in all_nodes]
end

function path(node::NodeT{S, A}) where {S, A}
    x, result = node, NodeT{S, A}[node]
    while x.parent !== nothing
        push!(result, x.parent)
        x = x.parent
    end
    return reverse!(result)
end

function get_path(node::NodeT)
    return reverse!(path(node))
end

function compare(node::NodeT)
    return node.path_cost
end

# assuming a positive cost function
function get_max_node(tree::Tree)
    sorted_leaves = sort(tree.leaves; by = compare, rev = true)
    return sorted_leaves[1]
end

function get_max_path_cost(tree::Tree)
    return get_max_node(tree).path_cost
end

function get_min_node(tree::Tree)
    return tree.root
end

function get_min_path_cost(tree::Tree)
    return get_min_node(tree).path_cost
end

####### tree with underlying distance between states #######

function find_k_min(tab, N)
    idx = sortperm(tab)
    Nidx = idx[1:min(N, length(tab))]
    return tab[Nidx], Nidx
end

# when you give a node, you return the k nearest neighbors, except those on the path from the node to the root
function k_nearest_neighbors(tree::Tree, node::NodeT, distance; k = 1)
    all_nodes = collect_nodes(tree)
    node_path = get_path(node)

    pertinent_nodes = filter(e -> !(e ∈ node_path), all_nodes)
    dists = map(e -> e === nothing ? Inf : distance(e.state, node.state), pertinent_nodes)

    d, idx = find_k_min(dists, k)
    return pertinent_nodes[idx], d
end

function k_nearest_neighbors(tree::Tree, state, distance; k = 1)
    all_nodes = collect_nodes(tree)
    dists = map(e -> e === nothing ? Inf : distance(e.state, state), all_nodes)
    d, idx = find_k_min(dists, k)
    return all_nodes[idx], d
end

# add a node in tree whose the parent'state is the closest to state
function add_closest_node!(tree::Tree, state, distance, get_action)
    closest_node, dists = k_nearest_neighbors(tree, state, distance)
    parent = closest_node[1]
    action, cost = get_action(state, parent.state)
    new_node = add_node!(tree, state, parent, action, cost)
    return new_node
end

function Base.show(io::IO, tree::Tree)
    println(io, "Number of nodes  : ", get_n_nodes(tree))
    println(io, "Number of leaves : ", get_n_leaves(tree))
    println(io, "Minimal value    : ", get_min_node(tree).path_cost)
    return println(io, "Maximum value    : ", get_max_node(tree).path_cost)
end

function cost_color(val, vmin, vmax)
    palette = [:lightblue, :deepskyblue, :dodgerblue, :blue, :darkblue]

    if vmax <= vmin + 1e-8
        return palette[cld(length(palette), 2)]
    end

    t = clamp((val - vmin) / (vmax - vmin), 0.0, 1.0)
    idx = clamp(round(Int, 1 + t * (length(palette) - 1)), 1, length(palette))
    return palette[idx]
end

@recipe function f(node::NodeT; pathB = false, cost = true)
    path = get_path(node)

    vmin = path[end].path_cost
    vmax = path[1].path_cost

    pathB ? path = sort(path; by = compare, rev = true) : path = [node]
    sortedPath = sort(path; by = compare, rev = true)

    for node in sortedPath
        @series begin
            color := cost ? cost_color(node.path_cost, vmin, vmax) : :yellow
            node.state
        end
    end

    for i in 1:(length(path) - 1)
        @series begin
            DrawArrow(path[i].state.c, path[i + 1].state.c)
        end
    end
end

@recipe function f(tree::Tree; with_arrows = true, cost = true)
    vmin = get_min_path_cost(tree)
    vmax = get_max_path_cost(tree)

    all_nodes = collect_nodes(tree)
    sort!(all_nodes; by = compare, rev = true)

    for node in all_nodes
        @series begin
            color := cost ? cost_color(node.path_cost, vmin, vmax) : :yellow
            node.state
        end
    end

    if with_arrows
        leaves = copy(tree.leaves)
        while !isempty(leaves)
            for leave in leaves
                if leave.parent !== nothing
                    @series begin
                        DrawArrow(leave.state.c, leave.parent.state.c)
                    end
                end
            end
            parents = filter(x -> x !== nothing, unique(map(x -> x.parent, leaves)))
            leaves = parents
        end
    end
end
