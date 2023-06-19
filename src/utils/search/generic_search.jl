"""
    SearchProblem
### Fields
    - `initial`  -- initial state of type `S` or a list of initial state of type `S`.
    - `goal`     -- possibly a goal state of type 'Union{Nothing,S}''.
### Example
    struct Problem{S} <: SearchProblem{S}
        initial::Union{S,Vector{S}}
        goal::Union{Nothing,S}
    end
    The constructor specifies the initial state, and possibly a goal
    state, if there is a unique goal.
    function Problem(initial::S; goal=nothing) where S
        return Problem{S}(initial,goal)
    end
"""
abstract type SearchProblem{S} end

"""
Given a state, return a sequence of (action, state) pairs reachable
from this state. If there are many successors, consider an iterator
that yields the successors one at a time, rather than building them
all at once.
"""
function successor(problem::SearchProblem{S}, state::S) where {S} end

"""
Return True if the state is a goal. The default method compares the
state to P.goal, as specified in the constructor. Implement this
method if checking against a single goal is not enough.
"""
function goal_test(problem::SearchProblem{S}, state::S) where {S}
    return state == problem.goal
end

"""
Return the cost of a solution path that arrives at state2 from
state1 via action, assuming cost c to get up to state1. If the problem
is such that the path doesn't matter, this function will only look at
state2.  If the path does matter, it will consider c and maybe state1
and action. The default method costs 1 for every step in the path.
"""
function path_cost(problem::SearchProblem{S}, c, state1::S, action, state2::S) where {S}
    return c + 1
end

##
"""
A node in a search tree. Contains a pointer to the parent (the node
that this is a successor of) and to the actual state for this node. Note
that if a state is arrived at by two paths, then there are two nodes with
the same state. Also includes the action that got us to this state, and
the total path_cost (also known as g) to reach the node. Other functions
may add an f and h value; see best_first_graph_search and astar_search for
an explanation of how the f and h values are handled.
"""
mutable struct Node{S}
    state::S
    parent::Union{Nothing, Node{S}}
    action::Any
    path_cost::Any
    depth::Int
end

"Create a search tree Node, derived from a parent by an action."
function Node(state; parent = nothing, action = nothing, path_cost = 0.0)
    depth = parent !== nothing ? parent.depth + 1 : 0
    return Node(state, parent, action, path_cost, depth)
end

"Create a list of nodes from the root to this node."
function path(node::Node)
    x, result = node, [node]
    while x.parent !== nothing
        push!(result, x.parent)
        x = x.parent
    end
    return reverse!(result)
end

"Yield the nodes reachable from this node."
function expand(node::Node, problem::SearchProblem)
    return [
        #yield() #not sure if usefull
        Node(
            next;
            parent = node,
            action = act,
            path_cost = path_cost(problem, node.path_cost, node.state, act, next),
        ) for (act, next) in successor(problem, node.state)
    ]
end

## Uninformed Search algorithms
"""
Search through the successors of a problem to find a goal.
The argument fringe should be an empty queue.
Don't worry about repeated paths to a state.
"""
function tree_search(problem::SearchProblem, fringe)
    append!(fringe, Node(problem.initial))
    n = 0
    while !isempty(fringe)
        node = pop!(fringe)
        n += 1
        if goal_test(problem, node.state)
            return node, n
        end
        extend!(fringe, expand(node, problem))
    end
    return nothing, n
end

"""
Search the shallowest nodes in the search tree first.
"""
function breadth_first_tree_search(problem::SearchProblem)
    return tree_search(problem, FIFOQueue{Node}())
end

"""
Search the deepest nodes in the search tree first.
"""
function depth_first_tree_search(problem::SearchProblem)
    return tree_search(problem, MyStack{Node}())
end

"""
Search through the successors of a problem to find a goal.
The argument fringe should be an empty queue.
If two paths reach a state, only use the best one.
"""
function graph_search(problem::SearchProblem{S}, fringe) where {S}
    problem.closed = Dict{S, Bool}()
    n = 0
    if typeof(problem.initial) == Vector{S}
        for init in problem.initial
            append!(fringe, Node(init))
        end
    else
        append!(fringe, Node(problem.initial))
    end
    while !isempty(fringe)
        node = pop!(fringe)
        n += 1
        if goal_test(problem, node.state)
            return node, n
        end
        if !haskey(problem.closed, node.state)
            problem.closed[node.state] = true
            extend!(fringe, expand(node, problem))
        end
    end
    return nothing, n
end

"""
Search the shallowest nodes in the search tree first.
"""
function breadth_first_graph_search(problem)
    return graph_search(problem, FIFOQueue{Node}())
end

"""
Search the deepest nodes in the search tree first.
"""
function depth_first_graph_search(problem::SearchProblem)
    return graph_search(problem, MyStack{Node}())
end

function depth_limited_search(problem; limit = 50)
    function recursive_dls(node, problem, limit)
        cutoff_occurred = false
        if goal_test(problem, node.state)
            return node
        elseif node.depth == limit
            return "cutoff"
        else
            for successor in expand(node, problem)
                result = recursive_dls(successor, problem, limit)
                if result == "cutoff"
                    cutoff_occurred = true
                elseif result != nothing
                    return result
                end
            end
        end
        if cutoff_occurred
            return "cutoff"
        else
            return nothing
        end
    end
    return recursive_dls(Node(problem.initial), problem, limit)
end

function iterative_deepening_search(problem; max_depth = 100)
    for depth in 1:max_depth
        result = depth_limited_search(problem; limit = depth)
        if result != "cutoff"
            return result
        end
    end
end

## Informed (Heuristic) Search

"""
Search the nodes with the lowest f scores first.
You specify the function f(node) that you want to minimize; for example,
if f is a heuristic estimate to the goal, then we have greedy best
first search; if f is node.depth then we have depth-first search.
"""
function best_first_graph_search(problem::SearchProblem, f)
    return graph_search(problem, MyPriorityQueue{Node, Float64}(f, problem))
end

"""
A* search is best-first graph search with f(n) = g(n)+h(n).
You need to specify the h function when you call astar_search.
"""
function astar_graph_search(problem::SearchProblem, h)
    f(n, problem) = n.path_cost + h(n, problem)
    return best_first_graph_search(problem, f)
end

"""
Search the nodes with the lowest f scores first.
You specify the function f(node) that you want to minimize; for example,
if f is a heuristic estimate to the goal, then we have greedy best
first search; if f is node.depth then we have depth-first search.
"""
function best_first_tree_search(problem::SearchProblem, f)
    return tree_search(problem, MyPriorityQueue{Node, Float64}(f, problem))
end

"""
A* search is best-first graph search with f(n) = g(n)+h(n).
You need to specify the h function when you call astar_search.
"""
function astar_tree_search(problem::SearchProblem, h)
    f(n, problem) = n.path_cost + h(n, problem)
    return best_first_tree_search(problem, f)
end
