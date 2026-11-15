module TestMain

using Test, StaticArrays
using Dionysos
const DI = Dionysos
const UT = DI.Utils

sleep(0.1) # used for good printing
println("Started test")

struct State
    x::Int
    y::Int
end

mutable struct PathProblem{T} <: UT.SearchProblem{T}
    initial::Any
    goal::Union{Nothing, T}
    map::Any
    closed::Any # to potentially change
end

function PathProblem(initial, map; goal = nothing)
    return PathProblem{State}(initial, goal, map, nothing)
end

Base.:(==)(s1::State, s2::State) = s1.x == s2.x && s1.y == s2.y

function UT.successor(problem::PathProblem, s::State)
    function is_admissible(map, x, y)
        dim = size(map)
        return x >= 1 && x <= dim[1] && y >= 1 && y <= dim[2] && map[x, y] == 1 ? true :
               false
    end
    function possible_move(map, x, y)
        L = []
        l1 = [-1, 1, 0, 0]
        l2 = [0, 0, -1, 1]
        for i in 1:4
            x_n = x + l1[i]
            y_n = y + l2[i]
            if is_admissible(map, x_n, y_n)
                push!(L, ((l1[i], l2[i]), State(x_n, y_n)))
            end
        end
        return L
    end
    return possible_move(problem.map, s.x, s.y)
end

function print_result(node; all = false)
    println("Number of moves: ", node.depth)
    println("Path cost: ", node.path_cost)
    path = UT.path(node)
    if all
        for n in path
            println(n.state)
        end
    end
end

# heurisic: Manhattan distance
function h1(node::UT.Node, problem::PathProblem)
    state = node.state
    tar = problem.goal
    return abs(state.x - tar.x) + abs(state.y - tar.y)
end

# no heurisic
function h2(node::UT.Node, problem::PathProblem)
    return 0.0
end

Map = transpose(
    SMatrix{7, 7}(
        [
            1 1 1 1 1 1 1
            1 0 0 1 1 1 0
            1 1 0 1 1 1 0
            1 1 0 1 1 1 0
            1 0 0 0 1 1 1
            1 1 1 1 0 0 1
            1 1 1 1 1 0 1
        ],
    ),
)

init = State(5, 7)
target = State(4, 4)
problem = PathProblem(init, Map; goal = target)
@testset "depth_first_graph_search" begin
    node, nb = UT.depth_first_graph_search(problem)
    @test node != nothing
    @test node.state == target
    path = UT.path(node)
    @test path[1].state == init
    @test path[end].state == target
    @test node.depth >= 16
    @test node.path_cost >= 16
end

@testset "breadth_first_graph_search" begin
    node, nb = UT.breadth_first_graph_search(problem)
    @test node != nothing
    @test node.state == target
    path = UT.path(node)
    @test path[1].state == init
    @test path[end].state == target
    @test node.path_cost == 16
    @test node.depth == 16
end
@testset "depth_limited_search" begin
    node = UT.depth_limited_search(problem; limit = 10)
    @test node == "cutoff"
    node = UT.depth_limited_search(problem; limit = 20)
    @test node != nothing
    @test node.state == target
    path = UT.path(node)
    @test path[1].state == init
    @test path[end].state == target
    @test node.path_cost >= 16
    @test node.depth >= 16
end
@testset "iterative_deepening_search" begin
    node = UT.iterative_deepening_search(problem)
    @test node != nothing
    @test node.state == target
    path = UT.path(node)
    @test path[1].state == init
    @test path[end].state == target
    @test node.path_cost == 16
    @test node.depth == 16
end
@testset "astar_graph_search" begin
    node, nb1 = UT.astar_graph_search(problem, h1)
    @test node != nothing
    @test node.state == target
    path = UT.path(node)
    @test path[1].state == init
    @test path[end].state == target
    @test node.path_cost == 16
    @test node.depth == 16

    node, nb2 = UT.astar_graph_search(problem, h2)
    @test node != nothing
    @test node.state == target
    path = UT.path(node)
    @test path[1].state == init
    @test path[end].state == target
    @test node.path_cost == 16
    @test node.depth == 16

    @test nb1 <= nb2
end

end # module TestMain
