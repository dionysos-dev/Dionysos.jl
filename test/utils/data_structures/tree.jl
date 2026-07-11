module TestMain

using Test
using Dionysos
import LinearAlgebra as LA
import LazySets

const DI = Dionysos
const UT = DI.Utils

distance(E1, E2) = UT.center_distance(E1, LazySets.center(E2))
get_action(E1, E2) = (1.0, 1.0)

make_ellipsoid(scale, center) =
    LazySets.Ellipsoid(center, Matrix{Float64}(LA.I, 2, 2) / scale)

@testset "Tree" begin
    ellipsoids = [
        make_ellipsoid(8.0, [-10.0; -10.0]),
        make_ellipsoid(5.0, [0.0; -10.0]),
        make_ellipsoid(1.0, [-10.0; 0.0]),
        make_ellipsoid(3.0, [20.0; -10.0]),
        make_ellipsoid(3.0, [-1.0; 0.0]),
        make_ellipsoid(3.0, [1.0; -8.0]),
        make_ellipsoid(3.0, [-1.0; 5.0]),
        make_ellipsoid(3.0, [3.0; 0.0]),
    ]

    tree = UT.Tree(ellipsoids[1])

    @test UT.get_n_leaves(tree) == 1
    @test UT.get_n_nodes(tree) == 1
    @test UT.is_leaf(tree.root)
    @test tree.root.depth == 0
    @test tree.root.path_cost == 0.0

    nearest_nodes, dists = UT.k_nearest_neighbors(tree, ellipsoids[2], distance; k = 1)
    @test length(nearest_nodes) == 1
    @test dists == [10.0]
    @test nearest_nodes[1] === tree.root

    action, cost = get_action(ellipsoids[2], nearest_nodes[1].state)
    node2 = UT.add_node!(tree, ellipsoids[2], nearest_nodes[1], action, cost)

    @test !UT.is_leaf(tree.root)
    @test UT.is_leaf(node2)
    @test node2.parent === tree.root
    @test node2.depth == 1
    @test node2.cost == 1.0
    @test node2.path_cost == 1.0

    @test UT.get_n_leaves(tree) == 1
    @test UT.get_n_nodes(tree) == 2

    nearest_nodes, dists = UT.k_nearest_neighbors(tree, ellipsoids[2], distance; k = 2)
    @test length(nearest_nodes) == 2
    @test dists[1] == 0.0
    @test dists[2] == 10.0

    node3 = UT.add_closest_node!(tree, ellipsoids[3], distance, get_action)
    node4 = UT.add_closest_node!(tree, ellipsoids[4], distance, get_action)
    node5 = UT.add_closest_node!(tree, ellipsoids[5], distance, get_action)
    node6 = UT.add_closest_node!(tree, ellipsoids[6], distance, get_action)
    node7 = UT.add_closest_node!(tree, ellipsoids[7], distance, get_action)
    node8 = UT.add_closest_node!(tree, ellipsoids[8], distance, get_action)

    @test UT.get_n_nodes(tree) == 8

    @test node3.path_cost == 1.0
    @test node5.path_cost == 2.0
    @test node7.path_cost == 3.0
    @test node8.path_cost == 3.0

    node3.path_cost = 3.0
    UT.propagate_cost_to_leaves(node3)

    @test node3.path_cost == 3.0
    @test node5.path_cost == 4.0
    @test node7.path_cost == 5.0
    @test node8.path_cost == 5.0

    old_parent = node5.parent
    @test old_parent === node3

    UT.rewire(tree, node5, node6, 1.0, 1.0)

    @test node5.parent === node6
    @test node5.action == 1.0
    @test node5.cost == 1.0
    @test node5.path_cost == 3.0
    @test node5.depth == node6.depth + 1

    @test node7.path_cost == 4.0
    @test node8.path_cost == 4.0
    @test node7.depth == node5.depth + 1
    @test node8.depth == node5.depth + 1

    @test UT.is_leaf(node3)
    @test !UT.is_leaf(node6)
end

@testset "Tree utility coverage" begin
    tree = UT.Tree(:root)
    root = tree.root

    a = UT.add_node!(tree, :a, root, :ua, 1.0)
    b = UT.add_node!(tree, :b, root, :ub, 2.0)
    c = UT.add_node!(tree, :c, a, :uc, 3.0)

    @test UT.get_state(c) == :c
    @test UT.get_parent(c) === a
    @test UT.get_action(c) == :uc
    @test UT.get_cost(c) == 3.0
    @test UT.get_path_cost(c) == 4.0

    @test getfield.(UT.path(c), :state) == [:root, :a, :c]
    @test getfield.(UT.get_path(c), :state) == [:c, :a, :root]

    @test Set(getfield.(UT.collect_children(root), :state)) == Set([:a, :b, :c])
    @test Set(getfield.(UT.collect_nodes(a), :state)) == Set([:a, :c])
    @test Set(getfield.(UT.collect_nodes(tree), :state)) == Set([:root, :a, :b, :c])

    @test Set(UT.collect_states(tree)) == Set([:root, :a, :b, :c])

    nodes = UT.get_nodes(tree, :a, ==)
    @test length(nodes) == 1
    @test nodes[1] === a

    @test UT.compare(c) == 4.0
    @test UT.get_min_node(tree) === root
    @test UT.get_min_path_cost(tree) == 0.0
    @test UT.get_max_node(tree) === c
    @test UT.get_max_path_cost(tree) == 4.0

    vals, idx = UT.find_k_min([3.0, 1.0, 2.0], 2)
    @test vals == [1.0, 2.0]
    @test idx == [2, 3]

    UT.delete_child!(tree, a, c)
    @test UT.is_leaf(a)
    @test !(c in a.children)

    io = IOBuffer()
    show(io, tree)
    s = String(take!(io))
    @test occursin("Number of nodes", s)
    @test occursin("Number of leaves", s)
end

@testset "kNearestNeighbors excludes path to root" begin
    tree = UT.Tree(0.0)
    root = tree.root
    a = UT.add_node!(tree, 1.0, root, :a, 1.0)
    b = UT.add_node!(tree, 2.0, a, :b, 1.0)
    c = UT.add_node!(tree, 10.0, root, :c, 1.0)

    nodes, _ = UT.k_nearest_neighbors(tree, b, (x, y) -> abs(x - y); k = 2)

    @test !(root in nodes)
    @test !(a in nodes)
    @test !(b in nodes)
    @test c in nodes
end

end # module TestMain
