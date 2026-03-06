module TestMain

using Test
using Dionysos
import LinearAlgebra as LA
const DI = Dionysos
const UT = DI.Utils

function distance(E1, E2)
    return UT.pointCenterDistance(E1, E2.c)
end
function get_action(E1, E2)
    return (1.0, 1.0)
end
    @testset "NodeT constructor" begin
        root = UT.NodeT(:root)

        @test UT.get_state(root) == :root
        @test UT.get_parent(root) === nothing
        @test UT.get_action(root) === nothing
        @test UT.get_cost(root) == 0.0
        @test UT.get_path_cost(root) == 0.0
        @test root.depth == 0
        @test isempty(root.children)

        child = UT.NodeT(:child; parent = root, action = :u, cost = 1.5, path_cost = 1.5)

        @test UT.get_state(child) == :child
        @test UT.get_parent(child) === root
        @test UT.get_action(child) == :u
        @test UT.get_cost(child) == 1.5
        @test UT.get_path_cost(child) == 1.5
        @test child.depth == 1
    end

    @testset "Tree constructor" begin
        tree = UT.Tree(:root)

        @test tree.root.state == :root
        @test UT.get_nLeaves(tree) == 1
        @test UT.get_nNodes(tree) == 1
        @test tree.leaves == [tree.root]
        @test UT.is_leave(tree.root)
    end

    @testset "add_node! updates leaves and counts" begin
        tree = UT.Tree(:root)
        root = tree.root

        n1 = UT.add_node!(tree, :n1, root, :a1, 2.0)
        n2 = UT.add_node!(tree, :n2, root, :a2, 3.0)

        @test UT.get_nNodes(tree) == 3
        @test UT.get_nLeaves(tree) == 2
        @test !UT.is_leave(root)
        @test UT.is_leave(n1)
        @test UT.is_leave(n2)
        @test Set(tree.leaves) == Set([n1, n2])

        @test n1.parent === root
        @test n1.cost == 2.0
        @test n1.path_cost == 2.0
        @test n1.depth == 1

        @test n2.path_cost == 3.0
        @test n2.depth == 1
    end

    @testset "delete_child! restores parent as leaf if needed" begin
        tree = UT.Tree(:root)
        root = tree.root
        n1 = UT.add_node!(tree, :n1, root, :a1, 1.0)

        UT.delete_child!(tree, root, n1)

        @test isempty(root.children)
        @test root in tree.leaves
        @test UT.is_leave(root)
    end

    @testset "propagate_cost_to_leaves" begin
        tree = UT.Tree(:root)
        root = tree.root
        n1 = UT.add_node!(tree, :n1, root, :a1, 1.0)
        n2 = UT.add_node!(tree, :n2, n1, :a2, 2.0)
        n3 = UT.add_node!(tree, :n3, n2, :a3, 3.0)

        n1.path_cost = 10.0
        n1.depth = 4
        UT.propagate_cost_to_leaves(n1)

        @test n2.path_cost == 12.0
        @test n2.depth == 5
        @test n3.path_cost == 15.0
        @test n3.depth == 6
    end

    @testset "rewire updates parent cost and depth recursively" begin
        tree = UT.Tree(:root)
        root = tree.root
        a = UT.add_node!(tree, :a, root, :ua, 1.0)
        b = UT.add_node!(tree, :b, root, :ub, 2.0)
        c = UT.add_node!(tree, :c, a, :uc, 3.0)
        d = UT.add_node!(tree, :d, c, :ud, 4.0)

        UT.rewire(tree, c, b, :newu, 5.0)

        @test c.parent === b
        @test c.action == :newu
        @test c.cost == 5.0
        @test c.path_cost == b.path_cost + 5.0
        @test c.depth == b.depth + 1

        @test d.path_cost == c.path_cost + d.cost
        @test d.depth == c.depth + 1
        @test c in b.children
        @test !(c in a.children)
    end

    @testset "collect helpers" begin
        tree = UT.Tree(:root)
        root = tree.root
        a = UT.add_node!(tree, :a, root, :ua, 1.0)
        b = UT.add_node!(tree, :b, root, :ub, 2.0)
        c = UT.add_node!(tree, :c, a, :uc, 3.0)

        children_root = UT.collect_children(root)
        @test Set(getfield.(children_root, :state)) == Set([:a, :b, :c])

        nodes_a = UT.collect_nodes(a)
        @test Set(getfield.(nodes_a, :state)) == Set([:a, :c])

        nodes_tree = UT.collect_nodes(tree)
        @test Set(getfield.(nodes_tree, :state)) == Set([:root, :a, :b, :c])

        states = UT.collect_states(tree)
        @test Set(states) == Set([:root, :a, :b, :c])
    end

    @testset "get_nodes with comparison" begin
        tree = UT.Tree((0, 0))
        root = tree.root
        UT.add_node!(tree, (1, 0), root, :a1, 1.0)
        UT.add_node!(tree, (1, 1), root, :a2, 2.0)
        UT.add_node!(tree, (2, 0), root, :a3, 3.0)

        nodes = UT.get_nodes(tree, 1, (x, y) -> x[1] == y)

        @test length(nodes) == 2
        @test Set(getfield.(nodes, :state)) == Set([(1, 0), (1, 1)])
    end

    @testset "path and get_path" begin
        tree = UT.Tree(:root)
        root = tree.root
        a = UT.add_node!(tree, :a, root, :ua, 1.0)
        b = UT.add_node!(tree, :b, a, :ub, 2.0)

        p1 = UT.path(b)
        p2 = UT.get_path(b)

        @test getfield.(p1, :state) == [:root, :a, :b]
        @test getfield.(p2, :state) == [:b, :a, :root]
    end

    @testset "min/max path cost" begin
        tree = UT.Tree(:root)
        root = tree.root
        a = UT.add_node!(tree, :a, root, :ua, 1.0)
        b = UT.add_node!(tree, :b, root, :ub, 5.0)
        c = UT.add_node!(tree, :c, a, :uc, 2.0)

        @test UT.get_min_Node(tree) === root
        @test UT.get_min_path_cost(tree) == 0.0
        @test UT.get_max_Node(tree) === c
        @test UT.get_max_path_cost(tree) == 3.0
    end

    @testset "findkmin" begin
        vals = [4.0, 1.0, 3.0, 2.0]
        d, idx = UT.findkmin(vals, 2)

        @test d == [1.0, 2.0]
        @test idx == [2, 4]
    end

    @testset "kNearestNeighbors with state" begin
        tree = UT.Tree(0.0)
        root = tree.root
        UT.add_node!(tree, 1.0, root, :a1, 1.0)
        UT.add_node!(tree, 3.0, root, :a2, 1.0)
        UT.add_node!(tree, 5.0, root, :a3, 1.0)

        nodes, dists = UT.kNearestNeighbors(tree, 2.2, (x, y) -> abs(x - y); k = 2)

        @test length(nodes) == 2
        @test dists == [0.8, 1.2]
        @test getfield.(nodes, :state) == [3.0, 1.0]
    end

    @testset "kNearestNeighbors with node excludes ancestors" begin
        tree = UT.Tree(0.0)
        root = tree.root
        a = UT.add_node!(tree, 1.0, root, :a, 1.0)
        b = UT.add_node!(tree, 2.0, a, :b, 1.0)
        c = UT.add_node!(tree, 10.0, root, :c, 1.0)

        nodes, dists = UT.kNearestNeighbors(tree, b, (x, y) -> abs(x - y); k = 2)

        @test !(root in nodes)
        @test !(a in nodes)
        @test !(b in nodes)
        @test c in nodes
        @test dists[1] == 8.0
    end

    @testset "add_closest_node!" begin
        tree = UT.Tree(0.0)
        root = tree.root
        a = UT.add_node!(tree, 2.0, root, :a, 2.0)

        newnode = UT.add_closest_node!(
            tree,
            2.5,
            (x, y) -> abs(x - y),
            (x, y) -> (:move, abs(x - y)),
        )

        @test newnode.parent === a
        @test newnode.action == :move
        @test newnode.cost == 0.5
        @test newnode.path_cost == a.path_cost + 0.5
        @test UT.get_nNodes(tree) == 3
    end

    @testset "show(Tree)" begin
        tree = UT.Tree(:root)
        io = IOBuffer()
        show(io, tree)
        str = String(take!(io))

        @test occursin("Number of nodes", str)
        @test occursin("Number of leaves", str)
        @test occursin("Minimal value", str)
        @test occursin("Maximum value", str)
    end
end
end
