module TestRRT

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

# Deterministic RRT* run that forces exactly one rewire, covering the `RRTstar` branch of
# `UT.RRT`. States are scalars on a line; the distance is the absolute difference.
#
# Iteration 1 adds A at 10.0 as a child of the root (path_cost 10). Iteration 2 adds B at
# 4.0 as a child of the root (path_cost 4). B is then a candidate parent for A: the
# transition B → A is feasible with cost 1, and 1 + 4 < 10, so A is rewired onto B.
@testset "RRT* rewire" begin
    distance(a, b) = abs(a - b)

    seq = [10.0, 4.0]          # the sequence of sampled states
    counter = Ref(0)
    function rand_state(tree, SI, SF, distance, data)
        counter[] += 1
        return seq[min(counter[], length(seq))]
    end

    # Reach the sampled state exactly; the edge cost is the distance travelled.
    new_conf(tree, Nnear, Srand, data) = (Srand, :move, distance(Nnear.state, Srand))
    keep(tree, LSACnew, SI, SF, distance, data) = LSACnew
    stop_crit(tree, LNnew, SI, SF, distance, data) = false
    # Any transition is feasible with a small fixed cost, enabling the cheaper reroute.
    compute_transition(from_state, to_state, data) = (true, :rewire, 1.0)

    tree = UT.RRT(
        0.0,
        100.0,
        distance,
        rand_state,
        new_conf,
        keep,
        stop_crit,
        nothing;
        maxIter = 2,
        RRTstar = true,
        compute_transition = compute_transition,
        k1 = 1,
        k2 = 1,
    )

    nodes = UT.collect_nodes(tree)
    @test UT.get_n_nodes(tree) == 3

    nodeA = only(filter(n -> n.state == 10.0, nodes))
    nodeB = only(filter(n -> n.state == 4.0, nodes))

    # A was rerouted from the root onto B, with the cheaper path cost 1 + 4.
    @test nodeA.parent === nodeB
    @test nodeA.action == :rewire
    @test nodeA.path_cost == 5.0
    @test nodeB.parent === tree.root   # B stays a direct child of the root
end

end # module
