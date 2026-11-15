module TestMain

using Test
using Dionysos
const DI = Dionysos
const UT = DI.Utils

sleep(0.1) # used for good printing
println("Started test")

@testset "Digraph1" begin
    testgraph = [("a", "b", 1), ("b", "e", 2), ("a", "e", 4)]
    g = UT.Digraph(testgraph)

    src, dst = "a", "e"
    path, cost = UT.dijkstrapath(g, src, dst)
    @test length(path) == 3
    @test path[1] == "a"
    @test path[2] == "b"
    @test path[3] == "e"
    @test cost[path[1]] == 0.0
    @test cost[path[2]] == 1.0
    @test cost[path[3]] == 3.0

    src, dst = "e", "a"
    path, cost = UT.dijkstrapath(g, src, dst)
    @test isempty(path) == true
end

@testset "Digraph2" begin
    testgraph = [("a", "b", 1), ("b", "e", 5), ("a", "e", 4)]
    g = UT.Digraph(testgraph)

    src, dst = "a", "e"
    path, cost = UT.dijkstrapath(g, src, dst)
    @test length(path) == 2
    @test path[1] == "a"
    @test path[2] == "e"
    @test cost[path[1]] == 0.0
    @test cost[path[2]] == 4.0
end

@testset "Digraph3" begin
    testgraph = [
        ("a", "b", 7.1),
        ("a", "c", 9.5),
        ("a", "f", 14.5),
        ("b", "c", 10.5),
        ("b", "d", 15.5),
        ("c", "d", 20.5),
        ("c", "f", 2.5),
        ("d", "e", 26.5),
        ("e", "f", 9.5),
        ("e", "g", 9.5),
    ]
    g = UT.Digraph(testgraph)

    src, dst = "a", "g"
    path, cost = UT.dijkstrapath(g, src, dst)
    @test length(path) == 5
    @test path[1] == "a"
    @test path[2] == "b"
    @test path[3] == "d"
    @test path[4] == "e"
    @test path[5] == "g"
    @test cost[path[1]] == 0.0
    @test cost[path[2]] == 7.1
    @test cost[path[3]] == 22.6
    @test cost[path[4]] == 49.1
    @test cost[path[5]] == 58.6

    src, dst = "f", "g"
    path, cost = UT.dijkstrapath(g, src, dst)
    @test isempty(path) == true

    src, dst = "f", "c"
    path, cost = UT.dijkstrapath(g, src, dst)
    @test isempty(path) == true

    src, dst = "f", "e"
    path, cost = UT.dijkstrapath(g, src, dst)
    @test isempty(path) == true

    src, dst = "f", "b"
    path, cost = UT.dijkstrapath(g, src, dst)
    @test isempty(path) == true
end

end  # module TestMain
