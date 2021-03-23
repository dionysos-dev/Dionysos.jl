include("../queue.jl")

module TestMain
using Test
using Main.MyQueue
const MQ = Main.MyQueue

sleep(0.1) # used for good printing
println("Started test")

@testset "MyStack" begin
    S = MQ.MyStack{Int64}()
    @test MQ.isempty(S)
    MQ.append!(S,1)
    MQ.append!(S,2)
    MQ.append!(S,3)
    MQ.append!(S,4)
    @test MQ.length(S) == 4
    MQ.append!(S,5)
    @test MQ.length(S) == 5

    @test MQ.pop!(S) == 5
    @test MQ.pop!(S) == 4
    @test MQ.length(S) == 3
    @test MQ.isempty(S) == false

    S2 = MQ.MyStack{Int64}()
    MQ.append!(S2,1)
    MQ.append!(S2,2)
    MQ.append!(S2,3)
    @test S==S2 ? true : false
    MQ.pop!(S2)
    @test S==S2 ? false : true

    empty!(S)
    @test MQ.isempty(S)
end

@testset "FIFOQueue" begin
    S = MQ.FIFOQueue{Int64}()
    @test MQ.isempty(S)
    MQ.append!(S,1)
    MQ.append!(S,2)
    MQ.append!(S,3)
    MQ.append!(S,4)
    @test MQ.length(S) == 4
    MQ.append!(S,5)
    @test MQ.length(S) == 5

    @test MQ.pop!(S) == 1
    @test MQ.pop!(S) == 2
    @test MQ.length(S) == 3
    @test MQ.isempty(S) == false

    S2 = MQ.FIFOQueue{Int64}()
    MQ.append!(S2,3)
    MQ.append!(S2,4)
    MQ.append!(S2,5)
    @test S==S2 ? true : false
    MQ.pop!(S2)
    @test S==S2 ? false : true

    empty!(S)
    @test MQ.isempty(S)
end

@testset "MyPriorityQueue" begin
    function f(x,ext)
     return x
    end

    S = MQ.MyPriorityQueue{Int64,Float64}(f,nothing)
    @test MQ.isempty(S)
    MQ.append!(S,4)
    MQ.append!(S,1)
    MQ.append!(S,3)
    MQ.append!(S,5)
    @test MQ.length(S) == 4
    MQ.append!(S,2)
    @test MQ.length(S) == 5

    @test MQ.pop!(S) == 1
    @test MQ.pop!(S) == 2
    @test MQ.length(S) == 3
    @test MQ.isempty(S) == false

    S2 = MQ.MyPriorityQueue{Int64,Float64}(f,nothing)
    MQ.append!(S2,5)
    MQ.append!(S2,3)
    MQ.append!(S2,4)
    @test S==S2 ? true : false
    MQ.pop!(S2)
    @test S==S2 ? false : true

    empty!(S)
    @test MQ.isempty(S)
end

end  # module TestMain
