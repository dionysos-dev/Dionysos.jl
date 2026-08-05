module TestMain

using Test
using Dionysos
const DI = Dionysos
const UT = DI.Utils

sleep(0.1) # used for good printing
println("Started test")

@testset "MyStack" begin
    S = UT.MyStack{Int64}()
    @test UT.isempty(S)
    UT.append!(S, 1)
    UT.append!(S, 2)
    UT.append!(S, 3)
    UT.append!(S, 4)
    @test UT.length(S) == 4
    UT.append!(S, 5)
    @test UT.length(S) == 5

    @test UT.pop!(S) == 5
    @test UT.pop!(S) == 4
    @test UT.length(S) == 3
    @test UT.isempty(S) == false

    S2 = UT.MyStack{Int64}()
    UT.append!(S2, 1)
    UT.append!(S2, 2)
    UT.append!(S2, 3)
    @test S == S2 ? true : false
    UT.pop!(S2)
    @test S == S2 ? false : true

    empty!(S)
    @test UT.isempty(S)
end

@testset "FIFOQueue" begin
    S = UT.FIFOQueue{Int64}()
    @test UT.isempty(S)
    UT.append!(S, 1)
    UT.append!(S, 2)
    UT.append!(S, 3)
    UT.append!(S, 4)
    @test UT.length(S) == 4
    UT.append!(S, 5)
    @test UT.length(S) == 5

    @test UT.pop!(S) == 1
    @test UT.pop!(S) == 2
    @test UT.length(S) == 3
    @test UT.isempty(S) == false

    S2 = UT.FIFOQueue{Int64}()
    UT.append!(S2, 3)
    UT.append!(S2, 4)
    UT.append!(S2, 5)
    @test S == S2 ? true : false
    UT.pop!(S2)
    @test S == S2 ? false : true

    empty!(S)
    @test UT.isempty(S)
end

@testset "MyPriorityQueue" begin
    function f(x, ext)
        return x
    end

    S = UT.MyPriorityQueue{Int64, Float64}(f, nothing)
    @test UT.isempty(S)
    UT.append!(S, 4)
    UT.append!(S, 1)
    UT.append!(S, 3)
    UT.append!(S, 5)
    @test UT.length(S) == 4
    UT.append!(S, 2)
    @test UT.length(S) == 5

    @test UT.pop!(S) == 1
    @test UT.pop!(S) == 2
    @test UT.length(S) == 3
    @test UT.isempty(S) == false

    S2 = UT.MyPriorityQueue{Int64, Float64}(f, nothing)
    UT.append!(S2, 5)
    UT.append!(S2, 3)
    UT.append!(S2, 4)
    @test S == S2 ? true : false
    UT.pop!(S2)
    @test S == S2 ? false : true

    empty!(S)
    @test UT.isempty(S)
end

end  # module TestMain
