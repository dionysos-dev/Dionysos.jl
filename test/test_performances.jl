using StaticArrays

struct T{F<:Function}
    func::F
end

function DoTest()
    sleep(0.1)
    println("Is (...) a perf-killer if we used with collect?")
    function foo1(args...)
        for x in args
            i = x
        end
    end
    function foo2(arg_coll)
        for x in arg_coll
            i = x
        end
    end
    arg_coll = 1:50
    @time foo1(arg_coll...)
    @time foo2(arg_coll) # much better
    println()

    println("Is (...) a perf-killer if used often with 1 arg?")
    function bar1(x...)
        x[1] + 1
    end

    function bar2(x)
        x + 1
    end
    x = 5
    @time bar1(x)
    @time bar2(x) # Better
    println()

    println("map vs generators with eachindex")
    n = 10
    x = Tuple(rand(n))
    y = Tuple(rand(n))
    v = Tuple(rand(n))
    @time Tuple(round.(Int, x .+ v.*y)) # better
    @time Tuple(round(Int, x[i] + v[i]*y[i]) for i in eachindex(x))
    println()

    println("Comparison UnitRange, expand..., and a:b")
    a = 20
    b = 100
    lb = collect(a:b)
    ub = collect(a+20:b+20)
    @time Tuple(UnitRange(x...) for x in zip(lb, ub))
    @time Tuple(lb[i]:ub[i] for i in eachindex(lb))
    @time Tuple(UnitRange(lb[i],ub[i]) for i in eachindex(lb)) # slightly better
    println()

    println("Check all")
    n = 100
    x = Tuple(rand(n))
    y = Tuple(rand(n))
    z = Tuple(rand(n))
    @time all(x .< y .<= z) # Much faster
    @time all(i -> (x[i] <= y[i] <= z[i]), eachindex(x))
    println()

    println("Delete value in vector")
    n = 50
    a1 = vcat(collect(1:2*n), collect(1:2*n))
    a2 = vcat(collect(1:2*n), collect(1:2*n))
    @time deleteat!(a1, findall(x -> x == n, a1))
    @time deleteat!(a2, a2 .== n) # faster
    println()

    println("Delete value in sorted array")
    function removeSorted!(V, x)
        idx = searchsorted(V, x)
        deleteat!(V, idx)
    end
    function removeSetdiff!(V, x)
        setdiff!(V, x)
    end
    n = 3000
    a1 = sort(vcat(collect(1:2*n), collect(1:2*n)))
    a2 = sort(vcat(collect(1:2*n), collect(1:2*n)))
    @time removeSorted!(a1, n)
    @time removeSetdiff!(a2, n) # Faster but make a2 unique...
    println()

    println("Push Set vs sorted vector")
    n = 3000000
    a1 = sort(collect(1:2*n))
    sizehint!(a1, n + 2)
    a2 = Set(1:2*n)
    sizehint!(a2, n + 2)
    @time push!(a1, n)
    @time push!(a2, n)
    println()

    println("Set vs sorted vector for membership")
    function inSortedVector(V, x)
        idx = searchsorted(V, x)
        return !isempty(idx)
    end
    function inVector(V, x)
        return in(V, x)
    end
    function inSet(S, x)
        in(S, x)
    end
    n = 3000000
    a1 = sort(collect(1:2*n))
    a2 = Set(1:2*n)
    @time inSortedVector(a1, n)
    @time inVector(a1, n)
    @time inSet(a2, n) # Faster !!!
    println()

    println("Create tuple with generator?")
    x = [1,2,3,3,3,3,6,5,9,4,5,6,9,9,8,9,7,10,2]
    sort!(x)
    idx = searchsorted(x, 3)
    @time Tuple(x[idx]) # Faster
    @time Tuple(x[i] for i in idx)
    println()

    println("Is sizehint preserved after empty!?")
    x = Int[]
    @time push!(x, 1)
    sizehint!(x, 10000)
    empty!(x)
    @time append!(x, 1:1001)
    # Yes
    println()

    println("Delete or not transitions?")
    function get(a, x)
        i = searchsorted(a, x)
        return i
    end
    function getremove(a, x)
        i = searchsorted(a, x)
        deleteat!(a, i)
        return i
    end
    n = Int(1e4)
    x1 = collect(1:n)
    x2 = collect(1:n)
    @time for x = 1:n/2 get(x1, x) end # Comparable but probably more deterministic
    @time for x = 1:n/2 getremove(x2, x) end
    println()

    println("Set or vec for transitions?")
    function filterVec(V, x)
        idx = searchsorted(V, x, by = x -> x[1])
        for i in idx
            V[i]
        end
    end
    function filterSet(S, x)
        for s in filter(e -> e[1] == x, S)
            s
        end
    end
    n = 200
    x1 = [(i÷100, i÷5, i) for i = 1:n]
    x2 = Set((i÷100, i÷5, i) for i = 1:n)
    @time for x = 1:n/2 filterVec(x1, x) end # Much much Faster
    @time for x = 1:n/2 filterSet(x2, x) end
    println()

    println("Set vs Vector for setdiff")
    n = 50000
    X1 = collect(1:n)
    X2 = Set(1:n)
    Z1 = collect(0.8*n:1.2*n)
    Z2 = Set(0.8*n:1.2*n)
    @time setdiff!(X1, Z1)
    @time setdiff!(X2, Z2)
    println()

    println("Tuple + SVector broadcast")
    function sum1(X::SVector{N}, Y::NTuple{N}) where N
        return Tuple(X .+ Y)
    end
    function sum2(X::SVector{N}, Y::NTuple{N}) where N
        return ntuple(i -> X[i] + Y[i], Val(N))
    end
    function sum3(X::SVector{N}, Y::NTuple{N}) where N
        return X .+ Y
    end
    function sum4(X::SVector{N}, Y::NTuple{N}) where N
        return X .+ SVector(Y)
    end
    X = @SVector [1,2,3,4,5]
    Y = (5,6,7,8,9)
    @time for i = 1:1000 x = sum1(X, Y) end
    @time for i = 1:1000 x = sum2(X, Y) end # Superfast both
    # @code_warntype sum3(X, Y)
    # @code_warntype sum4(X, Y)
    println()

    println("Float*Float or Float*Int?")
    N = 50000
    x = 1.0
    y = 2.0
    @time for i = 1:N x*2 + 5 + y/3 end
    @time for i = 1:N x*2 + 5.0 + y/3.0 end
    # Similar
    println()
end

println("---------------------------------------------------------------------")
DoTest()
