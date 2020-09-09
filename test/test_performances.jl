using StaticArrays
using LinearAlgebra

# Is (...) a perf-killer if we used with collection?
function TestDotsColl()
    function fooWithSpaltting(args...)
        for x in args
            i = x
        end
    end
    function fooWithCollection(arg_coll)
        for x in arg_coll
            i = x
        end
    end
    arg_coll = 1:50
    @time fooWithSpaltting(arg_coll...)
    @time fooWithCollection(arg_coll) # much better
end

# Is (...) a perf-killer if used often with 1 arg?
function TestDotsSingle()
    function barWithSplatting(x...)
        x[1] + 1
    end
    function barWithCollection(x)
        x + 1
    end
    x = 5
    @time barWithSplatting(x)
    @time barWithCollection(x) # Better
end

# map vs generators with eachindex
function TestMapGenerator()
    n = 10
    x = Tuple(rand(n))
    y = Tuple(rand(n))
    v = Tuple(rand(n))
    @time Tuple(round.(Int, x .+ v.*y)) # better
    @time Tuple(round(Int, x[i] + v[i]*y[i]) for i in eachindex(x))
end

# Comparison UnitRange, expand..., and a:b
function TestUnitRange()
    a = 20
    b = 100
    lb = collect(a:b)
    ub = collect(a+20:b+20)
    @time Tuple(UnitRange(x...) for x in zip(lb, ub))
    @time Tuple(lb[i]:ub[i] for i in eachindex(lb))
    @time Tuple(UnitRange(lb[i],ub[i]) for i in eachindex(lb)) # slightly better
end

# Check all vector
function TestCheckAllVectors()
    n = 100
    x = Tuple(rand(n))
    y = Tuple(rand(n))
    z = Tuple(rand(n))
    @time all(x .< y .<= z) # Much faster
    @time all(i -> (x[i] <= y[i] <= z[i]), eachindex(x))
end

# Check all matrix
function TestCheckAllSMatrices()
    A = SMatrix{6,6}(rand(6, 6))
    x = SVector{6}(rand(6))
    b = SVector{6}(rand(6))
    idx = SVector{6}(1:6)
    @time all(abs.(A*x) .<= b) # Much faster
    @time all(x -> x <= 0, abs.(A*x) - b)
    @time all(i -> abs(dot(A[i,idx],x)) <= b[i], eachindex(x))
end

# Delete value in vector
function TestDelete()
    n = 50
    a1 = vcat(collect(1:2*n), collect(1:2*n))
    a2 = vcat(collect(1:2*n), collect(1:2*n))
    @time deleteat!(a1, findall(x -> x == n, a1))
    @time deleteat!(a2, a2 .== n) # faster
end

# Delete value in sorted vector
function TestDeleteSorted()
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
end

# Set vs Sorted Vector for membership
function TestSetVsVectorMembership()
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
    @time inSet(a2, n) # Faster
end

# Create tuple with generator?
function TestVectorToTuple()
    x = [1,2,3,3,3,3,6,5,9,4,5,6,9,9,8,9,7,10,2]
    sort!(x)
    idx = searchsorted(x, 3)
    @time Tuple(x[idx]) # Faster
    @time Tuple(x[i] for i in idx)
end

# Is sizehint preserved after empty!?
function TestSizehintEmpty()
    x = Int[]
    @time push!(x, 1)
    sizehint!(x, 10000)
    empty!(x)
    @time append!(x, 1:1001)
    # Yes
end

# Delete or not transitions?
function TestTransitionsList()
    function get(a, x)
        i = searchsorted(a, x)
        return i
    end
    function getAndRemove!(a, x)
        i = searchsorted(a, x)
        deleteat!(a, i)
        return i
    end
    n = Int(1e4)
    x1 = collect(1:n)
    x2 = collect(1:n)
    @time for x in 1:n/2 get(x1, x) end
    @time for x in 1:n/2 getAndRemove!(x2, x) end
    # Comparable but get probably more regular in terms of time
end

# Set or vec for transitions?
function TestTransitionsSetVsVector()
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
    x1 = [(i÷100, i÷5, i) for i in 1:n]
    x2 = Set((i÷100, i÷5, i) for i in 1:n)
    @time for x in 1:n/2 filterVec(x1, x) end # Much much Faster
    @time for x in 1:n/2 filterSet(x2, x) end
end

# Set vs Vector for setdiff
function TestSetVsVectorDiff()
    n = 50000
    X1 = collect(1:n)
    X2 = Set(1:n)
    Z1 = collect(0.8*n:1.2*n)
    Z2 = Set(0.8*n:1.2*n)
    @time setdiff!(X1, Z1)
    @time setdiff!(X2, Z2) # Much faster
end

# Tuple + SVector broadcast
function TestTuplePlusSVector()
    function mySumTS1(X::SVector{N}, Y::NTuple{N}) where N
        return Tuple(X .+ Y)
    end # return Tuple
    function mySumTS2(X::SVector{N}, Y::NTuple{N}) where N
        return ntuple(i -> X[i] + Y[i], Val(N))
    end # return Tuple
    function mySumTS3(X::SVector{N}, Y::NTuple{N}) where N
        return X .+ Y
    end # return SVector
    function mySumTS4(X::SVector{N}, Y::NTuple{N}) where N
        return X .+ SVector(Y)
    end # return SVector
    X = @SVector [1,2,3,4,5]
    Y = (5,6,7,8,9)
    @time for i in 1:1000 x = mySumTS1(X, Y) end
    @time for i in 1:1000 x = mySumTS2(X, Y) end # Superfast both
    # @code_warntype mySumTS3(X, Y)
    # @code_warntype mySumTS4(X, Y)
end

# Float*Float or Float*Int?
function TestOperationFloatInt()
    N = 50000
    x = 1.0
    y = 2.0
    @time for i in 1:N x*2 + 5 + y/3 end
    @time for i in 1:N x*2 + 5.0 + y/3.0 end
    # Similar
end

# Tail-recursive hashing for tuples?
using Base.Cartesian

struct MyTuple{N,S<:NTuple{N}}
    T::S
end

Base.tail(x::MyTuple) = MyTuple(Base.tail(x.T))
Base.hash(::MyTuple{0,S}, h::UInt) where S = h + Base.tuplehash_seed
Base.hash(x::MyTuple, h::UInt) = Base.hash(Base.tail(x), hash(x.T[1], h))

struct YourTuple{N,S<:NTuple{N}}
    T::S
end

Base.hash(x::YourTuple{0,S}, h::UInt) where S = h + Base.tuplehash_seed
@generated function Base.hash(x::YourTuple{N}, h::UInt) where N
    quote
        h += Base.tuplehash_seed
        @nexprs $N i -> h = hash(x.T[$N-i+1], h)
    end
end

struct HisTuple{N,S<:NTuple{N}}
    T::S
end

Base.hash(x::HisTuple{N}, h::UInt) where N = _hash(x, Val(N), h + Base.tuplehash_seed)
_hash(x, ::Val{0}, h) = h
function _hash(x, ::Val{N}, h) where N
    return _hash(x, Val(N-1), hash(x.T[N], h))
end

function TestHashTuple(n)
    D1 = [(1,2,i) for i in 1:n]
    D2 = [MyTuple((1,2,i)) for i in 1:n]
    D3 = [YourTuple((1,2,i)) for i in 1:n]
    D4 = [HisTuple((1,2,i)) for i in 1:n]
    for (d1, d3, d4) in zip(D1, D3, D4)
        @assert hash(d1) == hash(d3) == hash(d4)
    end
    @time for d in D1
        x = hash(d)
    end
    @time for d in D2
        x = hash(d)
    end
    @time for d in D3
        x = hash(d)
    end
    @time for d in D4
        x = hash(d)
    end
end

function TestInTuple(n)
    D1 = Set((1,2,i) for i in 1:n)
    D2 = Set(MyTuple((1,2,i)) for i in 1:n)
    D3 = Set(YourTuple((1,2,i)) for i in 1:n)
    D4 = Set(HisTuple((1,2,i)) for i in 1:n)
    F1 = [(1,2,round(Int,i)) for i = 0.8*n:1.2*n]
    F2 = [MyTuple((1,2,round(Int,i))) for i = 0.8*n:1.2*n]
    F3 = [YourTuple((1,2,round(Int,i))) for i = 0.8*n:1.2*n]
    F4 = [HisTuple((1,2,round(Int,i))) for i = 0.8*n:1.2*n]
    @time for d in D1
        x = hash(d)
    end
    @time for d in D2
        x = hash(d)
    end
    @time for d in D3
        x = hash(d)
    end
    @time for d in D4
        x = hash(d)
    end
    println()
    @time for f in F1
        f ∈ D1
    end
    @time for f in F2
        f ∈ D2
    end
    @time for f in F3
        f ∈ D3
    end
    @time for f in F4
        f ∈ D4
    end
end

# Are function with condition with static values inlined smartly?
@inline foo(x, b) = b ? x - 4 : x + 7
testTrue(x) = foo(x, true)
testFalse(x) = foo(x, false)

function TestStaticCondition(n)
    @time for i in 1:n
        b = i != 0
        x = testTrue(i) end
    @time for i in 1:n
        b = i == 0
        x = testFalse(i)
    end
    @time for i in 1:n
        b = i != 0
        x = foo(i, b)
    end
    @time for i in 1:n
        b = i == 0
        x = foo(i, b)
    end
end

# Splatting in function
foo(x, y, z=1.0) = x .+ y .+ z
bar1(args...) = foo(ones(10), args...)
bar2(y, z=1.0) = foo(ones(10), y, z)

function TestSplatting(n)
    for i = 1:n res = bar1(8, 9.0) end
    for i = 1:n res = bar2(8, 9.0) end
    @time for i = 1:n res = bar1(8, 9.0) end
    @time for i = 1:n res = bar2(8, 9.0) end
    @time for i = 1:n res = bar1(8) end
    @time for i = 1:n res = bar2(8) end
end

# Anonymous function and type stability
foo(x::Int, y::Int) = 5*y
foo(x::Int, y::Float64) = "a"
foo(x::Float64, y) = 5.0*y
bar(u) = let x = u[1]
    y -> foo(x, y)
end
function vim(u)
    the_bar = bar(u)
    z = the_bar(8.0)
    z + 2
end
@code_warntype vim([6])

# From https://github.com/c42f/FastClosures.jl
# code_warntype problem
function f1()
    if true
    end
    r = 1
    cb = ()->r
    identity(cb)
end

# code_warntype clean
function f2()
    if true
    end
    r = 1
    cb = let r = r
        ()->r
    end
    identity(cb)
end

@code_warntype f1()
@code_warntype f2()

# Closure with multiple arguments length
foo(x, y) = x*y
foo(x, y, z) = x + y + z
bar(x) = let x = x
    (args...) -> foo(x[1], args...)
end
vim1(x) = let x = x
    y -> foo(x[1], y)
end
vim2(x) = let x = x
    (y, z) -> foo(x[1], y, z)
end

function TestClosure(n)
    _bar = bar([5])
    _vim1 = vim1([5])
    _vim2 = vim2([5])
    @time for i = 1:n display(_bar(i)) end
    @time for i = 1:n display(_bar(i, i + 1)) end
    @time for i = 1:n display(_vim1(i)) end
    @time for i = 1:n display(_vim2(i, i + 1)) end
end

@code_warntype bar([5])(6)
