function foo1(y...)
    i = 0
    for yy in y
        i += yy
    end
end

function foo2(y)
    i = 0
    for yy in y
        i += yy
    end
end

function bar(val)
    return val + 1
end

function predef()
    x = rand(1000)
    i = 10:400
    y = x[i]
    for j in 1:200
        z = y[j]
    end
end

function notpredef()
    x = rand(100)
    i = 10:40
    for j in 1:20
        z = x[i][j]
    end
end

function notpredef2()
    x = rand(100)
    i = 10:40
    im = minimum(i)
    iM = maximum(i)
    for j in 1:20
        z = x[im+j]
    end
end

function expand(x...)
end

function notexpand(x)
end

function DoTest()
    y = 1:5
    @time foo1(y...)
    @time foo2(y)
    D = Dict(5 => 6)
    v = get(D, 5, nothing)
    x = rand(2)
    y = rand(2)
    v = Tuple(rand(2))
    display(v[2])
    @time Tuple(round.(Int, x + v.*y))
    @time Tuple(round(Int, x[i] + v[i]*y[i]) for i in eachindex(x))
    lbI = collect(1:6)
    ubI = collect(6:11)
    display((lbI, ubI))
    @time Tuple(UnitRange(x...) for x in zip(lbI, ubI))
    @time Tuple(lbI[i]:ubI[i] for i in 1:6)
    @time Tuple(UnitRange(lbI[i],ubI[i]) for i in 1:6)
    x = rand(100)
    y = rand(100)
    display((x, y))
    @time all(x.<y.<=y)
    @time all(i -> (x[i]<=y[i]<=y[i]), eachindex(x))
    a1 = [1,2,3,3,3,5,3,5,5,6]
    a2 = [1,2,3,3,3,5,3,5,5,6]
    display(a1)
    @time deleteat!(a1, findall(x -> x == 5, a1))
    @time deleteat!(a2, a2 .== 5)
    display(a2)
    x = [1,2,3,3,3,3,6,5,9,4,5,6,9,9,8,9,7,10,2]
    sort!(x)
    idx = searchsorted(x, 3)
    @time Tuple(x[idx])
    @time Tuple(x[i] for i in idx)
end

# @code_warntype
DoTest()

@time predef()
@time notpredef()
@time notpredef2()

@time for i = 1:100000 expand(5) end
@time for i = 1:100000 notexpand(5) end

x = Int[]
@time push!(x, 1)
sizehint!(x, 10000)
empty!(x)
@time append!(x, 1:1001)

x = [1,2,3]
y = [7,8,9]

# @code_warntype Iterators.product((x[i]:y[i] for i in 1:length(x)))
