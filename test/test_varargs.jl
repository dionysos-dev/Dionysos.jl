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

function DoTest()::Int
    y = 1:5
    @time foo1(y...)
    @time foo2(y)
    D = Dict(5 => 6)
    v = get(D, 5, nothing)
end

@code_warntype DoTest()
