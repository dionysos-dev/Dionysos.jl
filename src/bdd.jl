using CUDD

struct IntTupleSet{N}
    manager::Ptr{Manager}
    variables::NTuple{N, Vector{Ptr{Nothing}}}
    size::NTuple{N, Int}
    function IntTupleSet(manager::Ptr{Manager}, size::NTuple{N, Int}) where N
        variables = create_variables.(manager, size)
        return new{N}(manager, variables, size)
    end
end

function nbits(n)
    n -= 1
    count = 0
    while n > 0
        count += 1
        n >>= 1
    end
    return count
end

function create_variables(manager, n)
    num = nbits(n)
    vars = [CUDD.Cudd_bddNewVar(manager) for i in 1:num]
end

manager = initialize_cudd()
IntTupleSet(manager, (2, 4))
