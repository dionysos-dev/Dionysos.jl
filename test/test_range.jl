using StaticArrays
using LinearAlgebra

struct SCont{N}
    elems::Set{NTuple{N,Int}}
end

function _ranges(ub::NTuple{N,Int}) where N
    return ntuple(i -> 1:ub[i], Val(N))
end

Base.in(x, S::SCont) = x ∈ S.elems

function DoTest()
    println("Any with Set, allocate?")
    S = SCont(Set(collect(Iterators.product(1:1000, 1:200))))
    display(S)
    @time (1, 2) ∈ S
    ub = (50, 100)
    I1 = Iterators.product(1:50, 1:100)
    I2 = Iterators.product(1:ub[1], 1:ub[2])
    I3 = Iterators.product(_ranges(ub)...)
    all(x -> x ∈ S, I1)
    all(x -> x ∈ S, I3)
    @time all(x -> x ∈ S, I1)
    @time all(x -> x ∈ S, I3)
    println()
end

println("---------------------------------------------------------------------")
DoTest()
