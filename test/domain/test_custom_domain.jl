module TestCustomDomain
using Test
using StaticArrays
using Dionysos
const DO = Dionysos.Domain
@testset "CustomList" begin
    # Create a small list of points
    elems = [SVector(0.0, 1.0), SVector(2.0, 3.0)]
    dom = DO.CustomList(elems)
    @test DO.get_ncells(dom) == 2
    @test DO.enum_elems(dom) == elems
    @test DO.get_elem_by_index(dom, 1) == SVector(0.0, 1.0)
    dom2 = DO.CustomList([SVector(2.0, 3.0), SVector(4.0, 5.0)])
    union!(dom, dom2)
    @test DO.get_ncells(dom) == 3
    setdiff!(dom, DO.CustomList([SVector(2.0, 3.0)]))
    @test DO.get_ncells(dom) == 2
    empty!(dom)
    @test isempty(dom)
end
end # module
