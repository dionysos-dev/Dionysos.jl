module TestCustomDomain
using Test
using StaticArrays
using Dionysos
const DO = Dionysos.Domain
@testset "CustomList" begin
    # Create a small list of points
    elems = [SVector(0.0, 1.0), SVector(2.0, 3.0), SVector(4.0, 5.0), SVector(20.0, 30.0)]
    dom = DO.CustomList(elems)
    @test DO.get_ncells(dom) == 4
    @test DO.enum_elems(dom) == elems
    @test DO.get_elem_by_index(dom, 4) == SVector(20.0, 30.0)
    dom2 = DO.CustomList([SVector(2.0, 3.0), SVector(40.0, 50.0)])
    union!(dom, dom2)
    @test DO.get_ncells(dom) == 5
    println(dom)
    setdiff!(dom, DO.CustomList([SVector(2.0, 3.0)]))
    @test DO.get_ncells(dom) == 4
    println(dom)
    empty!(dom)
    @test isempty(dom)
    println(dom)
end
end 
