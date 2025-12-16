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

@testset "Additional CustomList tests" begin
    # Test empty domain
    empty_dom = DO.CustomList(SVector{2, Float64}[])
    @test DO.get_ncells(empty_dom) == 0
    @test isempty(empty_dom)
    @test DO.enum_elems(empty_dom) == SVector{2, Float64}[]

    # Test type parameters
    elems = [SVector(1.0, 2.0), SVector(3.0, 4.0)]
    dom = DO.CustomList{2, Float64}(elems)
    @test dom isa DO.CustomList{2, Float64}
    @test dom isa DO.DomainType

    # Test convert_to_custom_domain returns the same object
    dom_conv = DO.convert_to_custom_domain(dom)
    @test dom_conv === dom

    # Test get_elem_by_index with out-of-bounds indices
    @test_throws BoundsError DO.get_elem_by_index(dom, 0)
    @test_throws BoundsError DO.get_elem_by_index(dom, 3)

    # Test union! with no new elements
    dom_before = copy(dom.elems)
    DO.union!(dom, DO.CustomList(copy(dom.elems)))
    @test dom.elems == dom_before

    # Test setdiff! with no matching elements
    dom2 = DO.CustomList([SVector(10.0, 20.0)])
    dom_before = copy(dom.elems)
    DO.setdiff!(dom, dom2)
    @test dom.elems == dom_before
end

end # module
