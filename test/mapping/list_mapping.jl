module TestListMapping

using Test
using StaticArrays
using Dionysos
const MP = Dionysos.Mapping

@testset "ListMapping" begin
    pts = [SVector(0.0, 1.0), SVector(2.0, 3.0)]
    m = MP.ListMapping(pts)

    @test MP.get_n_state(m) == 2
    @test collect(MP.enum_states(m)) == [1, 2]
    @test MP.get_coord_by_state(m, 1) == SVector(0.0, 1.0)
    @test MP.get_coord_by_state(m, 2) == SVector(2.0, 3.0)

    @test MP.get_state_by_coord(m, SVector(0.0, 1.0)) == 1
    @test MP.get_state_by_coord(m, [2.0, 3.0]) == 2   # accepts Vector -> SVector conversion

    # out of dictionary key => KeyError
    @test_throws KeyError MP.get_state_by_coord(m, SVector(9.0, 9.0))
end

@testset "ListMapping deduplication" begin
    pts = [SVector(1.0, 2.0), SVector(1.0, 2.0), SVector(3.0, 4.0)]
    m = MP.ListMapping(pts)

    # keeps first occurrence only
    @test MP.get_n_state(m) == 2
    @test MP.get_coord_by_state(m, 1) == SVector(1.0, 2.0)
    @test MP.get_coord_by_state(m, 2) == SVector(3.0, 4.0)
    @test MP.get_state_by_coord(m, SVector(1.0, 2.0)) == 1
end

@testset "ListMapping constructors" begin
    # NTuple constructor
    m1 = MP.ListMapping([(0.0, 1.0), (2.0, 3.0)])
    @test MP.get_n_state(m1) == 2
    @test MP.get_coord_by_state(m1, 2) == SVector(2.0, 3.0)

    # Vector-of-vectors constructor
    m2 = MP.ListMapping([[0.0, 1.0], [2.0, 3.0]])
    @test MP.get_n_state(m2) == 2
    @test MP.get_coord_by_state(m2, 1) == SVector(0.0, 1.0)

    # empty vector-of-vectors should assert
    @test_throws AssertionError MP.ListMapping(Vector{Vector{Float64}}())
end

@testset "ListMapping empty! and convert_to_list_mapping" begin
    m = MP.ListMapping([SVector(0.0, 0.0), SVector(1.0, 1.0)])

    @test MP.convert_to_list_mapping(m) === m

    empty!(m)
    @test MP.get_n_state(m) == 0
    @test collect(MP.enum_states(m)) == Int[]
    @test_throws BoundsError MP.get_coord_by_state(m, 1)
end

@testset "ListMapping union-like behavior via concatenation" begin
    mA = MP.ListMapping([SVector(0.0, 1.0), SVector(2.0, 3.0)])
    mB = MP.ListMapping([SVector(2.0, 3.0), SVector(4.0, 5.0)])

    # No union! method for mappings; emulate by building a new one from concatenated coords
    pts = vcat(
        [MP.get_coord_by_state(mA, q) for q in MP.enum_states(mA)],
        [MP.get_coord_by_state(mB, q) for q in MP.enum_states(mB)],
    )

    mU = MP.ListMapping(pts)

    @test MP.get_n_state(mU) == 3
    @test MP.get_state_by_coord(mU, SVector(0.0, 1.0)) == 1
    @test MP.get_state_by_coord(mU, SVector(2.0, 3.0)) in (2, 1, 3)  # id depends on insertion order; we only care it exists
    @test MP.get_state_by_coord(mU, SVector(4.0, 5.0)) >= 1
end

end # module
