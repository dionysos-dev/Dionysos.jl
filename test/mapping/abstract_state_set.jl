module TestMain

using Test
using StaticArrays
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const MP = DI.Mapping

@testset "AbstractStateSet" begin
    # ---- Make a tiny explicit grid mapping with a known universe ----
    grid = MP.GridFree(SVector(0.0, 0.0), SVector(1.0, 1.0))
    m = MP.ExplicitGridMapping(grid)

    # Add 4 positions => 4 states
    q11 = MP.add_pos!(m, (0, 0))  # state 1
    q12 = MP.add_pos!(m, (1, 0))  # state 2
    q21 = MP.add_pos!(m, (0, 1))  # state 3
    q22 = MP.add_pos!(m, (1, 1))  # state 4

    @test MP.get_n_state(m) == 4
    @test Set(MP.enum_states(m)) == Set(1:4)

    # ----------------------------
    # ExplicitIdSet
    # ----------------------------
    S = MP.ExplicitIdSet{2}()
    @test MP.get_n_state(S, m) == 0

    MP.add_state!(S, m, q11)
    MP.add_state!(S, m, q22)
    @test MP.contains_state(S, m, q11)
    @test !MP.contains_state(S, m, q12)
    @test MP.contains_state(S, m, q22)

    @test Set(MP.enum_states(S, m)) == Set([q11, q22])
    @test MP.get_n_state(S, m) == 2

    MP.remove_state!(S, m, q11)
    @test !MP.contains_state(S, m, q11)
    @test MP.get_n_state(S, m) == 1

    MP.empty_states!(S)
    @test MP.get_n_state(S, m) == 0

    # stateset_from_states default
    S2 = MP.stateset_from_states(m, [q12, q21])
    @test S2 isa MP.ExplicitIdSet{2}
    @test Set(MP.enum_states(S2, m)) == Set([q12, q21])

    # ----------------------------
    # MappingSet: read-only "all states"
    # ----------------------------
    All = MP.MappingSet{2}()
    @test MP.get_n_state(All, m) == 4
    @test MP.contains_state(All, m, 1)
    @test MP.contains_state(All, m, 4)
    @test !MP.contains_state(All, m, 999)

    @test_throws ErrorException MP.add_state!(All, m, 1)
    @test_throws ErrorException MP.remove_state!(All, m, 1)
    @test_throws ErrorException MP.empty_states!(All)

    # ----------------------------
    # add_set!/remove_set! on ExplicitIdSet through geometry + mapping
    # ----------------------------
    S3 = MP.ExplicitIdSet{2}()

    rect = UT.HyperRectangle(SVector(0.0, 0.0), SVector(1.0, 1.0))
    added = MP.add_set!(S3, m, rect, MP.OUTER)   # uses get_states_from_set(m, rect, ...)
    @test all(q -> MP.contains_state(S3, m, q), added)
    @test MP.get_n_state(S3, m) == length(Set(added))

    removed = MP.remove_set!(S3, m, rect, MP.OUTER)
    @test MP.get_n_state(S3, m) == 0
    @test all(q -> !MP.contains_state(S3, m, q), removed)

    # ----------------------------
    # UnionStateSet and SetMinusStateSet
    # ----------------------------
    A = MP.stateset_from_states(m, [q11, q12])
    B = MP.stateset_from_states(m, [q12, q21])

    U = MP.UnionStateSet{2, typeof(A), typeof(B)}(A, B)
    @test MP.contains_state(U, m, q11)
    @test MP.contains_state(U, m, q12)
    @test MP.contains_state(U, m, q21)
    @test !MP.contains_state(U, m, q22)
    @test Set(MP.enum_states(U, m)) == Set([q11, q12, q21])

    D = MP.SetMinusStateSet{2, typeof(A), typeof(B)}(A, B)  # A \ B = {q11}
    @test MP.contains_state(D, m, q11)
    @test !MP.contains_state(D, m, q12)
    @test !MP.contains_state(D, m, q21)
    @test Set(MP.enum_states(D, m)) == Set([q11])

    # ----------------------------
    # ImplicitStateSet (geometry-backed)
    # ----------------------------
    I = MP.ImplicitStateSet{2}()

    # Add a rectangle to A (allowed), no holes in B yet
    MP.add_set!(I, m, UT.HyperRectangle(SVector(-0.4, -0.4), SVector(0.4, 0.4)))
    # The cell at pos (0,0) has center (0,0), corners (+/-0.5)
    # CENTER should be true (center in rect), INNER should be false (corners not all inside)
    @test MP.contains_state(I, m, q11; incl_mode = MP.CENTER)
    @test !MP.contains_state(I, m, q11; incl_mode = MP.INNER)
    @test MP.contains_state(I, m, q11; incl_mode = MP.OUTER)

    # Now remove the center region as a hole: put same rect into B
    MP.remove_set!(I, m, UT.HyperRectangle(SVector(-0.2, -0.2), SVector(0.2, 0.2)))
    @test !MP.contains_state(I, m, q11; incl_mode = MP.CENTER)

    MP.empty_states!(I)
    # After empty, A and B are empty => no state contained
    @test !MP.contains_state(I, m, q11; incl_mode = MP.OUTER)
end

end # module
