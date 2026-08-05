module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

# A state set that overrides nothing, to check the abstract interface stubs error.
struct _UnimplementedStateSet <: MP.AbstractStateSet{2} end

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
    # FullStateSet: read-only "all states"
    # ----------------------------
    All = MP.FullStateSet{2}()
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

    rect = UT.box(SVector(0.0, 0.0), SVector(1.0, 1.0))
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
    MP.add_set!(I, m, UT.box(SVector(-0.4, -0.4), SVector(0.4, 0.4)), MP.OUTER)
    # The cell at pos (0,0) has center (0,0), corners (+/-0.5)
    # CENTER should be true (center in rect), INNER should be false (corners not all inside)
    @test MP.contains_state(I, m, q11)

    # Now remove the center region as a hole: put same rect into B
    MP.remove_set!(I, m, UT.box(SVector(-0.2, -0.2), SVector(0.2, 0.2)))
    @test !MP.contains_state(I, m, q11)

    MP.empty_states!(I)
    # After empty, A and B are empty => no state contained
    @test !MP.contains_state(I, m, q11)
end

@testset "ImplicitStateSet inclusion modes" begin
    # 1×1 cells: cell (0,0) has center (0,0) and corners (±0.5, ±0.5); cell (1,1) has
    # center (1,1) and corners in [0.5, 1.5]².
    grid = MP.GridFree(SVector(0.0, 0.0), SVector(1.0, 1.0))
    m = MP.ExplicitGridMapping(grid)
    q00 = MP.add_pos!(m, (0, 0))
    q11 = MP.add_pos!(m, (1, 1))

    # Invalid states are never contained (guard before any geometry).
    Sc = MP.ImplicitStateSet(UT.box(SVector(-0.4, -0.4), SVector(0.4, 0.4)), MP.CENTER)
    @test !MP.contains_state(Sc, m, 999)

    # CENTER: contained iff the cell center lies in the set.
    @test MP.contains_state(Sc, m, q00)     # center (0,0) ∈ [-0.4, 0.4]²
    @test !MP.contains_state(Sc, m, q11)    # center (1,1) ∉

    # INNER (conservative): contained iff *every* corner lies in the set.
    Si = MP.ImplicitStateSet(UT.box(SVector(-1.0, -1.0), SVector(1.0, 1.0)), MP.INNER)
    @test MP.contains_state(Si, m, q00)     # all corners (±0.5) inside [-1, 1]²
    @test !MP.contains_state(Si, m, q11)    # corner (1.5, 1.5) is outside

    # OUTER (sufficient): contained iff *some* sample lies in the set. Here the center
    # is outside but a corner is inside, so the corner scan decides membership.
    So = MP.ImplicitStateSet(UT.box(SVector(0.3, 0.3), SVector(1.0, 1.0)), MP.OUTER)
    @test MP.contains_state(So, m, q00)     # center (0,0) ∉, but corner (0.5, 0.5) ∈
    Sfar = MP.ImplicitStateSet(UT.box(SVector(5.0, 5.0), SVector(6.0, 6.0)), MP.OUTER)
    @test !MP.contains_state(Sfar, m, q00)  # no sample of cell (0,0) reaches the set
end

@testset "ImplicitStateSet periodic add_set!/remove_set!" begin
    # Period 4 in dim 1; the grid origin is aligned to start + h/2 = 0.5.
    periodic_dims = SVector(1)
    periods = SVector(4.0)
    start = SVector(0.0)
    h = SVector(1.0, 1.0)
    grid = MP.GridFree(SVector(0.5, 0.0), h)
    pm = MP.PeriodicGridMapping(periodic_dims, periods, start, MP.ExplicitGridMapping(grid))

    S = MP.ImplicitStateSet{2}()

    # A box outside the fundamental period is wrapped back into it: [4.5, 5.5] → [0.5, 1.5].
    MP.add_set!(S, pm, UT.box(SVector(4.5, 0.0), SVector(5.5, 1.0)), MP.OUTER)
    @test SVector(1.0, 0.5) ∈ S.set     # inside the wrapped region
    @test SVector(3.0, 0.5) ∉ S.set     # elsewhere in the period

    # Removal is wrapped the same way: carve [4.5, 5.0] → [0.5, 1.0] out of the set.
    MP.remove_set!(S, pm, UT.box(SVector(4.5, 0.0), SVector(5.0, 1.0)))
    @test SVector(0.7, 0.5) ∉ S.set     # in the carved-out region
    @test SVector(1.3, 0.5) ∈ S.set     # still inside the remainder
end

@testset "MappedStateSet (public surface)" begin
    grid = MP.GridFree(SVector(0.0, 0.0), SVector(1.0, 1.0))
    m = MP.ExplicitGridMapping(grid)
    q1 = MP.add_pos!(m, (0, 0))
    q2 = MP.add_pos!(m, (1, 0))
    q3 = MP.add_pos!(m, (0, 1))
    q4 = MP.add_pos!(m, (1, 1))

    # bundles a set with its mapping, so calls no longer thread (set, mapping)
    ms = MP.MappedStateSet(MP.ExplicitIdSet{2}(), m)
    @test MP.get_dim(ms) == 2
    @test MP.get_mapping(ms) === m
    @test MP.get_set(ms) isa MP.ExplicitIdSet{2}
    @test MP.get_n_state(ms) == 0

    MP.add_state!(ms, q1)
    @test MP.contains_state(ms, q1)
    @test !MP.contains_state(ms, q2)
    @test MP.get_n_state(ms) == 1

    MP.add_states!(ms, [q2, q3])
    @test Set(MP.enum_states(ms)) == Set([q1, q2, q3])
    @test MP.get_n_state(ms) == 3

    MP.remove_state!(ms, q1)
    @test !MP.contains_state(ms, q1)

    # add_set!/remove_set! forward through to the (set, mapping) methods
    rect = UT.box(SVector(0.0, 0.0), SVector(1.0, 1.0))
    MP.add_set!(ms, rect, MP.OUTER)
    @test all(q -> MP.contains_state(ms, q), (q1, q2, q3, q4))
    MP.remove_set!(ms, rect, MP.OUTER)
    @test MP.get_n_state(ms) == 0

    # copy is independent in its set but shares the (immutable) mapping
    MP.add_state!(ms, q1)
    ms2 = copy(ms)
    MP.add_state!(ms2, q4)
    @test MP.contains_state(ms2, q4)
    @test !MP.contains_state(ms, q4)     # mutating the copy leaves the original alone
    @test MP.get_mapping(ms2) === m
end

@testset "UniqueStates dedup (union get_n_state, unique_states)" begin
    grid = MP.GridFree(SVector(0.0, 0.0), SVector(1.0, 1.0))
    m = MP.ExplicitGridMapping(grid)
    q1 = MP.add_pos!(m, (0, 0))
    q2 = MP.add_pos!(m, (1, 0))
    q3 = MP.add_pos!(m, (0, 1))

    A = MP.stateset_from_states(m, [q1, q2])
    B = MP.stateset_from_states(m, [q2, q3])   # overlaps A on q2

    # get_n_state on the lazy union counts distinct states → exercises
    # UniqueStates length (it must not double-count the shared q2).
    U = MP.UnionStateSet{2, typeof(A), typeof(B)}(A, B)
    @test MP.get_n_state(U, m) == 3

    D = MP.SetMinusStateSet{2, typeof(A), typeof(B)}(A, B)  # A \ B = {q1}
    @test Set(MP.enum_states(D, m)) == Set([q1])
    @test MP.get_n_state(D, m) == 1   # SetMinus override (lazy filter has no length)

    # the deduplicating iterator directly
    dup = MP.unique_states([q1, q2, q2, q1, q3])
    @test eltype(dup) == Int
    @test length(dup) == 3
    @test Set(dup) == Set([q1, q2, q3])
end

@testset "AbstractStateSet interface stubs error" begin
    grid = MP.GridFree(SVector(0.0, 0.0), SVector(1.0, 1.0))
    m = MP.ExplicitGridMapping(grid)
    S = _UnimplementedStateSet()
    @test_throws ErrorException MP.contains_state(S, m, 1)
    @test_throws ErrorException MP.enum_states(S, m)
    @test_throws ErrorException MP.add_state!(S, m, 1)
    @test_throws ErrorException MP.remove_state!(S, m, 1)
    @test_throws ErrorException MP.empty_states!(S)
end

end # module
