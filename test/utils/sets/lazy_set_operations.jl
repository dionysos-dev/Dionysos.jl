module TestMain

using Test
using Dionysos
const DI = Dionysos
const UT = DI.Utils

@testset "LazySetOperations - broader coverage" begin
    # Common shapes
    A1 = UT.HyperRectangle([-1.0, -1.0], [1.0, 1.0])
    A2 = UT.HyperRectangle([2.0, 2.0], [3.0, 3.0])   # disjoint from A1
    A3 = UT.HyperRectangle([0.5, 0.5], [1.5, 1.5])   # overlaps A1

    c0 = [0.0, 0.0]
    c1 = [1.0, 1.0]
    P = [1.0 0.0; 0.0 1.0]
    B0 = UT.Ellipsoid(P, c0)
    B1 = UT.Ellipsoid(P, c1)

    # ----------------------------
    # LazySetUnion basic API: isempty, length, iterate
    # ----------------------------
    U0 = UT.LazySetUnion{2, Float64}()
    @test isempty(U0)
    @test length(U0) == 0
    @test collect(U0) == UT.AbstractSetNode{2, Float64}[]  # might be [] depending on exports

    Urect = UT.LazySetUnion([A1, A2])
    @test !isempty(Urect)
    @test length(Urect) == 2
    @test collect(Urect)[1] == A1
    @test collect(Urect)[2] == A2

    # ----------------------------
    # Membership: Union / Minus / Intersection
    # ----------------------------
    @test ([0.0, 0.0] ∈ Urect) == true     # inside A1
    @test ([2.5, 2.5] ∈ Urect) == true     # inside A2
    @test ([10.0, 10.0] ∈ Urect) == false

    M = UT.LazySetMinus(A1, B0)            # A1 \ B0
    @test ([0.0, 0.0] ∈ M) == false        # center likely inside ellipsoid
    @test ([1.0, 0.5] ∈ M) == true
    @test ([0.95, 0.5] ∈ M) == true

    I0 = UT.LazySetIntersection{2, Float64}()  # empty intersection: all(...) over empty vector is true
    @test ([123.0, 456.0] ∈ I0) == true

    I1 = intersect(I0, A1)                 # Intersection ∩ node
    @test ([0.0, 0.0] ∈ I1) == true
    @test ([10.0, 10.0] ∈ I1) == false

    # ----------------------------
    # add_set! / add_set: node + union
    # ----------------------------
    U = UT.LazySetUnion{2, Float64}()
    UT.add_set!(U, A1)
    @test length(U) == 1
    UT.add_set!(U, UT.LazySetUnion([A2]))
    @test length(U) == 2
    @test U.sets[1] == A1
    @test U.sets[2] == A2

    Ucopy = UT.add_set(U, A3)  # non-mutating
    @test length(U) == 2
    @test length(Ucopy) == 3

    # invalid type into union triggers error
    @test_throws ArgumentError UT.add_set!(UT.LazySetUnion{2, Float64}(), 123)

    # ----------------------------
    # add_set(Minus, s) and remove_set(Minus, s)
    # (forces _promote_copy paths for node child)
    # ----------------------------
    M0 = UT.LazySetMinus(A1, A2)        # B is node
    M_added = UT.add_set(M0, A3)        # promotes A to union and adds A3
    @test M_added.A isa UT.LazySetUnion{2, Float64}
    @test length(M_added) == 2          # length(Minus) proxies A: now union of 2 nodes

    M_removed = UT.remove_set(M0, A3)   # promotes B to union and adds A3 to holes
    @test M_removed.B isa UT.LazySetUnion{2, Float64}

    # ----------------------------
    # Intersection branches
    # ----------------------------

    # Union ∩ node already in your test, but also ensure empty intersections are filtered out
    Iu = intersect(Urect, A3)
    # A2 is disjoint from A3 -> should not appear if intersect returns empty rect filtered by _push_nonempty!
    @test length(Iu) == 1
    @test Iu.sets[1] == UT.HyperRectangle([0.5, 0.5], [1.0, 1.0])

    # Union ∩ Union (exercises nested loops + _push_nonempty!)
    Uleft = UT.LazySetUnion([A1, A2])
    Uright = UT.LazySetUnion([A3])
    Iuu = intersect(Uleft, Uright)
    @test length(Iuu) == 1
    @test Iuu.sets[1] == UT.HyperRectangle([0.5, 0.5], [1.0, 1.0])

    # Minus ∩ node : (A\B)∩S = (A∩S)\(B∩S)
    Mn = UT.LazySetMinus(A1, UT.LazySetUnion([A2])) # B as union (still allowed)
    MnI = intersect(Mn, A3)
    @test MnI isa UT.LazySetMinus
    @test !isempty(MnI.A)   # A∩S should be non-empty
    # B∩S should be empty union (since A2 disjoint from A3)
    @test isempty(MnI.B)

    # Union ∩ Minus and Minus ∩ Union
    Ux = UT.LazySetUnion([A1])
    My = UT.LazySetMinus(A3, A2)  # A3 \ A2 (A2 disjoint anyway)
    Um = intersect(Ux, My)
    @test Um isa UT.LazySetMinus
    @test !isempty(Um.A)
    @test isempty(Um.B)

    mu = intersect(My, Ux)
    @test mu isa UT.LazySetMinus

    # Minus ∩ Minus: (A\B)∩(C\D) = (A∩C)\(B∪D)
    M1 = UT.LazySetMinus(A1, A2)
    M2 = UT.LazySetMinus(A3, UT.LazySetUnion([A2]))
    M12 = intersect(M1, M2)
    @test M12 isa UT.LazySetMinus
    @test !isempty(M12.A)           # A1 ∩ A3 non-empty
    @test M12.B isa UT.LazySetUnion{2, Float64}  # Bunion created
    @test length(M12.B) >= 1

    # Intersection node accumulator (Intersection ∩ anything)
    Iacc = UT.LazySetIntersection{2, Float64}()

    # ----------------------------
    # set_in_period wrappers (Union and Minus)
    # (assumes you have set_in_period for HyperRectangle nodes as earlier)
    # ----------------------------
    periodic_dims = UT.SVector(1)
    periods = UT.SVector(1.0)
    start = UT.SVector(0.0)

    Rwrap = UT.HyperRectangle([0.8, 0.0], [1.2, 1.0])
    Uwrap = UT.LazySetUnion([Rwrap])
    WU = UT.set_in_period(Uwrap, periodic_dims, periods, start)
    @test WU isa UT.LazySetUnion{2, Float64}
    @test length(WU) == 2  # wrapped into two rectangles

    Mwrap = UT.LazySetMinus(Uwrap, UT.LazySetUnion([A2]))
    WM = UT.set_in_period(Mwrap, periodic_dims, periods, start)
    @test WM isa UT.LazySetMinus
    @test WM.A isa UT.LazySetUnion{2, Float64}
end

println("End test")
end # module
