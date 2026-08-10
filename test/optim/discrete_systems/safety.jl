module TestMain

import Dionysos
import LazySets
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

@testset "ControllerSafe" begin
    # ----------------------------
    # Build finite X mapping
    # ----------------------------
    lbX = SVector(-5.0, -5.0)
    ubX = SVector(5.0, 5.0)
    x0 = SVector(0.0, 0.0)
    hx = SVector(0.47, 0.23)

    Xgrid = MP.GridFree(x0, hx)
    Xmap = MP.ExplicitGridMapping(Xgrid)
    MP.cover!(Xmap, LazySets.Hyperrectangle(; low = lbX, high = ubX), MP.OUTER)

    # ----------------------------
    # Build finite U mapping
    # ----------------------------
    lbU = SVector(-4.0)
    ubU = SVector(4.0)
    u0 = SVector(0.0)
    hu = SVector(0.5)

    Ugrid = MP.GridFree(u0, hu)
    Umap = MP.ExplicitGridMapping(Ugrid)
    MP.cover!(Umap, LazySets.Hyperrectangle(; low = lbU, high = ubU), MP.OUTER)

    # ----------------------------
    # Concrete system + abstraction
    # ----------------------------
    tstep = 0.2
    F_sys(x, u) = SVector(u[1], -x[2] + u[1])

    jacobian_bound(u) = SMatrix{2, 2}(0.0, 0.0, 0.0, -1.0)

    concrete_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        F_sys,
        2,
        1,
        nothing,
        nothing,
    )

    continuous_approx =
        ST.ContinuousTimeGrowthBound(concrete_system; jacobian_bound = jacobian_bound)
    discrete_approx = ST.discretize(continuous_approx, tstep)

    symmodel = SY.SymbolicModelList(Xmap, Umap)
    SY.compute_abstract_system_from_concrete_system!(symmodel, discrete_approx)

    @test SY.ntransitions(symmodel.autom) == 60787

    # ----------------------------
    # Initial set -> initlist (states)
    # ----------------------------
    init_rect =
        LazySets.Hyperrectangle(; low = SVector(-3.0, -3.0), high = SVector(-2.9, -2.9))
    initlist = collect(MP.get_states_from_set(Xmap, init_rect, MP.OUTER))

    @test !isempty(initlist)

    # ----------------------------
    # Safe set: Xfull minus obstacle rectangle -> safelist
    # ----------------------------
    all_states = collect(MP.enum_states(Xmap))

    # The historical version had reversed x-bounds (an accidentally empty
    # obstacle under the old sentinel semantics); boxes now reject crossed
    # bounds, so the obstacle is a real thin strip at x ≈ -1.
    obstacle_rect =
        LazySets.Hyperrectangle(; low = SVector(-1.1, -2.0), high = SVector(-1.0, 4.0))
    bad_states = Set(MP.get_states_from_set(Xmap, obstacle_rect, MP.OUTER))

    safelist = [q for q in all_states if !(q in bad_states)]
    @test !isempty(safelist)

    # ----------------------------
    # Largest invariant set in safelist
    # ----------------------------
    contr, invariant_set_symbols, invariant_set_complement_symbols =
        OPDS.compute_largest_invariant_set(symmodel.autom, safelist)

    # Verify controller keeps you inside safelist
    safe_set = Set(safelist)
    correct = true

    for source in 1:SY.get_n_state(symmodel.autom)
        if !(source in ST.domain(contr))
            continue
        end

        targetlist = Int[]
        symbols = contr.controller_map(source)

        for symbol in symbols
            SY.compute_post!(targetlist, symmodel.autom, source, symbol)
        end

        if !(all(t -> t in safe_set, targetlist))
            correct = false
            break
        end
    end

    @test correct

    # Pick one initial state and recover a concrete coordinate
    q0 = first(initlist)
    x0_concrete = MP.get_coord_by_state(Xmap, q0)
    @test x0_concrete isa SVector{2, Float64}
end

end # module TestMain
