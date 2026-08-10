module TestMain

import Dionysos
import LazySets
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

@testset "ControllerReach" begin
    # ----------------------------
    # X mapping = full grid box minus obstacle
    # ----------------------------
    lbX = SVector(-5.0, -5.0)
    ubX = SVector(5.0, 5.0)
    x0 = SVector(0.0, 0.0)
    hx = SVector(0.47, 0.23)

    Xgrid = MP.GridFree(x0, hx)
    Xmap_full = MP.ExplicitGridMapping(Xgrid)
    MP.cover!(Xmap_full, LazySets.Hyperrectangle(; low = lbX, high = ubX), MP.OUTER)

    # The historical version had reversed x-bounds ([-1.0, -1.1]), which the old
    # sentinel semantics silently treated as an empty obstacle; boxes now reject
    # crossed bounds, so the obstacle is a real thin strip at x ≈ -1.
    obstacle =
        LazySets.Hyperrectangle(; low = SVector(-1.1, -2.0), high = SVector(-1.0, 4.0))
    bad = Set(MP.get_states_from_set(Xmap_full, obstacle, MP.OUTER))

    # Build a filtered explicit mapping containing only "safe" positions
    positions_ok = NTuple{2, Int}[]
    for q in MP.enum_states(Xmap_full)
        if q in bad
            continue
        end
        push!(positions_ok, MP.get_pos_by_state(Xmap_full, q))
    end

    Xmap = MP.ExplicitGridMapping{2, Float64}(Xgrid, positions_ok)

    # ----------------------------
    # U mapping
    # ----------------------------
    lbU = SVector(-2.0)
    ubU = SVector(2.0)
    u0 = SVector(0.0)
    hu = SVector(1.0)

    Ugrid = MP.GridFree(u0, hu)
    Umap = MP.ExplicitGridMapping(Ugrid)
    MP.cover!(Umap, LazySets.Hyperrectangle(; low = lbU, high = ubU), MP.OUTER)

    # ----------------------------
    # Concrete system + abstraction
    # ----------------------------
    tstep = 1.0

    F_sys(x, u) = SVector(1.0, u[1])
    jacobian_bound(u) = SMatrix{2, 2}(0.0, 0.0, 0.0, 0.0)

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

    # ----------------------------
    # Initial set -> initlist (states)
    # ----------------------------
    init_rect =
        LazySets.Hyperrectangle(; low = SVector(-3.0, -3.0), high = SVector(-2.9, -2.9))
    initlist = collect(MP.get_states_from_set(Xmap, init_rect, MP.OUTER))
    @test !isempty(initlist)

    # ----------------------------
    # Target set -> targetlist (states)
    # ----------------------------
    target_rect =
        LazySets.Hyperrectangle(; low = SVector(0.0, 0.0), high = SVector(4.0, 4.0))
    targetlist = collect(MP.get_states_from_set(Xmap, target_rect, MP.OUTER))
    @test !isempty(targetlist)

    contr, controllable_set, uncontrollable_set, value_fun_tab =
        OPDS.compute_worst_case_cost_controller(
            symmodel.autom,
            targetlist;
            initial_set = initlist,
        )

    @test !isempty(controllable_set)
    @test value_fun_tab !== nothing
end

end # module TestMain
