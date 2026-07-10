module TestMain

using Test
using StaticArrays
import LazySets
using MathematicalSystems
using Dionysos

const DI = Dionysos
const UT = DI.Utils
const MP = DI.Mapping
const ST = DI.System
const SY = DI.Symbolic
const OPDS = DI.Optim.DiscreteSystems

println("Started test")

@testset "Generic LazySets through the grid mapping" begin
    h = SVector(0.25, 0.25)
    grid = MP.GridFree(h; zero_origin = true)

    ball = LazySets.Ball2([0.0, 0.0], 1.0)
    outer = collect(MP.get_pos_from_set(grid, ball, MP.OUTER))
    inner = collect(MP.get_pos_from_set(grid, ball, MP.INNER))
    center = collect(MP.get_pos_from_set(grid, ball, MP.CENTER))

    @test !isempty(inner)
    @test issubset(Set(inner), Set(center))
    @test issubset(Set(center), Set(outer))
    @test (0, 0) in inner

    # INNER cells are fully inside the ball (corner certification is exact for
    # convex sets); CENTER cells have their center inside.
    for pos in inner
        @test all(v -> v ∈ ball, LazySets.vertices_list(MP.get_rec(grid, pos)))
    end
    for pos in center
        @test MP.get_coord_by_pos(grid, pos) ∈ ball
    end

    zono = LazySets.Zonotope([1.0, 1.0], [0.3 0.1; 0.0 0.2])
    zouter = collect(MP.get_pos_from_set(grid, zono, MP.OUTER))
    @test !isempty(zouter)
    @test MP.get_pos_by_coord(grid, SVector(1.0, 1.0)) in zouter
end

@testset "Reach-avoid with zonotope obstacle and ball target" begin
    Xgrid = MP.GridFree(SVector(0.0, 0.0), SVector(0.47, 0.23))
    Xmap_full = MP.ExplicitGridMapping(Xgrid)
    MP.add_set!(Xmap_full, UT.box(SVector(-5.0, -5.0), SVector(5.0, 5.0)), MP.OUTER)

    # thin vertical zonotope strip at x ≈ -1
    obstacle = LazySets.Zonotope([-1.0, 1.0], [0.1 0.0; 0.0 3.0])
    bad = Set(MP.get_states_from_set(Xmap_full, obstacle, MP.OUTER))
    @test !isempty(bad)

    positions_ok = [
        MP.get_pos_by_state(Xmap_full, q) for q in MP.enum_states(Xmap_full) if !(q in bad)
    ]
    Xmap = MP.ExplicitGridMapping{2, Float64}(Xgrid, positions_ok)

    Ugrid = MP.GridFree(SVector(0.0), SVector(1.0))
    Umap = MP.ExplicitGridMapping(Ugrid)
    MP.add_set!(Umap, UT.box(SVector(-2.0), SVector(2.0)), MP.OUTER)

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
        ST.ContinuousTimeGrowthBound_from_jacobian_bound(concrete_system, jacobian_bound)
    discrete_approx = ST.discretize(continuous_approx, 1.0)

    symmodel = SY.SymbolicModelList(Xmap, Umap)
    SY.compute_abstract_system_from_concrete_system!(symmodel, discrete_approx)

    init_rect = UT.box(SVector(-3.0, -3.0), SVector(-2.9, -2.9))
    initlist = collect(MP.get_states_from_set(Xmap, init_rect, MP.OUTER))
    @test !isempty(initlist)

    target = LazySets.Ball2([2.0, 2.0], 1.5)
    targetlist = collect(MP.get_states_from_set(Xmap, target, MP.INNER))
    @test !isempty(targetlist)

    contr, controllable_set, _, value_fun_tab = OPDS.compute_worst_case_cost_controller(
        symmodel.autom,
        targetlist;
        initial_set = initlist,
    )
    @test !isempty(controllable_set)
    @test value_fun_tab !== nothing
end

println("End test")

end # module TestMain
