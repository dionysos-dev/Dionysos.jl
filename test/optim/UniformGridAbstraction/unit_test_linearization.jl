module TestMain

using Test
using StaticArrays, MathematicalSystems
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic

sleep(0.1)
println("Started test")

@testset "FromControlSystem (Mapping-based)" begin
    # ----------------------------
    # Build XMapping: grid + implicit rectangular universe
    # ----------------------------
    lb = SVector(0.0, 0.0)
    ub = SVector(10.0, 11.0)

    x0 = SVector(0.0, 0.0)
    h = SVector(1.0, 2.0)

    Xgrid = MP.GridFree(x0, h)
    Xrect = UT.HyperRectangle(lb, ub)

    # Implicit universe = all grid cells covering Xrect (OUTER)
    Xmap = MP.ImplicitGridMapping(Xgrid, Xrect; incl_mode = MP.OUTER)

    # ----------------------------
    # Build UMapping: finite list of inputs
    # ----------------------------
    # Same as before: [-1, -0.5, 0, 0.5, 1]
    Uvals = [SVector(u) for u in (-1.0:0.5:1.0)]
    Umap = MP.ListMapping(Uvals)

    # ----------------------------
    # Build symbolic model with default domains = "all states of mapping"
    # (MappingSet = read-only set containing all mapping states)
    # ----------------------------
    symmodel = SY.SymbolicModelList(
        Xmap,
        Umap;
        Xset = MP.MappingSet{2}(),
        Rset = nothing,              # defaults to Xset
        Uset = MP.MappingSet{1}(),
        convert_U_to_list = false,   # Umap already is a ListMapping
    )

    # ----------------------------
    # Concrete system + approximation
    # ----------------------------
    tstep = 0.5
    F_sys(x, u) = SVector(u[1], -cos(x[1]))
    DF_sys(x, u) = SMatrix{2, 2}(0.0, sin(x[1]), 0.0, 0.0)
    bound_DF(u) = 1.0
    bound_DDF(u) = 1.0

    concrete_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        F_sys,
        2,
        1,
        nothing,
        nothing,
    )
    continuous_approx =
        ST.ContinuousTimeLinearized(concrete_system, DF_sys, bound_DF, bound_DDF)
    discrete_approx = ST.discretize(continuous_approx, tstep)

    SY.compute_abstract_system_from_concrete_system!(symmodel, discrete_approx)

    @test SY.ntransitions(symmodel.autom) == 2175

    # ----------------------------
    # Robust post test: find a state with at least one successor for some input symbol
    # ----------------------------
    inputs = collect(SY.enum_inputs(symmodel))
    @test !isempty(inputs)
    symbol = first(inputs)

    targetlist = Int[]
    source = nothing

    for q in SY.enum_states(symmodel)
        empty!(targetlist)
        SY.compute_post!(targetlist, symmodel.autom, q, symbol)
        if !isempty(targetlist)
            source = q
            break
        end
    end

    @test source !== nothing
    @test !isempty(targetlist)

    # ----------------------------
    # Sanity checks on mapping consistency
    # ----------------------------
    # concrete representative of the source state (center of cell, etc.)
    x = MP.get_coord_by_state(Xmap, source)
    u = SY.get_concrete_input(symmodel, symbol)
    @test x isa SVector{2, Float64}
    @test u isa SVector{1, Float64}

    # all targets should be valid mapping states
    for tgt in targetlist
        @test 1 <= tgt <= MP.get_n_state(Xmap)
    end
end

sleep(0.1)
println("End test")

end # module TestMain
