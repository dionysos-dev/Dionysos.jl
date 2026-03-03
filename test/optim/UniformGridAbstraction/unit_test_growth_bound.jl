module TestMain

using Test
using StaticArrays
using MathematicalSystems
using Dionysos

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic

sleep(0.1)
println("Started test")

@testset "FromControlSystem" begin
    # ----------------------------
    # X mapping
    # ----------------------------
    lbX = SVector(0.0, 0.0)
    ubX = SVector(10.0, 11.0)
    x0  = SVector(0.0, 0.0)
    hx  = SVector(1.0, 2.0)

    Xgrid = MP.GridFree(x0, hx)
    Xmap  = MP.ExplicitGridMapping(Xgrid)
    MP.add_set!(Xmap, UT.HyperRectangle(lbX, ubX), MP.OUTER)

    # ----------------------------
    # U mapping
    # ----------------------------
    lbU = SVector(-1.0)
    ubU = SVector( 1.0)
    u0  = SVector(0.0)
    hu  = SVector(0.5)

    Ugrid = MP.GridFree(u0, hu)
    Umap  = MP.ExplicitGridMapping(Ugrid)
    MP.add_set!(Umap, UT.HyperRectangle(lbU, ubU), MP.OUTER)

    # ----------------------------
    # Concrete system + abstraction
    # ----------------------------
    tstep = 5.0

    F_sys(x, u) = SVector(u[1], -cos(x[1]))
    jacobian_bound(u) = SMatrix{2,2}(0.0, 1.0, 0.0, 0.0)

    concrete_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        F_sys,
        2,
        1,
        nothing,
        nothing,
    )

    continuous_approx =
        ST.ContinuousTimeGrowthBound_from_jacobian_bound(concrete_system, jacobian_bound)
    discrete_approx = ST.discretize(continuous_approx, tstep)

    symmodel = SY.SymbolicModelList(Xmap, Umap; Xset=MP.MappingSet{2}(), Rset=MP.MappingSet{2}())
    SY.compute_abstract_system_from_concrete_system!(symmodel, discrete_approx)

    @test SY.ntransitions(symmodel.autom) == 1355

    # ----------------------------
    # One post computation sanity check (robust)
    # ----------------------------
    symbol = first(SY.enum_inputs(symmodel))

    source_found = nothing
    targetlist = Int[]

    for q in SY.enum_states(symmodel)
        empty!(targetlist)
        SY.compute_post!(targetlist, symmodel.autom, q, symbol)
        if !isempty(targetlist)
            source_found = q
            break
        end
    end

    @test source_found !== nothing   # there exists at least one enabled transition for this symbol

    Ysimple = MP.ExplicitGridMapping(Xgrid)
    for target in targetlist
        MP.add_pos!(Ysimple, MP.get_pos_by_state(Xmap, target))
    end

    @test !isempty(targetlist)
end

sleep(0.1)
println("End test")

end # module TestMain