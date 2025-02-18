module TestMain

using Test
using StaticArrays, MathematicalSystems
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

sleep(0.1) # used for good printing
println("Started test")

@testset "ControllerSafe" begin
    lb = SVector(-5.0, -5.0)
    ub = SVector(5.0, 5.0)
    x0 = SVector(0.0, 0.0)
    h = SVector(0.47, 0.23)
    Xgrid = DO.GridFree(x0, h)
    Xfull = DO.DomainList(Xgrid)
    DO.add_set!(Xfull, UT.HyperRectangle(lb, ub), DO.OUTER)

    lb = SVector(-4.0)
    ub = SVector(4.0)
    u0 = SVector(0.0)
    h = SVector(0.5)
    Ugrid = DO.GridFree(u0, h)
    Ufull = DO.DomainList(Ugrid)
    DO.add_set!(Ufull, UT.HyperRectangle(lb, ub), DO.OUTER)

    tstep = 0.2
    F_sys(x, u) = SVector(u[1], -x[2] + u[1])
    sysnoise = SVector(1.0, 1.0) * 0.001
    measnoise = SVector(1.0, 1.0) * 0.001
    jacobian_bound(u) = SMatrix{2, 2}(0.0, 0.0, 0.0, -1.0)

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

    symmodel = SY.NewSymbolicModelListList(Xfull, Ufull)
    SY.compute_abstract_system_from_concrete_system!(symmodel, discrete_approx)

    @test SY.ntransitions(symmodel.autom) == 60787

    Xinit = DO.DomainList(Xgrid)
    DO.add_set!(
        Xinit,
        UT.HyperRectangle(SVector(-3.0, -3.0), SVector(-2.9, -2.9)),
        DO.OUTER,
    )
    initlist = Int[]
    for pos in DO.enum_pos(Xinit)
        push!(initlist, SY.get_state_by_xpos(symmodel, pos))
    end
    Xsafe = DO.DomainList(Xgrid)
    union!(Xsafe, Xfull)
    DO.remove_set!(
        Xsafe,
        UT.HyperRectangle(SVector(-1.0, -2.0), SVector(-1.1, 4.0)),
        DO.OUTER,
    )
    safelist = Int[]
    for pos in DO.enum_pos(Xsafe)
        push!(safelist, SY.get_state_by_xpos(symmodel, pos))
    end

    contr, invariant_set_symbols, uninvariant_set_symbols =
        AB.UniformGridAbstraction.compute_largest_invariant_set(symmodel, safelist)
    @test length(contr) == 15045

    invlist = Int[]
    for source in 1:(symmodel.autom.nstates)
        if !isempty(UT.fix_and_eliminate_first(contr, source))
            push!(invlist, source)
        end
    end
    Xinv = DO.DomainList(Xgrid)
    Yinv = DO.DomainList(Xgrid)
    correct = true
    for source in invlist
        DO.add_pos!(Xinv, SY.get_xpos_by_state(symmodel, source))
        if !correct
            break
        end
        targetlist = Int[]
        for symbol in UT.fix_and_eliminate_first(contr, source)
            SY.compute_post!(targetlist, symmodel.autom, source, symbol)
        end
        for target in targetlist
            DO.add_pos!(Yinv, SY.get_xpos_by_state(symmodel, target))
        end
        correct = correct && targetlist ⊆ safelist
    end
    @test correct

    xpos = DO.get_somepos(Xinit)
    x0 = DO.get_coord_by_pos(Xgrid, xpos)
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
