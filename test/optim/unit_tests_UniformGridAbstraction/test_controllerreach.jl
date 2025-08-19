module Test

# using TestAbstraction
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

@testset "ControllerReach" begin
    lb = SVector(-5.0, -5.0)
    ub = SVector(5.0, 5.0)
    x0 = SVector(0.0, 0.0)
    h = SVector(0.47, 0.23)
    Xgrid = DO.GridFree(x0, h)
    Xfull = DO.DomainList(Xgrid)
    DO.add_set!(Xfull, UT.HyperRectangle(lb, ub), DO.OUTER)
    DO.remove_set!(
        Xfull,
        UT.HyperRectangle(SVector(-1.0, -2.0), SVector(-1.1, 4.0)),
        DO.OUTER,
    )

    lb = SVector(-2.0)
    ub = SVector(2.0)
    u0 = SVector(0.0)
    h = SVector(1.0)
    Ugrid = DO.GridFree(u0, h)
    Ufull = DO.DomainList(Ugrid)
    DO.add_set!(Ufull, UT.HyperRectangle(lb, ub), DO.OUTER)

    tstep = 1.0
    nsys = 3
    ngrowthbound = 3
    F_sys(x, u) = SVector(1.0, u[1])
    sysnoise = SVector(1.0, 1.0) * 0.001
    measnoise = SVector(1.0, 1.0) * 0.001
    jacobian_bound(u) = SMatrix{2, 2}(0.0, 0.0, 0.0, 0.0)

    concrete_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        F_sys,
        1,
        1,
        nothing,
        nothing,
    )
    continuous_approx =
        ST.ContinuousTimeGrowthBound_from_jacobian_bound(concrete_system, jacobian_bound)
    discrete_approx = ST.discretize(continuous_approx, tstep)

    symmodel = SY.NewSymbolicModelList(Xfull, Ufull)
    SY.compute_abstract_system_from_concrete_system!(symmodel, discrete_approx)

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
    Xtarget = DO.DomainList(Xgrid)
    DO.add_set!(Xtarget, UT.HyperRectangle(SVector(0.0, 0.0), SVector(4.0, 4.0)), DO.OUTER)
    targetlist = Int[]
    for pos in DO.enum_pos(Xtarget)
        push!(targetlist, SY.get_state_by_xpos(symmodel, pos))
    end

    contr, controllable_set, uncontrollable_set, value_fun_tab =
        SY.compute_worst_case_cost_controller(
            symmodel.autom,
            targetlist;
            initial_set = initlist,
        )

    @test length(ST.domain(contr)) == 412
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
