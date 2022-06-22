module DCDC
using Test     
# First, let us import [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl).

using StaticArrays
using MathematicalSystems

# At this point, we import the useful Dionysos sub-module for this problem: [Abstraction](@__REPO_ROOT_URL__/src/Abstraction/abstraction.jl).
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const CO = DI.Control
const SY = DI.Symbolic

function dynamicofsystem(vs = 1.0, rL = 0.05, xL = 3.0, rC = 0.005, xC = 70.0, r0 = 1.0,ngrowthbound = 5)
    # Definition of the dynamics functions $f_p$ of the system:
    b = SVector(vs/xL, 0.0);
    A1 = SMatrix{2,2}(-rL/xL, 0.0, 0.0, -1.0/xC/(r0+rC))
    A2 = SMatrix{2,2}(-(rL+r0*rC/(r0+rC))/xL, 5.0*r0/(r0+rC)/xC,
        -r0/(r0+rC)/xL/5.0, -1.0/xC/(r0+rC))
    F_sys = let b = b, A1 = A1, A2 = A2
        (x, u) -> u[1] == 1 ? A1*x + b : A2*x + b
    end
    A2_abs = SMatrix{2,2}(-(rL+r0*rC/(r0+rC))/xL, 5.0*r0/(r0+rC)/xC,
                        r0/(r0+rC)/xL/5.0, -1.0/xC/(r0+rC))
    L_growthbound = let A1 = A1, A2_abs = A2_abs
        u -> u[1] == 1 ? A1 : A2_abs
    end
    return F_sys, L_growthbound,ngrowthbound
end


function system(
    sysnoise = SVector(0.0, 0.0),
    measnoise = SVector(0.0, 0.0),
    tstep = 0.5,
    nsys = 5,
    _X_ = UT.HyperRectangle(SVector(1.15, 5.45), SVector(1.55, 5.85)),
    _U_ = UT.HyperRectangle(SVector(1), SVector(2)),
    xdim=2,
    udim=1
)
    # Definition of the dynamics functions $f_p$ of the system:
    F_sys,L_growthbound,ngrowthbound=dynamicofsystem();
    contsys = ST.NewControlSystemGrowthRK4(tstep, F_sys, L_growthbound, sysnoise,
                                        measnoise, nsys, ngrowthbound);
    sys = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(contsys, xdim, udim, _X_, _U_)
    return sys
end

# ### Definition of the control problem
function solveproblem(problem,x0 = SVector(0.0, 0.0), hx = SVector(2.0/4.0e3, 2.0/4.0e3), u0 = SVector(1), hu = SVector(1);)
    # Definition of the grid of the state-space on which the abstraction is based (origin `x0` and state-space discretization `h`):
    Xgrid = DO.GridFree(x0, hx);
    Xfull = DO.DomainList(Xgrid);
    # Construction of the struct `DomainList` containing the feasible cells of the state-space.
    # Note, we used `AB.INNER` to make sure to add cells entirely contained in the domain because we are working with a safety problem.
    DO.add_set!(Xfull, problem.X, DO.INNER);
    # Definition of the grid of the input-space on which the abstraction is based (origin `u0` and input-space discretization `h`):
    Ugrid = DO.GridFree(u0, hu);
    Ufull = DO.DomainList(Ugrid);
    # Construction of the struct `DomainList` containing the quantized inputs:
    DO.add_set!(Ufull, problem.U, DO.OUTER);
    # Construction of the abstraction:
    symmodel = SY.NewSymbolicModelListList(Xfull, Ufull);
    @time SY.compute_symmodel_from_controlsystem!(symmodel, problem.f)
    # ### Construction of the controller
    # In this problem, we consider both: the initial state-space and the safety state-space are equal to the entire state-space.
    #
    # Computation of the initial symbolic states:
    Xinit = DO.DomainList(Xgrid);
    union!(Xinit, Xfull)
    initlist = [SY.get_state_by_xpos(symmodel, pos) for pos in DO.enum_pos(Xinit)];
    # Computation of the safety symbolic states:
    Xsafe = DO.DomainList(Xgrid)
    union!(Xsafe, Xfull)
    safelist = [SY.get_state_by_xpos(symmodel, pos) for pos in DO.enum_pos(Xsafe)];
    # Construction of the controller:
    contr = CO.NewControllerList();
    @time CO.compute_controller_safe!(contr, symmodel.autom, initlist, safelist)
    print("finished")
   # return contr
end


end