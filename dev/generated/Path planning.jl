using StaticArrays

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const CO = DI.Control
const SY = DI.Symbolic

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "PathPlanning.jl"))

problem = PathPlanning.problem();

F_sys = problem.system.f;
_X_ = problem.system.X;
_U_ = problem.system.U;

ngrowthbound = 5;
function L_growthbound(u)
    β = abs(u[1]/cos(atan(tan(u[2])/2)))
    return SMatrix{3,3}(
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        β, β, 0.0)
end;

sysnoise = SVector(0.0, 0.0, 0.0);
measnoise = SVector(0.0, 0.0, 0.0);

tstep = 0.3;
nsys = 5;

contsys = ST.NewControlSystemGrowthRK4(tstep, F_sys, L_growthbound, sysnoise,
                                       measnoise, nsys, ngrowthbound);

x0 = SVector(0.0, 0.0, 0.0);
h = SVector(0.2, 0.2, 0.2);
Xgrid = DO.GridFree(x0, h);

Xfull = DO.DomainList(Xgrid);
DO.add_set!(Xfull, _X_, DO.OUTER);

u0 = SVector(0.0, 0.0);
h = SVector(0.3, 0.3);
Ugrid = DO.GridFree(u0, h);

Ufull = DO.DomainList(Ugrid);
DO.add_set!(Ufull, _U_, DO.OUTER);

symmodel = SY.NewSymbolicModelListList(Xfull, Ufull);
@time SY.compute_symmodel_from_controlsystem!(symmodel, contsys)

_I_ = problem.initial_set;
_T_ = problem.target_set;

Xinit = DO.DomainList(Xgrid);
DO.add_subset!(Xinit, Xfull, _I_, DO.OUTER)
initlist = [SY.get_state_by_xpos(symmodel, pos) for pos in DO.enum_pos(Xinit)];

Xtarget = DO.DomainList(Xgrid)
DO.add_subset!(Xtarget, Xfull, _T_, DO.OUTER)
targetlist = [SY.get_state_by_xpos(symmodel, pos) for pos in DO.enum_pos(Xtarget)];

contr = CO.NewControllerList();
@time CO.compute_controller_reach!(contr, symmodel.autom, initlist, targetlist)

nstep = 100;
x0 = SVector(0.4, 0.4, 0.0);

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

