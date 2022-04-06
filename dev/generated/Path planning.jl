using StaticArrays

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const CO = DI.Control
const SY = DI.Symbolic

function F_sys(x, u)
    α = atan(tan(u[2])/2)
    return SVector{3}(
        u[1]*cos(α + x[3])/cos(α),
        u[1]*sin(α + x[3])/cos(α),
        u[1]*tan(u[2]))
end;

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

_X_ = UT.HyperRectangle(SVector(0.0, 0.0, -pi - 0.4), SVector(4.0, 10.0, pi + 0.4));

X1_lb = [1.0, 2.2, 2.2, 3.4, 4.6, 5.8, 5.8, 7.0, 8.2, 8.4, 9.3, 8.4, 9.3, 8.4, 9.3];
X1_ub = [1.2, 2.4, 2.4, 3.6, 4.8, 6.0, 6.0, 7.2, 8.4, 9.3, 10.0, 9.3, 10.0, 9.3, 10.0];
X2_lb = [0.0, 0.0, 6.0, 0.0, 1.0, 0.0, 7.0, 1.0, 0.0, 8.2, 7.0, 5.8, 4.6, 3.4, 2.2];
X2_ub = [9.0, 5.0, 10.0, 9.0, 10.0, 6.0, 10.0, 10.0, 8.5, 8.6, 7.4, 6.2, 5.0, 3.8, 2.6];

_U_ = UT.HyperRectangle(SVector(-1.0, -1.0), SVector(1.0, 1.0));

_I_ = UT.HyperRectangle(SVector(0.4, 0.4, 0.0), SVector(0.4, 0.4, 0.0));

_T_ = UT.HyperRectangle(SVector(3.0, 0.5, -100.0), SVector(3.6, 0.8, 100.0));

x0 = SVector(0.0, 0.0, 0.0);
h = SVector(0.2, 0.2, 0.2);
Xgrid = DO.GridFree(x0, h);

Xfull = DO.DomainList(Xgrid);
DO.add_set!(Xfull, _X_, DO.OUTER)
for (x1lb, x2lb, x1ub, x2ub) in zip(X1_lb, X2_lb, X1_ub, X2_ub)
    box = UT.HyperRectangle(SVector(x1lb, x2lb, _X_.lb[3]), SVector(x1ub, x2ub, _X_.ub[3]))
    if box ⊆ _X_ && isempty(box ∩ _I_) && isempty(box ∩ _T_)
        DO.remove_set!(Xfull, box, DO.OUTER)
    end
end

u0 = SVector(0.0, 0.0);
h = SVector(0.3, 0.3);
Ugrid = DO.GridFree(u0, h);

Ufull = DO.DomainList(Ugrid);
DO.add_set!(Ufull, _U_, DO.OUTER)

symmodel = SY.NewSymbolicModelListList(Xfull, Ufull);
@time SY.compute_symmodel_from_controlsystem!(symmodel, contsys)

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

