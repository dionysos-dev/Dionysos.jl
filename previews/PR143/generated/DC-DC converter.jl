using StaticArrays

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const CO = DI.Control
const SY = DI.Symbolic

vs = 1.0; rL = 0.05; xL = 3.0; rC = 0.005; xC = 70.0; r0 = 1.0;

b = SVector(vs/xL, 0.0);
A1 = SMatrix{2,2}(-rL/xL, 0.0, 0.0, -1.0/xC/(r0+rC));
A2 = SMatrix{2,2}(-(rL+r0*rC/(r0+rC))/xL, 5.0*r0/(r0+rC)/xC,
    -r0/(r0+rC)/xL/5.0, -1.0/xC/(r0+rC));
F_sys = let b = b, A1 = A1, A2 = A2
    (x, u) -> u[1] == 1 ? A1*x + b : A2*x + b
end;

ngrowthbound = 5;
A2_abs = SMatrix{2,2}(-(rL+r0*rC/(r0+rC))/xL, 5.0*r0/(r0+rC)/xC,
                      r0/(r0+rC)/xL/5.0, -1.0/xC/(r0+rC));
L_growthbound = let A1 = A1, A2_abs = A2_abs
    u -> u[1] == 1 ? A1 : A2_abs
end;

sysnoise = SVector(0.0, 0.0);
measnoise = SVector(0.0, 0.0);

tstep = 0.5;
nsys = 5;

contsys = ST.NewControlSystemGrowthRK4(tstep, F_sys, L_growthbound, sysnoise,
                                       measnoise, nsys, ngrowthbound);

_X_ = UT.HyperRectangle(SVector(1.15, 5.45), SVector(1.55, 5.85));

_U_ = UT.HyperRectangle(SVector(1), SVector(2));

x0 = SVector(0.0, 0.0);
h = SVector(2.0/4.0e3, 2.0/4.0e3);
Xgrid = DO.GridFree(x0, h);

Xfull = DO.DomainList(Xgrid);
DO.add_set!(Xfull, _X_, DO.INNER)

u0 = SVector(1);
h = SVector(1);
Ugrid = DO.GridFree(u0, h);

Ufull = DO.DomainList(Ugrid);
DO.add_set!(Ufull, _U_, DO.OUTER);

symmodel = SY.NewSymbolicModelListList(Xfull, Ufull);
@time SY.compute_symmodel_from_controlsystem!(symmodel, contsys)

Xinit = DO.DomainList(Xgrid);
union!(Xinit, Xfull)
initlist = [SY.get_state_by_xpos(symmodel, pos) for pos in DO.enum_pos(Xinit)];

Xsafe = DO.DomainList(Xgrid)
union!(Xsafe, Xfull)
safelist = [SY.get_state_by_xpos(symmodel, pos) for pos in DO.enum_pos(Xsafe)];

contr = CO.NewControllerList();
@time CO.compute_controller_safe!(contr, symmodel.autom, initlist, safelist)

nstep = 300;
x0 = SVector(1.2, 5.6);

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

