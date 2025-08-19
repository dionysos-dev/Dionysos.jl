using Dionysos
using StaticArrays, MathematicalSystems
using LinearAlgebra
using Plots

const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic

rectX = UT.HyperRectangle(SVector(-2, -2), SVector(2, 2));
rectU = UT.HyperRectangle(SVector(-5), SVector(5));

x0 = SVector(0.0, 0.0);
h = SVector(1.0 / 5, 1.0 / 5);
Xgrid = DO.GridFree(x0, h);

domainX = DO.DomainList(Xgrid);
DO.add_set!(domainX, rectX, DO.INNER)

u0 = SVector(0.0);
h = SVector(1.0 / 5);
Ugrid = DO.GridFree(u0, h);
domainU = DO.DomainList(Ugrid);
DO.add_set!(domainU, rectU, DO.INNER);

tstep = 0.1;
nsys = 10; # Runge-Kutta pre-scaling

A = SMatrix{2, 2}(0.0, 1.0, -3.0, 1.0);
B = SMatrix{2, 1}(0.0, 1.0);

F_sys = let A = A
    (x, u) -> A * x + B * u
end;

ngrowthbound = 10; # Runge-Kutta pre-scaling
A_diag = diagm(diag(A));
A_abs = abs.(A) - abs.(A_diag) + A_diag
jacobian_bound = x -> abs.(A)

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

symmodel = SY.NewSymbolicModelList(domainX, domainU)
SY.compute_abstract_system_from_concrete_system!(symmodel, discrete_approx)

xpos = DO.get_pos_by_coord(Xgrid, SVector(1.1, 1.3))

x = DO.get_coord_by_pos(Xgrid, xpos)
abstract_input = 1
u = SY.get_concrete_input(symmodel, abstract_input)

post = Int[]
SY.compute_post!(post, symmodel.autom, symmodel.xpos2int[xpos], abstract_input)

domainPostx = DO.DomainList(Xgrid);
for pos in symmodel.xint2pos[post]
    DO.add_pos!(domainPostx, pos)
end

fig = plot(;
    aspect_ratio = :equal,
    xtickfontsize = 10,
    ytickfontsize = 10,
    guidefontsize = 16,
);
xlims!(-2, 2)
ylims!(-2, 2)
dims = [1, 2]

plot!(domainX; fc = "white", dims = dims);
domainx = DO.DomainList(Xgrid);
DO.add_pos!(domainx, xpos)
plot!(domainx; fc = "blue", dims = dims);
plot!(domainPostx; fc = "green", dims = dims)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
