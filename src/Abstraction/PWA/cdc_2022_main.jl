using Dionysos
using Polyhedra
using MathematicalSystems, HybridSystems
using CDDLib
using SemialgebraicSets
using StaticArrays
using LinearAlgebra

const AB = Dionysos.Abstraction

lib = CDDLib.Library()
eye(n) = diagm(ones(n))
# Define system
N_region = 3
n_sys = 2 
n_u = 1; 
Uaux = diagm(1:n_u)
Usz = 0.5
U = [(Uaux.==i)./Usz for i in 1:n_u];

repX1 = intersect(HalfSpace(SVector{2}([1, 0]), -1))
pX1 = polyhedron(repX1, lib);
repX2 = HalfSpace(SVector{2}([-1, 0]), 1) ∩ HalfSpace(SVector{2}([1, 0]), 1)
pX2 = polyhedron(repX2, lib);
repX3 = intersect(HalfSpace(SVector{2}([-1, 0]), -1))
pX3 = polyhedron(repX3, lib);
#repU = HalfSpace(SVector{1}([-1.0]), 1/U[1]) ∩ HalfSpace(SVector{1}([1.0]), 1/U[1]);
if n_u>1
      repU = intersect([HalfSpace(SVector{n_u}(-(1:n_u .==i)  ), Usz) ∩ HalfSpace(SVector{n_u}((1:n_u .==i)), Usz) for i in 1:n_u]...);
else
      repU = HalfSpace(SVector{n_u}(-[1.0]), Usz) ∩ HalfSpace(SVector{n_u}([1.0]), Usz) 
end
pU = polyhedron(repU, lib);
pX = [pX1 pX2 pX3];



W = [-1 -1  1 1;
     -1  1 -1 1]*0.001;

PWAdomains = [pX1, pX2, pX3];

A = Vector{SMatrix{n_sys,n_sys,Float64}}(undef, N_region)
B = Vector{SMatrix{(n_sys,1),Float64}}(undef, N_region)
H = Vector{SMatrix{(n_sys,2),Float64}}(undef, N_region)
g = Vector{SVector{n_sys,Float64}}(undef, N_region)
A[1] = SMatrix{2,2}(eye(n_sys)+[0.1 20 ;
    -20 0.1]*0.01);
A[2] = transpose(A[1]);
A[3] = A[1];
#    B = fill(SMatrix{2,1}([0.3; 0.4]*10), N_region)
#B = fill(SMatrix{2,2}(eye(n_u)), N_region)
B = fill(SMatrix{n_sys,n_u}([1.0 ,1.0]), N_region)
H = fill(SMatrix{n_sys,n_sys}([0.1 0.0; 0.0 0.1]), N_region)
g[1] = -SMatrix{2,1}([1.0; 1.0])*0.1;
g[2] = g[1]*0;
g[3] = -g[1];

a = LightAutomaton(N_region)
add_transition!(a, 1, 2, 2); 
add_transition!(a, 2, 1, 1);
add_transition!(a, 1, 1, 1);
add_transition!(a, 3, 2, 2);
add_transition!(a, 2, 2, 2);
add_transition!(a, 2, 3, 3);
add_transition!(a, 3, 3, 3);

systems = [ConstrainedAffineControlDiscreteSystem(A[i], B[i], g[i], pX[i], pU) for i in 1:N_region];



switching = AutonomousSwitching()
switchings = fill(switching, 1)


resetmaps = []


system = HybridSystem(a, systems, resetmaps, switchings)
# Create Abstraction



rectX = AB.HyperRectangle(SVector(-10.0, -10.0), SVector(10.0, 10.0));
rectU = rectX


step = 4
X_origin = SVector(0.0, 0.0);
X_step = SVector(1.0/step, 1.0/step);
Xgrid = AB.GridRectangular(X_origin, X_step, rectX);

domainX = AB.DomainList(Xgrid);
AB.add_set!(domainX, rectX, AB.OUTER)


domainU = domainX;
AB.add_set!(domainU, rectU, AB.INNER)


symmodel = AB.NewSymbolicModelListList(domainX, domainU);



L = [eye(n_sys+n_u) zeros(n_sys+n_u,1)];

transitionCost = Dict()
transitionKappa = Dict()
AB.compute_symmodel_from_hybridcontrolsystem!(symmodel,transitionCost, transitionKappa, system, W, L, U)

using PyPlot

PyPlot.pygui(true) #jl
fig = PyPlot.figure()

ax = PyPlot.axes(aspect = "equal")
ax.set_xlim(rectX.lb[1]-0.2,rectX.ub[1]+0.2)
ax.set_ylim(rectX.lb[2]-0.2,rectX.ub[2]+0.2)

vars = [1, 2];


include(dirname(pathof(Dionysos)) * "/plotting.jl")

Plot.domain!(ax, vars, domainX, fc = "white")


d=0.01;
for t in symmodel.autom.transitions.data
      dist = randn(1)*d
      arrow_x, arrow_y = AB.get_coord_by_pos(Xgrid,AB.get_xpos_by_state(symmodel,t[1])).+dist
      aux = (AB.get_coord_by_pos(Xgrid,AB.get_upos_by_symbol(symmodel,t[2]))-[arrow_x, arrow_y]).+dist
      arrow_dx, arrow_dy = aux*(norm(aux)-0.15)/norm(aux)
      if t[1]==t[2]
            PyPlot.scatter(arrow_x, arrow_y, ec="red")
      else
            PyPlot.arrow(arrow_x, arrow_y, arrow_dx, arrow_dy,  ec = "red",head_width=.1)

      end
      

end
# Define Specification

# Synthesis of the discrete controler

# Simulation


