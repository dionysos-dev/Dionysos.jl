using Dionysos
using Polyhedra
using MathematicalSystems, HybridSystems
using CDDLib
using SemialgebraicSets
using StaticArrays

lib = CDDLib.Library()
# Define system
N_region = 3
n_sys = 2 
X1 = intersect(HalfSpace(SVector{2}([1, 0]), -1))
pX1 = polyhedron(X1, lib);
X2 = HalfSpace(SVector{2}([-1, 0]), 1) ∩ HalfSpace(SVector{2}([1, 0]), 1)
pX2 = polyhedron(X2, lib);
X3 = intersect(HalfSpace(SVector{2}([-1, 0]), -1))
pX3 = polyhedron(X3, lib);
U = HalfSpace(SVector{3}([-1.0,0,0]), 10) ∩ HalfSpace(SVector{3}([1.0,0,0]), 10);
pU = polyhedron(U, lib);
pX = [pX1 pX2 pX3];

W = [ [1.0; 0.0];
      [-1.0; 0.0];
      [0.0; 1.0];
      [0.0; -1.0] ];

PWAdomains = [pX1, pX2, pX3];

A = Vector{SMatrix{n_sys,n_sys,Float64}}(undef, N_region)
B = Vector{SMatrix{(n_sys,1),Float64}}(undef, N_region)
H = Vector{SMatrix{(n_sys,2),Float64}}(undef, N_region)
g = Vector{SVector{n_sys,Float64}}(undef, N_region)
A[1] = SMatrix{2,2}([0.2 0.6 ;
       -0.6 0.2]);
A[2] = -A[1];
A[3] = 0.5*A[1];
B = fill(SMatrix{2,1}([0.1; 0.4]), N_region)
H = fill(SMatrix{2,2}([0.1 0.0; 0.0 0.1]), N_region)
g[1] = [0.0; 1.0];
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

systems = [ConstrainedAffineControlDiscreteSystem(A[i], [B[i] H[i]], g[i], pX[i], pU) for i in 1:N_region];



switching = AutonomousSwitching()
switchings = fill(switching, 1)


resetmaps = []


system = HybridSystem(a, systems, resetmaps, switchings)
# Create Abstraction
const AB = Dionysos.Abstraction;



rectX = AB.HyperRectangle(SVector(-2.0, -2), SVector(2.0, 2));
rectU = AB.HyperRectangle(SVector(-2.0, -2, -2), SVector(2.0, 2, 2));

X_origin = SVector(0.0, 0.0);
X_step = SVector(1.0/5, 1.0/5);
Xgrid = AB.GridRectangular(X_origin, X_step, rectX);

domainX = AB.DomainList(Xgrid);
AB.add_set!(domainX, rectX, AB.INNER)


U_origin = SVector(0.0, 0.0, 0.0);
U_step = SVector(1.0/5, 1.0/5, 10/5);
Ugrid = AB.GridFree(U_origin, U_step);

domainU = AB.DomainList(Ugrid);
AB.add_set!(domainU, rectU, AB.INNER)


symmodel = AB.NewSymbolicModelListList(domainX, domainU);


AB.compute_symmodel_from_hybridcontrolsystem!(symmodel, system, W)



# Define Specification

# Synthesis of the discrete controler

# Simulation


