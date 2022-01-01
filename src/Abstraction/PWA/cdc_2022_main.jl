using Dionysos
using Polyhedra
using MathematicalSystems, HybridSystems
using CDDLib
using SemialgebraicSets
using StaticArrays
using LinearAlgebra

const AB = Dionysos.Abstraction

lib = CDDLib.Library() #polyhedron lib
eye(n) = diagm(ones(n)) # I matrix
# Define system
N_region = 3
n_sys = 2 
n_u = 2; 
Uaux = diagm(1:n_u)
Usz = 30 # upper limit on |u|
U = [(Uaux.==i)./Usz for i in 1:n_u];

#PWA partitions
repX1 = intersect(HalfSpace(SVector{2}([1, 0]), -1))
pX1 = polyhedron(repX1, lib);
repX2 = HalfSpace(SVector{2}([-1, 0]), 1) ∩ HalfSpace(SVector{2}([1, 0]), 1)
pX2 = polyhedron(repX2, lib);
repX3 = intersect(HalfSpace(SVector{2}([-1, 0]), -1))
pX3 = polyhedron(repX3, lib);

#control input bounded region
if n_u>1 
      repU = intersect([HalfSpace(SVector{n_u}(-(1:n_u .==i)  ), Usz) ∩ HalfSpace(SVector{n_u}((1:n_u .==i)), Usz) for i in 1:n_u]...);
else
      repU = HalfSpace(SVector{n_u}(-[1.0]), Usz) ∩ HalfSpace(SVector{n_u}([1.0]), Usz) 
end
pU = polyhedron(repU, lib);
pX = [pX1 pX2 pX3];



dt = 0.01 

W = [-1 -1  1 1;
     -1  1 -1 1]*dt; # polytope of disturbances

PWAdomains = [pX1, pX2, pX3];

A = Vector{SMatrix{n_sys,n_sys,Float64}}(undef, N_region)
B = Vector{SMatrix{(n_sys,n_u),Float64}}(undef, N_region)
H = Vector{SMatrix{(n_sys,2),Float64}}(undef, N_region)
g = Vector{SVector{n_sys,Float64}}(undef, N_region)
A[1] = SMatrix{2,2}(eye(n_sys)+[1 30 ;
    -10 1]*dt);
A[2] = transpose(A[1]);
A[3] = A[1];
B = fill(SMatrix{2,2}(eye(n_u)*dt), N_region)
g[1] = -SMatrix{2,1}([10; 10])*0.01;
g[2] = g[1]*0;
g[3] = -g[1];

# automata for pwa switching between partitions (unecessary?)
a = LightAutomaton(N_region)
add_transition!(a, 1, 2, 2); 
add_transition!(a, 2, 1, 1);
add_transition!(a, 1, 1, 1);
add_transition!(a, 3, 2, 2);
add_transition!(a, 2, 2, 2);
add_transition!(a, 2, 3, 3);
add_transition!(a, 3, 3, 3);

# subsystems
systems = [ConstrainedAffineControlDiscreteSystem(A[i], B[i], g[i], pX[i], pU) for i in 1:N_region];

switching = AutonomousSwitching()
switchings = fill(switching, 1)
resetmaps = []


system = HybridSystem(a, systems, resetmaps, switchings) # pwa system
# Create Abstraction

max_x = 2 # bound on |x|
rectX = AB.HyperRectangle(SVector(-max_x, -max_x), SVector(max_x, max_x));
rectU = rectX


step = 4 # discretization of one unit of space
X_origin = SVector(0.0, 0.0);
X_step = SVector(1.0/step, 1.0/step);
Xgrid = AB.GridRectangular(X_origin, X_step, rectX);

domainX = AB.DomainList(Xgrid); # state space
AB.add_set!(domainX, rectX, AB.OUTER)


domainU = domainX; # input space
AB.add_set!(domainU, rectU, AB.INNER)

#symbolic model for tilde S
symmodel = AB.NewSymbolicModelListList(domainX, domainU);


# stage cost matrix (J = ||L*[x; u ; 1]||)

L = [eye(n_sys+n_u) zeros(n_sys+n_u,1)];

transitionCost = Dict()  #dictionary with cost of each transition
transitionKappa = Dict() #dictionary with controller associated each transition
empty!(symmodel.autom)

#building abstraction
@time AB.compute_symmodel_from_hybridcontrolsystem!(symmodel,transitionCost, transitionKappa, system, W, L, U)


# Define Specifications
x0 = SVector(2.0,-2.0); # initial condition
Xinit = AB.DomainList(Xgrid)
AB.add_coord!(Xinit, x0)
Xfinal = AB.DomainList(Xgrid) # goal set
#AB.add_coord!(Xfinal, SVector(0.0, 0.0))
AB.add_set!(Xfinal,AB.HyperRectangle(SVector(-1.1, 1.4), SVector(-1.1, 1.8)), AB.OUTER) 


initlist = [AB.get_state_by_xpos(symmodel, pos) for pos in AB.enum_pos(Xinit)]; 
finallist = [AB.get_state_by_xpos(symmodel, pos) for pos in AB.enum_pos(Xfinal)];

# Synthesis of the discrete controler

contr = AB.NewControllerList();
@time AB.compute_controller_reach!(contr, symmodel.autom, initlist, finallist)



# Simulation
function get_mode(x) # return pwa mode for a given x
      o = map(m-> (x ∈ m.X), system.modes).*(1:length(system.modes))
      return o[o.>0][1]
end

K = 100; #max num of steps
x_traj = zeros(n_sys,K+1);
u_traj = zeros(n_u,K+1);
x_traj[:,1] = x0;

k = 1; #iterator
costBound = 0;
costTrue = 0;
currState = AB.get_state_by_xpos(symmodel,AB.get_pos_by_coord(Xgrid,x_traj[:,k]));
println("started at: $(currState)")
while currState ∉ finallist && k ≤ K # While not at goal or not reached max iter 
      println("k: $(k)")
      next_action = nothing
      for action in contr.data
            if(action[1]==currState)
                  next_action = action
                  println("next action: $(next_action)")
                  break
            end
      end
      println("cm: $(AB.get_coord_by_pos(Xgrid, AB.get_xpos_by_state(symmodel, next_action[2])))")
      
      c = AB.get_coord_by_pos(Xgrid, AB.get_xpos_by_state(symmodel, currState))
      println("c: $(c)")


      u_traj[:,k] = transitionKappa[next_action]*vcat(x_traj[:,k]-c,1.0)


      global costBound = costBound + transitionCost[next_action]
      global costTrue += norm(L*vcat(x_traj[:,k], u_traj[:,k],1.0))^2


      m = get_mode(c)
      
      w = (2*rand(2).-1).*W[:,1]

      x_traj[:,k+1] = A[m]*x_traj[:,k]+B[m]*u_traj[:,k] + g[m] + w

      global k += 1;
      global currState = AB.get_state_by_xpos(symmodel, AB.get_pos_by_coord(Xgrid, x_traj[:,k]));
      println("Arrived in: $(currState): $(x_traj[:,k])")
end




#Plotting 

using PyPlot
PyPlot.pygui(true) 


fig = PyPlot.figure()

ax = PyPlot.axes(aspect = "equal")
ax.set_xlim(rectX.lb[1]-0.2,rectX.ub[1]+0.2)
ax.set_ylim(rectX.lb[2]-0.2,rectX.ub[2]+0.2)

vars = [1, 2];


include(dirname(pathof(Dionysos)) * "/plotting.jl")

Plot.domain!(ax, vars, domainX, fc = "white")


d=0.02;
for t in symmodel.autom.transitions.data
      offset = randn(1)[1]*d
      arrow_x, arrow_y = AB.get_coord_by_pos(Xgrid,AB.get_xpos_by_state(symmodel,t[2]))
      aux = (AB.get_coord_by_pos(Xgrid,AB.get_upos_by_symbol(symmodel,t[1]))-[arrow_x, arrow_y])
      arrow_dx, arrow_dy = aux*(norm(aux)-0.15)/norm(aux)
      if t[1]==t[2]
            PyPlot.scatter(arrow_x, arrow_y,s=30, ec="red")
      else
            PyPlot.arrow(arrow_x+offset, arrow_y+offset, arrow_dx, arrow_dy,  ec = "red",head_width=.06)

      end
      

end
gcf() 



fig = PyPlot.figure()

ax = PyPlot.axes(aspect = "equal")
ax.set_xlim(rectX.lb[1]-0.2,rectX.ub[1]+0.2)
ax.set_ylim(rectX.lb[2]-0.2,rectX.ub[2]+0.2)

vars = [1, 2];

include(dirname(pathof(Dionysos)) * "/plotting.jl")

Plot.domain!(ax, vars, domainX, fc = "white")
Plot.domain!(ax, vars, Xinit, fc = "gray")
Plot.domain!(ax, vars, Xfinal, fc = "green")

PyPlot.plot(x_traj[1,1:k],x_traj[2,1:k],"bo-")
gcf() 
