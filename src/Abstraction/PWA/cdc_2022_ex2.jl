using Dionysos
using Polyhedra
using MathematicalSystems, HybridSystems
using CDDLib
using SemialgebraicSets
using StaticArrays
using LinearAlgebra
using Mosek

const AB = Dionysos.Abstraction

lib = CDDLib.Library() #polyhedron lib
eye(n) = diagm(ones(n)) # I matrix
# Define system
N_region = 3
n_sys = 2 
n_u = 2; 
Uaux = diagm(1:n_u)
Usz = 50 # upper limit on |u|
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

W = 5*[-1 -1  1 1;
     -1  1 -1 1]*dt; # polytope of disturbances

# PWAdomains = [pX1, pX2, pX3];

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


n_step = 5 # discretization of one unit of space
X_origin = SVector(0.0, 0.0);
X_step = SVector(1.0/n_step, 1.0/n_step);
P = (1/n_sys)*diagm((X_step./2).^(-2))
Xgrid = AB.GridEllipsoidalRectangular(X_origin, X_step, P, rectX);

domainX = AB.DomainList(Xgrid); # state space
AB.add_set!(domainX, rectX, AB.OUTER)

domainU = domainX; # input space
AB.add_set!(domainU, rectU, AB.OUTER)

#symbolic model for tilde S
symmodel = AB.NewSymbolicModelListList(domainX, domainU);


# stage cost matrix (J = ||L*[x; u ; 1]||)

L = [eye(n_sys+n_u) zeros(n_sys+n_u,1)]*dt;

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
AB.add_coord!(Xfinal, SVector(-2.0, 1.0))
#AB.add_set!(Xfinal,AB.HyperRectangle(SVector(-2.0, 0.5), SVector(-2.0, 0.5)), AB.OUTER) 


Xobstacles = AB.DomainList(Xgrid) # goal set
AB.add_set!(Xobstacles,AB.HyperRectangle(SVector(0.0, -1.0), SVector(0.25, 1.5)), AB.OUTER) 
AB.add_set!(Xobstacles,AB.HyperRectangle(SVector(0.0, 1.25), SVector(1.0, 1.5)), AB.OUTER) 

initlist = [AB.get_state_by_xpos(symmodel, pos) for pos in AB.enum_pos(Xinit)]; 
finallist = [AB.get_state_by_xpos(symmodel, pos) for pos in AB.enum_pos(Xfinal)];
obstaclelist = [AB.get_state_by_xpos(symmodel, pos) for pos in AB.enum_pos(Xobstacles)];
# Synthesis of the discrete controler

contr = AB.NewControllerList();


# function to transform the transition cost from Dict to Vector{Tuple{String, String, Int64}}:
function transFdict(transitionCost)
    testgraph = Vector{Tuple{Int64, Int64, Float64}}()
    for (key, value) in transitionCost
        # append!(testgraph,(string(key(1)),string(key(2)),value))
        if key[2] ∉ obstaclelist
            push!(testgraph,((key[1]),(key[2]),value))
        end
    end
    return testgraph
end

testgraph = transFdict(transitionCost) # uses said function

# control design
src, dst = initlist[1], finallist[1] # initial and goal sets
rev_graph = [(t[2],t[1],t[3]) for t in testgraph]
gc = AB.Digraph(rev_graph) 
@time rev_path, cost = AB.dijkstrapath(gc, dst, src) # gets optimal path
path = reverse(rev_path)
println("Shortest path from $src to $dst: ", isempty(path) ? "no possible path" : join(path, " → "), " (cost $cost[dst])")

## converts path to Int64 and add to list of control actions
for l = 1:length(path)-1
    new_action = (path[l], path[l+1])
    AB.push_new!(contr, new_action)
end

# Simulation
function get_mode(x) # return pwa mode for a given x
      o = map(m-> (x ∈ m.X), system.modes).*(1:length(system.modes))
      return o[o.>0][1]
end

K = 100; #max num of steps
x_traj = zeros(n_sys,K+1);
u_traj = zeros(n_u,K+1);
x_traj[:,1] = x0;
state_traj = []

k = 1; #iterator
costBound = 0;
costTrue = 0;
currState = AB.get_all_states_by_xpos(symmodel,AB.crop_to_domain(domainX,AB.get_all_pos_by_coord(Xgrid,x_traj[:,k])))
push!(state_traj,currState)
println("started at: $(currState)")
while (currState ∩ finallist) == [] && k ≤ K # While not at goal or not reached max iter 
      println("k: $(k)")
      next_action = nothing
      for action in contr.data
            if (action[1] ∩ currState) ≠ []
                  next_action = action
                  println("next action: $(next_action)")
            end
      end
      println("cm: $(AB.get_coord_by_pos(Xgrid, AB.get_xpos_by_state(symmodel, next_action[2])))")
      
      c = AB.get_coord_by_pos(Xgrid, AB.get_xpos_by_state(symmodel, next_action[1]))
      println("c: $(c)")


      u_traj[:,k] = transitionKappa[next_action]*vcat(x_traj[:,k]-c,1.0)


      global costBound = costBound + transitionCost[next_action]
      global costTrue += norm(L*vcat(x_traj[:,k], u_traj[:,k],1.0))^2


      m = get_mode(c)
      
      w = (2*(rand(2).^(1/4)).-1).*W[:,1]

      x_traj[:,k+1] = A[m]*x_traj[:,k]+B[m]*u_traj[:,k] + g[m] + w

      global k += 1;
      global currState =  AB.get_all_states_by_xpos(symmodel,AB.crop_to_domain(domainX,AB.get_all_pos_by_coord(Xgrid,x_traj[:,k])));
      push!(state_traj,currState)
      println("Arrived in: $(currState): $(x_traj[:,k])")
end

println("Goal set reached")
println("Guaranteed cost:\t $(costBound)")
println("True cost:\t\t $(costTrue)")


#Plotting  #######################################################

using PyPlot
PyPlot.pygui(true) 


PyPlot.rc("text",usetex=true)
PyPlot.rc("font",family="serif")
PyPlot.rc("font",serif="Computer Modern Roman")
PyPlot.rc("text.latex",preamble="\\usepackage{amsfonts}")
##

fig = PyPlot.figure(tight_layout=true, figsize=(4,4))

ax = PyPlot.axes(aspect = "equal")
ax.set_xlim(rectX.lb[1]-0.2,rectX.ub[1]+0.2)
ax.set_ylim(rectX.lb[2]-0.2,rectX.ub[2]+0.2)

vars = [1, 2];


include(dirname(pathof(Dionysos)) * "/plotting.jl")

Plot.domain_ellips!(ax, vars, P, domainX, fc = "none", ew = 0.5)


d=0.02;
for t in symmodel.autom.transitions.data
      if isempty(t ∩ obstaclelist) 
            offset = randn(1)[1]*d
            arrow_x, arrow_y = AB.get_coord_by_pos(Xgrid,AB.get_xpos_by_state(symmodel,t[2]))
            aux = (AB.get_coord_by_pos(Xgrid,AB.get_upos_by_symbol(symmodel,t[1]))-[arrow_x, arrow_y])
            arrow_dx, arrow_dy = aux*(norm(aux)-0.15)/norm(aux)
            color = (abs(0.6*sin(t[1])), abs(0.6*sin(t[1]+2π/3)), abs(0.6*sin(t[1]-2π/3)));
            if t[1]==t[2]
                  PyPlot.scatter(arrow_x, arrow_y, s=10, fc=color, ec=color)
            else
                  PyPlot.arrow(arrow_x+offset, arrow_y+offset, arrow_dx, arrow_dy,fc=color, ec=color,width=0.01, head_width=.08)
            end
      end

end
PyPlot.xlabel("\$x_1\$", fontsize=14)
PyPlot.ylabel("\$x_2\$", fontsize=14)
PyPlot.title("Transitions", fontsize=14)
plt.savefig("ex2_trans.eps", format="eps")
#gcf() 



fig = PyPlot.figure(tight_layout=true, figsize=(4,4))

ax = PyPlot.axes(aspect = "equal")
ax.set_xlim(rectX.lb[1]-0.2,rectX.ub[1]+0.2)
ax.set_ylim(rectX.lb[2]-0.2,rectX.ub[2]+0.2)

vars = [1, 2];

include(dirname(pathof(Dionysos)) * "/plotting.jl")

#Plot.domain_ellips!(ax, vars, P, domainX, fc = "none")
# for i=2:(k-1)
#       lyap_val = cost[(state_traj[i] ∩ path)[1]]
#       aux_dom = AB.DomainList(Xgrid)
#       AB.add_coord!(aux_dom, x_traj[:,i])
#       Plot.domain_ellips!(ax, vars, P, aux_dom, fc = (0.8*(costBound-lyap_val)/costBound, 0.8*lyap_val/costBound, 0.5))
# end
lyap_max = max(filter(isfinite,getfield.([cost...],:second))...)

cost_ordered = reverse(sort(hcat([ (lyap,state) for (state,lyap) in cost]...),dims=2))
for (lyap,state) in cost_ordered
      aux_dom = AB.DomainList(Xgrid)
      x_aux = AB.get_coord_by_pos(Xgrid,AB.get_upos_by_symbol(symmodel,state))
      #println(state,": ",lyap)
      AB.add_coord!(aux_dom, x_aux)
      if (lyap ≠ Inf)
            Plot.domain_ellips!(ax, vars, P, aux_dom,  ew = 0, fc = (0.4*lyap/lyap_max+0.5, 0.4*(lyap_max-lyap)/lyap_max+0.5, 0.75), fa=1.0)
      end
      
end
for (lyap,state) in cost_ordered
      aux_dom = AB.DomainList(Xgrid)
      x_aux = AB.get_coord_by_pos(Xgrid,AB.get_upos_by_symbol(symmodel,state))
      #println(state,": ",lyap)
      AB.add_coord!(aux_dom, x_aux)
      if (lyap == Inf)
            Plot.domain_ellips!(ax, vars, P, aux_dom, fc = "none", ew = 0.2)
      else
            Plot.domain_ellips!(ax, vars, P, aux_dom, fc = "none", fa=1.0, ew = 0.2)
      end
      
end
Plot.domain_ellips!(ax, vars, P, Xinit, fc = "none", ew = 3)
Plot.domain_ellips!(ax, vars, P, Xfinal, fc = "none", ew = 3)

Plot.domain_ellips!(ax, vars, P, Xobstacles, fc = "black", ew = 0.5, fa=1.0)

PyPlot.plot(x_traj[1,1:k],x_traj[2,1:k],"bo-",markersize=4)
cmap = PyPlot.ColorMap("mycolor",hcat([0.0,0.8,0.5,0.5],[0.8,0.0,0.5,0.5])');
PyPlot.colorbar(PyPlot.ScalarMappable(norm=PyPlot.cm.colors.Normalize(vmin=0, vmax=lyap_max),cmap=cmap),shrink=0.7)

PyPlot.xlabel("\$x_1\$", fontsize=14)
PyPlot.ylabel("\$x_2\$", fontsize=14)
PyPlot.title("Trajectory and Lyapunov-like Fun.", fontsize=14)
plt.savefig("ex2_traj.pdf", format="pdf")
#gcf() 


