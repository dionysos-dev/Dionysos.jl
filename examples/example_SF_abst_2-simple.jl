using Polyhedra
using MathematicalSystems, HybridSystems
using CDDLib
using SemialgebraicSets
using StaticArrays
using LinearAlgebra
using SDPA, Ipopt, JuMP
using Plots, Colors
using Test

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const CO = DI.Control
const PR = DI.Problem

# Example of ellipsoidal based abstraction

if !isdefined(@__MODULE__, :Usz)
      Usz = 50 # upper limit on |u|
      Wsz = 5
      n_step = 3 # discretization of one unit of space
end

opt_sdp = optimizer_with_attributes(SDPA.Optimizer, MOI.Silent() => true)
opt_qp = optimizer_with_attributes(Ipopt.Optimizer, MOI.Silent() => true)

lib = CDDLib.Library() #polyhedron lib
include("../problems/PWAsys.jl")


dt = 0.01
x_0 = [-1.0,-0.6]
const problem = PWAsys.problem(lib, dt, Usz, x_0)
const system = problem.system

n_sys = size(system.resetmaps[1].A, 1);
n_u = size(system.resetmaps[1].B, 2);

W = Wsz*[-1 -1  1 1;
         -1  1 -1 1]*dt; # polytope of disturbances

Uaux = diagm(1:n_u)
U = [(Uaux.==i)./Usz for i in 1:n_u];

# Create Abstraction

rectX = UT.HyperRectangle(SVector(-2.0, -1.5), SVector(-0.5, 1.3));
rectU = rectX


X_origin = SVector(0.0, 0.0);
X_step = SVector(1.0/n_step, 1.0/n_step);
P = (1/n_sys)*diagm((X_step./2).^(-2))
Xgrid = DO.GridEllipsoidalRectangular(X_origin, X_step, P, rectX);

domainX = DO.DomainList(Xgrid); # state space
DO.add_set!(domainX, rectX, DO.OUTER)

domainU = domainX; # input space
DO.add_set!(domainU, rectU, DO.OUTER)

#symbolic model for tilde S
symmodel = SY.NewSymbolicModelListList(domainX, domainU);


# stage cost matrix (J = ||L*[x; u ; 1]||)
Q_aug = CO.get_full_psd_matrix(problem.transition_cost[1][1])
eigen_Q = eigen(Q_aug);
L = (sqrt.(eigen_Q.values).*(eigen_Q.vectors'))';

transitionCost = Dict()  #dictionary with cost of each transition
transitionKappa = Dict() #dictionary with controller associated each transition
empty!(symmodel.autom)

#building abstraction
@time SY.compute_symmodel_from_hybridcontrolsystem!(symmodel, transitionCost, transitionKappa, system, W, L, U, opt_sdp, opt_qp)

# Define Specifications
x0 = SVector{n_sys}(problem.initial_set) # initial condition
Xinit = DO.DomainList(Xgrid)
DO.add_coord!(Xinit, problem.initial_set)
Xfinal = DO.DomainList(Xgrid) # goal set
DO.add_coord!(Xfinal, SVector{n_sys}(problem.target_set))


Xobstacles = DO.DomainList(Xgrid) # obstacle set
for o in system.ext[:obstacles]
      DO.add_set!(Xobstacles, o, DO.OUTER) 
end

initlist = [SY.get_state_by_xpos(symmodel, pos) for pos in DO.enum_pos(Xinit)]; 
finallist = [SY.get_state_by_xpos(symmodel, pos) for pos in DO.enum_pos(Xfinal)];
filterTab = [pos for pos in DO.enum_pos(Xobstacles)]
filtered_array = filter(x->SY.is_in(symmodel, x), filterTab)
obstaclelist = [SY.get_state_by_xpos(symmodel, pos) for pos in filtered_array]; 
# Synthesis of the discrete controler

contr = CO.NewControllerList();


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
gc = UT.Digraph(rev_graph) 
UT.add_states!(gc, [src, dst])

@time rev_path, cost = UT.dijkstrapath(gc, dst, src) # gets optimal path
path = reverse(rev_path)
println("Shortest path from $src to $dst: ", isempty(path) ? "no possible path" : join(path, " → "), " (cost $cost[dst])")

## converts path to Int64 and add to list of control actions
for l = 1:length(path)-1
    new_action = (path[l], path[l+1])
    UT.push_new!(contr, new_action)
end

# Simulation

# return pwa mode for a given x
get_mode(x) = findfirst(m -> (x ∈ m.X), system.resetmaps)


K = typeof(problem.time) == PR.Infinity ? 100 : problem.time; #max num of steps
x_traj = zeros(n_sys,K+1);
u_traj = zeros(n_u,K+1);
x_traj[:,1] = x0;
state_traj = []

k = 1; #iterator
costBound = 0;
costTrue = 0;
currState = SY.get_all_states_by_xpos(symmodel,DO.crop_to_domain(domainX,DO.get_all_pos_by_coord(Xgrid,x_traj[:,k])))
push!(state_traj, currState)
println("started at: $(currState)")
while (currState ∩ finallist) == [] && k ≤ K && k≤ length(path)-1 # While not at goal or not reached max iter 
      println("k: $(k)")
      next_action = nothing
      for action in contr.data
            if (action[1] ∩ currState) ≠ []
                  next_action = action
                  println("next action: $(next_action)")
            end
      end
      println("cm: $(DO.get_coord_by_pos(Xgrid, SY.get_xpos_by_state(symmodel, next_action[2])))")
      
      c = DO.get_coord_by_pos(Xgrid, SY.get_xpos_by_state(symmodel, next_action[1]))
      println("c: $(c)")


      u_traj[:,k] = transitionKappa[next_action]*vcat(x_traj[:,k]-c,1.0)


      global costBound = costBound + transitionCost[next_action]
      global costTrue += norm(L*vcat(x_traj[:,k], u_traj[:,k],1.0))^2


      m = get_mode(c)
      
      w = (2*(rand(2).^(1/4)).-1).*W[:,1]

      x_traj[:,k+1] = system.resetmaps[m].A*x_traj[:,k]+system.resetmaps[m].B*u_traj[:,k] + system.resetmaps[m].c + w

      global k += 1;
      global currState = SY.get_all_states_by_xpos(symmodel,DO.crop_to_domain(domainX,DO.get_all_pos_by_coord(Xgrid,x_traj[:,k])));
      push!(state_traj,currState)
      println("Arrived in: $(currState): $(x_traj[:,k])")
end

println("Goal set reached")
println("Guaranteed cost:\t $(costBound)")
println("True cost:\t\t $(costTrue)")

@static if get(ENV, "CI", "false") == "false" && (isdefined(@__MODULE__, :no_plot) && no_plot==false)
      ############# Figure 1 #############

      fig = plot(aspect_ratio=:equal, xtickfontsize=10, ytickfontsize=10, guidefontsize=16, titlefontsize=14)
      xlims!(rectX.lb[1]-0.2,rectX.ub[1]+0.2)
      ylims!(rectX.lb[2]-0.2,rectX.ub[2]+0.2)
      vars = [1, 2]
      Plots.plot!(domainX; dims=[1,2], color=:yellow, opacity=0.2)
      for t in symmodel.autom.transitions.data
            if isempty(t ∩ obstaclelist) 
                  color = RGB(abs(0.6*sin(t[1])), abs(0.6*sin(t[1]+2π/3)), abs(0.6*sin(t[1]-2π/3)))
                  if t[1]==t[2]
                        p1 = DO.get_coord_by_pos(Xgrid, SY.get_xpos_by_state(symmodel, t[2]))
                        UT.plot_point!(p1; dims=vars, color=color)
                  else
                        p1 = DO.get_coord_by_pos(Xgrid, SY.get_xpos_by_state(symmodel, t[2]))
                        p2 = DO.get_coord_by_pos(Xgrid, SY.get_xpos_by_state(symmodel, t[1]))
                        UT.plot_arrow!(p1, p2; dims=vars, color=color)
                  end
            end
      end
      xlabel!("\$x_1\$")
      ylabel!("\$x_2\$")
      title!("Transitions")
      # savefig("ex2_trans.png")
      display(fig)

      ############# Figure 2 #############

      fig = plot(aspect_ratio=:equal, xtickfontsize=10, ytickfontsize=10, guidefontsize=16, titlefontsize=14)
      xlims!(rectX.lb[1]-0.2,rectX.ub[1]+0.2)
      ylims!(rectX.lb[2]-0.2,rectX.ub[2]+0.2)
      vars = [1, 2]

      LyapMax = max(filter(isfinite,getfield.([cost...],:second))...)
      colormap = Colors.colormap("Blues")
      mycolorMap = UT.Colormap([0.0,LyapMax], colormap)
      cost_ordered = reverse(sort(hcat([ (lyap,state) for (state,lyap) in cost]...), dims=2))
      for (lyap, state) in cost_ordered
            pos = SY.get_xpos_by_state(symmodel, state)
            elli = DO.get_elem_by_pos(Xgrid, pos)
            if (lyap ≠ Inf)
                  UT.plotE!(elli, color=UT.get_color(mycolorMap, lyap))
            else
                  UT.plotE!(elli, color=:yellow)
            end
      end
      Einit = UT.Ellipsoid(collect(P), collect(problem.initial_set))
      Etarget = UT.Ellipsoid(collect(P), collect(problem.target_set))
      UT.plotE!(Einit, color=:green)
      UT.plotE!(Etarget, color=:red)
      Plots.plot!(Xobstacles, color=:black, opacity=1.0)

      trajCoord = [[x_traj[1,i], x_traj[2,i]] for i in 1:k]
      UT.plot_traj!(trajCoord, color=:black)

      UT.plot_colorBar!(mycolorMap)
      xlabel!("\$x_1\$")
      ylabel!("\$x_2\$")
      title!("Trajectory and Lyapunov-like Fun.")
      # savefig("ex2_traj.png")
      display(fig)
end
