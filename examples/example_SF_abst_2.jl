using Dionysos
using Dionysos.Problem
using Polyhedra
using MathematicalSystems, HybridSystems
using CDDLib
using SemialgebraicSets
using StaticArrays
using LinearAlgebra
using SDPA, Ipopt, JuMP
using Plots, Colors
using Test

const UT = Dionysos.Utils


# Example of ellipsoidal based abstraction

if !isdefined(@__MODULE__, :Usz)
      Usz = 50 # upper limit on |u|
      Wsz = 5
      n_step = 5 # discretization of one unit of space
end

opt_sdp = optimizer_with_attributes(SDPA.Optimizer, MOI.Silent() => true)
opt_qp = optimizer_with_attributes(Ipopt.Optimizer, MOI.Silent() => true)

lib = CDDLib.Library() #polyhedron lib
include("../problems/PWAsys.jl")


dt = 0.01;

const problem = PWAsys.problem(lib, dt, Usz)
const system = problem.system

n_sys = size(system.resetmaps[1].A, 1);
n_u = size(system.resetmaps[1].B, 2);

W = Wsz*[-1 -1  1 1;
         -1  1 -1 1]*dt; # polytope of disturbances

Uaux = diagm(1:n_u)
U = [(Uaux.==i)./Usz for i in 1:n_u];

# Create Abstraction

max_x = 2 # bound on |x|_∞
rectX = Dionysos.Utils.HyperRectangle(SVector(-max_x, -max_x), SVector(max_x, max_x));
rectU = rectX


X_origin = SVector(0.0, 0.0);
X_step = SVector(1.0/n_step, 1.0/n_step);
P = (1/n_sys)*diagm((X_step./2).^(-2))
Xgrid = Dionysos.Domain.GridEllipsoidalRectangular(X_origin, X_step, P, rectX);

domainX = Dionysos.Domain.DomainList(Xgrid); # state space
Dionysos.Domain.add_set!(domainX, rectX, Dionysos.Domain.OUTER)

domainU = domainX; # input space
Dionysos.Domain.add_set!(domainU, rectU, Dionysos.Domain.OUTER)

#symbolic model for tilde S
symmodel = Dionysos.Symbolic.NewSymbolicModelListList(domainX, domainU);


# stage cost matrix (J = ||L*[x; u ; 1]||)
Q_aug = Dionysos.Control.get_full_psd_matrix(problem.transition_cost[1][1])
eigen_Q = eigen(Q_aug);
L = (sqrt.(eigen_Q.values).*(eigen_Q.vectors'))';

transitionCost = Dict()  #dictionary with cost of each transition
transitionKappa = Dict() #dictionary with controller associated each transition
empty!(symmodel.autom)

#building abstraction
@time Dionysos.Symbolic.compute_symmodel_from_hybridcontrolsystem!(symmodel, transitionCost, transitionKappa, system, W, L, U, opt_sdp, opt_qp)

# Define Specifications
x0 = SVector{n_sys}(problem.initial_set) # initial condition
Xinit = Dionysos.Domain.DomainList(Xgrid)
Dionysos.Domain.add_coord!(Xinit, problem.initial_set)
Xfinal = Dionysos.Domain.DomainList(Xgrid) # goal set
Dionysos.Domain.add_coord!(Xfinal, SVector{n_sys}(problem.target_set))
#Dionysos.Domain.add_set!(Xfinal,Dionysos.Domain.HyperRectangle(SVector(-2.0, 0.5), SVector(-2.0, 0.5)), Dionysos.Domain.OUTER) 


Xobstacles = Dionysos.Domain.DomainList(Xgrid) # obstacle set
for o in system.ext[:obstacles]
      Dionysos.Domain.add_set!(Xobstacles, o, Dionysos.Domain.OUTER) 
end

initlist = [Dionysos.Symbolic.get_state_by_xpos(symmodel, pos) for pos in Dionysos.Domain.enum_pos(Xinit)]; 
finallist = [Dionysos.Symbolic.get_state_by_xpos(symmodel, pos) for pos in Dionysos.Domain.enum_pos(Xfinal)];
obstaclelist = [Dionysos.Symbolic.get_state_by_xpos(symmodel, pos) for pos in Dionysos.Domain.enum_pos(Xobstacles)];
# Synthesis of the discrete controler

contr = Dionysos.Control.NewControllerList();


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


K = typeof(problem.time) == Infinity ? 100 : problem.time; #max num of steps
x_traj = zeros(n_sys,K+1);
u_traj = zeros(n_u,K+1);
x_traj[:,1] = x0;
state_traj = []

k = 1; #iterator
costBound = 0;
costTrue = 0;
currState = Dionysos.Symbolic.get_all_states_by_xpos(symmodel,Dionysos.Domain.crop_to_domain(domainX,Dionysos.Domain.get_all_pos_by_coord(Xgrid,x_traj[:,k])))
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
      println("cm: $(Dionysos.Domain.get_coord_by_pos(Xgrid, Dionysos.Symbolic.get_xpos_by_state(symmodel, next_action[2])))")
      
      c = Dionysos.Domain.get_coord_by_pos(Xgrid, Dionysos.Symbolic.get_xpos_by_state(symmodel, next_action[1]))
      println("c: $(c)")


      u_traj[:,k] = transitionKappa[next_action]*vcat(x_traj[:,k]-c,1.0)


      global costBound = costBound + transitionCost[next_action]
      global costTrue += norm(L*vcat(x_traj[:,k], u_traj[:,k],1.0))^2


      m = get_mode(c)
      
      w = (2*(rand(2).^(1/4)).-1).*W[:,1]

      x_traj[:,k+1] = system.resetmaps[m].A*x_traj[:,k]+system.resetmaps[m].B*u_traj[:,k] + system.resetmaps[m].c + w

      global k += 1;
      global currState =  Dionysos.Symbolic.get_all_states_by_xpos(symmodel,Dionysos.Domain.crop_to_domain(domainX,Dionysos.Domain.get_all_pos_by_coord(Xgrid,x_traj[:,k])));
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
                        p1 = Dionysos.Domain.get_coord_by_pos(Xgrid, Dionysos.Symbolic.get_xpos_by_state(symmodel, t[2]))
                        UT.plot_point!(p1; dims=vars, color=color)
                  else
                        p1 = Dionysos.Domain.get_coord_by_pos(Xgrid, Dionysos.Symbolic.get_xpos_by_state(symmodel, t[2]))
                        p2 = Dionysos.Domain.get_coord_by_pos(Xgrid, Dionysos.Symbolic.get_xpos_by_state(symmodel, t[1]))
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
            pos = Dionysos.Symbolic.get_xpos_by_state(symmodel, state)
            elli = Dionysos.Domain.get_elem_by_pos(Xgrid, pos)
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


# @static if get(ENV, "CI", "false") == "false" && (isdefined(@__MODULE__, :no_plot) && no_plot==false)
#       using PyPlot
#       include("../src/utils/plotting/plotting.jl")
#       PyPlot.pygui(true) 


#       PyPlot.rc("text",usetex=true)
#       PyPlot.rc("font",family="serif")
#       PyPlot.rc("font",serif="Computer Modern Roman")
#       PyPlot.rc("text.latex",preamble="\\usepackage{amsfonts}")
#       ##

#       fig = PyPlot.figure(tight_layout=true, figsize=(4,4))

#       ax = PyPlot.axes(aspect = "equal")
#       ax.set_xlim(rectX.lb[1]-0.2,rectX.ub[1]+0.2)
#       ax.set_ylim(rectX.lb[2]-0.2,rectX.ub[2]+0.2)

#       vars = [1, 2];



#       Plot.domain_ellips!(ax, vars, P, domainX, fc = "none", ew = 0.5)


#       d=0.02;
#       for t in symmodel.autom.transitions.data
#             if isempty(t ∩ obstaclelist) 
#                   offset = randn(1)[1]*d
#                   arrow_x, arrow_y = Dionysos.Domain.get_coord_by_pos(Xgrid,Dionysos.Symbolic.get_xpos_by_state(symmodel,t[2]))
#                   aux = (Dionysos.Domain.get_coord_by_pos(Xgrid,Dionysos.Symbolic.get_upos_by_symbol(symmodel,t[1]))-[arrow_x, arrow_y])
#                   arrow_dx, arrow_dy = aux*(norm(aux)-0.15)/norm(aux)
#                   color = (abs(0.6*sin(t[1])), abs(0.6*sin(t[1]+2π/3)), abs(0.6*sin(t[1]-2π/3)));
#                   if t[1]==t[2]
#                         PyPlot.scatter(arrow_x, arrow_y, s=10, fc=color, ec=color)
#                   else
#                         PyPlot.arrow(arrow_x+offset, arrow_y+offset, arrow_dx, arrow_dy,fc=color, ec=color,width=0.01, head_width=.08)
#                   end
#             end

#       end
#       PyPlot.xlabel("\$x_1\$", fontsize=14)
#       PyPlot.ylabel("\$x_2\$", fontsize=14)
#       PyPlot.title("Transitions", fontsize=14)
#       #plt.savefig("ex2_trans.eps", format="eps")
#       gcf() 



#       fig = PyPlot.figure(tight_layout=true, figsize=(4,4))

#       ax = PyPlot.axes(aspect = "equal")
#       ax.set_xlim(rectX.lb[1]-0.2,rectX.ub[1]+0.2)
#       ax.set_ylim(rectX.lb[2]-0.2,rectX.ub[2]+0.2)

#       vars = [1, 2];

#       #Plot.domain_ellips!(ax, vars, P, domainX, fc = "none")
#       # for i=2:(k-1)
#       #       lyap_val = cost[(state_traj[i] ∩ path)[1]]
#       #       aux_dom = AB.DomainList(Xgrid)
#       #       AB.add_coord!(aux_dom, x_traj[:,i])
#       #       Plot.domain_ellips!(ax, vars, P, aux_dom, fc = (0.8*(costBound-lyap_val)/costBound, 0.8*lyap_val/costBound, 0.5))
#       # end
#       lyap_max = max(filter(isfinite,getfield.([cost...],:second))...)

#       cost_ordered = reverse(sort(hcat([ (lyap,state) for (state,lyap) in cost]...),dims=2))
#       for (lyap,state) in cost_ordered
#             aux_dom = Dionysos.Domain.DomainList(Xgrid)
#             x_aux = Dionysos.Domain.get_coord_by_pos(Xgrid, Dionysos.Symbolic.get_upos_by_symbol(symmodel,state))
#             #println(state,": ",lyap)
#             Dionysos.Domain.add_coord!(aux_dom, x_aux)
#             if (lyap ≠ Inf)
#                   Plot.domain_ellips!(ax, vars, P, aux_dom,  ew = 0, fc = (0.4*lyap/lyap_max+0.5, 0.4*(lyap_max-lyap)/lyap_max+0.5, 0.75), fa=1.0)
#             end
            
#       end
#       for (lyap,state) in cost_ordered
#             aux_dom = Dionysos.Domain.DomainList(Xgrid)
#             x_aux = Dionysos.Domain.get_coord_by_pos(Xgrid, Dionysos.Symbolic.get_upos_by_symbol(symmodel,state))
#             #println(state,": ",lyap)
#             Dionysos.Domain.add_coord!(aux_dom, x_aux)
#             if (lyap == Inf)
#                   Plot.domain_ellips!(ax, vars, P, aux_dom, fc = "none", ew = 0.2)
#             else
#                   Plot.domain_ellips!(ax, vars, P, aux_dom, fc = "none", fa=1.0, ew = 0.2)
#             end
            
#       end
#       Plot.domain_ellips!(ax, vars, P, Xinit, fc = "none", ew = 3)
#       Plot.domain_ellips!(ax, vars, P, Xfinal, fc = "none", ew = 3)

#       Plot.domain_ellips!(ax, vars, P, Xobstacles, fc = "black", ew = 0.5, fa=1.0)

#       PyPlot.plot(x_traj[1,1:k],x_traj[2,1:k],"bo-",markersize=4)
#       cmap = PyPlot.ColorMap("mycolor",hcat([0.0,0.8,0.5,0.5],[0.8,0.0,0.5,0.5])');
#       PyPlot.colorbar(PyPlot.ScalarMappable(norm=PyPlot.cm.colors.Normalize(vmin=0, vmax=lyap_max),cmap=cmap),shrink=0.7)

#       PyPlot.xlabel("\$x_1\$", fontsize=14)
#       PyPlot.ylabel("\$x_2\$", fontsize=14)
#       PyPlot.title("Trajectory and Lyapunov-like Fun.", fontsize=14)
#       #plt.savefig("ex2_traj.pdf", format="pdf")
#       gcf() 

# end
