using Test     #src
# # Example: Optimal control of a PWA System by State-feedback Abstractions
#
#md # [![Binder](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/State-feedback Abstraction: PWA System.ipynb)
#md # [![nbviewer](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/State-feedback Abstraction: PWA System.ipynb)
# 
# This document reproduces [1, Example 2], which is a possible application of state-feedback transition system for the optimal control of piecewise-affine systems.
# Consider a system $\mathcal{S}:=(\mathcal{X}, \mathcal{U},\rightarrow_F)$ with the transition function  
# ```math
# \begin{equation}
#     \{A_{\psi(x)}x+B_{\psi(x)}u+g_{\psi(x)}\}\oplus\Omega_{\psi(x)}
# \end{equation}
# ```
# where $\psi:\mathcal{X}\rightarrow\{1,\dots,N_p\}$ selects one of the $N_p$ subsystems defined by the matrices 
# ```math
# \begin{equation}
# 	A_1=\begin{bmatrix}
#     1.01 & 0.3\\
#     -0.1 & 1.01
# \end{bmatrix}, ~B_1=\begin{bmatrix}
#     1&0\\ 0 & 1
# \end{bmatrix},~g_1=\begin{bmatrix}
# -0.1\\-0.1
# \end{bmatrix},
# \end{equation}
# ```
# $A_2=A_1^\top,~ A_3=A_1,~B_2=B_3=B_1,~g_2=0$ and $g_3=-g_1$. These systems are three spiral sources with unstable equilibria at $x_{e1}=[-0.9635~~0.3654]^\top,~x_{e2}=0,$ and $x_{e3}=-x_{e1}$. We also define the additive-noise sets $\Omega_1=\Omega_2=\Omega_3=[-0.05,0.05]^2$, the control-input set $\mathcal{U}=[-0.5,0.5]^2$ and the state space $\mathcal{X}=[-2,2]^2$. The $N_p=3$ partitions of $\mathcal{X}$ are $\mathcal{X}_1= \{x\in\mathcal{X}~:~x_1\leq-1 \},~\mathcal{X}_3= \{x\in\mathcal{X}~:~x_1>1 \},$ and $\mathcal{X}_2=\mathcal{X}\setminus(\mathcal{X}_1\cup\mathcal{X}_3)$. The goal is to bring the state $x$ from the initial set $\mathcal{X}_0$ to a final set $\mathcal{X}_*$, while avoiding the obstacle $\mathcal{O}$, which are to be defined. 

# First, let us import 
# [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl), 
# [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/), 
# [CDDLib](https://github.com/JuliaPolyhedra/CDDLib.jl), 
# [SDPA](https://github.com/jump-dev/SDPA.jl),
# [Ipopt](https://github.com/jump-dev/Ipopt.jl), and
# [JuMP](https://jump.dev/JuMP.jl/stable/). We also instantiate our optimizers and CDDLib.

using StaticArrays
using LinearAlgebra
using CDDLib
using SDPA, Ipopt, JuMP

opt_sdp = optimizer_with_attributes(SDPA.Optimizer, MOI.Silent() => true)
opt_qp = optimizer_with_attributes(Ipopt.Optimizer, MOI.Silent() => true)
lib = CDDLib.Library() #polyhedron lib

# We define now the of the polytope $\Omega$ with the disturbances and the input constraint set $\mathcal{U}$ and we imported the optimal control problem for this example, which is defined in the ```PWAsys.jl``` problem:
Wsz = 3 # `Wsz = 5` in [1]
Usz = 70 # upper limit on |u|, `Usz = 50` in [1]
dt = 0.01; # discretization step

include("../../../problems/PWAsys.jl")

const problem = PWAsys.problem(lib, dt, Usz)
const system = problem.system

n_x = size(system.resetmaps[1].A,1)
n_u = size(system.resetmaps[1].B,2) 

W = Wsz*[-1 -1  1 1;
         -1  1 -1 1]*dt; # polytope of disturbances 

Uaux = diagm(1:n_u)
U = [(Uaux.==i)./Usz for i in 1:n_u]; # matrices U_i

# notice that in [1] it was used `Wsz = 5` and `Usz = 50`. These, and other values were changed here to speed up the build time of the documentation.
# 



# Let us know get some more information about the problem as the dimension $n_x$ of the system


# At this point we'll import Dionysos in order to solve our optimal control problem
using Dionysos
using Dionysos.Problem

# Let us now define the state space $X$ within which we are searching for a optimal solution.
max_x = 2 # bound on |X|_∞
rectX = Dionysos.Utils.HyperRectangle(SVector(-max_x, -max_x), SVector(max_x, max_x));
rectU = rectX

# This is state-space is defined by the `HyperRectangle rectX`. We also define a control space with the same bounds. This is done because, for a state-feedback abstraction, selecting a controller out of the set of controllers is the same as selecting a destination state out of the set of cells $\mathcal{X}_d$, given it's determinism. 

# To build this deterministic state-feedback abstraction in alternating simulation relation  with the system as described in [1, Lemma 1], a set of balls of radius 0.2 covering the state space is adopted as cells $\xi\in\mathcal{X}_d$. We assume that inside cells intersecting the boundary of partitions of $\mathcal{X}$ the selected piecewise-affine mode is the same all over its interior and given by the mode defined at its center. An alternative to this are discussed in [1]. Let us define the corresponding grid:

n_step = 2 #  n_step = 5 in [1]
X_origin = SVector(0.0, 0.0);
X_step = SVector(1.0/n_step, 1.0/n_step);
P = (1/n_x)*diagm((X_step./2).^(-2))

Xgrid = Dionysos.Domain.GridEllipsoidalRectangular(X_origin, X_step, P, rectX) 

# and finally let us instantiate $\mathcal{X}_d$ as `domainX`

domainX = Dionysos.Domain.DomainList(Xgrid); # state space
Dionysos.Domain.add_set!(domainX, rectX, Dionysos.Domain.OUTER)

# the input space 

domainU = domainX; # input space
Dionysos.Domain.add_set!(domainU, rectU, Dionysos.Domain.OUTER)

#and the symbolic model for the state-feedback abstraction $\tilde{\mathcal{S}}$
symmodel = Dionysos.Symbolic.NewSymbolicModelListList(domainX, domainU);
empty!(symmodel.autom)


# Now let us define the L matrix defining the stage cost $\mathcal{J}(x,u) = ||L*[x; u ; 1]||^2_2$
Q_aug = Dionysos.Control.get_full_psd_matrix(problem.transition_cost[1][1])
eigen_Q = eigen(Q_aug);
L = (sqrt.(eigen_Q.values).*(eigen_Q.vectors'))'*dt

# ## Building the abstraction

# We then initialize the dictionaries for saving the cost and the controller associated with each transition in the abstraction

transitionCost = Dict()  #dictionary with cost of each transition
transitionKappa = Dict() #dictionary with controller associated each transition


# and finally build the state-feedback abstraction 
using Suppressor
@suppress  begin # this is a workaround to supress the undesired output of SDPA
    global t 
    t = @elapsed Dionysos.Symbolic.compute_symmodel_from_hybridcontrolsystem!(symmodel,transitionCost, transitionKappa, system, W, L, U, opt_sdp, opt_qp);
end
println("Abstraction created in $t seconds with $(length(transitionCost)) transitions")


# Now let us prepare to synthesize our controller. Define the specifications
x0 = SVector{n_x}(problem.initial_set); # initial condition

Xinit = Dionysos.Domain.DomainList(Xgrid) # set of initial cells
Dionysos.Domain.add_coord!(Xinit, problem.initial_set)

Xfinal = Dionysos.Domain.DomainList(Xgrid) # set of target cells
Dionysos.Domain.add_coord!(Xfinal, SVector{n_x}(problem.target_set))

Xobstacles = Dionysos.Domain.DomainList(Xgrid) # set of obstacle cells
for o in system.ext[:obstacles]
      Dionysos.Domain.add_set!(Xobstacles, o, Dionysos.Domain.OUTER) 
end

initlist = [Dionysos.Symbolic.get_state_by_xpos(symmodel, pos) for pos in Dionysos.Domain.enum_pos(Xinit)]; 
finallist = [Dionysos.Symbolic.get_state_by_xpos(symmodel, pos) for pos in Dionysos.Domain.enum_pos(Xfinal)];
obstaclelist = [Dionysos.Symbolic.get_state_by_xpos(symmodel, pos) for pos in Dionysos.Domain.enum_pos(Xobstacles)];


# ## Control synthesis

# Before synthesizing our controller let us define this auxiliary function to transform the transition cost from Dict to Vector{Tuple{String, String, Int64}} and use it. It also removes states associated with obstacles in $\mathcal{O}$

function transFdict(transitionCost)
    testgraph = Vector{Tuple{Int64, Int64, Float64}}()
    for (key, value) in transitionCost
        
        if key[2] ∉ obstaclelist
            push!(testgraph,((key[1]),(key[2]),value))
        end
    end
    return testgraph
end

testgraph = transFdict(transitionCost); # uses said function


# We perform the synthesis of the discrete controller through Dijkstra's algorithm [2, p.86] to the reversed graph associated to the abstraction $\tilde{\mathcal{S}}$

contr = Dionysos.Control.NewControllerList();

src, dst = initlist[1], finallist[1] # initial and goal sets

rev_graph = [(t[2],t[1],t[3]) for t in testgraph] # we applied dijkstra to the reversed graph

gc = Dionysos.Search.Digraph(rev_graph) 
t = @elapsed rev_path, cost = Dionysos.Search.dijkstrapath(gc, dst, src) # gets optimal path
path = reverse(rev_path)


println("Shortest path from $src to $dst found in $t seconds:\n ", isempty(path) ? "no possible path" : join(path, " → "), " (cost $(cost[dst]))")

# We create the list of control actions along the optimal path
for l = 1:length(path)-1 
    new_action = (path[l], path[l+1])
    Dionysos.Utils.push_new!(contr, new_action)
end

# and now we will simulate one trajectory of the system.

# ## Simulation
# Let us initialize variables used in our simulation
ϕ(x) = findfirst(m -> (x ∈ m.X), system.resetmaps) # returns pwa mode for a given x


K = typeof(problem.time) == Infinity ? 100 : problem.time; #max num of steps
x_traj = zeros(n_x,K+1);
u_traj = zeros(n_u,K+1);
x_traj[:,1] = x0;
state_traj = []

k = 1; #iterator
costBound = 0;
costTrue = 0;
currState = Dionysos.Symbolic.get_all_states_by_xpos(symmodel,Dionysos.Domain.crop_to_domain(domainX,Dionysos.Domain.get_all_pos_by_coord(Xgrid,x_traj[:,k])))
push!(state_traj,currState);

# and proceed with the simulation
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
      println("cm: $(Dionysos.Domain.get_coord_by_pos(Xgrid, Dionysos.Symbolic.get_xpos_by_state(symmodel, next_action[2])))")
      
      c = Dionysos.Domain.get_coord_by_pos(Xgrid, Dionysos.Symbolic.get_xpos_by_state(symmodel, next_action[1]))
      println("c: $(c)")


      u_traj[:,k] = transitionKappa[next_action]*vcat(x_traj[:,k]-c,1.0)
      println("x: $(x_traj[:,k])")
      println("u: $(u_traj[:,k])")

      global costBound = costBound + transitionCost[next_action]
      global costTrue += norm(L*vcat(x_traj[:,k], u_traj[:,k],1.0))^2


      m = ϕ(c)
      
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

# Finally let us visualize the results. Let us plot the transitions in the state-feedback abstraction

using PyPlot
include("../../../src/utils/plotting/plotting.jl")

PyPlot.pygui(true) #jl


PyPlot.rc("text",usetex=true)
PyPlot.rc("font",family="serif")
PyPlot.rc("font",serif="Computer Modern Roman")
PyPlot.rc("text.latex",preamble="\\usepackage{amsfonts}")


fig = PyPlot.figure(tight_layout=true, figsize=(4,4))

ax = PyPlot.axes(aspect = "equal")
ax.set_xlim(rectX.lb[1]-0.2,rectX.ub[1]+0.2)
ax.set_ylim(rectX.lb[2]-0.2,rectX.ub[2]+0.2)

vars = [1, 2];



Plot.domain_ellips!(ax, vars, P, domainX, fc = "none", ew = 0.5)


d=0.02;
for t in symmodel.autom.transitions.data
      if isempty(t ∩ obstaclelist) 
            offset = randn(1)[1]*d
            arrow_x, arrow_y = Dionysos.Domain.get_coord_by_pos(Xgrid,Dionysos.Symbolic.get_xpos_by_state(symmodel,t[2]))
            aux = (Dionysos.Domain.get_coord_by_pos(Xgrid,Dionysos.Symbolic.get_upos_by_symbol(symmodel,t[1]))-[arrow_x, arrow_y])
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
gcf() #md
gcf() #nb

# Then we plot a colormap with the Lyapunov-like function $v(x)$ and the simulated trajectory in blue.

fig = PyPlot.figure(tight_layout=true, figsize=(4,4))

ax = PyPlot.axes(aspect = "equal")
ax.set_xlim(rectX.lb[1]-0.2,rectX.ub[1]+0.2)
ax.set_ylim(rectX.lb[2]-0.2,rectX.ub[2]+0.2)

vars = [1, 2];

lyap_max = max(filter(isfinite,getfield.([cost...],:second))...)

cost_ordered = reverse(sort(hcat([ (lyap,state) for (state,lyap) in cost]...),dims=2))
for (lyap,state) in cost_ordered
      aux_dom = Dionysos.Domain.DomainList(Xgrid)
      x_aux = Dionysos.Domain.get_coord_by_pos(Xgrid, Dionysos.Symbolic.get_upos_by_symbol(symmodel,state))
      Dionysos.Domain.add_coord!(aux_dom, x_aux)
      if (lyap ≠ Inf)
            Plot.domain_ellips!(ax, vars, P, aux_dom,  ew = 0, fc = (0.4*lyap/lyap_max+0.5, 0.4*(lyap_max-lyap)/lyap_max+0.5, 0.75), fa=1.0)
      end
      
end
for (lyap,state) in cost_ordered
      aux_dom = Dionysos.Domain.DomainList(Xgrid)
      x_aux = Dionysos.Domain.get_coord_by_pos(Xgrid, Dionysos.Symbolic.get_upos_by_symbol(symmodel,state))
      Dionysos.Domain.add_coord!(aux_dom, x_aux)
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
gcf() #md
gcf() #nb

# We recall that, to speed up the build time of this documentation, some values were modified in comparison with [1, Example 2]. To obtain the same figures use `Usz = 50`, `Wsz = 5` and `n_step = 5`.

@test contr.data[1] == (45,21)                  #src
@test costBound ≈ 3.230370549 rtol=1e-3         #src
@test costTrue ≈ 2.365335178 rtol=1e-3          #src

# ## References
#
# 1. L. N. Egidio, T. Alves Lima, R. M. Jungers, "State-feedback Abstractions for Optimal Control of Piecewise-affine Systems", IEEE 61st Conference on Decision and Control (CDC), 2022, accepted.
# 1. D. Bertsekas, "Dynamic programming and optimal control". Volume I, Athena scientific, 2012.
