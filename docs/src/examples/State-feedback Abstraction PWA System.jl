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
opt_ip = optimizer_with_attributes(Ipopt.Optimizer, MOI.Silent() => true)
lib = CDDLib.Library() #polyhedron lib

# We define now the of the polytope $\Omega$ with the disturbances and the input constraint set $\mathcal{U}$ and we imported the optimal control problem for this example, which is defined in the ```PWAsys.jl``` problem:
Wsz = 3 # `Wsz = 5` in [1]
Usz = 70 # upper limit on |u|, `Usz = 50` in [1]
dt = 0.01; # discretization step

include("../../../problems/PWAsys.jl")

problem = PWAsys.problem(lib, dt, Usz)
system = problem.system

n_x = size(system.resetmaps[1].A,1)
n_u = size(system.resetmaps[1].B,2) 

W = Wsz*[-1 -1  1 1;
         -1  1 -1 1]*dt; # polytope of disturbances 

Uaux = diagm(1:n_u)
U = [(Uaux.==i)./Usz for i in 1:n_u]; # matrices U_i

system.ext[:W] = W
system.ext[:U] = U

# notice that in [1] it was used `Wsz = 5` and `Usz = 50`. These, and other values were changed here to speed up the build time of the documentation.



# Let us now get some more information about the problem as the dimension $n_x$ of the system


# At this point we'll import Dionysos in order to solve our optimal control problem
using Dionysos

const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const SY = DI.Symbolic
const CO = DI.Control
const PR = DI.Problem
const OP = DI.Optim
const AB = OP.Abstraction

# Let us now define the state space $X$ within which we are searching for a optimal solution.
max_x = 2 # bound on |X|_∞
rectX = UT.HyperRectangle(SVector(-max_x, -max_x), SVector(max_x, max_x));
rectU = rectX

# This is state-space is defined by the `HyperRectangle rectX`. We also define a control space with the same bounds. This is done because, for a state-feedback abstraction, selecting a controller out of the set of controllers is the same as selecting a destination state out of the set of cells $\mathcal{X}_d$, given it's determinism. 

# To build this deterministic state-feedback abstraction in alternating simulation relation  with the system as described in [1, Lemma 1], a set of balls of radius 0.2 covering the state space is adopted as cells $\xi\in\mathcal{X}_d$. We assume that inside cells intersecting the boundary of partitions of $\mathcal{X}$ the selected piecewise-affine mode is the same all over its interior and given by the mode defined at its center. An alternative to this are discussed in [1]. Let us define the corresponding grid:

n_step = 2 #  n_step = 5 in [1]
X_origin = SVector(0.0, 0.0);
X_step = SVector(1.0/n_step, 1.0/n_step);
P = (1/n_x)*diagm((X_step./2).^(-2))

state_grid = DO.Grid
soidalRectangular(X_origin, X_step, P, rectX) 

# At this point, we instantiate the optimizer provided in Dionysos that creates ellipsoidal-based 
# abstractions `EllipsoidsAbstraction.Optimizer`

using JuMP
optimizer = MOI.instantiate(AB.EllipsoidsAbstraction.Optimizer) # 

MOI.set(optimizer, MOI.RawOptimizerAttribute("problem"), problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("sdp_solver"), opt_sdp)
MOI.set(optimizer, MOI.RawOptimizerAttribute("ip_solver"), opt_ip)


# and run it to build the state-feedback abstraction and solve the optimal control problem by through Dijkstra's algorithm [2, p.86].
using Suppressor
@suppress  begin # this is a workaround to supress the undesired output of SDPA
      MOI.optimize!(optimizer)
end
contr = MOI.get(optimizer, MOI.RawOptimizerAttribute("controller"))
symmodel = MOI.get(optimizer, MOI.RawOptimizerAttribute("symmodel"))
transitionCost = MOI.get(optimizer, MOI.RawOptimizerAttribute("transitionCost"))
transitionKappa = MOI.get(optimizer, MOI.RawOptimizerAttribute("transitionKappa"))
lyap_fun = MOI.get(optimizer, MOI.RawOptimizerAttribute("lyap_fun"));


# ## Simulation
# Now let us get the Q_aug matrix defining the stage cost $\mathcal{J}(x,u) = ||Q^{1/2} \cdot [x; u ; 1]||^2_2$
Q_aug = CO.get_full_psd_matrix(problem.transition_cost[1][1])
# Let us initialize variables used in our simulation and now we will simulate one trajectory of the system.
x0 = Vector(problem.initial_set); # initial condition

domainX = symmodel.Xdom

ϕ(x) = findfirst(m -> (x ∈ m.X), system.resetmaps) # returns pwa mode for a given x

K = typeof(problem.time) == PR.Infinity ? 100 : problem.time; #max num of steps
x_traj = zeros(n_x,K+1);
u_traj = zeros(n_u,K+1);
x_traj[:,1] = x0;
state_traj = []

k = 1; #iterator
costBound = 0;
costTrue = 0;
currState = SY.get_all_states_by_xpos(symmodel,DO.crop_to_domain(domainX,DO.get_all_pos_by_coord(state_grid,x_traj[:,k])))
push!(state_traj,currState);

Xinit = DO.DomainList(state_grid) # set of initial cells
DO.add_coord!(Xinit, problem.initial_set)

Xfinal = DO.DomainList(state_grid) # set of target cells
DO.add_coord!(Xfinal, Vector(problem.target_set))

Xobstacles = DO.DomainList(state_grid) # set of obstacle cells
for o in system.ext[:obstacles]
      DO.add_set!(Xobstacles, o, DO.OUTER) 
end

initlist = [SY.get_state_by_xpos(symmodel, pos) for pos in DO.enum_pos(Xinit)]; 
finallist = [SY.get_state_by_xpos(symmodel, pos) for pos in DO.enum_pos(Xfinal)];
obstaclelist = [SY.get_state_by_xpos(symmodel, pos) for pos in DO.enum_pos(Xobstacles)];

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
      println("cm: $(DO.get_coord_by_pos(state_grid, SY.get_xpos_by_state(symmodel, next_action[2])))")
      
      c = DO.get_coord_by_pos(state_grid, SY.get_xpos_by_state(symmodel, next_action[1]))
      println("c: $(c)")


      u_traj[:,k] = transitionKappa[next_action]*vcat(x_traj[:,k]-c,1.0)
      println("x: $(x_traj[:,k])")
      println("u: $(u_traj[:,k])")

      global costBound = costBound + transitionCost[next_action]
      xk_aug = vcat(x_traj[:,k], u_traj[:,k],1.0)
      global costTrue += xk_aug'Q_aug*xk_aug


      m = ϕ(c)
      
      w = (2*(rand(2).^(1/4)).-1).*W[:,1]

      x_traj[:,k+1] = system.resetmaps[m].A*x_traj[:,k]+system.resetmaps[m].B*u_traj[:,k] + system.resetmaps[m].c + w

      global k += 1;
      global currState =  SY.get_all_states_by_xpos(symmodel,DO.crop_to_domain(domainX,DO.get_all_pos_by_coord(state_grid,x_traj[:,k])));
      push!(state_traj,currState)
      println("Arrived in: $(currState): $(x_traj[:,k])")
end

println("Goal set reached")
println("Guaranteed cost:\t $(costBound)")
println("True cost:\t\t $(costTrue)")

# Finally let us visualize the results. Let us plot the transitions in the state-feedback abstraction
using Plots

fig = plot(aspect_ratio=:equal, xtickfontsize=10, ytickfontsize=10, guidefontsize=16, titlefontsize=14);
xlims!(rectX.lb[1]-0.2,rectX.ub[1]+0.2)
ylims!(rectX.lb[2]-0.2,rectX.ub[2]+0.2)
dims = [1, 2]
plot!(domainX; dims=dims, color=:yellow, opacity=0.2);
for t in symmodel.autom.transitions.data
      if isempty(t ∩ obstaclelist) 
            color = RGB(abs(0.6*sin(t[1])), abs(0.6*sin(t[1]+2π/3)), abs(0.6*sin(t[1]-2π/3)))
            if t[1]==t[2]
                  p1 = DO.get_coord_by_pos(state_grid, SY.get_xpos_by_state(symmodel, t[2]))
                  plot!(fig, UT.DrawPoint(p1), color = color);
            else
                  p1 = DO.get_coord_by_pos(state_grid, SY.get_xpos_by_state(symmodel, t[2]))
                  p2 = DO.get_coord_by_pos(state_grid, SY.get_xpos_by_state(symmodel, t[1]))
                  plot!(fig, UT.DrawArrow(p1, p2), color = color);
            end
      end
end
xlabel!("\$x_1\$")
ylabel!("\$x_2\$")
title!("Transitions")

# Then we plot a colormap with the Lyapunov-like function $v(x)$ and the simulated trajectory in blue.

fig = plot(aspect_ratio=:equal, xtickfontsize=10, ytickfontsize=10, guidefontsize=16, titlefontsize=14);
xlims!(rectX.lb[1]-0.2,rectX.ub[1]+0.2)
ylims!(rectX.lb[2]-0.2,rectX.ub[2]+0.2)
dims = [1, 2]

LyapMax = max(filter(isfinite,getfield.([lyap_fun...],:second))...)
colormap = Colors.colormap("Blues")
mycolorMap = UT.Colormap([0.0,LyapMax], colormap)
cost_ordered = reverse(sort(hcat([ (lyap,state) for (state,lyap) in lyap_fun]...), dims=2))
for (lyap, state) in cost_ordered
      pos = SY.get_xpos_by_state(symmodel, state)
      elli = DO.get_elem_by_pos(state_grid, pos)
      if (lyap ≠ Inf)
            plot!(fig, elli, color = UT.get_color(mycolorMap, lyap));
      else
            plot!(fig, elli, color = :yellow);
      end
end
Einit = UT.Ellipsoid(collect(P), collect(problem.initial_set))
Etarget = UT.Ellipsoid(collect(P), collect(problem.target_set))

plot!(fig, Einit, color=:green);
plot!(fig, Etarget, color=:red);
plot!(Xobstacles, color=:black, opacity=1.0);

trajCoord = [[x_traj[1,i], x_traj[2,i]] for i in 1:k]
plot!(fig, UT.DrawTrajectory(trajCoord));
plot!(mycolorMap);
xlabel!("\$x_1\$")
ylabel!("\$x_2\$")
title!("Trajectory and Lyapunov-like Fun.")

# We recall that, to speed up the build time of this documentation, some values were modified in comparison with [1, Example 2]. To obtain the same figures use `Usz = 50`, `Wsz = 5` and `n_step = 5`.

@test contr.data[1] == (45,21)            #src
@test costBound ≈ 3.230 rtol=1e-3         #src
@test costTrue <= costBound               #src

# ## References
#
# 1. L. N. Egidio, T. Alves Lima, R. M. Jungers, "State-feedback Abstractions for Optimal Control of Piecewise-affine Systems", IEEE 61st Conference on Decision and Control (CDC), 2022, accepted.
# 1. D. Bertsekas, "Dynamic programming and optimal control". Volume I, Athena scientific, 2012.
