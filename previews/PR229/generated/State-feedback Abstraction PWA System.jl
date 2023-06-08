using StaticArrays
using LinearAlgebra
using CDDLib
using SDPA, Ipopt, JuMP

opt_sdp = optimizer_with_attributes(SDPA.Optimizer, MOI.Silent() => true)
opt_ip = optimizer_with_attributes(Ipopt.Optimizer, MOI.Silent() => true)
lib = CDDLib.Library() #polyhedron lib

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

using Dionysos

const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const SY = DI.Symbolic
const CO = DI.Control
const PR = DI.Problem
const OP = DI.Optim
const AB = OP.Abstraction

max_x = 2 # bound on |X|_∞
rectX = UT.HyperRectangle(SVector(-max_x, -max_x), SVector(max_x, max_x));
rectU = rectX

n_step = 2 #  n_step = 5 in [1]
X_origin = SVector(0.0, 0.0);
X_step = SVector(1.0/n_step, 1.0/n_step);
P = (1/n_x)*diagm((X_step./2).^(-2))

state_grid = DO.GridEllipsoidalRectangular(X_origin, X_step, P, rectX)

using JuMP
optimizer = MOI.instantiate(AB.EllipsoidsAbstractions.Optimizer) #

MOI.set(optimizer, MOI.RawOptimizerAttribute("problem"), problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("sdp_solver"), opt_sdp)
MOI.set(optimizer, MOI.RawOptimizerAttribute("ip_solver"), opt_ip)

using Suppressor
@suppress  begin # this is a workaround to supress the undesired output of SDPA
      MOI.optimize!(optimizer)
end
contr = MOI.get(optimizer, MOI.RawOptimizerAttribute("controller"))
symmodel = MOI.get(optimizer, MOI.RawOptimizerAttribute("symmodel"))
transitionCost = MOI.get(optimizer, MOI.RawOptimizerAttribute("transitionCost"))
transitionKappa = MOI.get(optimizer, MOI.RawOptimizerAttribute("transitionKappa"))
lyap_fun = MOI.get(optimizer, MOI.RawOptimizerAttribute("lyap_fun"));

Q_aug = CO.get_full_psd_matrix(problem.transition_cost[1][1])

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

using Plots

fig = plot(aspect_ratio=:equal, xtickfontsize=10, ytickfontsize=10, guidefontsize=16, titlefontsize=14)
xlims!(rectX.lb[1]-0.2,rectX.ub[1]+0.2)
ylims!(rectX.lb[2]-0.2,rectX.ub[2]+0.2)
dims = [1, 2]
plot!(domainX; dims=dims, color=:yellow, opacity=0.2)
for t in symmodel.autom.transitions.data
      if isempty(t ∩ obstaclelist)
            color = RGB(abs(0.6*sin(t[1])), abs(0.6*sin(t[1]+2π/3)), abs(0.6*sin(t[1]-2π/3)))
            if t[1]==t[2]
                  p1 = DO.get_coord_by_pos(state_grid, SY.get_xpos_by_state(symmodel, t[2]))
                  plot!(fig, UT.DrawPoint(p1), color = color)
            else
                  p1 = DO.get_coord_by_pos(state_grid, SY.get_xpos_by_state(symmodel, t[2]))
                  p2 = DO.get_coord_by_pos(state_grid, SY.get_xpos_by_state(symmodel, t[1]))
                  plot!(fig, UT.DrawArrow(p1, p2), color = color)
            end
      end
end
xlabel!("\$x_1\$")
ylabel!("\$x_2\$")
title!("Transitions")
display(fig)

fig = plot(aspect_ratio=:equal, xtickfontsize=10, ytickfontsize=10, guidefontsize=16, titlefontsize=14)
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
            plot!(fig, elli, color = UT.get_color(mycolorMap, lyap))
      else
            plot!(fig, elli, color = :yellow)
      end
end
Einit = UT.Ellipsoid(collect(P), collect(problem.initial_set))
Etarget = UT.Ellipsoid(collect(P), collect(problem.target_set))

plot!(fig, Einit, color=:green)
plot!(fig, Etarget, color=:red)
plot!(Xobstacles, color=:black, opacity=1.0)

trajCoord = [[x_traj[1,i], x_traj[2,i]] for i in 1:k]
plot!(fig, UT.DrawTrajectory(trajCoord))
plot!(mycolorMap)
xlabel!("\$x_1\$")
ylabel!("\$x_2\$")
title!("Trajectory and Lyapunov-like Fun.")
display(fig)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

