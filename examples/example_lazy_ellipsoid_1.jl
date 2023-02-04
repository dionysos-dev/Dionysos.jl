# include("../src/Dionysos.jl")
# using .Dionysos
using Dionysos
UT = Dionysos.Utils
SY = Dionysos.System
SC = Dionysos.Symbolic
PR = Dionysos.Problem

using Plots, Colors, LinearAlgebra
using Symbolics
using IntervalArithmetic
using LinearAlgebra
using Mosek
using MosekTools
using JuMP
import Random
Random.seed!(0)

#######################################################
################## Problem setting ####################
#######################################################

# System
function unstableSimple()
    Symbolics.@variables px py vx vy wx wy T

    f = [1.1*px-0.2*py-0.00005*py^3+T*vx;
         1.1*py+0.2*px+0.00005*px^3+T*vy];
    
    x = [px; py] # state
    u = [vx; vy] # control
    w = [wx; wy]
    return f, x, u, w, T
end
f, x, u, w, T = unstableSimple()
n_x = length(x)
n_u = length(u)
n_w = length(w)

# Bounds on u
Usz = 10
Uaux = diagm(1:n_u)
Ub = [(Uaux.==i)./Usz for i in 1:n_u];
# Box bounding x, u, w
X = IntervalBox(-15..15,2);
U = IntervalBox(-Usz..Usz,n_u)

##########################
Ts = 1
f_eval = eval(build_function(f,x,u,w,T)[1])
# augmented argument
xi = [x;u;w]
fT = Symbolics.substitute(f,Dict([T => Ts]))
# Boxes on which J must be bounded
maxStep = 5

maxRadius = 1
maxΔu = Usz*2
ΔX = IntervalBox(-maxRadius..maxRadius,2) #× (-0.01..0.01);
ΔU = IntervalBox(-maxΔu..maxΔu,2)
ΔW = IntervalBox(-0.0..0.0,1)
##########################


# Cost function
S = Matrix{Float64}(I(n_x+n_u+1)) #TODO

# Control problem
xinit = [-10.0;-10.0]
Einit = UT.Ellipsoid(Matrix{Float64}(I(n_x))*10.0, xinit)
xtarget = [10.0;10.0]
Etarget = UT.Ellipsoid(Matrix{Float64}(I(n_x))*1.0, xtarget)
obstacles = [UT.Ellipsoid(Matrix{Float64}(I(n_x))*1/50, [0.0;0.0])]

system = SY.EllipsoidalLazySystem(f_eval, X, U, Ub, Ts, fT, x, u, w, maxRadius, maxΔu, ΔX, ΔU, ΔW, obstacles)
problem = PR.OptimalControlProblem(system, Einit, Etarget, S, nothing, 0.0)

#######################################################
## Concrete implementation of generic ellipsoial RRT ##
#######################################################

# SI, SF
function distance(E1,E2)
    return UT.pointCenterDistance(E1, E2.c)
end

function get_candidate(tree, X::IntervalBox, E0; probSkew=0.0, probE0=0.05, intialDist=1)
    guess = SC.sample_box(X)
    randVal = rand()
    if randVal>probSkew+probE0
        return guess
    elseif randVal>probSkew
        return E0.c
    else 
        closestNode, dist  = UT.findNClosestNode(tree, UT.Ellipsoid(Matrix{Float64}(I(n_x)), E0))        
        l = randVal/probSkew
        r = dist/intialDist
        return (E0.c*l + closestNode.state.c*(1-l))*(1-0.3*r) +(0.3*r)*guess #heuristic bias
    end
end

# tree, SI, SF, distance, data
function rand_state(tree, EF, EI, distance, problem)
    X = problem.system.X
    xrand = get_candidate(tree, X, EI)
    return UT.Ellipsoid(Matrix{Float64}(I(length(xrand))), xrand)
end

# data-driven technique on nominal system (without noise)
function get_closest_reachable_point(f_eval, Ts, xinit, xtarget, U; nSamples=500)
    wnew = zeros(n_w)

    unew = SC.sample_box(U) #(0.8+0.2*norm(xPar-X0.c)/intialDist)
    xnew = f_eval(xinit, unew, wnew, Ts)
    uBestDist = norm(xnew-xtarget)
    for i in 1:nSamples
        ucandnew = SC.sample_box(U)*0.002*i
        xcandnew = f_eval(xinit, ucandnew, wnew, Ts)
        if norm(xcandnew-xtarget)< uBestDist
            uBestDist = norm(xcandnew-xtarget)
            xnew = xcandnew
            unew = ucandnew
        end
    end
    return (unew, xnew, uBestDist)
end

# tree, Nnear, Srand, data
function new_conf(tree, Nnear, Erand, problem)
    sys = problem.system
    (unew, xnew, uBestDist) = get_closest_reachable_point(sys.f_eval, -sys.Ts, Nnear.state.c, Erand.c, sys.U)
    n_w = 2
    wnew = zeros(n_w)
    X̄ = IntervalBox(xnew .+ ΔX)
    Ū = IntervalBox(unew .+ ΔU)
    W̄ = IntervalBox(wnew .+ ΔW)
    (affineSys, L) = Dionysos.System.buildAffineApproximation(sys.fT,sys.x,sys.u,sys.w,xnew,unew,wnew,X̄,Ū,W̄)
    S = problem.state_cost
    return Dionysos.Symbolic.hasTransition(xnew, unew, Nnear.state, affineSys, L, S, sys.Ub, sys.maxRadius, sys.maxΔu, sdp_opt; λ=0.015) 
end
#c, u, Ep::UT.Ellipsoid,subsys::AffineSys,L, S, U, maxRadius, maxΔu, optimizer; λ=0.01
# heuristic: keep only the closest ellipsoid to the initial ellipsoid
# tree, LSACnew, SI, SF, distance, data
function keep(tree, LSACnew, EF, EI, distance, problem; scale_for_obstacle=true)
    obstacles = problem.system.obstacles
    minDist = Inf
    iMin = 0
    for (i,data) in enumerate(LSACnew)
        Enew, kappa, cost, Nnear = data
        if Enew===nothing             
            print("\tInfeasible")
        elseif EI ∈ Enew
            iMin = i
            break
        elseif minDist > norm(EI.c-Enew.c) # minPathCost > cost + Eclosest.path_cost
            if Nnear==tree.root || eigmin(EI.P*0.5-Enew.P)>0 # E ⊂ E0 => P-P0>0
                iMin = i
                minDist = norm(EI.c-Enew.c)
            else
            end
        else
        end 
    end
    if iMin == 0
        return [] 
    end
    ElMin, kappaMin, costMin, NnearMin = LSACnew[iMin]
    if ElMin !== nothing 
        if all(O -> !(ElMin ∩ O), obstacles) #all(O -> !(ElMin.c ∈ O), obstacles)
            return [LSACnew[iMin]]
        elseif scale_for_obstacle
            for O in obstacles
                ElMin = UT.compress_if_intersection(ElMin, O)
                if ElMin === nothing
                    return []
                end
            end
            return [(ElMin, kappaMin, costMin, NnearMin)]
        else
            return []
        end
    else
        return []
    end
end

#state1, state2, data
function compute_transition(E1, E2, problem)
    E1, E2, U, Ub, S, sdp_opt, fT, x, u, w
    sys = problem.system
    xnew = E1.c
    unew = [0.0;0.0] #SC.sample_box(U)
    wnew = zeros(n_w)
    X̄ = IntervalBox(xnew .+ ΔX)
    Ū = IntervalBox(unew .+ ΔU)
    W̄ = IntervalBox(wnew .+ ΔW)
    (affineSys, L) = Dionysos.System.buildAffineApproximation(sys.fT,sys.x,sys.u,sys.w,xnew,unew,wnew,X̄,Ū,W̄)
    W = 0.0*[-1 -1  1 1;
         -1  1 -1 1]
    S = problem.state_cost
    ans, cost, kappa = Dionysos.Symbolic.my_has_transition(affineSys.A, affineSys.B, affineSys.c, sys.W, sys.Ub, S, E1, E2, sdp_opt)
    return ans, cost, kappa
end

# tree, LNnew, SI, SF, distance, data
function stop_crit(tree, LNnew, EF, EI, distance, problem; continues=false)
    newEllipsoids = [newNode.state for newNode in LNnew]
    return any(map(E->(EI ∈ E), newEllipsoids)) || continues
end

sdp_opt =  optimizer_with_attributes(Mosek.Optimizer, MOI.Silent() => true)
maxIter = 100
RRTstar = false
optimizer = PR.Abstraction.build_OptimizerLazyEllipsoids(problem, distance, rand_state, new_conf, keep, stop_crit, RRTstar, compute_transition, maxIter)

MOI.optimize!(optimizer)

tree = optimizer.tree

p = plot(aspect_ratio=:equal)
for obs in obstacles
    UT.plotE!(obs,color=:black)
end

UT.plot_Tree!(tree)
UT.plotE!(Einit, color=:green)
UT.plotE!(Etarget, color=:red)

# x = Einit.c #[3.0;-8.0] #E0.c
# trajx, trajE = SC.simulate(tree, f_eval, Ts, x) 
# SC.plot_traj!(trajx, trajE, color=:green)
display(p)

# c = [-10.0;-10.0]
# P = [2.0 6.0; 6.0 20.0]
# E = UT.Ellipsoid(P, c)
# box = UT.get_min_bounding_box(E)
# p = plot(aspect_ratio=:equal)
# SC.plot_box!(box)
# UT.plotE!(E)
# display(p)




