include("../src/Dionysos.jl")
using .Dionysos
UT = Dionysos.Utils
SY = Dionysos.System
SC = Dionysos.Symbolic
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
x0 = [-10.0;-10.0]
E0 = UT.Ellipsoid(Matrix{Float64}(I(n_x))*10.0, x0)
xf = [10.0;10.0]
EF = UT.Ellipsoid(Matrix{Float64}(I(n_x))*1.0, xf)
obstacles = [UT.Ellipsoid(Matrix{Float64}(I(n_x))*1/50, [0.0;0.0])]

#######################################################
## Concrete implementation of generic ellipsoial RRT ##
#######################################################

function distance(E1,E2)
    return UT.pointCenterDistance(E1, E2.c)
end

function get_candidate(rrt, X::IntervalBox, E0; probSkew=0.0, probE0=0.05, intialDist=1)
    guess = SC.sample_box(X)
    randVal = rand()
    if randVal>probSkew+probE0
        return guess
    elseif randVal>probSkew
        return E0.c
    else 
        closestNode, dist  = UT.findNClosestNode(rrt, UT.Ellipsoid(Matrix{Float64}(I(n_x)), E0))        
        l = randVal/probSkew
        r = dist/intialDist
        return (E0.c*l + closestNode.state.c*(1-l))*(1-0.3*r) +(0.3*r)*guess #heuristic bias
    end
end

function get_random_state(rrt,X,E0)
    xrand = get_candidate(rrt,X,E0)
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

function get_new_state(f_eval, Ts, Eclosest, Erand, U, S, Ub, maxRadius, maxΔu, sdp_opt, fT, x, u, w)
    (unew, xnew, uBestDist) = get_closest_reachable_point(f_eval, -Ts, Eclosest.state.c, Erand.c, U)
    wnew = zeros(n_w)
    X̄ = IntervalBox(xnew .+ ΔX)
    Ū = IntervalBox(unew .+ ΔU)
    W̄ = IntervalBox(wnew .+ ΔW)
    (sys, L) = Dionysos.System.buildAffineApproximation(fT,x,u,w,xnew,unew,wnew,X̄,Ū,W̄)
    return Dionysos.Symbolic.hasTransition(xnew, unew, Eclosest.state, sys, L, S, Ub, maxRadius, maxΔu, sdp_opt; λ=0.015) 
end

# heuristic: keep only the closest ellipsoid to the initial ellipsoid
function keep(rrt, newStates, E0, obstacles;scale_for_obstacle=true) 
    minDist = Inf
    iMin = 0
    for (i,data) in enumerate(newStates)
        Enew, kappa, cost, Eclosest = data
        if Enew===nothing             
            print("\tInfeasible")
        elseif E0 ∈ Enew
            iMin = i
            break
        elseif minDist > norm(E0.c-Enew.c) # minPathCost > cost + Eclosest.path_cost
            if Eclosest==rrt.treeRoot || eigmin(E0.P*0.5-Enew.P)>0 # E ⊂ E0 => P-P0>0
                iMin = i
                minDist = norm(E0.c-Enew.c)
            else
            end
        else
        end 
    end
    if iMin == 0
        return [] 
    end
    ElMin, kappaMin, costMin, EclosestMin = newStates[iMin]
    if ElMin !== nothing 
        if all(O -> !(ElMin ∩ O), obstacles) #all(O -> !(ElMin.c ∈ O), obstacles)
            return [newStates[iMin]]
        elseif scale_for_obstacle
            for O in obstacles
                ElMin = UT.compress_if_intersection(ElMin, O)
                if ElMin === nothing
                    return []
                end
            end
            return [(ElMin, kappaMin, costMin, EclosestMin)]
        else
            return []
        end
    else
        return []
    end
end


function compute_transition(E1, E2, U, Ub, S, sdp_opt, fT, x, u, w)
    xnew = E1.c
    unew = [0.0;0.0] #SC.sample_box(U)
    wnew = zeros(n_w)
    X̄ = IntervalBox(xnew .+ ΔX)
    Ū = IntervalBox(unew .+ ΔU)
    W̄ = IntervalBox(wnew .+ ΔW)
    (sys, L) = Dionysos.System.buildAffineApproximation(fT,x,u,w,xnew,unew,wnew,X̄,Ū,W̄)
    W = 0.0*[-1 -1  1 1;
         -1  1 -1 1]
    ans, cost, kappa = Dionysos.Symbolic.my_has_transition(sys.A, sys.B, sys.c, W, Ub, S, E1, E2, sdp_opt)
    return ans, cost, kappa
end


sdp_opt =  optimizer_with_attributes(Mosek.Optimizer, MOI.Silent() => true)
maxIter = 100
RRTstar = false
continues = true
rrt = SC.build_lazy_ellipsoidal_abstraction(E0, EF, obstacles, S, f_eval, X, U, Ub, Ts, fT, x, u, w, maxRadius, maxΔu, ΔX, ΔU, ΔW, sdp_opt, distance, get_random_state, get_new_state, keep;maxIter=maxIter,RRTstar=RRTstar,compute_transition=compute_transition,continues=continues)

p = plot(aspect_ratio=:equal)
for obs in obstacles
    UT.plotE!(obs,color=:black)
end
UT.plot_RRT!(rrt)
UT.plotE!(E0,color=:green)
UT.plotE!(EF,color=:red)

# x = E0.c #[3.0;-8.0] #E0.c
# trajx, trajE = SC.simulate(rrt, f_eval, Ts, x) 
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




