using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const CO = DI.Control
const PR = DI.Problem
const OP = DI.Optim
const AB = OP.Abstraction

const LEA = AB.LazyEllipsoidsAbstraction

using LinearAlgebra, Random, Symbolics, IntervalArithmetic
using JuMP, Mosek, MosekTools
using Plots, Colors
Random.seed!(0)

#######################################################
################## Problem setting ####################
#######################################################

# convert square box into intersection of degenerated ellipsoids
# square box [-x,x]
function convert_H_into_E(Ub::IntervalBox)
    nu = length(Ub)
    Uaux = diagm(1:nu)
    U = [(Uaux .== i) ./ Ub[i].hi for i in 1:nu]
    return U
end

function unstableSimple()
    Symbolics.@variables px py vx vy wx wy T

    f = [
        1.1 * px - 0.2 * py - 0.00005 * py^3 + T * vx
        1.1 * py + 0.2 * px + 0.00005 * px^3 + T * vy
    ]

    x = [px; py] # state
    u = [vx; vy] # control
    w = [wx; wy]
    return f, x, u, w, T
end

########## Function description #########
f, x, u, w, T = unstableSimple()
nx = length(x)
nu = length(u)
nw = length(w)
Ts = 1
f_eval = eval(build_function(f, x, u, w, T)[1])

#### PWA approximation description #####
fsymbolic = Symbolics.substitute(f, Dict([T => Ts]))
ΔX = IntervalBox(-1.0 .. 1.0, 2)
ΔU = IntervalBox(-10 * 2 .. 10 * 2, 2)
ΔW = IntervalBox(-0.0 .. 0.0, 1)

########## Inputs description ##########
Usz = 10
Ub = IntervalBox(-Usz .. Usz, nu)
U = convert_H_into_E(Ub)

########## Noise description ##########
W = 0.0 * [
    -1 -1 1 1
    -1 1 -1 1
]

########### Cost description ############
S = Matrix{Float64}(I(nx + nu + 1)) #TODO

########## Control description ##########
xinit = [-10.0; -10.0]
Einit = UT.Ellipsoid(Matrix{Float64}(I(nx)) * 10.0, xinit)
xtarget = [10.0; 10.0]
Etarget = UT.Ellipsoid(Matrix{Float64}(I(nx)) * 1.0, xtarget)
X = IntervalBox(-20 .. 20, 2);
obstacles = [UT.Ellipsoid(Matrix{Float64}(I(nx)) * 1 / 50, [0.0; 0.0])]

#########################################
system = ST.EllipsoidalLazySystem(
    f_eval,
    Ts,
    nx,
    nu,
    nw,
    U,
    Ub,
    W,
    X,
    obstacles,
    fsymbolic,
    x,
    u,
    w,
    ΔX,
    ΔU,
    ΔW,
)
problem = PR.OptimalControlProblem(system, Einit, Etarget, S, nothing, 0.0)

#######################################################
## Concrete implementation of generic ellipsoial RRT ##
#######################################################

# SI, SF
function distance(E1, E2)
    return UT.centerDistance(E1, E2)
end

function get_candidate(
    tree,
    X::IntervalBox,
    E0;
    probSkew = 0.0,
    probE0 = 0.05,
    intialDist = 1,
)
    guess = UT.sample_box(X)
    randVal = rand()
    if randVal > probSkew + probE0
        return guess
    elseif randVal > probSkew
        return E0.c
    else
        closestNode, dist =
            UT.findNClosestNode(tree, UT.Ellipsoid(Matrix{Float64}(I(nx)), E0))
        l = randVal / probSkew
        r = dist / intialDist
        return (E0.c * l + closestNode.state.c * (1 - l)) * (1 - 0.3 * r) +
               (0.3 * r) * guess #heuristic bias
    end
end

# tree, SI, SF, distance, data
function rand_state(tree, EF, EI, distance, optimizer)
    problem = optimizer.problem
    xrand = get_candidate(tree, problem.system.X, EI)
    return UT.Ellipsoid(Matrix{Float64}(I(length(xrand))), xrand)
end

# data-driven technique on nominal system (without noise)
function get_closest_reachable_point(sys, xinit, xtarget, U, Ub; nSamples = 500)
    wnew = zeros(sys.nw)
    unew = UT.sample_box(Ub) #(0.8+0.2*norm(xPar-X0.c)/intialDist)
    xnew = sys.f_eval(xinit, unew, wnew, -sys.Ts)
    uBestDist = norm(xnew - xtarget)
    for i in 1:nSamples
        ucandnew = UT.sample_box(Ub) * 0.002 * i
        xcandnew = sys.f_eval(xinit, ucandnew, wnew, -sys.Ts)
        if norm(xcandnew - xtarget) < uBestDist
            uBestDist = norm(xcandnew - xtarget)
            xnew = xcandnew
            unew = ucandnew
        end
    end
    return (unew, xnew, uBestDist)
end

# tree, Nnear, Srand, data
function new_conf(tree, Nnear, Erand, optimizer)
    problem = optimizer.problem
    sys = problem.system
    (unew, xnew, uBestDist) =
        get_closest_reachable_point(sys, Nnear.state.c, Erand.c, sys.U, sys.Ub)
    wnew = zeros(sys.nw)
    X̄ = IntervalBox(xnew .+ sys.ΔX)
    Ū = IntervalBox(unew .+ sys.ΔU)
    W̄ = IntervalBox(wnew .+ sys.ΔW)
    (affineSys, L) = ST.buildAffineApproximation(
        sys.fsymbolic,
        sys.x,
        sys.u,
        sys.w,
        xnew,
        unew,
        wnew,
        X̄,
        Ū,
        W̄,
    )
    S = problem.state_cost
    return SY.transition_backward(
        affineSys,
        Nnear.state,
        xnew,
        unew,
        sys.U,
        S,
        L,
        optimizer.sdp_opt;
        λ = optimizer.λ,
        maxδx = optimizer.maxδx,
        maxδu = optimizer.maxδu,
    )
end

# heuristic: keep only the closest ellipsoid to the initial ellipsoid
# tree, LSACnew, SI, SF, distance, data
function keep(tree, LSACnew, EF, EI, distance, optimizer; scale_for_obstacle = true)
    problem = optimizer.problem
    obstacles = problem.system.obstacles
    minDist = Inf
    iMin = 0
    for (i, data) in enumerate(LSACnew)
        Enew, cont, cost, Nnear = data
        if Enew === nothing
            print("\tInfeasible")
        elseif EI ∈ Enew
            iMin = i
            break
        elseif minDist > norm(EI.c - Enew.c) # minPathCost > cost + Eclosest.path_cost
            if Nnear == tree.root || eigmin(EI.P * 0.5 - Enew.P) > 0 # E ⊂ E0 => P-P0>0
                iMin = i
                minDist = norm(EI.c - Enew.c)
            else
            end
        else
        end
    end
    if iMin == 0
        return []
    end
    ElMin, contMin, costMin, NnearMin = LSACnew[iMin]
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
            return [(ElMin, contMin, costMin, NnearMin)]
        else
            return []
        end
    else
        return []
    end
end

#state1, state2, data
function compute_transition(E1, E2, optimizer)
    problem = optimizer.problem
    sys = problem.system
    xnew = E1.c
    unew = zeros(sys.nu)
    wnew = zeros(sys.nw)
    X̄ = IntervalBox(xnew .+ sys.ΔX)
    Ū = IntervalBox(unew .+ sys.ΔU)
    W̄ = IntervalBox(wnew .+ sys.ΔW)
    (affineSys, L) = ST.buildAffineApproximation(
        sys.fsymbolic,
        sys.x,
        sys.u,
        sys.w,
        xnew,
        unew,
        wnew,
        X̄,
        Ū,
        W̄,
    )
    S = problem.state_cost
    ans, cont, cost =
        SY.transition_fixed(affineSys, E1, E2, sys.U, sys.W, S, optimizer.sdp_opt)
    return ans, cont, cost
end

global myBool = true
global myBool2 = false
global NI = nothing
# tree, LNnew, SI, SF, distance, data
function stop_crit(tree, LNnew, EF, EI, distance, optimizer; continues = false)
    problem = optimizer.problem
    minDist = 10.0
    for Nnew in LNnew
        E = Nnew.state
        if distance(EI, E) <= minDist
            ans, cont, cost = compute_transition(EI, E, optimizer)
            if ans
                if myBool
                    global NI = UT.add_node!(tree, EI, Nnew, cont, cost)
                    println("Path cost from EI : ", NI.path_cost)
                    global myBool = false
                end
                if !continues
                    return true
                else
                    if myBool2 && cost + Nnew.path_cost < NI.path_cost
                        UT.rewire(tree, NI, Nnew, cont, cost)
                        println("Path cost from EI : ", NI.path_cost)
                    end
                end
                global myBool2 = true
            end
        end
    end
    newEllipsoids = [newNode.state for newNode in LNnew]
    return any(map(E -> (EI ∈ E), newEllipsoids)) && !continues
end

sdp_opt = optimizer_with_attributes(Mosek.Optimizer, MOI.Silent() => true)
maxδx = 100 # 100
maxδu = Usz * 2
maxIter = 100 # 100
RRTstar = false
k1 = 1
k2 = 1
λ = 0.01 # 0.01
optimizer = LEA.build_OptimizerLazyEllipsoids(
    problem,
    distance,
    rand_state,
    new_conf,
    keep,
    stop_crit,
    RRTstar,
    compute_transition,
    maxIter,
    maxδx,
    maxδu,
    λ,
    sdp_opt,
    k1,
    k2,
)

MOI.optimize!(optimizer)

tree = optimizer.tree
println("Path cost from EI : ", NI.path_cost)
fig = plot(; aspect_ratio = :equal)
for obs in obstacles
    plot!(obs; color = :black)
end

plot!(tree; arrowsB = true, cost = true)
# plot!(NI; pathB=true, cost=true)
# plot!(Einit, color = :green)
# plot!(Etarget, color = :red)
display(fig)
