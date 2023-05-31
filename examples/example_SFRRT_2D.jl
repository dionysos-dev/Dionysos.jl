# Warning : deprecated example, see https://github.com/dionysos-dev/Dionysos.jl/issues/221

using Dionysos

const SY = Dionysos.System
const UT = Dionysos.Utils

using Symbolics
using IntervalArithmetic
using LinearAlgebra
using Mosek
using MosekTools

import Random
Random.seed!(0)

using Plots
#Def sym variables


# discrete-time dynamic function
#T = 0.1
# f = [px+2*v/ω*sin(T*ω/2)*cos(th+T*ω/2);
#      py+2*v/ω*sin(T*ω/2)*sin(th+T*ω/2);
#      th+T*ω
#     ];
function unicycleRobot()
    # due to "Global Observability Analysis of a Nonholonomic Robot using Range Sensors"
    Symbolics.@variables px py th v ω w1 T

    # sinc not implemented in symbolic and division by `th` makes IntervalArithmetic bug
    mysinc(th) = sin(sqrt(th^2+1e-3))/sqrt(th^2+1e-3)
    f = [px+T*v*mysinc(T*ω/2)*cos(th+T*ω/2);
    py+T*v*mysinc(T*ω/2)*sin(th+T*ω/2);
    th+T*ω
    ];
    
    x = [px; py; th] # state
    u = [v; ω] # control
    w = [w1]
    
    return f, x, u, w, T    
end


function pathPlaning()
    Symbolics.@variables px py vx vy wx wy T

    f = [px+T*vx;
         py+T*vy];
    
    x = [px; py] # state
    u = [vx; vy] # control
    w = [wx; wy]
    return f, x, u, w, T
end

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

# sys dimensions
n_x = length(x)
n_u = length(u)
n_w = length(w)

Ts = 1
# sys eval function
# function f_eval(x̄,ū,w̄,T̄)
#     rules = Dict([x[i] => x̄[i] for i=1:n_x] ∪ 
#     [u[i] => ū[i] for i=1:n_u] ∪ 
#     [w[i] => w̄[i] for i=1:n_w] ∪
#     [T => T̄])
    
#     ff = Symbolics.substitute(f,rules)
#     Base.invokelatest(eval(eval(build_function(ff)[1])))
# end
f_eval = eval(build_function(f,x,u,w,T)[1])


# augmented argument
xi = [x;u;w]

fT = Symbolics.substitute(f,Dict([T => Ts]))

# Sym Jacobian
#Jxi_s = Symbolics.jacobian(fT,xi)



# Bounds on u
Usz = 10
Uaux = diagm(1:n_u)
Ub = [(Uaux.==i)./Usz for i in 1:n_u];
# Box bounding x, u, w
X = IntervalBox(-15..15,2);
U = IntervalBox(-Usz..Usz,n_u)

# Boxes on which J must be bounded
maxStep = 5

maxRadius = 1
maxΔu = Usz*2
ΔX = IntervalBox(-maxRadius..maxRadius,2) #× (-0.01..0.01);
ΔU = IntervalBox(-maxΔu..maxΔu,2)

ΔW = IntervalBox(-0.0..0.0,1)


# Cost
S = Matrix{Float64}(I(n_x+n_u+1)) #TODO


using JuMP
sdp_opt =  optimizer_with_attributes(Mosek.Optimizer, MOI.Silent() => true)



x0 = [-10.0;-10.0]
X0 = UT.Ellipsoid(Matrix{Float64}(I(n_x))*10.0, x0)
xf = [10.0;10.0]
XF = UT.Ellipsoid(Matrix{Float64}(I(n_x))*1.0, xf)
obstacles = [UT.Ellipsoid(Matrix{Float64}(I(n_x))*1/50, [0.0;0.0])]
intialDist = norm(x0-xf) 
treeRoot = UT.Node(XF)
treeLeaves = [treeRoot]

function findNClosestNode(nodeList, x; N=1) #TODO N>1
    dists = map(e-> e===nothing ? Inf : UT.pointCenterDistance(e.state, x), nodeList)
    d, idx = findmin(dists)
    parents = filter(x -> x!==nothing, unique(map(x-> x.parent,nodeList)))
    if !isempty(parents)
        bestPar, dBestPar = findNClosestNode(parents, x) 
        if d<dBestPar
            return nodeList[idx], d
        else
            return bestPar, dBestPar
        end
    else
        return nodeList[idx], d
    end
end


function findCloseNodes(nodeList, x; d=1) 
    filterFun = e-> e===nothing ? false : UT.pointCenterDistance(e.state, x)<=d
    closeNodes = filter(filterFun, nodeList)
    parents = filter(x -> x!==nothing, unique(map(x-> x.parent,nodeList)))
    if !isempty(parents)
        parCloseNodes = findCloseNodes(parents, x; d=d) 
        return unique(parCloseNodes ∪ closeNodes)

    else
        return closeNodes
    end
end
function sample_x(;probSkew=0.0, probX0=0.05)
    guess = map(x-> x.lo + (x.hi-x.lo)*rand(), X.v)
    randVal = rand()
    if randVal>probSkew+probX0
        # println("random guess")
        return guess
    elseif randVal>probSkew
        # println("X0.c guess")
        return X0.c
    else 
        closestNode, dist  = findNClosestNode(treeLeaves,X0.c)
        
        l = randVal/probSkew
        # println("skewed guess")
        r = dist/intialDist
        return (X0.c*l + closestNode.state.c*(1-l))*(1-0.3*r) +(0.3*r)*guess #heuristic bias
    end
end
sample_u() = map(x-> (x.lo + (x.hi-x.lo)*rand()), U.v)



global maxIter = 10000
global Xnew = XF
global bestDist = UT.centerDistance(X0,Xnew)
while !(X0 ∈ Xnew) && maxIter>0
    xsample = Vector(sample_x())
    #Xclosest_, _ = findNClosestNode(treeLeaves,xsample)
    #closeNodes = [Xclosest_]
    
    closeNodes= first(sort(findCloseNodes(treeLeaves,xsample; d=Inf); lt=((e1,e2)-> UT.pointCenterDistance(e1.state, xsample)<UT.pointCenterDistance(e2.state, xsample))), 5)
    # print(">>>>>>>>to: \t")
    # println(xAim)
    if isempty(closeNodes)
        continue
    end
    print("Iterations2Go:\t")
    println(maxIter)
    minPathCost = Inf
    minDist = Inf
    ElMin = nothing
    kappaMin = nothing
    parentMin = nothing
    for Xclosest in closeNodes
        xPar = UT.get_center(Xclosest.state)
        if norm(xsample-xPar)>maxStep
            xAim = maxStep*(xsample-xPar)./norm(xsample-xPar)+xPar
        else
            xAim = xsample
        end
        # println(xsample)
    
        wnew = zeros(n_w)
        # unew = finv(xAim,xPar)
        # xnew = f_eval(xPar, unew, wnew,-Ts)
        
        # print(xPar)
        # print(xAim)
        unew = sample_u()*1#(0.8+0.2*norm(xPar-X0.c)/intialDist)
        xnew = f_eval(xPar, unew, wnew, -Ts)
        
        # println(xnew)
        uBestDist = norm(xnew-xAim)
        for i in 1:500
            ucandnew = sample_u()*0.002*i
            xcandnew = f_eval(xPar, ucandnew, wnew, -Ts)
            if norm(xcandnew-xAim)< uBestDist
                uBestDist = norm(xcandnew-xAim)
                xnew = xcandnew
                unew = ucandnew
            end
        end
        # println(xnew)
        # println(xnew)
        # println(unew)
        # println()
        # print(">>>>>> from: \t")
        # println(xPar)
        # print(">>>>>>>> to: \t")
        # println(xnew)
        #xaux = f_eval(xPar, unew, wnew, -Ts)
        #unew = zeros(n_u)
        #wnew = zeros(n_w)
        #xaux = xsample-xPar
        #xaux = 2*maxRadius*xaux./norm(xaux) 
        #xnew = xPar + xaux
        # println(xnew)

        X̄ = IntervalBox(xnew .+ ΔX)
        Ū = IntervalBox(unew .+ ΔU)
        W̄ = IntervalBox(wnew .+ ΔW)
        (sys, L) = Dionysos.System.buildAffineApproximation(fT,x,u,w,xnew,unew,wnew,X̄,Ū,W̄)
        # println(L)
        El, kappa, cost = Dionysos.Symbolic.hasTransition(xnew, unew, Xclosest.state, sys, L, S, Ub, maxRadius, maxΔu, sdp_opt; λ=0.015)
        if El===nothing             
            # print("\tInfeasible")
        elseif X0 ∈ El
            ElMin = El
            kappaMin = kappa
            minDist = norm(X0.c-El.c)
            minPathCost = Xclosest.path_cost+cost
            parentMin = Xclosest
            break
        elseif minDist > norm(X0.c-El.c) # minPathCost > cost + Xclosest.path_cost
            # print("\tFeasible")
            if Xclosest==treeRoot || eigmin(X0.P*0.5-El.P)>0 # E ⊂ E0 => P-P0>0
                ElMin = El
                kappaMin = kappa
                minDist = norm(X0.c-El.c)
                minPathCost = Xclosest.path_cost+cost
                parentMin = Xclosest
            else
                # print("\tRejecting...")
            end
        else
            # print("\tNotTheBest")
        end 
        # println()
    end

    if ElMin !== nothing && all(O -> !(ElMin.c ∈ O), obstacles)
        # if sqrt(eigmin(ElMin.P))<1/(maxRadius)
        #     println("STOP!")
        #     break
        # end
        global Xnew = ElMin
        push!(treeLeaves, UT.Node(Xnew; parent=parentMin, action=kappaMin, path_cost=minPathCost) )
        setdiff!(treeLeaves, [parentMin])
        dAux = UT.centerDistance(X0,Xnew)
        print("\t ")
        print(norm(X0.c-ElMin.c))
        print("\t ")
        println(UT.get_volume(ElMin))
        if bestDist > dAux
            global bestDist = dAux
            print("******")
        end
    end
    print("\tClosest Dist: ")
    println(bestDist)
    global maxIter-=1
end
global xSpan = [x0]
global uSpan = []
global xk = x0
global k = 1
if X0 ∈ Xnew
    global currNode = last(treeLeaves)
    while !(xk ∈ XF)
        global xk
        println(k)
        println(xk)
        while (xk ∈ currNode.parent.state)
            global currNode = currNode.parent
        end
        if !(xk ∈ currNode.state)
            println("ERROR")
            break
        end
        uk = currNode.action*[xk-currNode.state.c;1]
        xk = f_eval(xk,uk,zeros(n_w),Ts)
        push!(xSpan, xk)
        push!(uSpan, uk)
        global k+=1
    end
    println("OK!")
else
    println("KO")
    return
end


if true
    @static if get(ENV, "CI", "false") == "false" && (isdefined(@__MODULE__, :no_plot) && no_plot==false)
        fig = plot(aspect_ratio=:equal, xtickfontsize=10, ytickfontsize=10, guidefontsize=16, titlefontsize=14)

        xlims!(X[1].lo-0.2, X[1].hi+0.2)
        ylims!(X[2].lo-0.2, X[2].hi+0.2)

        LyapMax = (max(map(x-> x.path_cost, treeLeaves)...))
        colormap = Colors.colormap("Blues")
        mycolorMap = UT.Colormap([0.0,LyapMax], colormap)

        plotted = []
        currNode = last(treeLeaves)

        for n in treeLeaves
            push!(plotted, n)
            lyap = (n.path_cost)
            plot!(fig, n.state, color=UT.get_color(mycolorMap, lyap))
            nodePar = n.parent
            aTail = n.state.c
            while nodePar !== nothing 
                aHead = nodePar.state.c 
                if !(nodePar ∈ plotted)
                    lyap = (nodePar.path_cost)
                    plot!(fig, nodePar.state, color=UT.get_color(mycolorMap, lyap))
                    push!(plotted, nodePar)
                    plot!(fig, UT.DrawArrow(aTail, aHead))
                    aTail = nodePar.state.c
                    nodePar = nodePar.parent
                else
                    plot!(fig, UT.DrawArrow(aTail, aHead))
                    nodePar = nothing
                end
            end
        end
        xlabel!("\$x_1\$")
        ylabel!("\$x_2\$")
        title!("Trajectory and Lyapunov-like Fun.")

        while currNode !== nothing
            lyap = (currNode.path_cost)
            plot!(fig, currNode.state, color=UT.get_color(mycolorMap, lyap))
            global currNode = currNode.parent
        end
        plot!(fig, X0, color = :green)
        plot!(fig, obstacles[1], color = :black)
        plot!(fig, XF, color=:red)
        UT.plot_colorBar!(mycolorMap)

        # savefig("ex2_traj.png")          
        display(fig)
    end
 end #if 