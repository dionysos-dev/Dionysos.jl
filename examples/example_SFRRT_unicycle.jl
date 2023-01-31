using Dionysos

# const PR = Dionysos.Problem
const SY = Dionysos.System
const UT = Dionysos.Utils

using Symbolics
using IntervalArithmetic
using LinearAlgebra
using Mosek
using MosekTools

import Random
Random.seed!(0)
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

    f = [1.3*px+2*py+0.005*py^3+T*vx;
         1.1*py-px-0.005*px^3+T*vy];
    
    x = [px; py] # state
    u = [vx; vy] # control
    w = [wx; wy]
    return f, x, u, w, T
end


f, x, u, w, T = unicycleRobot()
finv(x,x₊) =begin
    ω = (x₊[3]-x[3])/Ts
    [0.5*(x₊[1]-x[1])/(Ts*sinc(Ts*ω/2/pi)*cos(x[3]+Ts*ω/2))+0.5*(x₊[2]-x[2])/(Ts*sinc(Ts*ω/2/pi)*sin(x[3]+Ts*ω/2))
    ω]
end 
# sys dimensions
n_x = length(x)
n_u = length(u)
n_w = length(w)

Ts = 0.05
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
Usz = 100
Uaux = diagm(1:n_u)
Ub = [(Uaux.==i)./Usz for i in 1:n_u];
# Box bounding x, u, w
X = IntervalBox(-15..15,2) × (-pi..pi);
U = IntervalBox(-Usz..Usz,n_u)

# Boxes on which J must be bounded
maxStep = 2

maxRadius = 2
maxΔu = 2
ΔX = IntervalBox(-maxRadius..maxRadius,3) #× (-0.01..0.01);
ΔU = IntervalBox(-maxΔu..maxΔu,2)

ΔW = IntervalBox(-0.0..0.0,1)


# Cost
S = Matrix{Float64}(I(n_x+n_u+1)) #TODO


using JuMP
sdp_opt =  optimizer_with_attributes(Mosek.Optimizer, MOI.Silent() => true)



x0 = [-7.0;0.0;pi/2]
X0 = UT.Ellipsoid(Matrix{Float64}(I(n_x))*10000, x0)
xf = [7.0;0.0;-pi/2]
XF = UT.Ellipsoid(Matrix{Float64}(I(n_x))*1, xf)

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



global maxIter = 300
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
        El, kappa, cost = Dionysos.Symbolic.hasTransition(xnew, unew, Xclosest.state, sys, L, S, Ub, maxRadius, maxΔu, sdp_opt; λ=0.000001)
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
    if ElMin !== nothing
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
        println(UT.volume(ElMin))
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
# @static if get(ENV, "CI", "false") == "false" && (isdefined(@__MODULE__, :no_plot) && no_plot==false)
    using PyPlot
    include("../src/utils/plotting/plotting.jl")
    PyPlot.pygui(true) 


    PyPlot.rc("text",usetex=true)
    PyPlot.rc("font",family="serif")
    PyPlot.rc("font",serif="Computer Modern Roman")
    PyPlot.rc("text.latex",preamble="\\usepackage{amsfonts}")
    ##

    fig = PyPlot.figure(tight_layout=true, figsize=(4,4))

    ax = PyPlot.axes(aspect = "equal")
    ax.set_xlim(X[1].lo-0.2, X[1].hi+0.2)
    ax.set_ylim(X[2].lo-0.2, X[2].hi+0.2)

    vars = [1, 2];

    function projectEllipsoid(elli::UT.Ellipsoid)
        D = elli.P[1:2,3:3]
        r = elli.P[3,3]
        P = elli.P[1:2,1:2]-D*D'./r 
        UT.Ellipsoid(P,elli.c[1:2])
    end


    function plotEllipsoid!(ax, elli; ngrid=25, alpha=1.0, color="blue")
        # Radii corresponding to the coefficients:
        #rx, ry, rz = 1 ./ sqrt.(coefs)
        L = cholesky((elli.P+elli.P')/2).U;
    
        center = elli.c
        # Set of all spherical angles:
        u = range(0, 2pi, length=ngrid)
        v = range(0, pi, length=ngrid)
    
        # Cartesian coordinates that correspond to the spherical angles:
        # (this is the equation of an ellipsoid):
        x = [ x * y for (x, y) in Iterators.product(cos.(u), sin.(v))]
        y = [ x * y for (x, y) in Iterators.product(sin.(u), sin.(v))]
        z = [ x * y for (x, y) in Iterators.product(ones(length(u)), cos.(v))]
        Lm = inv(L)
        xR = [ [l[1]*x+l[2]*y+l[3]*z + c  for (x,y,z) in Iterators.zip(x,y,z)] for (l,c) in Iterators.zip(eachrow(Lm),center)]
        ax.plot_surface(xR[1], xR[2], xR[3], alpha=alpha, color=color)

    end

    function plotEllipse!(ax,elli::UT.Ellipsoid;
            fc = "red", fa = 0.5, ec = "black", ea = 1.0, ew = 0.5, n_points=30) 
        @assert length(vars) == 2 
        fca = Plot.FC(fc, fa)
        eca = Plot.FC(ec, ea)
        L = cholesky((elli.P+elli.P')/2).U;
        theta = range(0,2π,length=n_points);
        x = L\hcat(sin.(theta),cos.(theta))'


        vertslist = NTuple{n_points,Vector}[]

        push!(vertslist, tuple(Vector.(eachcol(x.+elli.c))...))
        
        polylist = matplotlib.collections.PolyCollection(vertslist)
        polylist.set_facecolor(fca)
        polylist.set_edgecolor(eca)
        polylist.set_linewidth(ew)
        ax.add_collection(polylist)

    end
    plotted = []
    currNode = last(treeLeaves)
    LyapMax = (max(map(x-> x.path_cost, treeLeaves)...))
    for n in treeLeaves
        push!(plotted, n)
        lyap = (n.path_cost)
        plotEllipse!(ax, projectEllipsoid(n.state), fc =  (0.4*lyap/LyapMax+0.5, 0.4*(LyapMax-lyap)/LyapMax+0.5, 0.75), ew = 0.5);
        nodePar = n.parent
        aTail = n.state.c
        while nodePar !== nothing 
            aDir = (nodePar.state.c-aTail)*0.8
            if !(nodePar ∈ plotted)
                lyap = (nodePar.path_cost)
                plotEllipse!(ax, projectEllipsoid(nodePar.state), fc =  (0.4*lyap/LyapMax+0.5, 0.4*(LyapMax-lyap)/LyapMax+0.5, 0.75), ew = 0.5);
                push!(plotted, nodePar)
                PyPlot.arrow(aTail[1], aTail[2], aDir[1], aDir[2], fc=(0,0,0), ec=(0,0,0),width=0.01, head_width=.18)
                aTail = nodePar.state.c
                nodePar = nodePar.parent
            else
                PyPlot.arrow(aTail[1], aTail[2], aDir[1], aDir[2], fc=(0,0,0), ec=(0,0,0),width=0.01, head_width=.18)
                nodePar = nothing
            end

        end
    end
    currNode = last(treeLeaves)
    while currNode !== nothing
        lyap = (currNode.path_cost)
        plotEllipse!(ax,projectEllipsoid(currNode.state), fc =  (0.4*lyap/LyapMax+0.5, 0.4*(LyapMax-lyap)/LyapMax+0.5, 0.75), ew = 2);
        global currNode = currNode.parent
    end

    cmap = PyPlot.ColorMap("mycolor",hcat([0.0,0.8,0.5,0.5],[0.8,0.0,0.5,0.5])');
    
    PyPlot.colorbar(PyPlot.ScalarMappable(norm=PyPlot.cm.colors.Normalize(vmin=0, vmax=LyapMax),cmap=cmap),shrink=0.7)

    PyPlot.xlabel("\$x_1\$", fontsize=14)
    PyPlot.ylabel("\$x_2\$", fontsize=14)
    xSpanp = hcat(xSpan...)
    PyPlot.plot(xSpanp[1,1:k],xSpanp[2,1:k],"bo-",markersize=4)
      
    PyPlot.title("Trajectory and Lyapunov-like Fun.", fontsize=14)
    plt.savefig("unicycle_2.pdf", format="pdf")
    gcf() 
    ################################################
    fig = PyPlot.figure(tight_layout=true, figsize=(4,4))

    ax = PyPlot.axes(projection="3d")
    ax.set_xlim(X[1].lo-0.2, X[1].hi+0.2)
    ax.set_ylim(X[2].lo-0.2, X[2].hi+0.2)
    ax.set_zlim(X[3].lo-0.2, X[3].hi+0.2)
    currNode = last(treeLeaves)  
    while currNode !== nothing
        lyap = (currNode.path_cost)
        plotEllipsoid!(ax ,currNode.state, color =  (0.4*lyap/LyapMax+0.5, 0.4*(LyapMax-lyap)/LyapMax+0.5, 0.75),alpha=0.4);
        global currNode = currNode.parent
    end 
    
    PyPlot.xlabel("\$x\$", fontsize=14)
    PyPlot.ylabel("\$y\$", fontsize=14)
    PyPlot.zlabel("\$\\theta\$", fontsize=14)
    PyPlot.plot3D(xSpanp[1,1:k],xSpanp[2,1:k],xSpanp[3,1:k],"bo-",markersize=4)
      
    PyPlot.title("Trajectory and Lyapunov-like Fun.", fontsize=14)
    plt.savefig("unicycle_3d.pdf", format="pdf")
    gcf() 
end #if 

