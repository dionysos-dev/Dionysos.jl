using Dionysos

# const PR = Dionysos.Problem
const SY = Dionysos.System
const UT = Dionysos.Utils

using Symbolics
using IntervalArithmetic
using LinearAlgebra

#Def sym variables
@variables px py th v ω w1 T
x = [px; py; th] # state
u = [v; ω] # control
w = [w1]

Ts = 0.1
# discrete-time dynamic function
# due to "Global Observability Analysis of a Nonholonomic Robot using Range Sensors"
#T = 0.1
# f = [px+2*v/ω*sin(T*ω/2)*cos(th+T*ω/2);
#      py+2*v/ω*sin(T*ω/2)*sin(th+T*ω/2);
#      th+T*ω
#     ];

# sinc not implemented in symbolic and division by `th` makes IntervalArithmetic bug
mysinc(th) = sin(sqrt(th^2+1e-14))/sqrt(th^2+1e-14)
f = [px+T*v*mysinc(T*ω/2)*cos(th+T*ω/2);
     py+T*v*mysinc(T*ω/2)*sin(th+T*ω/2);
     th+T*ω
    ];

# sys dimensions
n = length(x)
m = length(u)
p = length(w)

# sys eval function
function f_eval(x̄,ū,w̄,T̄)
    rules = Dict([x[i] => x̄[i] for i=1:n] ∪ 
            [u[i] => ū[i] for i=1:m] ∪ 
            [w[i] => w̄[i] for i=1:p] ∪
            [T => T̄])

    ff = Symbolics.substitute(f,rules)
    Base.invokelatest(eval(eval(build_function(ff)[1])))
end

# augmented argument
xi = [x;u;w]

fT = Symbolics.substitute(f,Dict([T => Ts]))

# Sym Jacobian
#Jxi_s = Symbolics.jacobian(fT,xi)


# Point of approximation
x̄ = [2.0;2.0;2.0]
ū = [2.0;2.0];
w̄ = 0.0;

# Box bounding x, u, w
X = IntervalBox(-10..10,2) × (-pi..pi);
U = IntervalBox(-10.0..10.0,2)

# Boxes on which J must be bounded
ΔX = IntervalBox(-0.1..0.1,2) × (-0.1..0.1);

ΔU = IntervalBox(-1.0..1.0,2)

ΔW = -1..1



x0 = [10.0;10.0;0.0]
X0 = UT.Ellipsoid(Matrix{Float64}(I(n)),x0)

XF = UT.Ellipsoid(Matrix{Float64}(I(n)),zeros(n))

treeRoot = UT.Node(XF)
treeLeaves = [treeRoot]
Xnew = XF
sample_x() = map(x-> x.lo + (x.hi-x.lo)*rand(), X.v)
sample_u() = map(x-> x.lo + (x.hi-x.lo)*rand(), U.v)

function findClosestNode(nodeList, x)
    dists = map(e-> e===nothing ? Inf : UT.pointCenterDistance(e.state, x), nodeList)
    d, idx = findmin(dists)
    parents = filter(x -> x!==nothing, unique(map(x-> x.parent,nodeList)))
    if !isempty(parents)
        n, dpar = findClosestNode(parents, x) 
        if d<dpar
            return nodeList[idx], d
        else
            return n, dpar
        end
    else
        return nodeList[idx], d
    end
end
maxIter = 100
while !(X0 ∈ Xnew) && maxIter>0
    xsample = sample_x()
    Xclosest, _ = findClosestNode(treeLeaves,xsample)
    xPar = UT.get_center(Xclosest.state)
    unew = sample_u()
    xnew = f_eval(xPar,unew,0.0,-Ts)

    X̄ = IntervalBox(xnew .+ ΔX)
    Ū = IntervalBox(unew .+ ΔU)
    W̄ = IntervalBox(0.0 .+ ΔW)
    (sys, L) = Dionysos.System.buildAffineApproximation(fT,x,u,w,xnew,unew,0.0,X̄,Ū,W̄)
    El, K, cost = SY.hasTransition(xnew,unew,Xclosest,sys,L)
    if El!==nothing
        UT.Node(El; parent=Xclosest, action=K, path_cost=cost)
        push!(treeLeaves, Xclosest)
        setdiff!(treeleaves, Xclosest)
    end 
    
    maxIter-=1
end


# fix system to be discrete-time (paper ECC20)
# calculate xcenter from the backwards integration
# 