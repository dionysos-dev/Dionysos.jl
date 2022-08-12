using Dionysos
using Dionysos.Problem
using Dionysos.System
using Dionysos.Utils

const PR = Dionysos.Problem
const SY = Dionysos.System
const UT = Dionysos.Utils

using Symbolics
using IntervalArithmetic
using LinearAlgebra

#Def sym variables
@variables px py th dpx dpy dth F τ w1
x = [px; py; th; dpx; dpy; dth] # state
u = [F; τ] # control
w = [w1]

# dynamic function
f = [dpx;
     dpy;
     dth;
     cos(th)*F;
     sin(th)*F;
     τ
    ];

# sys dimensions
n = length(x)
m = length(u)

# augmented argument
xi = [x;u;w]

# Sym Jacobian
Jxi_s = Symbolics.jacobian(f,xi)


# Point of approximation
x̄ = [2.0;2.0;2.0;2.0;2.0;2.0]
ū = [2.0;2.0];
w̄ = 0.0;

# Box bounding x, u, w
X = IntervalBox(-10..10,2) × (-pi..pi) × IntervalBox(-5..5,3);

# Boxes on which J must be bounded
ΔX = IntervalBox(-0.1..0.1,2) × (-0.1..0.1) × IntervalBox(-0.1..0.1,3);
X̄ = IntervalBox(x̄ .+ ΔX)

ΔU = IntervalBox(-1.0..1.0,2)
Ū = IntervalBox(ū .+ ΔU)

ΔW = -∞..∞
W̄ = IntervalBox(w̄ .+ ΔW)


(sys, L) = Dionysos.System.buildAffineApproximation(f,x,u,w,x̄,ū,w̄,X̄,Ū,W̄)

x0 = [10.0;10.0;zeros(n-2)]
X0 = UT.Ellipsoid(I(n),x0)

XF = UT.Ellipsoid(I(n),zeros(n))

TreeRoot = UT.Search.Node(XF)
