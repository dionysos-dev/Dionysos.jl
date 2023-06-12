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

using LinearAlgebra, StaticArrays, Random, Symbolics, IntervalArithmetic
using MathematicalSystems, HybridSystems
using JuMP, Mosek, MosekTools, SDPA
using Plots, Colors
using CDDLib, SemialgebraicSets, Polyhedra
Random.seed!(0)

lib = CDDLib.Library() #polyhedron lib
# aux functions
sm(M) = SMatrix{size(M,1),size(M,2)}(M)
sv(M) = SVector{size(M,1)}(M)

optimizer = optimizer_with_attributes(SDPA.Optimizer, MOI.Silent() => true)

# convert square box into intersection of degenerated ellipsoids
# square box [-x,x]
function convert_H_into_E(Ub::IntervalBox)
    nu = length(Ub)
    Uaux = diagm(1:nu)
    U = [(Uaux.==i)./Ub[i].hi for i in 1:nu]

    return U
end


function get_cost_eval()
    return cost_eval(x, u) = u'*u
end

# for a two-dimensional state space system
function plot_results_2D(sys, cont, E1, E2)
    f_eval = ST.get_f_eval(sys)
    aff_sys = sys.constrainedAffineSys
    
    fig = plot(aspect_ratio=:equal)
    plot!(E1, color = :green)
    plot!(E2, color = :blue)
    Ef = UT.affine_transformation(E1, aff_sys.A+aff_sys.B*cont.K, aff_sys.B*cont.ℓ+aff_sys.c)
    plot!(Ef, color = :red)
    ST.plot_transitions!(E1, f_eval, cont.c_eval, 1; N=100)
    display(fig)
end

function trial_jc(dt, Ubound, Wmax, contraction, initial_vol)
    ########## PWA approximation description ##########
    Ac = sm([0.0  1.0;;
             1.0 -1.0]);

    Bc = sm([1.0; 1.0]);

    gc = sv(zeros(2,1))

    E = sm([1.0; 1.0]);
    nx = 2
    nu = 1
    nw = 1 

    xbar = zeros(nx)
    ubar = zeros(nu)
    wbar = zeros(nw)
    ΔX = IntervalBox(-1.0..1.0, nx) 
    ΔU = IntervalBox(-10*2..10*2, nu)
    ΔW = IntervalBox(-0.0..0.0, nw)
    X̄ =  IntervalBox(xbar .+ ΔX)
    Ū =  IntervalBox(ubar .+ ΔU)
    W̄ =  IntervalBox(wbar .+ ΔW)
    L = [0.0;0.0;0.0]

    ########## Inputs description ##########
    #Ubound = 10
    Ub = IntervalBox(-Ubound..Ubound, nu)
    U = convert_H_into_E(Ub)
    println(U)
    ########## Noise description ##########
    W = Wmax*[-1 1]
    ########## Cost description ##########
    S = Matrix{Float64}(I(nx+nu+1))

    ########## Ellipsoids ##########
    E1 = UT.Ellipsoid(Matrix{Float64}(I(nx)), [0.0;0.0])
    E2 = UT.Ellipsoid(Matrix{Float64}(I(nx))*(1/15), [3.0;3.0])

    ################################################################################
    sys = ST.build_AffineApproximationDiscreteSystem(Ac, Bc, gc, E, X̄, Ū, W̄, L)
    has_transition, cont, cost = SY.transition_fixed(sys.constrainedAffineSys, E1, E2, U, W, S, optimizer)
    sr = 1.0
    println("Has transition: $(has_transition)")
    if has_transition
        #println("K:\t $(K)\nell:\t $(ell)")
        println("cost:\t $(cost)")
        println("s.r.:\t $(sr)")
    end

    f_eval = ST.get_f_eval(sys)
    c_eval = ST.get_c_eval(cont)

    cost_eval = get_cost_eval()
    bool = ST.check_feasibility(E1, E2, f_eval, c_eval, nw, Ub; N=100) 
    println(bool)

    fig = plot(aspect_ratio=:equal)
    ST.plot_check_feasibility!(E1, E2, f_eval, c_eval, nw; N=100)
    display(fig)

    fig = plot(aspect_ratio=:equal)
    ST.plot_controller_cost!(E1, c_eval, cost_eval; N=4000, scale=0.004)
    display(fig)

    plot_results_2D(sys, cont, E1, E2)


    # U = sm([1/10.0])
    # println("JULIEN CALBERT")
    # println(size(U))
    # c = [0.0]
    # UEllipsoid = UT.DegenerateEllipsoid(U, c)
    # plot!(UEllipsoid)
    
    # plot!(Ub)


    println()
end


#to vary 
# - initial volume
# - dt
# - contraction 

dt = 0.5
Ubound = 10.0 # upper limit on |u|
Wmax = 0.0
initial_vol = 10
contraction = 0.8 #1.0

trial_jc(dt, Ubound, Wmax, contraction, initial_vol)





# ########## PWA approximation description ##########
# Ac = sm([0.0  1.0  0.0;
# 0.0  0.0  1.0;
# 1.0 -1.0 -1.0]);

# Bc = sm([0.0; 0.0; 1.0]);

# gc = sv(zeros(3,1))

# E = sm([1.0; 1.0; 1.0]);
# nx = 3
# nu = 1
# nw = 1 

# xbar = zeros(nx)
# ubar = zeros(nu)
# wbar = zeros(nw)
# ΔX = IntervalBox(-1.0..1.0, nx) 
# ΔU = IntervalBox(-10*2..10*2, nu)
# ΔW = IntervalBox(-0.0..0.0, nw)
# X̄ =  IntervalBox(xbar .+ ΔX)
# Ū =  IntervalBox(ubar .+ ΔU)
# W̄ =  IntervalBox(wbar .+ ΔW)
# L = [0.0;0.0;0.0]

