# Warning : deprecated example, see https://github.com/dionysos-dev/Dionysos.jl/issues/221

using .Dionysos
using Polyhedra
using MathematicalSystems, HybridSystems
using CDDLib
using SemialgebraicSets
using StaticArrays
using LinearAlgebra
using Plots
using SDPA, JuMP
using IntervalArithmetic
UT = Dionysos.Utils
SI = Dionysos.Symbolic


lib = CDDLib.Library() #polyhedron lib
# aux functions
eye(n) = diagm(ones(n)) # I matrix
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

function get_controllers_matrices(m)
    return m[:, 1:end-1], m[:, end]
end


function get_f_eval(sys)
    return f_eval(x, u, w) = sys.A*x+sys.B*u+sys.D*w+sys.c
end
function get_c_eval(kappa, c)
    K, ℓ = get_controllers_matrices(kappa)
    return c_eval(x) = K*(x-c)+ℓ
end

# for a two-dimensional state space system
function plot_results_2D(sys, kappa, E1, E2)
    aff_sys = sys.constrainedAffineSys
    f_eval = get_f_eval(aff_sys)
    c_eval = get_c_eval(kappa, E1.c)
    
    p = plot(aspect_ratio=:equal)
    plot!(p, E1, color = :green)
    plot!(p, E2, color = :blue)

    K, ℓ = get_controllers_matrices(kappa)
    Ef = UT.affine_transformation(E1, aff_sys.A+aff_sys.B*K, aff_sys.B*ℓ+aff_sys.c)
    plot!(p, Ef, color = :red)
    SI.plot_transitions!(E1, f_eval, c_eval, 1; N=100)
    display(p)
end

function trial(dt, Ubound, Wmax, contraction, initial_vol)
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
    sys = Dionysos.System.build_AffineApproximationDiscreteSystem(Ac, Bc, gc, E, X̄, Ū, W̄, L)
    has_transition, kappa, cost = Dionysos.Symbolic.transition_fixed(sys.constrainedAffineSys, E1, E2, U, W, S, optimizer)
    sr = 1.0
    println("Has transition: $(has_transition)")
    if has_transition
        #println("K:\t $(K)\nell:\t $(ell)")
        println("cost:\t $(cost)")
        println("s.r.:\t $(sr)")
    end

    f_eval = get_f_eval(sys.constrainedAffineSys)
    c_eval = get_c_eval(kappa, E1.c)

    bool = SI.check_controller(E1, E2, f_eval, c_eval, nw, Ub; N=500) 
    println(bool)
    println(kappa)


    
    plot_results_2D(sys, kappa, E1, E2)
    K, ℓ = get_controllers_matrices(kappa)
    EU = UT.affine_transformation(E1, K, ℓ)
    println(EU)
    println(Ub)
    # p = plot(aspect_ratio=:equal)
    # #UT.plotE!(EU)

    # U = sm([1/10.0])
    # println("JULIEN CALBERT")
    # println(size(U))
    # c = [0.0]
    # UEllipsoid = UT.DegenerateEllipsoid(U, c)
    # UT.plotE!(UEllipsoid)
    
    # #plot!(Ub)
    # display(p)



    println()

    

    #SI.plot_check(E1, f_eval, c_eval, 1, E2; N=100)
    #return has_transition, cost, sr
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

trial(dt, Ubound, Wmax, contraction, initial_vol)





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

