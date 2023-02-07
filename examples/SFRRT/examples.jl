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


function example_box_ellipsoid()
    c = [-10.0;-10.0]
    P = [2.0 6.0; 6.0 20.0]
    E = UT.Ellipsoid(P, c)
    box = UT.get_min_bounding_box(E)
    p = plot(aspect_ratio=:equal)
    UT.plot_box!(box)
    UT.plotE!(E)
    display(p)
end



function create_sys()
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
    system = SY.EllipsoidalLazySystem(f_eval, X, U, Ub, Ts, fT, x, u, w, maxRadius, maxΔu, ΔX, ΔU, ΔW, [])
    return system 
end
function example_test_controller() # could be added as a test for ellipsoidal_transition
    E1 = UT.Ellipsoid([1.0 0.5 ; 0.5 1.0], [0.0 ; 0.0])
    E2 = UT.Ellipsoid([2.0 0.2 ; 0.2 0.5], [3.0 ; 3.0])

    ##### get the system #####
    sys = create_sys()
    xnew = E1.c
    unew = [0.0;0.0] 
    wnew = zeros(2)
    X̄ = IntervalBox(xnew .+ ΔX)
    Ū = IntervalBox(unew .+ ΔU)
    W̄ = IntervalBox(wnew .+ ΔW)
    (affineSys, L) = Dionysos.System.buildAffineApproximation(sys.fT,sys.x,sys.u,sys.w,xnew,unew,wnew,X̄,Ū,W̄)
    W = 0.0*[-1 -1  1 1;
             -1  1 -1 1]
    S = Matrix{Float64}(I(n_x+n_u+1))
    sdp_opt =  optimizer_with_attributes(Mosek.Optimizer, MOI.Silent() => true)

    A = affineSys.A
    B = affineSys.B
    c = affineSys.c
    ##########################
    ans, cost, kappa = Dionysos.Symbolic.my_has_transition(A, B, c, W, sys.Ub, S, E1, E2, sdp_opt)
    println("succeed to compute the controller : ", ans)
    function f_eval(x, u, w, Ts)        
        return A*x + B*u + c
    end
    function f_cost(x, u)
        return norm(x) + norm(u)# 1.0# norm(x) + norm(u)
    end

    #could be good to print the cost in the original ellipsoid from te transition with a color bar.
    SC.check_controller(E1, kappa, E2, f_eval, sys.Ts; N=500)

    # SC.plot_controller(E1, kappa, E2, f_eval, sys.Ts; N=100)

    SC.plot_controller_cost(E1, kappa, E2, f_eval, sys.Ts, f_cost; N=2000)

end


# example_box_ellipsoid()
example_test_controller()