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
using JuMP, Mosek, MosekTools, SCS
using Plots, Colors
Random.seed!(0)


function example_box_ellipsoid()
    c = [-10.0;-10.0]
    P = [2.0 6.0; 6.0 20.0]
    E = UT.Ellipsoid(P, c)
    box = UT.get_min_bounding_box(E)
    p = plot(aspect_ratio=:equal)
    UT.plot_box!(box)
    plot!(p, E)
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
    
    ########## Function description #########
    f, x, u, w, T = unstableSimple()
    nx = length(x)
    nu = length(u)
    nw = length(w)
    Ts = 1
    f_eval = eval(build_function(f,x,u,w,T)[1])
    
    #### PWA approximation description #####
    fsymbolic = Symbolics.substitute(f,Dict([T => Ts]))
    ΔX = IntervalBox(-1.0..1.0,2) 
    ΔU = IntervalBox(-10*2..10*2,2)
    ΔW = IntervalBox(-0.0..0.0,1)
    
    ########## Inputs description ##########
    Usz = 10
    Ub = IntervalBox(-Usz..Usz,nu)
    U = convert_H_into_E(Ub)
    
    ########## Noise description ##########
    W = 0.0*[-1 -1  1 1;
             -1  1 -1 1]
    
    system = ST.EllipsoidalLazySystem(f_eval, Ts, nx, nu, nw, U, Ub, W, X, obstacles, fsymbolic, x, u, w, ΔX, ΔU, ΔW)
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
    X̄ = IntervalBox(xnew .+ sys.ΔX)
    Ū = IntervalBox(unew .+ sys.ΔU)
    W̄ = IntervalBox(wnew .+ sys.ΔW)
    (affineSys, L) = ST.buildAffineApproximation(sys.fT, sys.x, sys.u, sys.w, xnew, unew, wnew, X̄, Ū, W̄)
    W = 0.0*[-1 -1  1 1;
             -1  1 -1 1]
    n_x = 2
    n_u = 2
    n_w = 2
    ########### Cost description ############
    S = Matrix{Float64}(I(nx+nu+1)) #TODO
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


function test_backward_transition()
    c = [0.0 ; 0.0]
    E2 = UT.Ellipsoid([2.0 0.2 ; 0.2 0.5], [3.0 ; 3.0])

    ##### get the system #####
    sys = create_sys()
    xnew = c
    unew = [0.0;0.0] 
    wnew = zeros(sys.nw)
    X̄ = IntervalBox(xnew .+ sys.ΔX)
    Ū = IntervalBox(unew .+ sys.ΔU)
    W̄ = IntervalBox(wnew .+ sys.ΔW)
    (affineSys, L) = ST.buildAffineApproximation(sys.fsymbolic, sys.x, sys.u, sys.w, xnew, unew, wnew, X̄, Ū, W̄)
    S = problem.state_cost
    E1, kappa, cost = SY.transition_backward(affineSys, E2, xnew, unew, sys.U, S, L, optimizer.sdp_opt ; λ=optimizer.λ, maxδx=optimizer.maxδx, maxδu=optimizer.maxδu)

    # W = 0.0*[-1 -1  1 1;
    #          -1  1 -1 1]
    
    # n_x = 2
    # n_u = 2
    # n_w = 2
    # S = Matrix{Float64}(I(n_x+n_u+1))
    # sdp_opt =  optimizer_with_attributes(SCS.Optimizer, MOI.Silent() => true)

    # A = affineSys.A
    # B = affineSys.B
    # g = affineSys.c
    # Lip = L
    # U = sys.Ub
    # u = zeros(2)
    # ##########################
    # #A, B, g, ct, Pt, c, U, W, S, Lip, optimizer;  maxRadius=Inf, maxΔu=Inf, λ=0.01

    # maxδx = 100  #
    # maxδu = 10.0*2
    # #E1, kappa, cost = Dionysos.Symbolic.transition_backward(A, B, g, E2.c, E2.P, c, u, U, W, S, Lip, sdp_opt, affineSys; maxδx=maxδx, maxδu=maxδu, λ=0.08) #0.01
    # E1, kappa, cost = SY.transition_backward(affineSys, E2, c, u, U, S, Lip, sdp_opt; maxδx=maxδx, maxδu=maxδu, λ=0.01)

    # E1 = UT.Ellipsoid([2.0 0.2 ; 0.2 0.5], [0.0 ; 0.0])
    # ans, kappa, cost = Dionysos.Symbolic.transition_fixed(affineSys, E1, E2, sys.Ub, W, S, sdp_opt)

    
    p = plot(aspect_ratio=:equal)
    plot!(p, E1, color = :green)
    plot!(p, E2, color = :red)
    display(p)
    println(kappa)
    println(cost)

    function f_eval(x, u, w, Ts)        
        return A*x + B*u + g
    end
    function c_eval(x)
        return zeros(2)
    end
    println(sys.U)
    Uset = UT.Ellipsoid([2.0 0.2 ; 0.2 0.5], [0.0 ; 0.0])
    println(SY.check_controller(E1, E2, f_eval, c_eval, 2, Uset; N=500)) #check_controller(E1, kappa, E2, f_eval, sys.Ts; N=500))
    #SY.plot_controller(E1, kappa, E2, f_eval, sys.Ts; N=100)
    #println(vol) # on a vol <= 4/3 π δx^3


end

using JuMP, Mosek, MosekTools, SCS
function test_log_det()

    # sdp_opt =  optimizer_with_attributes(Mosek.Optimizer, MOI.Silent() => true)

    # model = Model(sdp_opt)
    # @variable(model, X[1:3, 1:3])
    # @variable(model, t)
    # @constraint(model, [t; 1; vec(X)] in MOI.LogDetConeSquare(3))
    # @objective(model, Max, t)
    # optimize!(model)

    sdp_opt = optimizer_with_attributes(SCS.Optimizer, MOI.Silent() => true) 
    model = Model(sdp_opt)
    @variable(model, Q[1:3, 1:3] in PSDCone())
    @variable(model, t)
    u_q = [Q[i, j] for j in 1:3 for i in 1:j]
    @constraint(model, vcat(t, 1, u_q) in MOI.LogDetConeTriangle(3))
    @objective(model, Max, t)
    optimize!(model)
    println(value(t))
end

# example_box_ellipsoid()
# example_test_controller()
test_backward_transition()
# test_log_det()