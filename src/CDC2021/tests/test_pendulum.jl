include(joinpath("..", "..", "Abstraction", "abstraction.jl"))
include("../general_domain.jl")
include("../utils.jl")
include("../alternating_simulation.jl")
include("../partition.jl")
include("../lazy_abstraction.jl")
include("../branch_and_bound.jl")
include("../optimal_control.jl")

module TestDoublePendulum
using RigidBodyDynamics, MeshCat, MeshCatMechanisms, SymPy, ForwardDiff, IntervalArithmetic
using LinearAlgebra,StaticArrays,Random
using Plots, StaticArrays, JuMP
using ReachabilityAnalysis

using ..Abstraction
const AB = Abstraction

using ..DomainList
const D = DomainList

using ..Utils
const U = Utils

using ..BranchAndBound
const BB = BranchAndBound

using ..OptimalControl
const OC = OptimalControl

using ..Lazy_abstraction ##
const LA = Lazy_abstraction

using ..AlternatingSimulation ##
const AS = AlternatingSimulation

## data of the problem

#later, I think that L_growthbound could be computed in this constructor
function NewControlSystemGrowthUrdf(urdf,tstep,sysnoise::SVector{N,T},measnoise::SVector{N,T},L_growthbound,ngrowthbound,X) where {N,T}

    function sys_map(x, u, tstep)
        return AB.RungeKutta4(f1, x, u, tstep, ngrowthbound)
    end
    #=mechanism = parse_urdf(urdf)

    # define the system function
    function sys_map(x, u, tstep)
        q = [x[1],x[2]] #SVector(x[1],x[2]) #[x[1],x[2]]  SVector ne fonctionne pas
        v = [x[3],x[4]] #SVector(x[4],x[3]) #[x[3],x[4]]
        function control!(torques::AbstractVector, t, state::MechanismState)
            torques[1] = u[1]
            torques[2] = 0.0
        end
        state = MechanismState(mechanism, q, v)
        t, q, v = simulate(state, tstep,control!; Δt = tstep)
        return SVector(q[2]...,v[2]...)
    end=#
    # define the growthbound_map function
    function growthbound_map(r::SVector{N,T}, u, tstep, x)
        K = compute_K(X,r,u,tstep,x)
        L = L_growthbound(u,K)
        function F_growthbound(r, u)
            return L*r + sysnoise
        end
        return AB.RungeKutta4(F_growthbound, r, u, tstep, ngrowthbound)
    end
    return AB.ControlSystemGrowth(tstep, sysnoise, measnoise, sys_map, growthbound_map)
end

function f1(x,u)
    lc1 = 0.5; l1 = 1.; m1 = 1.; I1 = 0.083; lc2 = 1.; l2 = 2.; m2 = 1.; I2 = 0.33
    g = 2.0 #9.81

    p1 =  m1*lc1^2+m2*l1^2+I1; p2 = m2*lc2^2+I2; p3 =  m2*l1*lc2; p4 = m1*lc1+m2*l1; p5 =  m2*lc2
    q1 = x[1]; q2 = x[2]; v1 = x[3]; v2 = x[4]
    F1 = x[3]
    F2 = x[4]

    s1 = sin(q1); s2 = sin(q2); s12 = sin(q1-q2); c12 = cos(q1-q2)
    coef1 = 1/((p3^2)*c12^2-p1*p2)
    coef2 = (p3^2/2)*sin(2*(q1-q2))
    coef3 = p2*p3*s12
    coef4 = -p3*p5*g*s2*c12+p2*p4*g*s1-p2*u[1]
    F3 = coef1*(coef2*v1^2+coef3*v2^2 + coef4)

    coef1 = 1/((p3^2)*c12^2-p1*p2)
    coef2 = -p1*p3*s12
    coef3 = -(p3^2/2)*sin(2*(q1-q2))
    coef4 = -p3*p4*g*s1*c12+p1*p5*g*s2+p3*c12*u[1]
    F4 = coef1*(coef2*v1^2+coef3*v2^2 + coef4)

    F = SVector(F1,F2,F3,F4)
    return F
end
function f(x,u)
    lc1 = 0.5; l1 = 1.; m1 = 1.; I1 = 0.083; lc2 = 1.; l2 = 2.; m2 = 1.; I2 = 0.33
    g = 2.0#9.81

    p1 =  m1*lc1^2+m2*l1^2+I1; p2 = m2*lc2^2+I2; p3 =  m2*l1*lc2; p4 = m1*lc1+m2*l1; p5 =  m2*lc2
    q1 = x[1]; q2 = x[2]; v1 = x[3]; v2 = x[4]
    F1 = x[3]
    F2 = x[4]

    s1 = sin(q1); s2 = sin(q2); s12 = sin(q1-q2); c12 = cos(q1-q2)
    coef1 = 1/((p3^2)*c12^2-p1*p2)
    coef2 = (p3^2/2)*sin(2*(q1-q2))
    coef3 = p2*p3*s12
    coef4 = -p3*p5*g*s2*c12+p2*p4*g*s1-p2*u[1]
    F3 = coef1*(coef2*v1^2+coef3*v2^2 + coef4)

    coef1 = 1/((p3^2)*c12^2-p1*p2)
    coef2 = -p1*p3*s12
    coef3 = -(p3^2/2)*sin(2*(q1-q2))
    coef4 = -p3*p4*g*s1*c12+p1*p5*g*s2+p3*c12*u[1]
    F4 = coef1*(coef2*v1^2+coef3*v2^2 + coef4)

    F = IntervalBox(F1,F2,F3,F4)
    return F
end

function compute_K(X,r,u,tstep,x)
    X = IntervalBox(X.lb, X.ub)
    B = f(X,u)
    R = IntervalBox(-r,r)
    K = x + R + B*tstep
    return K
end

function L_growthbound(u,K)
    function f(x)
        lc1 = 0.5; l1 = 1.; m1 = 1.; I1 = 0.083; lc2 = 1.; l2 = 2.; m2 = 1.; I2 = 0.33
        g = 2.0#9.81

        p1 =  m1*lc1^2+m2*l1^2+I1; p2 = m2*lc2^2+I2; p3 =  m2*l1*lc2; p4 = m1*lc1+m2*l1; p5 =  m2*lc2
        q1 = x[1]; q2 = x[2]; v1 = x[3]; v2 = x[4]
        F1 = x[3]
        F2 = x[4]

        s1 = sin(q1); s2 = sin(q2); s12 = sin(q1-q2); c12 = cos(q1-q2)
        coef1 = 1/((p3^2)*c12^2-p1*p2)
        coef2 = (p3^2/2)*sin(2*(q1-q2))
        coef3 = p2*p3*s12
        coef4 = -p3*p5*g*s2*c12+p2*p4*g*s1-p2*u[1]
        F3 = coef1*(coef2*v1^2+coef3*v2^2 + coef4)

        coef1 = 1/((p3^2)*c12^2-p1*p2)
        coef2 = -p1*p3*s12
        coef3 = -(p3^2/2)*sin(2*(q1-q2))
        coef4 = -p3*p4*g*s1*c12+p1*p5*g*s2+p3*c12*u[1]
        F4 = coef1*(coef2*v1^2+coef3*v2^2 + coef4)

        F = similar(x)
        F[1] = F1 #note, we trivially know the two first lines of the jacobian matrix
        F[2] = F2
        F[3] = F3
        F[4] = F4
        return F
    end
    L_interval = ForwardDiff.jacobian(f, K)
    dims = size(L_interval)
    L = zeros(dims)
    for i=1:dims[1]
        for j=1:dims[2]
            L[i,j] =  max(abs(L_interval[i,j].lo),abs(L_interval[i,j].hi))
        end
    end
    return L
end
function build_system(X)
    urdf = "Acrobot.urdf"
    tstep = 0.1 #0.15#0.004
    sysnoise = SVector(0.0, 0.0, 0.0, 0.0)
    measnoise = SVector(0.0, 0.0, 0.0, 0.0)
    ngrowthbound = 4

    return NewControlSystemGrowthUrdf(urdf,tstep,sysnoise,measnoise,L_growthbound,ngrowthbound,X)
end
function build_input()
    U = AB.HyperRectangle(SVector(-5.0), SVector(5.0))
    x0 = SVector(0.0); hu = SVector(0.5)
    Ugrid = AB.GridFree(x0,hu)
    Udom = AB.DomainList(Ugrid)
    AB.add_set!(Udom, U, AB.OUTER)
    return Udom
end
function transition_cost(x,u)
    return 1.0
end

## problem specific-functions required
function minimum_transition_cost(symmodel,contsys,source,target)
    return 1.0
end

function compute_reachable_set(rect::AB.HyperRectangle,contsys,Udom)
    println(rect)
    tstep = contsys.tstep
    r = (rect.ub-rect.lb)/2.0 + contsys.measnoise
    x = U.center(rect)
    n =  U.dims(rect)
    lb = fill(Inf,n)
    ub = fill(-Inf,n)
    for upos in AB.enum_pos(Udom)
        u = AB.get_coord_by_pos(Udom.grid,upos)
        Fx = contsys.sys_map(x, u, tstep)
        Fr = contsys.growthbound_map(r, u, tstep, x)
        lb = min.(lb,Fx .- Fr)
        ub = max.(ub,Fx .+ Fr)
    end
    println(AB.HyperRectangle(lb,ub))
    return AB.HyperRectangle(lb,ub)
end


function post_image(symmodel,contsys,xpos,u)
    Xdom = symmodel.Xdom
    x = AB.get_coord_by_pos(Xdom.grid, xpos)
    tstep = contsys.tstep
    Fx = contsys.sys_map(x, u, tstep)
    r = Xdom.grid.h/2.0 + contsys.measnoise
    Fr = contsys.growthbound_map(r, u, tstep, x)
    post_rec = AB.HyperRectangle(Fx .- Fr, Fx .+ Fr)
    rectI = AB.get_pos_lims_outer(Xdom.grid, AB.HyperRectangle(Fx .- Fr, Fx .+ Fr))
    ypos_iter = Iterators.product(AB._ranges(rectI)...)
    over_approx = []
    allin = true
    for ypos in ypos_iter
        ypos = D.set_in_period_pos(Xdom,ypos)
        if !(ypos in Xdom)
            allin = false
            break
        end
        target = AB.get_state_by_xpos(symmodel, ypos)
        push!(over_approx, target)
    end
    return allin ? over_approx : []
end


function pre_image(symmodel,contsys,xpos,u)
    Xdom = symmodel.Xdom
    x = AB.get_coord_by_pos(Xdom.grid, xpos)
    tstep = contsys.tstep
    X = AB.HyperRectangle(SVector(0.0,0.0,-5.0,-5.0),SVector(2*π,2*π,5.0,5.0))
    r = Xdom.grid.h/2.0 + contsys.measnoise
    K = compute_K(X,r,u,-tstep,x)
    n = length(xpos); lb = zeros(n); ub = zeros(n)
    for i=1:n
        lb[i] = K[i].lo
        ub[i] = K[i].hi
    end
    R = AB.HyperRectangle(lb,ub)
    println(x)
    println(R)
    rectI = AB.get_pos_lims_outer(Xdom.grid, R)
    ypos_iter = Iterators.product(AB._ranges(rectI)...)
    over_approx = []
    allin = true
    for ypos in ypos_iter
        ypos = D.set_in_period_pos(Xdom,ypos)
        if ypos in Xdom
            target = AB.get_state_by_xpos(symmodel, ypos)
            push!(over_approx, target)
        end
    end
    return over_approx
end

##

function test2()
    println("start")
    tol = 10^(-8)
    X = AB.HyperRectangle(SVector(0.0,0.0,-5.0,-5.0),SVector(2*π+tol,2*π+tol,5.0,5.0))
    contsys = build_system(X)
    Udom = build_input()
    u = [5.0]
    r = SVector(0.01,0.01,0.1,0.1)
    x = SVector(0.1,0.1,0.1,0.1)
    K = compute_K(X,r, u, contsys.tstep, x)
    println(K)
    Fr = contsys.growthbound_map(r, u, contsys.tstep, x)
    post_rec = AB.HyperRectangle(x .- Fr, x .+ Fr)
    println(post_rec)
end



function test3()
    println()
    X = AB.HyperRectangle(SVector(0.0,0.0,-0.5,-0.5),SVector(0.1,0.1,0.5,0.5))
    contsys = build_system(X)
    Udom = build_input()
    hx = [0.01, 0.01, 0.05, 0.05]*2.0#[0.01, 0.01, 0.02, 0.02]
    periodic = [1,2]
    periods = [2*π,2*π]
    Xdom = D.GeneralDomainList(hx;periodic=periodic,periods=periods)
    AB.add_set!(Xdom, X, AB.INNER)

    pos = AB.get_pos_by_coord(Xdom.grid,SVector(0.01,0.01,0.0,0.0))
    x = AB.get_coord_by_pos(Xdom.grid, pos)
    tstep = contsys.tstep
    println("x:  ",x)
    u = [4.0]
    Fx = contsys.sys_map(x, u, tstep)
    println("Fx:  ",Fx)
    r = Xdom.grid.h/2.0 + contsys.measnoise
    rec = AB.HyperRectangle(x .- r, x .+ r)
    println("rec:  ",rec)
    Fr = contsys.growthbound_map(r, u, tstep, x)
    R = AB.HyperRectangle(Fx .- Fr, Fx .+ Fr)
    println("R:  ",R)
    fig = plot(aspect_ratio = 1,legend = false)
    plot!(U.rectangle(rec.lb[1:2],rec.ub[1:2]),color=:blue)
    plot!(U.rectangle(R.lb[1:2],R.ub[1:2]),color=:red)
    display(fig)
    fig = plot(aspect_ratio = 1,legend = false)
    plot!(U.rectangle(rec.lb[3:4],rec.ub[3:4]),color=:blue)
    plot!(U.rectangle(R.lb[3:4],R.ub[3:4]),color=:red)
    display(fig)
end

function test()
    println("start")
    X = AB.HyperRectangle(SVector(0.0,0.0,-5.0,-5.0),SVector(2*π,2*π,5.0,5.0))
    contsys = build_system(X)
    Udom = build_input()

    # control problem
    #_T_ = AB.HyperRectangle(SVector(π, 0.1,0.0,0.0)-SVector(0.1,0.0,0.1,0.1), SVector(π, 0.1,0.0,0.0)+SVector(0.1,0.1,0.1,0.1))
    _T_ = AB.HyperRectangle(SVector(0.1, 0.1,-0.8,-0.8), SVector(0.15, 0.15,0.8,0.8))
    x0 = SVector(0.2,0.2,0.0,0.0)#SVector(π,π,0.0,0.0)
    RI = SVector(0.01, 0.01, 0.2, 0.2)
    _I_ = AB.HyperRectangle(x0-RI, x0+RI)

    # problem-specific functions
    functions = [compute_reachable_set,minimum_transition_cost,post_image,pre_image]

    hx_coarse = [0.8, 0.8, 1.0, 1.0]
    hx_medium = [0.2, 0.2, 0.4, 0.4]
    hx_fine = [0.01, 0.01, 0.05, 0.05]*10
    periodic = [1,2]
    periods = [2*π,2*π]

    fig = plot(aspect_ratio = 1,legend = false)
    plot!(U.rectangle(_I_.lb,_I_.ub), opacity=.8,color=:green)
    plot!(U.rectangle(_T_.lb,_T_.ub), opacity=.8,color=:red)
    optimal_control_prob = OC.OptimalControlProblem(x0,_I_,_T_,contsys,periodic,periods,Udom,transition_cost,(X,hx_coarse),hx_medium,hx_fine,functions)

    max_iter = 1
    max_time = 1000
    optimizer = BB.Optimizer(optimal_control_prob,max_iter,max_time,log_level=2)
    println("optimize")
    @time MOI.optimize!(optimizer)
    println(optimizer.status)
    println(optimizer.best_sol)

    display(fig)
end




function h2(node::LA.S.Node,problem::LA.LazyAbstraction)
    source = node.state.source
    symmodel = problem.symmodel
    xpos = AB.get_xpos_by_state(symmodel, source)
    x = AB.get_coord_by_pos(symmodel.Xdom.grid, xpos)

    heuristic = problem.heuristic_data
    symmodel2 = heuristic.symmodel
    xpos2 = AB.get_pos_by_coord(symmodel2.Xdom.grid, x)
    source2 = AB.get_state_by_xpos(symmodel2, xpos2)
    return heuristic.dists[source2]
end

function build_heuristic_data(X,contsys,Udom,_I_)
    # build the alternating simulation
    hx = [0.1, 0.1, 0.5, 0.5]
    periodic = [1,2]
    periods = [2*π,2*π]
    Xdom = D.GeneralDomainList(hx;periodic=periodic,periods=periods)
    AB.add_set!(Xdom, X , AB.INNER)
    symmodel = AB.NewSymbolicModelListList(Xdom, Udom)
    problem = AS.symmodelProblem(symmodel,contsys,compute_reachable_set,minimum_transition_cost,AS.get_possible_transitions_2)
    symmodel.autom = AS.build_alternating_simulation(problem)
    # build the heurisitic
    initlist = U.get_symbol(symmodel,_I_,AB.OUTER)
    heuristic_data = AS.build_heuristic(symmodel,initlist)
    println("Heuristic ended")
    fig = plot(aspect_ratio = 1,legend = false)
    U.plot_domain!(heuristic_data.symmodel.Xdom)
    AS.plot_heuristic!(heuristic_data)
    display(fig)
    return heuristic_data
end
function test4()
    println("start")
    X = AB.HyperRectangle(SVector(0.0,0.0,-3.0,-3.0),SVector(0.3,0.3,3.0,3.0))
    contsys = build_system(X)
    Udom = build_input()

    # control problem
    #_T_ = AB.HyperRectangle(SVector(π, 0.1,0.0,0.0)-SVector(0.1,0.0,0.1,0.1), SVector(π, 0.1,0.0,0.0)+SVector(0.1,0.1,0.1,0.1))
    _T_ = AB.HyperRectangle(SVector(0.1, 0.1,-0.4,-0.4), SVector(0.15, 0.15,0.4,0.4))
    x0 = SVector(0.2,0.2,0.0,0.0)#SVector(π,π,0.0,0.0)
    RI = SVector(0.01, 0.01, 0.2, 0.2)
    _I_ = AB.HyperRectangle(x0-RI, x0+RI)

    # problem-specific functions
    functions = [compute_reachable_set,minimum_transition_cost,post_image,pre_image]
    hx = [0.05, 0.05, 0.1, 0.1]
    periodic = [1,2]
    periods = [2*π,2*π]

    Xdom = D.GeneralDomainList(hx;periodic=periodic,periods=periods)
    AB.add_set!(Xdom, X , AB.INNER)
    symmodel = AB.NewSymbolicModelListList(Xdom, Udom)
    initlist = U.get_symbol(symmodel,_I_,AB.OUTER)
    targetlist = U.get_symbol(symmodel,_T_,AB.OUTER)

    # Heuristic data
    fig = plot(aspect_ratio = 1,legend = false)
    heuristic_data = build_heuristic_data(X,contsys,Udom,_I_)
    println("initset :")
    println(initlist)
    println("targetset :")
    println(targetlist)
    # Lazy Abstraction implementation
    time = @elapsed begin
    #problem = LA.compute_controller(symmodel, contsys, initlist, targetlist, pre_image, post_image, h1)
    problem,sucess = LA.compute_controller(symmodel, contsys, initlist, targetlist, transition_cost, pre_image, post_image, h2, heuristic_data=heuristic_data)
    contr = problem.contr
    end
    println("total time: lazy abstraction + controller: ", time)
    #fig = plot(aspect_ratio = 1,legend = false)
    x0 = SVector(0.2,0.2,0.0,0.0)
    LA.plot_result!(problem,x0=x0)
    display(fig)
end
println()
#test2()
#test()
test()
end # end module





#=
function reachable_set(X0,u,tstep)
    function dynamics!(dx, x, p, t)
        # write the nonlinear differential equations defining the model
        lc1 = -0.5
        l1 = -1.
        m1 = 1.
        I1 = 0.333 # about joint instead of CoM in URDF
        lc2 = -1.
        l2 = -2.
        m2 = 1.
        I2 = 1.33 # about joint instead of CoM in URDF
        g = -9.81

        q1 = x[1]; q2 = x[2]; v1 = x[3]; v2 = x[4]
        c1 = cos(q1); c2 = cos(q2); s1 = sin(q1); s2 = sin(q2); s12 = sin(q1 + q2)

        M11 = I1 + I2 + m2 * l1^2 + 2 * m2 * l1 * lc2 * c2
        M12 = I2 + m2 * l1 * lc2 * c2
        M22 = I2
        M = [M11 M12; M12 M22]

        C11 = -2 * m2 * l1 * lc2 * s2 * v2
        C12 = -m2 * l1 * lc2 * s2 * v2
        C21 = m2 * l1 * lc2 * s2 * v1
        C22 = 0
        C = [C11 C12; C21 C22]

        G = [m1 * g * lc1 * s1 + m2 * g * (l1 * s1 + lc2 * s12); m2 * g * lc2 * s12]

        v = [v1;v2]
        τ = [u[1];0.0]
        r = inv(M)*(-C*v-G+τ)
        dx[1] = x[3]
        dx[2] = x[4]
        dx[3] = r[1]
        dx[4] = r[2]
    end

    # formulate the initial-value problem
    prob = @ivp(x' = dynamics!(x), x(0) ∈ X0, dim=4)
    # solve using a Taylor model set representation
    sol_forward = ReachabilityAnalysis.solve(prob, T=tstep, alg=TMJets())
    #sol_forward = ReachabilityAnalysis.solve(prob, T=tstep, alg=GLGM06(δ=tstep))#alg=TMJets())

    #ReachabilityAnalysis.post(TMJets(), prob, tspan=(0.0, tstep))


    #solz1 = overapproximate(sol_forward, Zonotope);
    #println(sol_forward)
    # formulate the initial-value problem
    #prob = @ivp(x' = duffing!(x), x(0) ∈ X0, dim=2)
    #sol_backward = solve(prob, tspan=(0.0, -10.0), alg=TMJets())


    # plot the flowpipe in state-space
    #println("plot")
    #=fig = plot(sol_forward, vars=(1, 2), xlab="x", ylab="y", lw=0.5, color=:red)
    display(fig)
    fig2 = plot(sol_forward, vars=(3, 4), xlab="v1", ylab="v2", lw=0.5, color=:blue)
    display(fig2)=#
    println(4)
end

function test_4()
    function f(x)
        f1 = 4*x[1]+x[2]
        f2 = x[2]^2
        return SVector(f1,f2)
    end
    function f2(x)
        F = similar(x)
        F[1] = 4*x[1]+x[2]
        F[2] = x[2]^2
        return F
    end

    l = @vars q1 q2 v1 v2
    println(l)
    l = SVector(l)
    println(l)
    println(l[1])
    ForwardDiff.jacobian(f2, l)

    #X = IntervalBox(1..2, 2..4)
    X = IntervalBox(q1..q2, q3..q4)
    result = ForwardDiff.jacobian(f2, X)
    println(result)
    println(X)
    L = [X...]
    println(L[1].lo)
    println(L[1].hi)
end

=#

#=

function L_growthbound(u,K)
    function f(x)
        # the nonlinear differential equations defining the model
        lc1 = -0.5
        l1 = -1.
        m1 = 1.
        I1 = 0.333 # about joint instead of CoM in URDF
        lc2 = -1.
        l2 = -2.
        m2 = 1.
        I2 = 1.33 # about joint instead of CoM in URDF
        g = -9.81

        q1 = x[1]; q2 = x[2]; v1 = x[3]; v2 = x[4]
        c1 = cos(q1); c2 = cos(q2); s1 = sin(q1); s2 = sin(q2); s12 = sin(q1 + q2)

        M11 = I1 + I2 + m2 * l1^2 + 2 * m2 * l1 * lc2 * c2
        M12 = I2 + m2 * l1 * lc2 * c2
        M22 = I2
        M = @SMatrix [M11 M12; M12 M22]

        C11 = -2 * m2 * l1 * lc2 * s2 * v2
        C12 = -m2 * l1 * lc2 * s2 * v2
        C21 = m2 * l1 * lc2 * s2 * v1
        C22 = 0
        C = @SMatrix [C11 C12; C21 C22]

        G = @SMatrix [m1 * g * lc1 * s1 + m2 * g * (l1 * s1 + lc2 * s12); m2 * g * lc2 * s12]

        v = @SMatrix [v1;v2]
        τ = @SMatrix [u[1];0.0]
        r = inv(M)*(-C*v-G+τ)
        F = similar(x)
        F[1] = x[3] #note, we trivially know the two first lines of the jacobian matrix
        F[2] = x[4]
        F[3] = r[1]
        F[4] = r[2]
        return F
    end
    L_interval = ForwardDiff.jacobian(f, K)
    dims = size(L_interval)
    L = zeros(dims)
    for i=1:dims[1]
        for j=1:dims[2]
            L[i,j] =  max(abs(L_interval[i,j].lo),abs(L_interval[i,j].hi))
        end
    end
    return L
end

=#
