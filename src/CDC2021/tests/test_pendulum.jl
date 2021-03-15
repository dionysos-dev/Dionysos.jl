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
AB = Abstraction

using ..DomainList
D = DomainList

using ..Utils
U = Utils

using ..BranchAndBound
BB = BranchAndBound

using ..OptimalControl
OC = OptimalControl

## data of the problem

#later, I think that L_growthbound could be computed in this constructor
function NewControlSystemGrowthUrdf(urdf,tstep,sysnoise::SVector{N,T},measnoise::SVector{N,T},L_growthbound,ngrowthbound) where {N,T}
    mechanism = parse_urdf(urdf)

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
    end
    # define the growthbound_map function
    function growthbound_map(r::SVector{N,T}, u, tstep, x)
        l = 0.1 #0.1
        K = IntervalBox(x[1]-l..x[1]+l,x[2]-l..x[2]+l, x[3]-l..x[3]+l, x[4]-l..x[4]+l)
        L = L_growthbound(u,K)
        function F_growthbound(r, u)
            return L*r + sysnoise
        end
        return AB.RungeKutta4(F_growthbound, r, u, tstep, ngrowthbound)+SVector(0.4,0.4,0.0,0.0)
    end
    return AB.ControlSystemGrowth(tstep, sysnoise, measnoise, sys_map, growthbound_map)
end


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
function build_system()
    urdf = "Acrobot.urdf"
    tstep = 0.04#0.004
    sysnoise = SVector(0.0, 0.0, 0.0, 0.0)
    measnoise = SVector(0.0, 0.0, 0.0, 0.0)
    ngrowthbound = 4

    return NewControlSystemGrowthUrdf(urdf,tstep,sysnoise,measnoise,L_growthbound,ngrowthbound)
end
function build_input()
    U = AB.HyperRectangle(SVector(-2.0), SVector(2.0))
    x0 = SVector(0.0); hu = SVector(0.5)
    Ugrid = AB.GridFree(x0,hu)
    Udom = AB.DomainList(Ugrid)
    AB.add_set!(Udom, U, AB.OUTER)
    box = AB.HyperRectangle(SVector(-0.1), SVector(0.1))
    AB.remove_set!(Udom, box, AB.OUTER)
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
    #println(rect)
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
    #println(AB.HyperRectangle(lb,ub))
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
    Fx = contsys.sys_map(x, -u, tstep)
    r = Xdom.grid.h/2.0 + contsys.measnoise
    Fr = contsys.growthbound_map(r, -u, tstep, x)

    rectI = AB.get_pos_lims_outer(Xdom.grid, AB.HyperRectangle(Fx .- Fr, Fx .+ Fr))
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

function test()
    println("start")
    X = AB.HyperRectangle(SVector(0.0,0.0,-5.0,-5.0),SVector(2*π,2*π,5.0,5.0))
    contsys = build_system()
    Udom = build_input()

    # control problem
    _T_ = AB.HyperRectangle(SVector(π, 0.1,0.0,0.0)-SVector(0.1,0.0,0.1,0.1), SVector(π, 0.1,0.0,0.0)+SVector(0.1,0.1,0.1,0.1))
    x0 = SVector(π,π,0.0,0.0)
    RI = SVector(0.01, 0.01, 0.01, 0.01)
    _I_ = AB.HyperRectangle(x0-RI, x0+RI)

    # problem-specific functions
    functions = [compute_reachable_set,minimum_transition_cost,post_image,pre_image]

    hx_coarse = [1.0, 1.0, 2.0, 2.0]
    hx_medium = [0.3, 0.3, 2.0, 2.0]
    hx_fine = [0.1, 0.1, 1.0, 1.0]
    periodic = [1,2]
    periods = [2*π,2*π]
    fig = plot(aspect_ratio = 1,legend = false)
    plot!(U.rectangle(_I_.lb,_I_.ub), opacity=.8,color=:green)
    plot!(U.rectangle(_T_.lb,_T_.ub), opacity=.8,color=:red)
    functions = [compute_reachable_set,minimum_transition_cost,post_image,pre_image]
    optimal_control_prob = OC.OptimalControlProblem(x0,_I_,_T_,contsys,periodic,periods,Udom,transition_cost,(X,hx_coarse),hx_medium,hx_fine,functions)

    max_iter = 5
    max_time = 1000
    optimizer = BB.Optimizer(optimal_control_prob,max_iter,max_time,log_level=2)
    println("optimize")
    @time MOI.optimize!(optimizer)
    println(optimizer.status)
    println(optimizer.best_sol)

    display(fig)
end

println()
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
