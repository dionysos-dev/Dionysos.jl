# A comparative study of different solutions to the path-tracking problem for an articulated vehicle
# Paolo Bolzern, Arturo Locatelli

include("control_system.jl")

function f(x, u)
    v, tanδ = u
    L1 = 1.0
    L2 = 0.5
    c = 0.5
    # dx = v1 * cos(θ)    (1a)
    # dy = v1 * sin(θ)    (1b)
    # dθ = v1/L1 * tan(δ) (1c)
    # dφ = -v1/(L1*L2) * (L1 * sin(φ) + (c * cos(φ) + L2) * tan(φ))
    sθ, cθ = sincos(x[3])
    dx = v * sθ
    dy = v * cθ
    dθ = v * tanδ / L1
    sφ, cφ = sincos(x[4])
    dφ = -v * (sφ + (c * cφ + L2) * tanδ) / (L1 * L2)
    return IntervalBox(dx, dy, dθ, dφ)
end

function test()
    println("start")
    X = AB.HyperRectangle(SVector(0.0,0.0,-π,-π),SVector(5.0,5.0,-π,π))
    tstep = 0.1
    sysnoise = SVector(0.0, 0.0, 0.0, 0.0)
    measnoise = SVector(0.0, 0.0, 0.0, 0.0)
    nsys = 5
    ngrowthbound = 5
    contsys = NewControlSystemGrowthLx(tstep, f, sysnoise, measnoise, nsys, ngrowthbound, X)
    return contsys
end

function build_input()
    U = AB.HyperRectangle(SVector(-1.0,-80/180*π), SVector(1.0,80/180*π))
    x0 = SVector(-1.0,0.0); hu = SVector(2.0,0.2)
    Ugrid = AB.GridFree(x0,hu)
    Udom = AB.DomainList(Ugrid)
    AB.add_set!(Udom, U, AB.OUTER)
    println(AB.get_ncells(Udom))
    #=for pos in AB.enum_pos(Udom)
        u = AB.get_coord_by_pos(Udom.grid,pos)
        println(u)
    end
    println(AB.get_ncells(Udom))
    println(Udom)=#
    return Udom
end
function transition_cost(x,u)
    return 1.0
end

function test2()
    println()
    X,O = build_domain()
    contsys = build_system(X)
    Udom = build_input()
    periodic = [3,4]
    periods = [2*π,2*π]
    T0 = [-π,-π]
    x = SVector(0.01,0.01,π/2.0+0.2,0.0)
    tstep = 1.8
    u = [-1.0,80/180*π]#[3.0,π/3.0]
    h = SVector(20.0, 20.0, 1.2, 2*π)
    h = SVector(1.0, 1.0, 0.05, 0.05)
    r = h/2
    rec = AB.HyperRectangle(x .- r, x .+ r)
    # post-image
    Fx = contsys.sys_map(x, u, tstep)
    Fr = contsys.growthbound_map(r, u, tstep, x; nstep=5)
    R = AB.HyperRectangle(Fx .- Fr, Fx .+ Fr)
    #post-K
    K = compute_K(X,r,u,tstep,x,nstep=6)
    n = 4; lb = zeros(n); ub = zeros(n)
    for i=1:n
        lb[i] = K[i].lo
        ub[i] = K[i].hi
    end
    K = AB.HyperRectangle(lb,ub)

    #pre-K
    K2 = compute_K(X,r,u,-tstep,x,nstep=5)
    n = 4; lb = zeros(n); ub = zeros(n)
    for i=1:n
        lb[i] = K2[i].lo
        ub[i] = K2[i].hi
    end
    K2 = AB.HyperRectangle(lb,ub)
    #reachable set
    RS = compute_reachable_set(rec,contsys,Udom)
    # pre-image
    u2 = SVector(-u[1],u[2])
    Fx2 = contsys.sys_map(x, u2, tstep)
    Fr2 = contsys.growthbound_map(r, u2, tstep, x; nstep=5)
    R2 = AB.HyperRectangle(Fx2 .- Fr2, Fx2 .+ Fr2)

    fig = plot(aspect_ratio = 1,legend = false)
    plot!(U.rectangle(RS.lb[1:2],RS.ub[1:2]),color=:black)
    plot!(U.rectangle(K.lb[1:2],K.ub[1:2]),color=:green)
    plot!(U.rectangle(K2.lb[1:2],K2.ub[1:2]),color=:yellow)
    plot!(U.rectangle(rec.lb[1:2],rec.ub[1:2]),color=:blue)
    plot!(U.rectangle(R.lb[1:2],R.ub[1:2]),color=:red)
    plot!(U.rectangle(R2.lb[1:2],R2.ub[1:2]),color=:pink)
    display(fig)
    fig = plot(aspect_ratio = 1,legend = false)
    plot!(U.rectangle(RS.lb[3:4],RS.ub[3:4]),color=:black)
    plot!(U.rectangle(K.lb[3:4],K.ub[3:4]),color=:green)
    plot!(U.rectangle(K2.lb[3:4],K2.ub[3:4]),color=:yellow)
    plot!(U.rectangle(rec.lb[3:4],rec.ub[3:4]),color=:blue)
    plot!(U.rectangle(R.lb[3:4],R.ub[3:4]),color=:red)
    plot!(U.rectangle(R2.lb[3:4],R2.ub[3:4]),color=:pink)
    display(fig)
end

function test()

    println("start")
    X,O = build_domain()
    contsys = build_system(X)

    Udom = build_input()

    # control problem
    x0 = SVector(30.0, 10.0,π/2.0,0.0)
    RI = SVector(1.0, 1.0, 0.05, 0.05)
    _I_ = AB.HyperRectangle(SVector(30.0, 5.0,π/2.0,0.0)-RI, SVector(30.0, 5.0,π/2.0,0.0)+RI)
    _T_ = AB.HyperRectangle(SVector(25.0, 12.0,π/2.0-0.1,-0.1), SVector(35.0, 19.0,π/2.0+0.1,0.1))
    # problem-specific functions
    functions = (compute_reachable_set,minimum_transition_cost,post_image,pre_image)

    hx_coarse = [20.0, 20.0, 1.2, 2*π]
    hx_medium = [5.0, 5.0, 0.4, 0.4]
    hx_fine = [1.0, 1.0, 0.05, 0.05]
    periodic = [3,4]
    periods = [2*π,2*π]
    T0 = [-π,-π]

    fig = plot(aspect_ratio = 1,legend = false)
    plot!(U.rectangle(X.lb,X.ub), opacity=.8,color=:blue)
    plot!(U.rectangle(_I_.lb,_I_.ub), opacity=.8,color=:green)
    plot!(U.rectangle(_T_.lb,_T_.ub), opacity=.8,color=:red)
    for o in O
        plot!(U.rectangle(o.lb,o.ub), opacity=.8,color=:black)
    end
    display(fig)

    fig = plot(aspect_ratio = 1,legend = false)
    plot!(U.rectangle(X.lb[3:4],X.ub[3:4]), opacity=.8,color=:blue)
    plot!(U.rectangle(_I_.lb[3:4],_I_.ub[3:4]), opacity=.8,color=:green)
    plot!(U.rectangle(_T_.lb[3:4],_T_.ub[3:4]), opacity=.8,color=:red)
    display(fig)
    optimal_control_prob = OC.OptimalControlProblem(x0,_I_,_T_,contsys,periodic,periods,T0,Udom,transition_cost,(X,hx_coarse,O),hx_medium,hx_fine,functions)
    fig = plot(aspect_ratio = 1,legend = false)
    max_iter = 1
    max_time = 10000
    optimizer = BB.Optimizer(optimal_control_prob,max_iter,max_time,log_level=2)
    println("optimize")
    MOI.optimize!(optimizer)
    println(optimizer.status)
    println(optimizer.best_sol)
    #display(fig)
    (traj,cost,sucess) = OC.simulate_trajectory(optimal_control_prob, optimizer.best_sol)
    println(sucess)
    #Xdom = optimal_control_prob.coarse_abstraction.Xdom
    #plot_articulted_vehicle_traj(traj;Xdom=Xdom)
end


test2()

end # end module



#=
function pre_image(symmodel,contsys,xpos,u)
    Xdom = symmodel.Xdom
    x = AB.get_coord_by_pos(Xdom.grid, xpos)
    tstep = contsys.tstep
    X = AB.HyperRectangle(SVector(0.0,0.0,-π,-π/3.0),SVector(6.0,6.0,π,π/3.0))
    r = Xdom.grid.h/2.0 + contsys.measnoise
    K = compute_K(X,r,u,-tstep,x,nstep=5)
    n = length(xpos); lb = zeros(n); ub = zeros(n)
    for i=1:n
        lb[i] = K[i].lo
        ub[i] = K[i].hi
    end
    R = AB.HyperRectangle(lb,ub)
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
end=#
