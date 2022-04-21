
function Van_der_pol_oscillator(μ;k=1)
    function f(x::SVector{N,T}, u) where {N,T}
        x = (1/k)*x
        y1 = x[2]
        y2 = μ*(1-x[1]*x[1])*x[2]-x[1]
        y = SVector(y1,y2)
        return y
    end
    return f
end

function reverse_Van_der_pol_oscillator(μ;k=1)
    function f(x::SVector{N,T}, u) where {N,T}
        x = (1/k)*x
        y1 = x[2]
        y2 = μ*(1-x[1]*x[1])*x[2]-x[1]
        y = -SVector(y1,y2)
        return y
    end
    return f
end

function Van_der_pol_oscillator_Jacobian(μ)
    function Jacobian(x, u) where {N,T}
        M = @SMatrix [ 0 1 ;
                      -2*μ*x[1]*x[2]-1 μ*(1-x[1]*x[1])]
        return -M
    end
    return Jacobian
end

function system_with_boundary(f1,radius,a)
    function f(x::SVector{N,T}, u) where {N,T}
        y = f1(x,u)
        if norm(x,2)>radius
            A = @SMatrix [ a 0.0 ;
                           0.0 a ]
            y = A*x
        end
        return y
    end
    return f
end


function Jacobian_with_boundary(J1,radius,a)
    function Jacobian(x, u) where {N,T}
        J = J1(x,u)
        # if norm(x,2)>radius
        #     J = @SMatrix [ a 0.0 ;
        #                    0.0 a ]
        # end
        return J
    end
    return Jacobian
end

struct SimpleSystem{N,T,F<:Function,F2<:Function} <: ST.ControlSystem{N,T}
    tstep::Float64
    measnoise::SVector{N,T}
    sys_map::F
    f::F2
end

function RungeKutta4(F, x, u, tstep, nsub::Int)
    τ = tstep/nsub
    for i in 1:nsub
        Fx1 = F(x, u)
        xrk = x + Fx1*(τ/2.0)
        Fx2 = F(xrk, u)
        xrk = x + Fx2*(τ/2.0)
        Fx3 = F(xrk, u)
        xrk = x + Fx3*τ
        Fx4 = F(xrk, u)
        x = x + (Fx1 + Fx2*2.0 + Fx3*2.0 + Fx4)*(τ/6.0)
    end
    return x
end

function NewControlSystemGrowthRK4(tstep, F_sys, measnoise::SVector{N,T}, nsys) where {N,T}
    sys_map = let nsys = nsys
        (x::SVector{N,T}, u, tstep) ->
            RungeKutta4(F_sys, x, u, tstep, nsys)::SVector{N,T}
    end
    return SimpleSystem(tstep, measnoise, sys_map, F_sys)
end

function build_system_Van_der_pol_oscillator(μ = 2.0,tstep = 0.08,measnoise = SVector(0.0, 0.0),nsys = 4)
    function F_sys(x::SVector{N,T}, u) where {N,T}
        y1 = x[2]
        y2 = μ*(1-x[1]*x[1])*x[2]-x[1]
        y = SVector(y1,y2)
        if norm(x,2)>8.0
            a = -0.8
            b = 0.0
            A = @SMatrix [ a -b ;
                           b a ]
            y = A*x
        end
        return y
    end
    return NewControlSystemGrowthRK4(tstep,F_sys,measnoise,nsys)
end

struct ControlSystemGrowth{N,T,F1<:Function,F2<:Function,F3<:Function,F4<:Function} <: ST.ControlSystem{N,T}
    tstep::Float64
    sysnoise::SVector{N,T}
    measnoise::SVector{N,T}
    sys_map::F1
    reverse_sys_map::F2
    growthbound_map::F3
    reverse_growthbound_map::F4
end

function NewControlSystemGrowthLx(tstep, f, jacobian,f2, measnoise::SVector{N,T}, sysnoise, nsys, ngrowthbound, X) where {N,T}
    X = UT.HyperRectangle(SVector(-7.0, -9.0), SVector(7.0, 9.0))
    function L_growthbound(u,K)
        #L = ForwardDiff.jacobian(x->f(x,u),K)
        L = jacobian(K,u)
        M = @SMatrix [i == j ? max(L[i,j].lo,L[i,j].hi) : max(abs(L[i,j].lo), abs(L[i,j].hi)) for i in 1:2, j in 1:2]
        return M
    end
    sys_map = let nsys = nsys
        (x::SVector{N,T}, u, tstep) ->
            ST.RungeKutta4(f, x, u, tstep, nsys)::SVector{N,T}
    end
    reverse_sys_map = let nsys = nsys
        (x::SVector{N,T}, u, tstep) ->
            ST.RungeKutta4(f2, x, u, tstep, nsys)::SVector{N,T}
    end
    function growthbound_map(r::SVector{N,T}, u, tstep, x; nstep=1)
        K = compute_K(X,r,u,tstep,x,f, nstep=nstep)
        L = L_growthbound(u, K)
        function F_growthbound(r, u)
            return L*r + sysnoise
        end
        return ST.RungeKutta4(F_growthbound, r, u, tstep, ngrowthbound)
    end
    function reverse_growthbound_map(r::SVector{N,T}, u, tstep, x; nstep=1)
        K = compute_K(X,r,u,tstep,x,f2, nstep=nstep)
        L = L_growthbound(u, K)
        function F_growthbound(r, u)
            return L*r + sysnoise
        end
        return ST.RungeKutta4(F_growthbound, r, u, tstep, ngrowthbound)
    end
    return ControlSystemGrowth(tstep, sysnoise, measnoise, sys_map, reverse_sys_map, growthbound_map, reverse_growthbound_map)
end


function compute_K(X,r,u,tstep,x,f;nstep=1)
    R = IntervalBox(-r,r)
    D = x + R + R*3.0
    K = x + R + R*3.0
    for k=1:nstep
        # println(K)
        B = f(K,u)
        # println(B)
        K = x + R + B*tstep
        K = K ∩ D
        K = K ∪ (x+R)
        # println(K)
        # println()
    end
    return K

    # D = IntervalBox(X.lb, X.ub)
    # K = IntervalBox(X.lb, X.ub)
    # R = IntervalBox(-r,r)
    # for k=1:nstep
    #     B = f(K,u)
    #     K = x + R + B*tstep
    #     K = K ∩ D
    #     K = K ∪ (x+R)
    # end
    # return K
end


###############################################################################

function f1(x)
    return x
end
function fi1(x)
    return x
end

#counter-clockwise
function rotate(x,θ)
    R = @SMatrix [ cos(θ) -sin(θ) ;
                   sin(θ)  cos(θ)]
    return R*x
end
function f2(x)
    c =  SVector(0.0,0.0) #SVector(2.0,3.0)
    θ = π/3.0
    return rotate(x-c,θ)+c
end
function fi2(x)
    c =  SVector(0.0,0.0) #SVector(2.0,3.0)
    θ = π/3.0
    return rotate(x-c,-θ)+c
end

function f3(x)
    return SVector(x[2]+sin(x[1]),x[1])
end
function fi3(x)
    return SVector(x[2],x[1]-sin(x[2]))
end

function f4(x)
    return SVector(x[1]*cos(x[2]),x[1]*sin(x[2]))
end
function fi4(x)
    return SVector(sqrt(x[1]*x[1] + x[2]*x[2]),atan(x[2],x[1]))
end

function build_f_rotation(θ,c =  SVector(0.0,0.0))
    function f(x)
        return rotate(x-c,θ)+c
    end
    function fi(x)
        return rotate(x-c,-θ)+c
    end
    return f,fi
end

function build_circular_grid(origin,h)
    free_grid = DO.GridFree(origin+h/2.0,h)
    #θ = h[2]*2.0
    return DO.DeformedGrid(free_grid,f4,fi4)
end


#######################################################
function get_entropy(transitions)
    h = 0.0
    for (t,s,u,p) in transitions
        h = h - p*log2(p)
    end
    return h
end

function plot_limit_cycle(symmodel,sys)
    l = 12
    fig = plot(aspect_ratio = 1,legend = false,title="steady-state probability", xlim=(-l, l), ylim=(-l,l))
    plot_steady_state!(symmodel,fact=0.0,color=:blue)
    x0 = SVector(2.0,2.0)
    plot_trajectory!(sys,x0,200)
    display(fig)
end
