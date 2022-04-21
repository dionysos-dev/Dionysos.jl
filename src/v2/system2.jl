using LinearAlgebra
function HD_f()
    function f(x::SVector{N,T}, u) where {N,T}
        #x = inv(R)*x
        μ = 2.0
        y1 = x[2]
        y2 = μ*(1-x[1]*x[1])*x[2]-x[1]
        y3 = -2*x[3]
        y4 = -4*x[4]
        y5 = -x[5]
        y = SVector(y1,y2,y3,y4,y5)
        #R*y
        return y
    end
    return f
end

function Givens_transformation(n,d1,d2,θ)
    A =  Matrix(1.0I, n, n)
    A[d1,d1] = cos(θ)
    A[d1,d2] = -sin(θ)
    A[d2,d1] = sin(θ)
    A[d2,d2] = cos(θ)
    return SMatrix{n,n}(A)
end

function HD_f_with_boundary(f1,radius,a)
    θ = 0.8 #1.2
    R1 = Givens_transformation(5,1,2,θ)
    R2 = Givens_transformation(5,1,3,0.2)
    R3 = Givens_transformation(5,3,4,0.2)
    R4 = Givens_transformation(5,4,5,0.2)
    R = R4*R3*R2*R1
    Rinv = inv(R)
    function f(x::SVector{N,T}, u) where {N,T}
        x = Rinv*x
        y = f1(x,u)
        r = sqrt(x[1]*x[1]+x[2]*x[2])
        if r > radius
            y1 = a*x[1]
            y2 = a*x[2]
            y3 = y[3]
            y4 = y[4]
            y5 = y[5]
            y = SVector(y1,y2,y3,y4,y5)
        end
        y = R*y
        # r = sqrt(x[1]*x[1]+x[2]*x[2])
        # if r > radius
        #     y1 = a*x[1]
        #     y2 = a*x[2]
        #     y3 = y[3]
        #     y4 = y[4]
        #     y5 = y[5]
        #     y = SVector(y1,y2,y3,y4,y5)
        # end
        return y
    end
    return f
end

struct SimpleSystem{N,T,F<:Function,F2<:Function} <: AB.ControlSystem{N,T}
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
