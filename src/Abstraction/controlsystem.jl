abstract type ControlSystem{N,T} end

struct ControlSystemGrowth{N,T,F1<:Function,F2<:Function,F3<:Function} <: ControlSystem{N,T}
    tstep::Float64
    sysnoise::SVector{N,T}
    measnoise::SVector{N,T}
    sys_map::F1
    growthbound_map::F2
    sys_inv_map::F3
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

function NewControlSystemGrowthRK4(tstep, F_sys, L_growthbound, sysnoise::SVector{N,T},
        measnoise::SVector{N,T}, nsys, ngrowthbound) where {N,T}
    sys_map = let nsys = nsys
        (x::SVector{N,T}, u, tstep) ->
            RungeKutta4(F_sys, x, u, tstep, nsys)::SVector{N,T}
    end
    sys_inv_map = let nsys = nsys
        (x::SVector{N,T}, u, tstep) ->
            RungeKutta4(F_sys, x, u, -tstep, nsys)::SVector{N,T}
    end
    F_growthbound = let sysnoise = sysnoise
        (r, u) -> L_growthbound(u)*r + sysnoise
    end
    growthbound_map = let ngrowthbound = ngrowthbound
        (r::SVector{N,T}, u, tstep) ->
            RungeKutta4(F_growthbound, r, u, tstep, ngrowthbound)::SVector{N,T}
    end
    return ControlSystemGrowth(tstep, sysnoise, measnoise, sys_map, growthbound_map,sys_inv_map)
end

struct ControlSystemLinearized{N,T,F1<:Function,F2<:Function,F3<:Function,F4<:Function} <: ControlSystem{N,T}
    tstep::Float64
    measnoise::SVector{N,T}
    sys_map::F1
    linsys_map::F2
    error_map::F3
    sys_inv_map::F4
end

function RungeKutta4Linearized(F, DF, x, dx, u, tstep, nsub::Int)
    τ = tstep/nsub
    for i in 1:nsub
        Fx1 = F(x, u)
        DFx1 = DF(x, u)*dx
        xrk = x + Fx1*(τ/2.0)
        dxrk = dx + DFx1*(τ/2.0)
        Fx2 = F(xrk, u)
        DFx2 = DF(xrk, u)*dxrk
        xrk = x + Fx2*(τ/2.0)
        dxrk = dx + DFx2*(τ/2.0)
        Fx3 = F(xrk, u)
        DFx3 = DF(xrk, u)*dxrk
        xrk = x + Fx3*τ
        dxrk = dx + DFx3*τ
        Fx4 = F(xrk, u)
        DFx4 = DF(xrk, u)*dxrk
        x += (Fx1 + Fx2*2.0 + Fx3*2.0 + Fx4)*(τ/6.0)
        dx += (DFx1 + DFx2*2.0 + DFx3*2.0 + DFx4)*(τ/6.0)
    end
    return (x, dx)
end

# Give an upper-bound on x(t) staisfying x'(t) ≦ a*x(t) + b*exp(2at)
function BoundSecondOrder(a, b, tstep)
    if a ≈ 0.0
        return b*tstep
    else
        ρ = exp(a*tstep)
        return (b/a)*ρ*(ρ - 1.0)/2.0
    end
end

function NewControlSystemLinearizedRK4(tstep, F_sys, DF_sys, bound_DF, bound_DDF,
        measnoise::SVector{N,T}, nsys) where {N,T}
    sys_map = let nsys = nsys
        (x::SVector{N,T}, u, tstep) ->
            RungeKutta4(F_sys, x, u, tstep, nsys)::SVector{N,T}
    end
    linsys_map = let nsys = nsys
        (x::SVector{N,T}, dx::SMatrix{N,N,T}, u, tstep) ->
            RungeKutta4Linearized(
                F_sys, DF_sys, x, dx, u, tstep, nsys)::Tuple{SVector{N,T},SMatrix{N,N,T}}
    end
    sys_inv_map = let nsys = nsys
        (x::SVector{N,T}, u, tstep) ->
            RungeKutta4(F_sys, x, u, -tstep, nsys)::SVector{N,T}
    end
    error_map = (r::T, u, tstep) ->
        BoundSecondOrder(bound_DF(u), bound_DDF(u), tstep)*(r^2)::T
    return ControlSystemLinearized(tstep, measnoise, sys_map, linsys_map, error_map,sys_inv_map)
end

