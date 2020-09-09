"""
    ControlSystem{N,T}

Abstract type describing a control system on the space of N-dimensional vectors with
elements of type `T`.

Typical methods include:

* Compute the image of a point for a given input and after a given time step
* Compute an over-approximation of the reachable set of a given region.
   For the moment, over-approximations methods can be computed:
   + Via a growth bound function as in https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=7519063
   + Via a first-order approximation as in https://drive.google.com/drive/folders/1Czz9RnDp8r8DRRJHIMOt3jvGIMwLpsem
"""
abstract type ControlSystem{N,T} end

xspacetype(::ControlSystem{N,T}) where {N,T} = (N, T)

struct ControlSystemGrowth{N,T,F1<:Function,F2<:Function} <: ControlSystem{N,T}
    tstep::Float64
    sysnoise::SVector{N,T}
    measnoise::SVector{N,T}
    sys_map::F1
    growthbound_map::F2
end

# Use integer litteral; see https://docs.julialang.org/en/v1/manual/style-guide/#Avoid-using-floats-for-numeric-literals-in-generic-code-when-possible-1
function _RK4(F, x, u, tstep, nsub::Int)
    τ = tstep/nsub
    for i in 1:nsub
        Fx1 = F(x, u)
        xrk = x + Fx1*(τ/2)
        Fx2 = F(xrk, u)
        xrk = x + Fx2*(τ/2)
        Fx3 = F(xrk, u)
        xrk = x + Fx3*τ
        Fx4 = F(xrk, u)
        x = x + (Fx1 + Fx2*2 + Fx3*2 + Fx4)*(τ/6)
    end
    return x
end

function ControlSystemGrowthRK4(
        tstep,
        F_sys, L_growthbound,
        sysnoise::SVector{N,T}, measnoise::SVector{N,T},
        nsys, ngrowthbound
    ) where {N,T}
    sys_map = let nsys = nsys
        (x::SVector{N,T}, u, tstep) ->
            _RK4(F_sys, x, u, tstep, nsys)::SVector{N,T}
    end
    F_growthbound = let sysnoise = sysnoise
        (r, u) -> L_growthbound(u)*r + sysnoise
    end
    growthbound_map = let ngrowthbound = ngrowthbound
        (r::SVector{N,T}, u, tstep) ->
            _RK4(F_growthbound, r, u, tstep, ngrowthbound)::SVector{N,T}
    end
    return ControlSystemGrowth(tstep, sysnoise, measnoise, sys_map, growthbound_map)
end

struct ControlSystemLinearized{
        N,T,F1<:Function,F2<:Function,F3<:Function} <: ControlSystem{N,T}
    tstep::Float64
    measnoise::SVector{N,T}
    sys_map::F1
    linsys_map::F2
    error_map::F3
end

function _RK4Linearized(F, DF, x, dx, u, tstep, nsub::Int)
    τ = tstep/nsub
    for i in 1:nsub
        Fx1 = F(x, u)
        DFx1 = DF(x, u)*dx
        xrk = x + Fx1*(τ/2)
        dxrk = dx + DFx1*(τ/2)
        Fx2 = F(xrk, u)
        DFx2 = DF(xrk, u)*dxrk
        xrk = x + Fx2*(τ/2)
        dxrk = dx + DFx2*(τ/2)
        Fx3 = F(xrk, u)
        DFx3 = DF(xrk, u)*dxrk
        xrk = x + Fx3*τ
        dxrk = dx + DFx3*τ
        Fx4 = F(xrk, u)
        DFx4 = DF(xrk, u)*dxrk
        x += (Fx1 + Fx2*2 + Fx3*2 + Fx4)*(τ/6)
        dx += (DFx1 + DFx2*2 + DFx3*2 + DFx4)*(τ/6)
    end
    return (x, dx)
end

# Give an upper-bound on x(t) staisfying x'(t) ≦ a*x(t) + b*exp(2at)
function _bound_2e_order(a, b, tstep)
    if a ≈ 0.0
        return b*tstep
    else
        ρ = exp(a*tstep)
        return (b/a)*ρ*(ρ - 1)/2
    end
end

function ControlSystemLinearizedRK4(
        tstep,
        F_sys, DF_sys, bound_DF, bound_DDF,
        measnoise::SVector{N,T},
        nsys
        ) where {N,T}
    sys_map = let nsys = nsys
        (x::SVector{N,T}, u, tstep) ->
            _RK4(F_sys, x, u, tstep, nsys)::SVector{N,T}
    end
    linsys_map = let nsys = nsys
        (x::SVector{N,T}, dx::SMatrix{N,N,T}, u, tstep) ->
            _RK4Linearized(
                F_sys, DF_sys, x, dx, u, tstep, nsys)::Tuple{SVector{N,T},SMatrix{N,N,T}}
    end
    error_map = (r::T, u, tstep) ->
        _bound_2e_order(bound_DF(u), bound_DDF(u), tstep)*(r*r)::T
    return ControlSystemLinearized(tstep, measnoise, sys_map, linsys_map, error_map)
end
