struct ControlSystem{N,T,F1<:Function,F2<:Function}
    tstep::Float64
    sysnoise::SVector{N,T}
    measnoise::SVector{N,T}
    sys_map::F1
    bound_map::F2
end

function RungeKutta4(F, x, u, tstep, nsub::Int)
    τ = tstep/nsub
    for i = 1:nsub
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

function NewControlSystemRK4(tstep, F_sys, L_bound, sysnoise::SVector{N,T},
        measnoise::SVector{N,T}, nsys, nbound) where {N,T}
    sys_map = let nsys = nsys
        (x, u, tstep) -> RungeKutta4(F_sys, x, u, tstep, nsys)
    end
    F_bound = let sysnoise = sysnoise
        (r, u) -> L_bound(u)*r + sysnoise
    end
    bound_map = let nbound = nbound
        (r, u, tstep) -> RungeKutta4(F_bound, r, u, tstep, nbound)
    end
    return ControlSystem(tstep, sysnoise, measnoise, sys_map, bound_map)
end

struct ControlSystemLin{N,T,F1<:Function,F2<:Function,F3<:Function}
    tstep::Float64
    sys_map::F1
    measnoise::SVector{N,T}
    linsys_map::F2
    bound_map::F3
end

function RungeKutta4Lin(F, DF, x, dx, u, tstep, nsub::Int)
    τ = tstep/nsub
    for i = 1:nsub
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
function BoundExp(a, b, tstep)
    ρ = exp(a*tstep)
    return b/a*ρ*(ρ - 1.0)
end

function NewControlSystemRK4(tstep, F_sys, DF_sys, DF_bound, DDF_bound,
        measenoise::SVector{N,T}, nsys, nbound) where {N,T}
    sys_map = let nsys = nsys
        (x, u, tstep) -> RungeKutta4(F_sys, x, u, tstep, nsys)
    end
    linsys_map = let nsys = nsys
        (x, dx, u, tstep) -> RungeKutta4Lin(F_sys, x, dx, u, tstep, nsys)
    end
    bound_map = (r, u, tstep) -> BoundExp(DDF_bound(u), DF_bound(u), tstep)*r
    return ControlSystemLin(tstep, sysnoise, measnoise, sys_map, bound_map)
end
