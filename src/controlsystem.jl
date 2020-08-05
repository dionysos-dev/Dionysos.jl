struct ControlSystem{N,T,Fsys<:Function,Fbound<:Function}
    tstep::Float64
    sysnoise::SVector{N,T}
    measnoise::SVector{N,T}
    sys_map::Fsys
    bound_map::Fbound
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
    sys_map = (x, u, tstep) -> RungeKutta4(F_sys, x, u, tstep, nsys)
    F_bound = (r, u) -> L_bound(u)*r + sysnoise
    bound_map = (r, u, tstep) -> RungeKutta4(F_bound, r, u, tstep, nbound)
    return ControlSystem(tstep, sysnoise, measnoise, sys_map, bound_map)
end
