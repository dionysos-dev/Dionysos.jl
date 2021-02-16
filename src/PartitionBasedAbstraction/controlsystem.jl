import .Abstraction.ControlSystem

struct ControlSystemLipschitz{N,T,F<:Function} <: ControlSystem{N,T}
    tstep::Float64
    measnoise::SVector{N,T}
    sys_map::F
    L::T
end

function NewControlSystemLipschitzRK4(tstep, F_sys, L,measnoise::SVector{N,T}, nsys) where {N,T}
    sys_map = let nsys = nsys
        (x::SVector{N,T}, u, tstep) ->
            AB.RungeKutta4(F_sys, x, u, tstep, nsys)::SVector{N,T}
    end
    return ControlSystemLipschitz(tstep, measnoise, sys_map, L)
end
