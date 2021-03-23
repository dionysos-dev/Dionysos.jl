function NewControlSystemGrowthLx(tstep, F_sys, sysnoise::SVector{N,T}, measnoise::SVector{N,T}, L_growthbound, ngrowthbound, X) where {N,T}
    # Use `let` following the advice of
    # https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-captured
    L_growthbound = let f = f
        (u, K) -> map(ForwardDiff.jacobian(f, K)) do L
            max(abs(L.lo), abs(L.hi))
        end
    end
    sys_map = let nsys = nsys
        (x, u, tstep) ->
            AB.RungeKutta4(F_sys, x, u, tstep, nsys)
    end
    growthbound_map = let ngrowthbound = ngrowthbound
        let X = X
        (r::SVector{N,T}, u, tstep, x) ->
            K = compute_K(X,r,u,tstep,x)
            L = L_growthbound(u, K)
            function F_growthbound(r, u)
                return L*r + sysnoise
            end
            return AB.RungeKutta4(F_growthbound, r, u, tstep, ngrowthbound)
        end
    end
    return AB.ControlSystemGrowth(tstep, sysnoise, measnoise, sys_map, growthbound_map)
end
