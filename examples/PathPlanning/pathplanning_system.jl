tstep = 0.3
nsys = 5
ngrowthbound = 5
function F_sys(x, u)
    α = atan(tan(u[2])/2)
    return SVector{3}(
        u[1]*cos(α + x[3])/cos(α),
        u[1]*sin(α + x[3])/cos(α),
        u[1]*tan(u[2]))
end
function L_growthbound(u)
    β = abs(u[1]/cos(atan(tan(u[2])/2)))
    # Both have the same speed
    # return @SMatrix [
    #     0.0 0.0 β;
    #     0.0 0.0 β;
    #     0.0 0.0 0.0]
    return SMatrix{3,3}(
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        β, β, 0.0)
end
function DF_sys(x, u)
    α = atan(tan(u[2])/2)
    β = u[1]/cos(α)
    return SMatrix{3,3}(
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            -β*sin(α + x[3]), β*cos(α + x[3]), 0.0)
end
bound_DF(u) = abs(u[1]/cos(atan(tan(u[2])/2)))
# DDF_1 = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 -β*cos(α + x[3])]
# DDF_2 = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 -β*sin(α + x[3])]
# DDF_3 = 0.0
bound_DDF(u) = abs(u[1]/cos(atan(tan(u[2])/2)))
sysnoise = SVector(0.0, 0.0, 0.0)
measnoise = SVector(0.0, 0.0, 0.0)
