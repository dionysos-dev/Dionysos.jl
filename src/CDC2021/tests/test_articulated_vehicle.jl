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
    X = AB.HyperRectangle(SVector(0.0,0.0,-π,-π),SVector(5.0,5.0,π,π))
    tstep = 0.1
    sysnoise = SVector(0.0, 0.0, 0.0, 0.0)
    measnoise = SVector(0.0, 0.0, 0.0, 0.0)
    nsys = 5
    ngrowthbound = 5
    contsys = NewControlSystemGrowthLx(tstep, f, sysnoise, measnoise, nsys, ngrowthbound)
end

test()
