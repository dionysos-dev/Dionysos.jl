module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

const PROBLEMS = joinpath(dirname(dirname(pathof(Dionysos))), "problems")
include(joinpath(PROBLEMS, "CartPendulum", "point_mass_cart_pendulum.jl"))
const PM = PointMassCartPendulum
const P = PM.Params()

@testset "the states of the two models correspond" begin
    for (θ, p, v) in ((0.3, -0.5, 1.2), (2.9, 1.1, -0.7), (-1.0, 0.0, 0.0))
        y = PM.force_state(P, SVector(θ, p), v)
        @test y[2] == v
        @test PM.velocity_state(P, y) ≈ SVector(θ, p)
    end
end

# The servo that the velocity model idealizes: from rest, hanging down, a step of the setpoint is
# reached within a few periods, the force clamped to its limit on the way.
@testset "the velocity loop tracks its setpoint under the force limit" begin
    f = PM.force_dynamics(P)
    dt = 1e-3
    loop = PM.velocity_loop(; Kp = 100.0, Ki = 500.0, F_max = 20.0, dt)
    mem = ST.initial_state(loop)
    y = SVector(π, 0.0, 0.0)
    peak = 0.0
    for _ in 1:300
        input = SVector(y[1], y[2], y[3], 0.5)
        F = ST.output_control(loop, mem, input)
        mem = ST.update_state(loop, mem, input)
        peak = max(peak, abs(F))
        y = ST.runge_kutta4(f, y, F, dt, 1)
    end
    @test peak ≈ 20.0
    @test abs(y[2] - 0.5) < 0.01
end

end # module
