module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import LazySets
import Random

# A growth bound is the one input whose being wrong invalidates everything downstream without
# ever failing loudly: the abstraction simply stops covering the real dynamics, and the
# "guarantee" it produces is worthless. The bounds in `problems/` are written by hand, and
# `SMatrix` fills column by column — so an off-diagonal entry landing in the transposed slot
# type-checks, runs, and silently under-bounds a whole direction.
#
# The contract (`src/system/approximation/growth_bound.jl`) is `ṙ = L(u) · r`, which needs
#
#     L[i,j] >= |J[i,j]|   off the diagonal
#     L[i,i] >=  J[i,i]    on it (a contracting direction may legitimately be negative)
#
# so this samples the state/input domain, takes the Jacobian by central differences, and
# checks exactly that. Sampling can refute a bound but never certify one; the point here is
# refutation, which is what catches a transposition.

const PROBLEMS = joinpath(dirname(dirname(pathof(Dionysos))), "problems")

function numerical_jacobian(f, x, u; h = 1e-6)
    n = length(x)
    J = zeros(length(f(x, u)), n)
    for j in 1:n
        e = zeros(n)
        e[j] = h
        xp = SVector{n}(collect(x) .+ e)
        xm = SVector{n}(collect(x) .- e)
        J[:, j] = (collect(f(xp, u)) .- collect(f(xm, u))) ./ (2h)
    end
    return J
end

"Largest amount by which the sampled Jacobian exceeds `L`. Positive means unsound."
function worst_violation(f, L, X, U; nsample = 2000, seed = 42)
    lo, hi = LazySets.low(X), LazySets.high(X)
    ulo, uhi = LazySets.low(U), LazySets.high(U)
    n = length(lo)
    rng = Random.MersenneTwister(seed)

    worst = -Inf
    for _ in 1:nsample
        x = SVector{n}(lo .+ rand(rng, n) .* (hi .- lo))
        u = ulo .+ rand(rng, length(ulo)) .* (uhi .- ulo)
        Lm = Matrix(L(u))
        J = numerical_jacobian(f, x, u)
        for i in 1:size(J, 1), j in 1:size(J, 2)
            need = i == j ? J[i, j] : abs(J[i, j])
            worst = max(worst, need - Lm[i, j])
        end
    end
    return worst
end

# Finite differences carry a little noise of their own, so allow a hair of slack.
const TOL = 1e-5

@testset "problem growth bounds dominate their dynamics" begin
    @testset "Integrator" begin
        include(joinpath(PROBLEMS, "Integrator", "integrator.jl"))
        s = Integrator.system()
        @test worst_violation(
            Integrator.dynamic(),
            Integrator.jacobian_bound(),
            s.X,
            s.U,
        ) <= TOL
    end

    @testset "Pendulum/simple" begin
        include(joinpath(PROBLEMS, "Pendulum", "simple_pendulum.jl"))
        s = SimplePendulum.system()
        @test worst_violation(
            SimplePendulum.dynamic(),
            SimplePendulum.jacobian_bound(),
            s.X,
            s.U,
        ) <= TOL
    end

    @testset "Pendulum/double" begin
        include(joinpath(PROBLEMS, "Pendulum", "double_pendulum.jl"))
        s = DoublePendulum.system()
        @test worst_violation(
            DoublePendulum.dynamic(),
            DoublePendulum.jacobian_bound(),
            s.X,
            s.U,
        ) <= TOL
    end

    @testset "CartPendulum" begin
        include(joinpath(PROBLEMS, "CartPendulum", "cart_pendulum.jl"))
        s = CartPendulum.system()
        @test worst_violation(
            CartPendulum.dynamic(),
            CartPendulum.jacobian_bound(),
            s.X,
            s.U,
        ) <= TOL
    end

    @testset "PathPlanning" begin
        include(joinpath(PROBLEMS, "PathPlanning", "path_planning.jl"))
        s = PathPlanning.system(
            LazySets.Hyperrectangle(;
                low = SVector(0.0, 0.0, -pi - 0.4),
                high = SVector(4.0, 10.0, pi + 0.4),
            ),
        )
        @test worst_violation(
            PathPlanning.dynamic(),
            PathPlanning.jacobian_bound(),
            s.X,
            s.U,
        ) <= TOL
    end

    @testset "VectoredThrustAircraft" begin
        include(joinpath(PROBLEMS, "VectoredThrustAircraft", "vectored_thrust_aircraft.jl"))
        s = VectoredThrustAircraft.system()
        @test worst_violation(
            VectoredThrustAircraft.dynamic(),
            VectoredThrustAircraft.jacobian_bound(),
            s.X,
            s.U,
        ) <= TOL
    end

    @testset "ArticulatedVehicle" begin
        include(joinpath(PROBLEMS, "ArticulatedVehicle", "articulated_vehicle.jl"))
        s = ArticulatedVehicle.system(
            LazySets.Hyperrectangle(;
                low = SVector(-5.0, -5.0, -pi, -pi),
                high = SVector(5.0, 5.0, pi, pi),
            ),
        )
        @test worst_violation(
            ArticulatedVehicle.dynamic(),
            ArticulatedVehicle.jacobian_bound(),
            s.X,
            s.U,
        ) <= TOL
        # `jacobian()` feeds the LINEARIZED approximation mode and had the same transposition,
        # so pin it against the numerical Jacobian directly.
        J = ArticulatedVehicle.jacobian()
        x = SVector(1.0, -2.0, 0.7, -0.4)
        u = SVector(0.8, 0.3)
        @test Matrix(J(x, u)) ≈ numerical_jacobian(ArticulatedVehicle.dynamic(), x, u) atol =
            1e-5
    end

    @testset "DCDC" begin
        include(joinpath(PROBLEMS, "DCDC", "dcdc_converter.jl"))
        s = DCDC.system()
        @test worst_violation(DCDC.dynamic(), DCDC.jacobian_bound(), s.X, s.U) <= TOL
    end

    @testset "Thermostat" begin
        include(joinpath(PROBLEMS, "Thermostat", "thermostat_system.jl"))
        s = ThermostatSystem.system()
        @test worst_violation(
            ThermostatSystem.dynamic(),
            ThermostatSystem.jacobian_bound(),
            s.X,
            s.U,
        ) <= TOL

        include(joinpath(PROBLEMS, "Thermostat", "thermostat_time_system.jl"))
        st = ThermostatTimeSystem.system()
        @test worst_violation(
            ThermostatTimeSystem.dynamic(),
            ThermostatTimeSystem.jacobian_bound(),
            st.X,
            st.U,
        ) <= TOL
    end
end

# The double pendulum is conservative with no input, so its energy is a second, independent
# check on the *dynamics* rather than the bound — the missing parentheses that used to divide
# only the gravity term made it gain 30x its own initial energy out of nothing.
@testset "double pendulum conserves energy" begin
    include(joinpath(PROBLEMS, "Pendulum", "double_pendulum.jl"))
    P = DoublePendulum.Params()
    f = DoublePendulum.dynamic(P)

    function energy(x)
        θ1, θ2, ω1, ω2 = x
        T =
            0.5 * P.m1 * P.l1^2 * ω1^2 +
            0.5 *
            P.m2 *
            (P.l1^2 * ω1^2 + P.l2^2 * ω2^2 + 2 * P.l1 * P.l2 * ω1 * ω2 * cos(θ1 - θ2))
        V = -(P.m1 + P.m2) * P.g * P.l1 * cos(θ1) - P.m2 * P.g * P.l2 * cos(θ2)
        return T + V
    end

    function step(x, dt)
        z = SVector(0.0)
        k1 = f(x, z)
        k2 = f(x + dt / 2 * k1, z)
        k3 = f(x + dt / 2 * k2, z)
        k4 = f(x + dt * k3, z)
        return x + dt / 6 * (k1 + 2k2 + 2k3 + k4)
    end

    for x0 in (SVector(1.0, 0.5, 0.0, 0.0), SVector(2.0, -1.0, 0.5, -0.3))
        x = x0
        E0 = energy(x)
        drift = 0.0
        for _ in 1:20_000
            x = step(x, 1e-3)
            drift = max(drift, abs(energy(x) - E0))
        end
        @test drift < 1e-6
    end
end

end # module TestMain
