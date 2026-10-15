# Tests for BipedFourLink.jl (simulation and visualization with RigidBodyDynamics). Run with an
# environment that has RigidBodyDynamics, MeshCatMechanisms, MeshCat, Plots, Dionysos, Symbolics and
# the packages of robot_problem.jl:
#   julia --project=<env> tests/test_biped_four_link.jl
using Test
using StaticArrays
using Symbolics
import Dionysos
include(joinpath(@__DIR__, "..", "BipedFourLink.jl"))
using .BipedFourLink
const RBD = BipedFourLink.RBD
include(joinpath(@__DIR__, "..", "robot_problem.jl"))
import .UnderactuatedRobotProblem as RP
const MB = RP.Multibody2D

model = build_model(joinpath(@__DIR__, "..", "four_link_biped.urdf"))
zero_torque!(u, t, q, q̇) = (u .= 0)

@testset "BipedFourLink" begin
    @testset "state conversion round trip, for both mechanisms" begin
        q, q̇ = [0.1, -0.2, 0.3, -0.4], [0.5, 0.1, -0.2, 0.3]
        qr, q̇r = coordinates_from_state(model, state_from_coordinates(model, q, q̇))
        @test qr ≈ q && q̇r ≈ q̇
        compliant = BipedFourLink._compliant_mechanism(model, 0.8)
        state = state_from_coordinates(model, q, q̇; mechanism = compliant)
        qc, q̇c = coordinates_from_state(model, state)
        @test qc ≈ q && q̇c ≈ q̇
        @test RBD.num_positions(compliant) == 6                  # [x, z, q₁…q₄]
        @test_throws DimensionMismatch state_from_coordinates(model, [0.0, 0.0, 0.0])
    end

    # The URDF angles are clockwise, those of the design model counter-clockwise: the same posture
    # has `q_design = -q_urdf`, and the world positions are the same.
    @testset "forward kinematics agree with the design model" begin
        geometry = RP.default_geometry()
        for k in 1:5
            q = [0.7 * sin(12.9898 * k + 78.233 * j) for j in 1:4]
            state = state_from_coordinates(model, q)
            design = RP.cartesian_coordinates(SVector{4}(-q), geometry)
            xz(p) = [p[1], p[3]]                                 # (x, z) of a world position
            @test xz(foot_position(model, state; foot = :grounded)) ≈ [0.0, 0.0] atol =
                1e-12
            @test xz(world_position(model, state, "grounded_knee")) ≈ design[2] atol = 1e-12
            @test xz(world_position(model, state, "hip")) ≈ design[3] atol = 1e-12
            @test xz(world_position(model, state, "swing_knee")) ≈ design[4] atol = 1e-12
            @test xz(foot_position(model, state; foot = :swing)) ≈ design[5] atol = 1e-12
        end
        @test_throws ArgumentError foot_position(
            model,
            state_from_coordinates(model, zeros(4));
            foot = :left,
        )
    end

    # The whole-robot simulation (RigidBodyDynamics, torques as inputs, URDF, angles clockwise) and
    # the symbolic design model (Multibody2D, counter-clockwise) describe the same mechanism: the
    # torque-level accelerations agree once `q`, `q̇`, `τ` and `q̈` change sign.
    @testset "torque-level dynamics agree with Multibody2D" begin
        bm = RP.symbolic_model()
        Dsym = RP.substitute_parameters(bm, MB.mass_matrix(bm.sys))
        hsym = RP.substitute_parameters(bm, MB.bias_vector(bm.sys))
        ev(e, d) = Float64(Symbolics.value(Symbolics.substitute(e, d; fold = Val(true))))
        B = [0 0 0; 1 0 0; 0 1 0; 0 0 1]
        for k in 1:10
            q = [0.5 * sin(1.3k + j) for j in 1:4]
            q̇ = [cos(0.7k + j) for j in 1:4]
            u = [0.3 * sin(2k + j) for j in 1:3]
            d = merge(Dict(zip(bm.sys.q, -q)), Dict(zip(bm.sys.qd, -q̇)))
            D = reshape([ev(e, d) for e in Dsym], 4, 4)
            h = [ev(e, d) for e in hsym]
            q̈_multibody = -(D \ (h + B * (-u)))
            @test accelerations(model, q, q̇, u) ≈ q̈_multibody rtol = 1e-7 atol = 1e-8
        end
    end

    @testset "locked-ankle simulation conserves energy without torques" begin
        q0, q̇0 = [0.1, -0.2, 0.3, -0.4], [0.5, 0.1, -0.2, 0.3]
        ts, qh, q̇h, qs, vs = simulate_trajectory(
            model,
            q0,
            q̇0,
            0.2,
            zero_torque!;
            foot_model = :locked,
            Δt = 1e-3,
        )
        @test length(ts) == length(qh) == length(q̇h) == length(qs) > 100
        @test qh[1] ≈ q0 && q̇h[1] ≈ q̇0
        function energy(q, v)
            s = RBD.MechanismState(model.mechanism)
            RBD.set_configuration!(s, q)
            RBD.set_velocity!(s, v)
            return RBD.kinetic_energy(s) + RBD.gravitational_potential_energy(s)
        end
        E = [energy(q, v) for (q, v) in zip(qs, vs)]
        @test maximum(abs.(E .- E[1])) < 1e-8 * abs(E[1])
        # A constant torque on the hip changes the motion.
        _, qh_τ, = simulate_trajectory(
            model,
            q0,
            q̇0,
            0.2,
            (u, t, q, q̇) -> (u .= (0.0, 0.5, 0.0));
            Δt = 1e-3,
        )
        @test qh_τ[end] != qh[end]
    end

    @testset "compliant-foot simulation and contact status" begin
        ts, qh, q̇h, qs, vs, status = simulate_trajectory(
            model,
            zeros(4),
            zeros(4),
            0.05,
            zero_torque!;
            foot_model = :compliant,
            Δt = 1e-4,
            return_contact_status = true,
            contact_foot = :both,
        )
        @test length(qh[1]) == 4 && length(qs[1]) == 6        # [x, z, q₁…q₄] in the full state
        @test keys(status) == (:grounded, :swing)
        @test status.grounded.was_in_contact
        @test contact_slip_status(model, qs; contact_foot = :grounded).was_in_contact
    end

    @testset "reaction force, friction cone, plots, MeshCat" begin
        q, q̇ = [0.1, -0.2, 0.3, -0.4], [0.5, 0.1, -0.2, 0.3]
        F = ankle_reaction_force(model, q, q̇, accelerations(model, q, q̇, zeros(3)))
        @test length(F) == 3 && F[3] > 0                       # the ground pushes up
        @test in_friction_cone([0.1, 0.0, 1.0], 0.8) &&
              !in_friction_cone([1.0, 0.0, 1.0], 0.8)
        @test !in_friction_cone([0.0, 0.0, -1.0], 0.8)
        ts = collect(0.0:0.01:0.1)
        qh = [q for _ in ts]
        q̇h = [q̇ for _ in ts]
        plots = plot_trajectory(ts, qh, q̇h)
        @test plots.forces === nothing && plots.positions !== nothing
        qs = [RBD.configuration(state_from_coordinates(model, q, q̇)) for _ in ts]
        @test animate_trajectory(model, ts, qs) !== nothing     # builds the MeshCat visualizer
    end
end
